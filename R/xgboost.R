#' Compute Smoother Matrix for a Single Tree
#' 
#' @description Computes the smoother matrix \eqn{S} for a single decision tree.
#'
#' @param current_tree_indices_train Numeric vector. Leaf indices for the training data.
#' @param current_tree_indices_test Numeric vector. Leaf indices for the test data.
#' @param lambda Numeric. The L2 regularization parameter.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{S_train}: The smoother matrix for the training set (n_train x n_train).
#'   \item \code{S_test}: The smoother matrix for the test set (n_test x n_train).
#' }
#' 
#' @keywords internal
create_S_from_tree <- function(current_tree_indices_train,
                               current_tree_indices_test,
                               lambda) {
  S_train <- outer(current_tree_indices_train, current_tree_indices_train, FUN = "==")
  S_train <- S_train / (rowSums(S_train) + lambda)
  
  S_test <- outer(current_tree_indices_test, current_tree_indices_train, FUN = "==")
  S_test <- S_test / (rowSums(S_test) + lambda)
  
  return(list(S_train = S_train, S_test = S_test))
}


#' Create Smoother Matrix for a Single Boosted Tree
#' 
#' @description Computes the smoother matrix contributions for a specific boosted tree step,
#' applying corrections based on the previous ensemble state (\code{S_gb_prev}).
#'
#' @param current_tree_indices_train Numeric vector. Leaf indices for the training data.
#' @param current_tree_indices_test Numeric vector. Leaf indices for the test data.
#' @param S_gb_prev Matrix or NULL. The accumulated smoother matrix from previous trees.
#'                  If NULL, this is treated as the first tree.
#' @param lambda Numeric. The L2 regularization parameter.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{S_train}: The smoother matrix for the training set (n_train x n_train).
#'   \item \code{S_test}: The smoother matrix for the test set (n_test x n_train).
#' }
#' @keywords internal
create_S_from_single_boosted_tree <- function(current_tree_indices_train,
                                              current_tree_indices_test,
                                              S_gb_prev,
                                              lambda) {
  S <- create_S_from_tree(current_tree_indices_train, current_tree_indices_test, lambda)
  
  if (is.null(S_gb_prev)) {
    # first tree: just normal tree
    return(list(S_train = S$S_train, S_test = S$S_test))
  }
  
  n_train <- length(current_tree_indices_train)
  n_test <- length(current_tree_indices_test)
  
  all_nodes <- unique(current_tree_indices_train)
  n_nodes <- length(all_nodes)
  
  node_corrections <- matrix(0, nrow = n_nodes, ncol = n_train)
  S_train_correction <- matrix(0, nrow = n_train, ncol = n_train)
  S_test_correction <- matrix(0, nrow = n_test, ncol = n_train)
  
  for (i in 1:length(all_nodes)) {
    n <- all_nodes[i]
    
    # Create correction matrix
    leaf_id_train <- current_tree_indices_train == n
    node_corrections[i, ] <- colSums(S_gb_prev[leaf_id_train, , drop = FALSE]) / (sum(leaf_id_train) + lambda)
    
    S_train_correction[leaf_id_train, ] <- matrix(rep(node_corrections[i, ], sum(leaf_id_train)), nrow = sum(leaf_id_train), byrow = TRUE)
    
    leaf_id_test <- current_tree_indices_test == n
    S_test_correction[leaf_id_test, ] <- matrix(rep(node_corrections[i, ], sum(leaf_id_test)), nrow = sum(leaf_id_test), byrow = TRUE)
  }
  
  S_train <- S$S_train - S_train_correction
  S_test <- S$S_test - S_test_correction
  
  return(list(S_train = S_train, S_test = S_test))
}


#' Iterative Smoother Matrix Construction from XGBoost Regressor
#' 
#' @description Iterates through all trees in an XGBoost model to 
#' construct the final smoother matrix.
#'
#' @param model An object of class \code{xgb.Booster}.
#' @param leaf_indices_train Matrix. Leaf indices for the training data (rows=samples, cols=trees).
#' @param leaf_indices_test Matrix. Leaf indices for the test data.
#' @param base_score Numeric or character. The model’s initial prediction (0 or "mean").
#' @param output_dir Character. Directory path to save intermediate RDS files.
#' @param save_output Logical. If TRUE, saves intermediate matrices to disk.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{S_train}: The smoother matrix for the training set (n_train x n_train).
#'   \item \code{S_test}: The smoother matrix for the test set (n_test x n_train).
#' }
#' @keywords internal
create_S_from_gbtregressor <- function(model,
                                       leaf_indices_train,
                                       leaf_indices_test,
                                       base_score,
                                       save_output,
                                       output_dir) {
  if (save_output) {
    # Check if the directory exists, and create it if it doesn't
    if (!file.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }
  
  n_train <- nrow(leaf_indices_train)
  n_test <- nrow(leaf_indices_test)
  
  # RR: Check for aliases and apply a default [xgboost specificity]
  get_param <- function(params, aliases, default) {
    for (alias in aliases) {
      if (!is.null(params[[alias]])) { return(params[[alias]]) }}
    return(default)
  }
  
  # Fetch eta and lambda from model
  params <- attr(model, "params")
  eta <- get_param(params, c("eta", "learning_rate"), 0.3)
  lambda <- get_param(params, c("lambda", "reg_lambda"), 1)
  
  if (base_score == 0) {
    S_train_curr <- matrix(0, nrow = n_train, ncol = n_train)
    S_test_curr <- matrix(0, nrow = n_test, ncol = n_train)
  } else if (base_score == "mean") {
    S_train_curr <- matrix(1 / n_train, nrow = n_train, ncol = n_train)
    S_test_curr <- matrix(1 / n_train, nrow = n_test, ncol = n_train)
  }
  
  number_of_trees <- ncol(leaf_indices_train)
  
  for (col in 1:number_of_trees) {
    current_tree_indices_train <- leaf_indices_train[, col]
    current_tree_indices_test <- leaf_indices_test[, col]
    
    if (col == 1 && base_score == 0) {
      S <- create_S_from_single_boosted_tree(
        current_tree_indices_train,
        current_tree_indices_test,
        NULL,
        lambda
      )
    } else {
      S <- create_S_from_single_boosted_tree(
        current_tree_indices_train,
        current_tree_indices_test,
        S_train_curr,
        lambda
      )
    }
    
    S_train_curr <- S_train_curr + eta * S$S_train
    S_test_curr <- S_test_curr + eta * S$S_test
    
    if (save_output) {
      # Save the current matrix to an RDS file
      output_filename_train <- paste(output_dir, "/S_curr_iteration_train_", col, ".rds", sep = "")
      output_filename_test <- paste(output_dir, "/S_curr_iteration_test_", col, ".rds", sep = "")
      
      saveRDS(S_train_curr, file = output_filename_train)
      saveRDS(S_test_curr, file = output_filename_test)
      
      # Save the current tree matrix
      output_filename_train <- paste(output_dir, "/S_iteration_train_", col, ".rds", sep = "")
      output_filename_test <- paste(output_dir, "/S_iteration_test_", col, ".rds", sep = "")
      
      saveRDS(S_train_curr, file = output_filename_train)
      saveRDS(S_test_curr, file = output_filename_test)
    }
  }
  
  return(list(S_train = S_train_curr, S_test = S_test_curr))
}


#' Extract XGBoost Smoother Weights
#' 
#' @description Main function to extract the smoother matrix from a trained 
#' XGBoost model. It validates model parameters and computes 
#' the smoother matrices for training and test sets.
#'
#' @param model An object of class \code{xgb.Booster}. The trained model.
#' @param dtrain An \code{xgb.DMatrix} object. The training data used to build the model.
#' @param dtest An \code{xgb.DMatrix} object. The test data.
#' @param save_output Logical. Whether to save intermediate matrices to disk. Defaults to FALSE.
#' @param output_dir Character. Directory for saving outputs if \code{save_output} is TRUE.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{S_train}: The smoother matrix for the training set (n_train x n_train).
#'   \item \code{S_test}: The smoother matrix for the test set (n_test x n_train).
#' }
#' @keywords internal
#' 
get_xgboost_weights <- function(model,
                                dtrain,
                                dtest,
                                save_output = FALSE,
                                output_dir = NULL) {
  params <- attr(model, "params")
  
  label_mean <- mean(xgboost::getinfo(dtrain, "label"))
  base_score <- NULL

  if (!is.null(params$base_score)) {
    if (isTRUE(all.equal(params$base_score, label_mean))) {
      base_score <- "mean"
    } else if (identical(params$base_score, 0)) {
      base_score <- 0
    }
  }
  
  if (is.null(base_score)) stop("base_score must be 0 or mean(Y) for xgboost call; DoubleML call requires 0.")
  
  leaf_indices_train <- predict(model, dtrain, predleaf = TRUE)
  leaf_indices_test <- predict(model, dtest, predleaf = TRUE)
  
  smoothers <- create_S_from_gbtregressor(
    model, leaf_indices_train, leaf_indices_test, base_score, save_output, output_dir
  )
  
  return(list(S_train = smoothers$S_train, S_test = smoothers$S_test))
}
