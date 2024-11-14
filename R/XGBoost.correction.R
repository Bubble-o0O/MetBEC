#' @title QC-based BEC algorithm: eXtreme Gradient Boost (XGBoost) regression
#'
#' @description Briefly, the intra-BEC principle employs the regression model based on QC samples batchwise, namely:\cr
#' \deqn{y_{i}\sim f(x,\boldsymbol{y}_{(i)}\mid \boldsymbol{\theta}),i=1,2,…,p,}\cr
#' where \eqn{y_{i}} denotes the intensity of the \eqn{i}th metabolite;\cr
#' \eqn{x} denotes the injection order;\cr
#' \eqn{\boldsymbol{y}_{(i)}} denotes the several (i.e., \code{cor_variable_num}) variables with the highest correlations to the \eqn{i}th variable;\cr
#' \eqn{\boldsymbol{\theta}} denotes the optimal hyperparameters;\cr
#' \eqn{p} denotes the variable (metabolite) number.
#'
#' @param data A dataframe.\cr
#' Execute \code{data(Dataset_I)} and \code{View(Dataset_I)} for formats.
#' @param batch_ratio Default is \code{ratio-A} when the batch number \eqn{B>1}, or else is \code{NULL} when \eqn{B=1}.
#' @param cor_variable_num Default is \code{NULL}, which equals 10 when the variable number \eqn{p>10} or else equals \eqn{p-1}. Otherwise, it should be a numeric scalar. The hyperparameter usually has slight influence on correction.
#' @param nrounds A numeric scalar or vector. Default is 50.
#' @param eta A numeric scalar or vector. Default is \code{c(0.05, 0.1)}.
#' @param gamma A numeric scalar or vector. Default is 1. The hyperparameter usually has slight influence on correction.
#' @param min_child_weight A numeric scalar or vector. Default is 2.
#' @param max_depth A numeric scalar or vector. Default is 4. The hyperparameter is similar to \code{nodesize} or \code{maxnodes} in \code{randomForest::\link{randomForest}}.
#' @param subsample A numeric scalar or vector. Default is 0.7. The hyperparameter is applied without replacement, contrary to bagging (\code{replace = TRUE}) in \code{randomForest::\link{randomForest}}.
#' @param colsample_bytree A numeric scalar or vector. Default is 0.8. The hyperparameter is similar to \code{mtry} in \code{randomForest::\link{randomForest}}.
#' @param cl Default is \code{NULL}, which uses all the CPU cores for parallel computing. Otherwise, it should be a numeric scalar.
#'
#' @return A dataframe. The corrected data, whose format is the same as \code{data}.
#'
#' @export
#' @import stats parallel caret
#'
#' @note
#' \itemize{
#'  \item{XGBoost has more hyperparameters when establishing the regression model, such that it is more sophisticated and powerful to handle different scales of datasets.}
#'  \item{5-fold-cross-validation is used for hyperparameter optimization.}
#'  \item{If \code{batch_ratio = NULL}, only intra-BEC will be implemented. If \code{batch_ratio = ratio-A} or \code{mean}, it can ensure the consistency among QC samples' mean vectors across different batches after correction.}
#' }
#' @author Zhendong Guo (\email{guozhendong19@mails.ucas.ac.cn}).
#' @references
#' \itemize{
#'  \item{Chen, T.; Guestrin, C.; Machinery, A. C. XGBoost: A Scalable Tree Boosting System. In \emph{KDD'16: Proceedings of the 22Nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining}, 2016; pp 785-794. DOI: 10.1145/2939672.2939785.}
#'  \item{Bentéjac, C.; Csrg, A.; Martínez-Muoz, G. A Comparative Analysis of XGBoost. \emph{arXiv} \strong{2019}. DOI: arXiv:1911.01914.}
#' }
#' @seealso \code{\link{batch_ratio.correction}}, \code{caret::\link[caret]{train}}, \code{xgboost::xgboost}.
#' @examples
#' data(Dataset_I)
#' data <- Dataset_I
#' data.XGBoost <- XGBoost.correction(data)

XGBoost.correction <- function(data, batch_ratio = c(NULL, 'ratio-A', 'median', 'mean'),
                               cor_variable_num = NULL,
                               nrounds = 50, eta = c(0.05, 0.1),
                               gamma = 1, min_child_weight = 2, max_depth = 4,
                               subsample = 0.7, colsample_bytree = 0.8,
                               cl = NULL){
  info <- data[, 1:4]

  batch.level <- unique(data$batch)
  B <- length(batch.level) # 批次数

  QC <- subset(data, data$type == 'qc')

  p <- ncol(data[, -1:-4]) # 变量数

  if (is.null(cor_variable_num) == TRUE){
    if (p > 10){
      cor_variable_num <- 10
    }else {
      cor_variable_num <- p - 1
    }
  }else if (cor_variable_num < 0 || cor_variable_num > p - 1){
    stop("'cor_variable_num' must be in [0,p-1],
         where p denotes the variable number.")
  }

  if (any(nrounds < 0) == TRUE){ # 该超参数越大，拟合效果越好（但太大将导致过拟合）
    # 经测试，该超参数对校正效果的影响较大
    stop("'nrounds' must be non-negative.")
  }

  if (any(eta < 0) == TRUE){ # 该超参数越大，拟合效果越好（但太大将导致过拟合）
    # 经测试，该超参数对校正效果的影响较大
    stop("'eta' must be non-negative.")
  }

  if (any(gamma < 0) == TRUE){ # 该超参数越小，拟合效果越好
    # 经测试，该超参数对校正效果几乎无影响
    stop("'gamma' must be non-negative.")
  }

  if (any(min_child_weight < 0) == TRUE){ # 该超参数越小，拟合效果越好（但太小将导致过拟合）
    stop("'min_child_weight' must be non-negative.")
  }

  if (any(max_depth <= 0) == TRUE){ # 该超参数越大，拟合效果越好（但太大将导致过拟合）
    # 经测试，该超参数对校正效果的影响不大
    # 该超参数类似于随机森林的nodesize或maxnodes
    stop("'max_depth' must be positive.")
  }

  if (any(subsample <= 0) == TRUE || any(subsample > 1) == TRUE){
    # 该超参数执行不放回抽样，这不同于随机森林的bagging
    stop("'subsample' must be in (0,1].")
  }

  if (any(colsample_bytree <= 0) == TRUE || any(colsample_bytree > 1) == TRUE){
    # 经测试，该超参数对校正效果的影响不大
    # 该超参数类似于随机森林的mtry（每次分裂只纳入部分变量），每棵树只纳入部分变量
    stop("'colsample_bytree' must be in (0,1].")
  }

  grid_search <- expand.grid(nrounds = unique(nrounds), # 迭代次数
                             eta = unique(eta), # 默认：0.3
                             gamma = unique(gamma), # 默认：0
                             min_child_weight = unique(min_child_weight), # 默认：1
                             max_depth = unique(max_depth), # 默认：6
                             subsample = unique(subsample), # 默认：1
                             colsample_bytree = unique(colsample_bytree)) # 默认：1

  if (is.null(cl) == TRUE){
    cl <- parallel::makeCluster(parallel::detectCores())
  }else {
    cl <- parallel::makeCluster(cl)
  }

  all <- do.call(rbind, lapply(1:B, function(b){
    data_b <- subset(data, data$batch == batch.level[b])
    order_b <- data_b$order # 数值型变量（进样顺序）
    Data_b <- data_b[, -1:-4]

    QC_b.data <- subset(QC, QC$batch == batch.level[b])
    QC.order_b <- QC_b.data$order
    QC_b.Data <- QC_b.data[, -1:-4]

    QC_b.var <- apply(QC_b.Data, 2, var)
    if (any(QC_b.var < 1e-6) == TRUE){
      var_zero <- colnames(QC_b.Data)[which(QC_b.var < 1e-6)]
      stop(sprintf("For QC samples in Batch %s, %d variables' variances are 0 (or close to 0), including: %s.",
                   batch.level[b], length(var_zero), paste(var_zero, collapse = ", ")))
    }

    # 备注：XGBoost模型不需要对数据进行标准化
    all_b <- parallel::parSapply(cl, 1:p, function(j){
      set.seed(j)
      correlation <- abs(cor(QC_b.Data[, j], QC_b.Data)[1, ]) # 对象范围：每个批次的QC样本；参考Han (2022)
      cor_index <- match(names(sort(correlation, decreasing = TRUE)[1:(cor_variable_num + 1)][-1]),
                         names(correlation)) # 得到与第i个变量相关性（绝对值）最强的若干个变量的索引；[-1]表示除去第i个变量自身
      QC.Data <- data.frame(y = QC_b.Data[, j], QC_b.Data[, cor_index], order = QC.order_b)
      colnames(QC.Data) <- c('y', names(correlation)[cor_index], 'order')

      XGBoost_model <- caret::train(method = 'xgbTree', y ~ ., data = QC.Data, tuneGrid = grid_search,
                                    trControl = caret::trainControl(method = 'cv', number = 5)) # 交叉验证

      # best_para_b <- XGBoost_model$bestTune
      all.Data <- data.frame(Data_b[, cor_index], order = order_b) # 根据order排序
      colnames(all.Data) <- c(names(correlation)[cor_index], 'order')

      all.correction <- Data_b[, j] / predict(XGBoost_model, all.Data) * median(QC_b.Data[, j]) # 除法校正

      return(all.correction)
    })

    # 替换非正值
    if (any(all_b <= 0) == TRUE){
      # warning(sprintf("%.d samples' intensity in Batch %s is non-positive and has been replaced.",
      #                 sum(all_b <= 0), batch.level[b]))
      all_b[all_b <= 0] <- Data_b[all_b <= 0]
    }

    return(all_b)
  }))

  # 关闭集群
  parallel::stopCluster(cl)

  colnames(all) <- colnames(data[, -1:-4])
  all <- cbind(info, all)

  if (length(batch_ratio) > 1){
    batch_ratio <- 'ratio-A'
  }

  if (B == 1){
    batch_ratio <- NULL
  }

  if (is.null(batch_ratio) == FALSE){
    all <- batch_ratio.correction(all, method = batch_ratio)
  }

  return(all)
}
