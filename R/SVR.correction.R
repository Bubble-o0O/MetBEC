#' @title QC-based BEC algorithm: Support Vector Regression (SVR)
#'
#' @description Briefly, the intra-BEC principle employs the regression model based on QC samples batchwise, namely:\cr
#' \deqn{y_{i}\sim f(x,\boldsymbol{y}_{(i)}\mid \boldsymbol{\theta}),i=1,2,…,p,}\cr
#' where \eqn{y_{i}} denotes the intensity of the \eqn{i}th metabolite;\cr
#' \eqn{x} denotes the injection order of the corresponding batch;\cr
#' \eqn{\boldsymbol{y}_{(i)}} denotes the several (i.e., \code{cor_variable_num}) variables with the highest correlations to \eqn{y_{i}};\cr
#' \eqn{\boldsymbol{\theta}} denotes the optimal hyperparameters;\cr
#' \eqn{p} denotes the variable (metabolite) number.
#'
#' @param data A dataframe. \strong{Use \code{data(Dataset_I)} for formats.}
#' @param batch_ratio Default is \code{ratio-A} when the batch number \eqn{B>1}, otherwise default is \code{NULL} when \eqn{B=1}. Noted that if \code{batch_ratio = NULL}, only intra-BEC will be implemented.
#' @param cor_variable_num Default is \code{NULL}, which equals to 10 when the variable number \eqn{p>10} or else equals to \eqn{p-1}. Otherwise, it should be a numeric scalar. The hyperparameter usually has slight influence on correction. Denoted as \eqn{p'}.
#' @param gamma Default is \code{NULL}, which is \eqn{\frac{1}{p'+1}·} \code{c(1,2,4)} for hyperparameter optimization. Otherwise, it should be a numeric scalar or vector.
#' @param cost A numeric scalar or vector. Default is 1.
#' @param epsilon A numeric scalar or vector. Default is 0.1. The hyperparameter usually has slight influence on correction.
#' @param cl Default is \code{NULL}, which uses all the CPU cores for parallel computing. Otherwise, it should be a numeric scalar.
#'
#' @return A dataframe. The corrected data, whose format is the same as \code{data}.
#'
#' @export
#' @import stats parallel e1071
#'
#' @note
#' \itemize{
#'  \item{5-fold-cross-validation (CV) is used for hyperparameter optimization.}
#'  \item{When \code{batch_ratio = ratio-A} or \code{mean}, batch-ratio can ensure the consistency among QC samples' mean vectors across different batches after correction.}
#' }
#' @author Zhendong Guo (\email{guozhendong19@mails.ucas.ac.cn}).
#' @references
#' \itemize{
#'  \item{Kuligowski, J.; Sanchez-Illana, A.; Sanjuan-Herraez, D.; Vento, M.; Quintas, G. Intra-batch effect correction in liquid chromatography-mass spectrometry using quality control samples and support vector regression (QC-SVRC). \emph{Analyst} \strong{2015}, 140 (22), 7810-7817. DOI: 10.1039/c5an01638j.}
#'  \item{Shen, X.; Gong, X.; Cai, Y.; Guo, Y.; Tu, J.; Li, H.; Zhang, T.; Wang, J.; Xue, F.; Zhu, Z.-J. Normalization and integration of large-scale metabolomics data using support vector regression. \emph{Metabolomics} \strong{2016}, 12 (5). DOI: 10.1007/s11306-016-1026-5.}
#'  \item{Ding, X.; Yang, F.; Chen, Y.; Xu, J.; He, J.; Zhang, R.; Abliz, Z. Norm ISWSVR: A Data Integration and Normalization Approach for Large-Scale Metabolomics. \emph{Analytical Chemistry} \strong{2022}, 94 (21), 7500-7509. DOI: 10.1021/acs.analchem.1c05502.}
#' }
#' @seealso \code{\link{batch_ratio.correction}}, \code{e1071::\link[e1071]{tune}}, \code{e1071::\link[e1071]{svm}}.
#' @examples
#' data(Dataset_I)
#' data <- Dataset_I
#' data.SVR <- SVR.correction(data)

SVR.correction <- function(data, batch_ratio = c(NULL, 'ratio-A', 'median', 'mean'),
                           cor_variable_num = NULL,
                           gamma = NULL, cost = 1, epsilon = 0.1, # cost与epsilon均为svm函数的默认值
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

  if (is.null(gamma) == TRUE){
    gamma_ref <- 1 / (cor_variable_num + 1) # 该超参数的默认值为自变量个数的倒数
    gamma <- gamma_ref * 2^(0:2)
  }else if (any(gamma <= 0) == TRUE){ # 该超参数越大，拟合效果越好（但太大将导致过拟合）
    stop("'gamma' must be positive.")
  }

  if (any(cost <= 0) == TRUE){ # 该超参数越大，拟合效果越好（但太大将导致过拟合）
    # 经测试，该超参数对校正效果的影响较大
    stop("'cost' must be positive.")
  }

  if (any(epsilon <= 0) == TRUE){ # 该超参数越小，拟合效果越好（但太小将导致过拟合）
    # 经测试，该超参数对校正效果的影响很小
    stop("'epsilon' must be positive.")
  }

  grid_search <- expand.grid(gamma = unique(gamma),
                             cost = unique(cost),
                             epsilon = unique(epsilon))

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

    if (any(apply(QC_b.Data, 2, var) == 0) == TRUE){
      stop(sprintf("QC samples' some variances in Batch %s are 0.",
                   batch.level[b]))
    }

    all_b <- parallel::parSapply(cl, 1:p, function(j){
      set.seed(j)
      correlation <- abs(cor(QC_b.Data[, j], QC_b.Data)[1, ]) # 对象范围：每个批次的QC样本；参考Han (2022)
      cor_index <- match(names(sort(correlation, decreasing = TRUE)[1:(cor_variable_num + 1)][-1]),
                         names(correlation)) # 得到与第i个变量相关性（绝对值）最强的若干个变量的索引；[-1]表示除去第i个变量自身
      QC.Data <- data.frame(y = QC_b.Data[, j], QC_b.Data[, cor_index], order = QC.order_b)
      colnames(QC.Data) <- c('y', names(correlation)[cor_index], 'order')

      SVR_model <- e1071::tune(e1071::svm, y ~ ., data = QC.Data, ranges = grid_search,
                               tunecontrol = e1071::tune.control(cross = 5)) # 交叉验证
      # tune与svm函数中默认包含scale = TRUE；
      # 如果不进行归一化，对x进行缩放后，模型将发生变化，但这其实是不合理的

      # best_para_b <- SVR_model$best.parameters
      all.Data <- data.frame(Data_b[, cor_index], order = order_b) # 根据order排序
      colnames(all.Data) <- c(names(correlation)[cor_index], 'order')

      all.correction <- Data_b[, j] / predict(SVR_model$best.model, all.Data) * median(QC_b.Data[, j]) # 除法校正

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
