#' @title QC-based BEC algorithm: Random Forest (RF) regression
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
#' @param cor_variable_num Default is \code{NULL}, which equals to 10 when the variable number \eqn{p>10} or else equals to \eqn{p-1}. Otherwise, it should be a numeric scalar. The hyperparameter usually has slight influence on correction.
#' @param mtry Default is \code{NULL}, which calls \code{randomForest::\link{tuneRF}} to automatically obtain a numeric vector for hyperparameter optimization. Otherwise, it should be a numeric scalar.
#' @param nodesize_percent A numeric scalar or vector. Default is 0.2.
#' @param ntree A numeric scalar. Default is 500. The hyperparameter usually has slight influence on correction.
#' @param cl Default is \code{NULL}, which uses all the CPU cores for parallel computing. Otherwise, it should be a numeric scalar.
#'
#' @return A dataframe. The corrected data, whose format is the same as \code{data}.
#'
#' @export
#' @import stats parallel
#' @importFrom randomForest tuneRF randomForest
#'
#' @note
#' \itemize{
#'  \item{Bootstrap is used for hyperparameter optimization instead of cross-validation (CV).}
#'  \item{When \code{batch_ratio = ratio-A} or \code{mean}, batch-ratio can ensure the consistency among QC samples' mean vectors across different batches after correction.}
#' }
#' @author Zhendong Guo (\email{guozhendong19@mails.ucas.ac.cn}).
#' @references
#' \itemize{
#'  \item{Luan, H.; Ji, F.; Chen, Y.; Cai, Z. statTarget: A streamlined tool for signal drift correction and interpretations of quantitative mass spectrometry-based omics data. \emph{Analytica Chimica Acta} \strong{2018}, 1036, 66-72. DOI: 10.1016/j.aca.2018.08.002.}
#'  \item{Fan, S.; Kind, T.; Cajka, T.; Hazen, S. L.; Tang, W. H. W.; Kaddurah-Daouk, R.; Irvin, M. R.; Arnett, D. K.; Barupal, D. K.; Fiehn, O. Systematic Error Removal Using Random Forest for Normalizing Large-Scale Untargeted Lipidomics Data. \emph{Analytical Chemistry} \strong{2019}, 91 (5), 3590-3596. DOI: 10.1021/acs.analchem.8b05592.}
#' }
#' @seealso \code{\link{batch_ratio.correction}}, \code{randomForest::\link[randomForest]{tuneRF}}, \code{randomForest::\link[randomForest]{randomForest}}.
#' @examples
#' data(Dataset_I)
#' data <- Dataset_I
#' data.RF <- RF.correction(data)

RF.correction <- function(data, batch_ratio = c(NULL, 'ratio-A', 'median', 'mean'),
                         cor_variable_num = NULL,
                         mtry = NULL, nodesize_percent = 0.2, ntree = 500,
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

  if (is.null(mtry) == TRUE){
    mtryStart <- max(1, floor((cor_variable_num + 1)/3)) # mtry的默认值为自变量个数的1/3
  }else if (mtry < 1 || mtry > (cor_variable_num + 1)){
    stop("'mtry' must be in [1,cor_variable_num + 1].")
  }else {
    mtryStart <- mtry
  }

  if (any(nodesize_percent < 0) == TRUE || any(nodesize_percent > 1) == TRUE){  # 该超参数越小，拟合效果越好（但太小将导致过拟合）
    stop("'nodesize_percent' must be in [0,1].")
  }
  nodesize_percent <- unique(nodesize_percent)

  if (ntree <= 0){
    # 经测试，该超参数对校正效果的影响不大
    stop("'ntree' must be positive.")
  }

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

    # 备注：随机森林模型不需要对数据进行标准化
    all_b <- parallel::parSapply(cl, 1:p, function(j){
      set.seed(j)
      correlation <- abs(cor(QC_b.Data[, j], QC_b.Data)[1, ]) # 对象范围：每个批次的QC样本；参考Han (2022)
      cor_index <- match(names(sort(correlation, decreasing = TRUE)[1:(cor_variable_num + 1)][-1]),
                         names(correlation)) # 得到与第i个变量相关性（绝对值）最强的若干个变量的索引；[-1]表示除去第i个变量自身
      QC.Data <- data.frame(y = QC_b.Data[, j], QC_b.Data[, cor_index], order = QC.order_b)
      colnames(QC.Data) <- c('y', names(correlation)[cor_index], 'order')

      nodesize <- sapply(1:length(nodesize_percent), function(k)
                         max(1, floor(nodesize_percent[k] * nrow(QC.Data))))

      if (is.null(mtry) == TRUE && cor_variable_num > 0){ # 如果cor_variable_num = 0，则mtry只能取1，因此tuneRF函数将报错
        OOB_error <- do.call(rbind, lapply(1:length(nodesize), function(k){
          tune_RF <- data.frame(randomForest::tuneRF(x = QC.Data[, -1], y = QC.Data[, 1], nodesize = nodesize[k],
                                                     mtryStart = mtryStart, ntreeTry = ntree,
                                                     trace = FALSE, plot = FALSE)) # 该函数没有subset的功能
          # 使用默认参数stepFactor = 2；mtry的优化范围: 将mtry分别不断乘以和除以stepFactor
          best_k <- which.min(tune_RF[, 2])
          return(tune_RF[best_k, ]) # 选择OOB error最小的行
        }))

        best <- which.min(OOB_error[, 2])
        mtry <- OOB_error[best, 1] # 选择OOB error最小的mtry及其对应的nodesize
        nodesize <- nodesize[best]
      }else {
        mtry <- mtryStart
        if (length(nodesize) > 1){
          OOB_error <- sapply(1:length(nodesize), function(k){
            return(randomForest::randomForest(y ~ ., data = QC.Data, ntree = ntree,
                                              mtry = mtry, nodesize = nodesize[k])$mse[ntree])
            # mse为ntree维向量，且呈递减趋势；最后一个数值（mse[ntree]）是该模型最终的MSE（参考官方tuneRF函数的代码）
          })

          best <- which.min(OOB_error)
          nodesize <- nodesize[best] # 选择OOB error最小的nodesize
        }
      }

      RF_model <- randomForest::randomForest(y ~ ., data = QC.Data, ntree = ntree,
                                             mtry = mtry, nodesize = nodesize)
      # 因变量是y；自变量是data中去除y的所有列；subset可以选择data的子集作为训练集

      all.Data <- data.frame(Data_b[, cor_index], order = order_b)
      colnames(all.Data) <- c(names(correlation)[cor_index], 'order')

      all.correction <- Data_b[, j] / predict(RF_model, all.Data) * median(QC_b.Data[, j]) # 除法校正

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
