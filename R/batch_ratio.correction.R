#' @title QC-based BEC algorithm: batch-ratio
#'
#' @description Batch-ratio is used for inter-BEC, which has several versions. Here, we provide three versions and recommend to use \code{ratio-A}.
#'
#' @param data A dataframe. \strong{Use \code{data(Dataset_I)} for formats.}
#' @param method Default is \code{ratio-A}, which combines median values and mean values, and is unnecessary to be tuned for most of the time.
#'
#' @return A dataframe. The corrected data, whose format is the same as \code{data}.
#'
#' @export
#' @import stats
#'
#' @note
#' \itemize{
#'  \item{When \code{method = ratio-A} or \code{mean}, batch-ratio can ensure the consistency among QC samples' mean vectors across different batches after correction.}
#' }
#' @author Zhendong Guo (\email{guozhendong19@mails.ucas.ac.cn}).
#' @references
#' \itemize{
#'  \item{van der Kloet, F.; Bobeldijk, I.; Verheij, E.; Jellema, R. Analytical Error Reduction Using Single Point Calibration for Accurate and Precise Metabolomic Phenotyping. \emph{Journal of Proteome Research} \strong{2009}, 8, 5132-5141. DOI: 10.1021/pr900499r.}
#'  \item{Kamleh, M. A.; Ebbels, T. M. D.; Spagou, K.; Masson, P.; Want, E. J. Optimizing the Use of Quality Control Samples for Signal Drift Correction in Large-Scale Urine Metabolic Profiling Studies. \emph{Analytical Chemistry} \strong{2012}, 84 (6), 2670-2677. DOI: 10.1021/ac202733q.}
#'  \item{Wang, S.-Y.; Kuo, C.-H.; Tseng, Y. J. Batch Normalizer: A Fast Total Abundance Regression Calibration Method to Simultaneously Adjust Batch and Injection Order Effects in Liquid Chromatography/Time-of-Flight Mass Spectrometry-Based Metabolomics Data and Comparison with Current Calibration Methods. \emph{Analytical Chemistry} \strong{2012}, 85 (2), 1037-1046. DOI: 10.1021/ac302877x.}
#' }
#' @examples
#' data(Dataset_I)
#' data <- Dataset_I
#' data.batch_ratio <- batch_ratio.correction(data)

batch_ratio.correction <- function(data,
                                   method = c('ratio-A', 'median', 'mean')){
  # 如果method = 'ratio-A'，则校正因子为median(QC.all)/mean(QC.b)，参考Wang (2013)；
  # 如果method = 'median'，则校正因子为median(QC.all)/median(QC.b)；
  # 如果method = 'mean'，则校正因子为mean(QC.all)/mean(QC.b)；
  # 校正因子的分母部分若为mean(QC.b)，可使每个批次的QC样本的校正后均值恰为校正因子的分子部分，
  # 因此method = 'median'不具备该性质；
  # 又因为中位数不受离群点的影响，而ratio-A综合了median与mean，所以推荐使用method = 'ratio-A'
  batch.level <- unique(data$batch)
  B <- length(batch.level) # 批次数
  if (B == 1){
    stop("There is only one batch.")
  }

  if (length(method) > 1){
    method <- 'ratio-A'
  }

  Data <- data[, -1:-4]
  p <- ncol(Data) # 变量数

  QC <- subset(data, data$type == 'qc')
  QC.Data <- subset(data, data$type == 'qc')[, -1:-4] # 所有的QC样本



  for (b in 1:B){
    Data_b <- subset(data, data$batch == batch.level[b])[, -1:-4]

    QC_b <- subset(QC, QC$batch == batch.level[b])
    QC_b.Data <- QC_b[, -1:-4]

    Data_b.scale <- sapply(1:p, function(j) { # 获取校正因子（向量形式）
      if (method == 'ratio-A'){
        median(QC.Data[, j]) / mean(QC_b.Data[, j])
      }else if (method == 'median'){
        median(QC.Data[, j]) / median(QC_b.Data[, j])
      }else if (method == 'mean'){
        mean(QC.Data[, j]) / mean(QC_b.Data[, j])
      }
    })

    Data_b.correction <- t(Data_b.scale * t(Data_b))
    data[data$batch == batch.level[b], -1:-4] <- Data_b.correction
  }

  return(data)
}
