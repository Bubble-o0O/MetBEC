#' @title QC-based BEC algorithm: Technical variation elImination with ensemble learninG architEctuRe (TIGER)
#'
#' @description Briefly, the intra-BEC principle is divided into two parts based on QC samples batchwise: the base model and the meta model. The former employs Random Forest (RF) to establish the regression model similar to \code{\link{SVR.correction}}, \code{\link{RF.correction}}, and \code{\link{XGBoost.correction}}, and the latter integrates several base models with different hyperparameters weighted by the loss function.
#'
#' @param data A dataframe. \strong{Use \code{data(Dataset_I)} for formats.}
#' @param batch_ratio Default is \code{ratio-A} when the batch number \eqn{B>1}, otherwise default is \code{NULL} when \eqn{B=1}. Noted that if \code{batch_ratio = NULL}, only intra-BEC will be implemented.
#' @param cor_variable_num A numeric scalar. Default is 10. The hyperparameter usually has slight influence on correction.
#' @param mtry_percent A numeric scalar or vector. Default is \code{seq(0.2, 0.8, 0.2)}.
#' @param nodesize_percent A numeric scalar or vector. Default is 0.2.
#' @param ntree A numeric scalar or vector. Default is 500. The hyperparameter usually has slight influence on correction.
#' @param cl Default is \code{NULL}, which uses all the CPU cores for parallel computing. Otherwise, it should be a numeric scalar.
#'
#' @return A dataframe. The corrected data, whose format is the same as \code{data}.
#'
#' @export
#' @import stats parallel TIGERr
#'
#' @note
#' \itemize{
#'  \item{TIGER has a more complex correction principle than \code{\link{SVR.correction}}, \code{\link{RF.correction}}, and \code{\link{XGBoost.correction}} (see the reference below for details).}
#'  \item{5-fold-cross-validation (CV) is used for hyperparameter optimization.}
#'  \item{When \code{batch_ratio = ratio-A} or \code{mean}, batch-ratio can ensure the consistency among QC samples' mean vectors across different batches after correction.}
#' }
#' @author Zhendong Guo (\email{guozhendong19@mails.ucas.ac.cn}).
#' @references
#' \itemize{
#'  \item{Han, S.; Huang, J.; Foppiano, F.; Prehn, C.; Adamski, J.; Suhre, K.; Li, Y.; Matullo, G.; Schliess, F.; Gieger, C.; et al. TIGER: technical variation elimination for metabolomics data using ensemble learning architecture. \emph{Brief Bioinform} \strong{2022}, 23 (2). DOI: 10.1093/bib/bbab535.}
#' }
#' @seealso \code{\link{batch_ratio.correction}}, \code{TIGERr::\link[TIGERr]{run_TIGER}}, \code{randomForest::\link[randomForest]{randomForest}}.
#' @examples
#' data(Dataset_I)
#' data <- Dataset_I
#' data.TIGER <- TIGER.correction(data)

TIGER.correction <- function(data, batch_ratio = c(NULL, 'ratio-A', 'median', 'mean'),
                             cor_variable_num = 10,
                             mtry_percent = seq(0.2, 0.8, 0.2),
                             nodesize_percent = 0.2, ntree = 500,
                             cl = NULL){
  batch.level <- unique(data$batch)
  B <- length(batch.level) # 批次数

  p <- ncol(data[, -1:-4]) # 变量数

  if (cor_variable_num < 1 || cor_variable_num > p - 1){
    stop("'cor_variable_num' must be in [1,p-1],
         where p denotes the variable number.")
  }

  if (any(mtry_percent < 0) == TRUE || any(mtry_percent > 1) == TRUE){
    stop("'mtry_percent' must be in [0,1].")
  }

  if (any(nodesize_percent < 0) == TRUE || any(nodesize_percent > 1) == TRUE){  # 该超参数越小，拟合效果越好（但太小将导致过拟合）
    stop("'nodesize_percent' must be in [0,1].")
  }

  if (any(ntree <= 0) == TRUE){
    stop("'ntree' must be positive.")
  }

  if (is.null(cl) == TRUE){
    cl <- parallel::detectCores()
  }

  # QC samples as training samples; subject samples as test samples:
  QC <- subset(data, data$type == "qc")
  samples  <- subset(data, data$type == "sample")

  # Only numeric data of metabolite variables are allowed:
  QC.Data <- QC[, -1:-4]
  sample.Data <- samples[, -1:-4]

  for (b in 1:B){
    QC_b.data <- subset(QC, QC$batch == batch.level[b])
    QC_b.Data <- QC_b.data[, -1:-4]

    if (any(apply(QC_b.Data, 2, var) == 0) == TRUE){
      stop(sprintf("QC samples' some variances in Batch %s are 0.",
                   batch.level[b]))
    }
  }

  #### compute_targetVal ####
  tarVal <- TIGERr::compute_targetVal(QC_num = QC.Data,
                                      sampleType = QC$type,
                                      batchID = QC$batch,
                                      targetVal_method = "median",
                                      targetVal_batchWise = TRUE,
                                      # 默认是FALSE；经测试，FALSE的结果与TRUE + batch ratio的结果无明显差异；但参考该文献，还是建议选择TRUE
                                      targetVal_removeOutlier = FALSE)

  #### select_variable ####
  selected_var <- TIGERr::select_variable(train_num = QC.Data,
                                          test_num  = NULL,
                                          train_batchID = QC$batch,
                                          test_batchID  = NULL,
                                          selectVar_batchWise = TRUE,
                                          # 如果选择TRUE，就对训练集的每个批次进行操作；如果选择FALSE，就对全部训练集进行操作
                                          selectVar_corMethod = "pearson",
                                          selectVar_minNum = cor_variable_num,
                                          selectVar_maxNum = cor_variable_num)
  # Shen (2016)选择了QC样本中相关性最强的5个变量作为自变量，但自变量未纳入进样顺序

  #### run_TIGER ####
  all <- TIGERr::run_TIGER(test_samples = data, # 应为全部的原始数据，否则输出结果将没有训练集的校正结果
                           train_samples = QC,
                           col_sampleID  = "name",    # input column name
                           col_sampleType = "type",   # input column name
                           col_batchID = "batch",     # input column name
                           col_order = "order",       # input column name
                           targetVal_external = tarVal,
                           selectVar_external = selected_var,
                           mtry_percent = mtry_percent, # 备注：randomForest函数的mtry_percent默认值是1/3
                           nodesize_percent = nodesize_percent, # 备注：randomForest函数的nodesize默认值是5
                           ntree = ntree, # 备注：randomForest函数的ntree默认值是500
                           parallel.cores = cl)

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
