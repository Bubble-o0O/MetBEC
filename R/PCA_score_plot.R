#' @title Principal Component Analysis (PCA) score plot
#'
#' @description A two-dimensional PCA score plot for multivariate visualization.
#'
#' @param data A dataframe.\cr
#' Execute \code{data(Dataset_I)} and \code{View(Dataset_I)} for formats.
#' @param sample_type \code{all} denotes all the samples; \code{qc} denotes QC samples; \code{sample} denotes subject samples.
#' @param PCi A numeric scalar. Default is 1. The \eqn{i}th principal component of PCA.
#' @param PCj A numeric scalar. Default is 2. The \eqn{j}th principal component of PCA.
#' @param PCA_scale Logical. Default is \code{TRUE}, which determines whether to autoscale (standardize) \code{data} and is unnecessary to be tuned for most of the time.
#' @param point_size A numeric scalar. Default is 1. Determines the point size in the scatter plot.
#' @param point_alpha A numeric scalar. Default is 0.6. Determines the point opacity in the scatter plot.
#'
#' @return A two-dimensional PCA score plot.
#'
#' @export
#' @import stats ggplot2
#'
#' @note
#' \itemize{
#'  \item{PCA is a qualitative method for correction assessment, whose two-dimensional score plot might be limited and even misleading. Hence, pay attention to \eqn{R^2}.}
#' }
#' @author Zhendong Guo (\email{guozhendong19@mails.ucas.ac.cn}).
#' @seealso \code{\link{compare_PCA_plot}}, \code{stats::\link[stats]{prcomp}}.
#' @examples
#' data(Dataset_I)
#' data <- Dataset_I
#' PCA_plot(data, sample_type = "qc")

PCA_plot <- function(data, sample_type = c("all", "qc", "sample"),
                     PCi = 1, PCj = 2, PCA_scale = TRUE,
                     point_size = 1, point_alpha = 0.6){
  batch.level <- unique(data$batch)
  B <- length(batch.level) # 批次数

  QC <- subset(data, data$type == 'qc')
  sample <- subset(data, data$type == 'sample')

  data$type[data$type == 'qc'] <- 'QC'
  data$type[data$type == 'sample'] <- 'subject'

  if (sample_type == 'all'){
    batch <- data$batch
    X <- data[, -1:-4]
    title <- 'QC samples + Subject samples'
  }else if (sample_type == 'qc'){
    batch <- QC$batch
    X <- QC[, -1:-4]
    title <- 'QC samples'
  }else if (sample_type == 'sample'){
    batch <- sample$batch
    X <- sample[, -1:-4]
    title <- 'Subject samples'
  }

  if (PCA_scale == TRUE){
    pca <- prcomp(X, scale. = TRUE)
  }else {
    pca <- prcomp(X, scale. = FALSE) # 该函数默认为中心化
  }

  scores <- data.frame(pca$x[, c(PCi, PCj)])
  scores$Batch <- factor(batch, levels = batch.level)

  R2 <- pca$sdev^2/sum(pca$sdev^2) # 每个主成分的方差贡献率

  PCi_text <- sprintf("PC%.d (%.1f%%)", PCi, 100 * R2[PCi]) # 数值型变量的格式化输出
  PCj_text <- sprintf("PC%.d (%.1f%%)", PCj, 100 * R2[PCj])

  if (sample_type == 'all'){
    scores$Type <- data$type

    QC.scores <- subset(scores, scores$Type == 'QC')

    ggplot2::ggplot(scores, aes(x = scores[, 1], y = scores[, 2],
                       color = Batch, shape = Type)) +
      geom_point(size = point_size, alpha = point_alpha) +
      scale_shape_manual(values = c(19, 1)) +
      stat_ellipse(data = QC.scores,
                   mapping = aes(x = QC.scores[, 1], y = QC.scores[, 2]),
                   level = 0.95, type = 'norm', alpha = 0.8) +
      labs(x = PCi_text, y = PCj_text, title = title) +
      theme_minimal()
  }else {
    ggplot2::ggplot(scores, aes(x = scores[, 1], y = scores[, 2],
                       color = Batch)) +
      geom_point(size = point_size, alpha = point_alpha) +
      stat_ellipse(level = 0.95, type = 'norm', alpha = 0.8) +
      labs(x = PCi_text, y = PCj_text, title = title) +
      theme_minimal()
  }
}


#' @title Principal Component Analysis (PCA) score plot of different data for comparison
#'
#' @description A two-dimensional PCA score plot for multivariate visualization.
#'
#' @param data1 A dataframe.\cr
#' Execute \code{data(Dataset_I)} and \code{View(Dataset_I)} for formats.
#' @param ... Alternatives. More data (e.g., \code{data2}, \code{data3}) can be input.
#' @param BEC_method A character vector. The length must be the same as \code{c(data1, ...)}. For example, \code{c("raw", "RF", "XGBoost")}.
#' @param batch_index A numeric scalar.
#' @param sample_type \code{all} denotes all the samples; \code{qc} denotes QC samples; \code{sample} denotes subject samples.
#' @param PCi A numeric scalar. Default is 1. The \eqn{i}th principal component of PCA.
#' @param PCj A numeric scalar. Default is 2. The \eqn{j}th principal component of PCA.
#' @param PCA_scale Logical. Default is \code{TRUE}, which determines whether to autoscale (standardize) \code{c(data, ...)} and is unnecessary to be tuned for most of the time.
#' @param point_size A numeric scalar. Default is 1. Determines the point size in the scatter plot.
#' @param point_alpha A numeric scalar. Default is 0.6. Determines the point opacity in the scatter plot.
#'
#' @return A two-dimensional PCA score plot.
#'
#' @export
#' @import stats ggplot2
#'
#' @author Zhendong Guo (\email{guozhendong19@mails.ucas.ac.cn}).
#' @note
#' \itemize{
#'  \item{PCA is a qualitative method for correction assessment, whose two-dimensional score plot might be limited and even misleading. Hence, pay attention to \eqn{R^2}.}
#' }
#' @seealso \code{\link{PCA_plot}}, \code{stats::\link[stats]{prcomp}}.
#' @examples
#' data(Dataset_I)
#' data <- Dataset_I
#' data.RF <- RF.correction(data)
#' compare_PCA_plot(data, data.RF, BEC_method = c("raw", "RF"),
#'                  batch_index = 2, sample_type = "qc")

compare_PCA_plot <- function(data1, ...,
                             BEC_method, batch_index,
                             sample_type = c('all', 'qc', 'sample'),
                             PCi = 1, PCj = 2, PCA_scale = TRUE,
                             point_size = 1, point_alpha = 0.6){
  data <- list(data1, ...)

  if (length(BEC_method) != length(data)){
    stop("'BEC_method' should be correted.")
  }

  batch.level <- unique(data1$batch)
  b <- batch.level[batch_index] # 选择作图的批次

  data_b <- do.call(rbind, lapply(1:length(data), function(x){
    Data_b <- subset(data[[x]], data[[x]]$batch == b)

    Data_b$type[Data_b$type == 'qc'] <- 'QC'
    Data_b$type[Data_b$type == 'sample'] <- 'subject'

    colnames(Data_b)[colnames(Data_b) == 'batch'] <- 'Method'
    Data_b$Method <- rep(BEC_method[x], nrow(Data_b))
    return(Data_b)
  }))

  QC_b <- subset(data_b, data_b$type == 'QC')
  sample_b <- subset(data_b, data_b$type == 'subject')

  if (sample_type == 'all'){
    Method <- data_b$Method
    X <- data_b[, -1:-4]
    title <- sprintf('QC samples + Subject samples in Batch %s', b) # 字符串型变量的格式化输出
  }else if (sample_type == 'qc'){
    Method <- QC_b$Method
    X <- QC_b[, -1:-4]
    title <- sprintf('QC samples in Batch %s', b)
  }else if (sample_type == 'sample'){
    Method <- sample_b$Method
    X <- sample_b[, -1:-4]
    title <- sprintf('Subject samples in Batch %s', b)
  }

  if (PCA_scale == TRUE){
    pca <- prcomp(X, scale. = TRUE)
  }else {
    pca <- prcomp(X, scale. = FALSE) # 该函数默认为中心化
  }

  scores <- data.frame(pca$x[, c(PCi, PCj)])
  scores$Method <- factor(Method, levels = BEC_method)

  R2 <- pca$sdev^2/sum(pca$sdev^2) # 每个主成分的方差贡献率

  PCi_text <- sprintf("PC%.d (%.1f%%)", PCi, 100 * R2[PCi]) # 数值型变量的格式化输出
  PCj_text <- sprintf("PC%.d (%.1f%%)", PCj, 100 * R2[PCj])

  if (sample_type == 'all'){
    scores$Type <- data_b$type

    QC.scores <- subset(scores, scores$Type == 'QC')

    ggplot2::ggplot(scores, aes(x = scores[, 1], y = scores[, 2],
                                color = Method, shape = Type)) +
      geom_point(size = point_size, alpha = point_alpha) +
      scale_shape_manual(values = c(19, 1)) +
      stat_ellipse(data = QC.scores,
                   mapping = aes(x = QC.scores[, 1], y = QC.scores[, 2]),
                   level = 0.95, type = 'norm', alpha = 0.8) +
      labs(x = PCi_text, y = PCj_text, title = title) +
      theme_minimal()
  }else {
    ggplot2::ggplot(scores, aes(x = scores[, 1], y = scores[, 2],
                                color = Method)) +
      geom_point(size = point_size, alpha = point_alpha) +
      stat_ellipse(level = 0.95, type = 'norm', alpha = 0.8) +
      labs(x = PCi_text, y = PCj_text, title = title) +
      theme_minimal()
  }
}
