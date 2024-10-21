#' @title D-ratio distribution of data
#'
#' @description Dispersion-ratio (D-ratio). The metric examines both QC samples and subject samples, whose formula is:\cr
#' \deqn{\mathrm{D\text{-}ratio}_{i}=\sqrt{\frac{\mathrm{Var}(y^{\mathrm{(QC)}}_{i})}{\mathrm{Var}(y^{\mathrm{(QC)}}_{i}) + \mathrm{Var}(y^{\mathrm{(ss)}}_{i})}},i=1,2,...,p,}\cr
#' where \eqn{y^{\mathrm{(QC)}}_{i}} denotes the intensity of QC samples' \eqn{i}th metabolite;\cr
#' \eqn{y^{\mathrm{(ss)}}_{i}} denotes the intensity of subject samples' \eqn{i}th metabolite;\cr
#' \eqn{p} denotes the variable (metabolite) number.
#'
#' @param data A dataframe. \strong{Use \code{data(Dataset_I)} for formats.}
#' @param print_plot Logical. Default is \code{TRUE}. Determines whether to plot the histogram and the cumulative distribution plot.
#' @param D_ratio_value Logical. Default is \code{FALSE}. Determines whether to output the D-ratio values.
#'
#' @return
#' \item{`D-ratio: 0 ~ 50\%`}{A numeric scalar. The cumulative frequency of D-ratio within 50\%.}
#' \item{D_ratio}{A numeric vector. Output only when \code{D_ratio_value = TRUE}.}
#' \item{histogram}{Plotted only when \code{print_plot = TRUE}.}
#' \item{cumulative distribution plot}{Plotted only when \code{print_plot = TRUE}.}
#'
#' @export
#' @import stats ggplot2
#' @importFrom scales label_number
#'
#' @note
#' \itemize{
#'  \item{The acceptance criterion is 50\%.}
#' }
#' @author Zhendong Guo (\email{guozhendong19@mails.ucas.ac.cn}).
#' @references
#' \itemize{
#'  \item{Broadhurst, D.; Goodacre, R.; Reinke, S. N.; Kuligowski, J.; Wilson, I. D.; Lewis, M. R.; Dunn, W. B. Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies. \emph{Metabolomics} \strong{2018}, 14 (6), 72. DOI: 10.1007/s11306-018-1367-3.}
#' }
#' @seealso \code{\link{compare_D_ratio}}.
#' @examples
#' data(Dataset_I)
#' data <- Dataset_I
#' D_ratio(data)

D_ratio <- function(data, print_plot = TRUE, D_ratio_value = FALSE){
  QC <- subset(data, data$type == 'qc')
  sample <- subset(data, data$type == 'sample')

  QC.Data <- QC[ ,-1:-4]
  sample.Data <- sample[ ,-1:-4]
  p <- ncol(QC.Data) # 变量数

  QC.var <- apply(QC.Data, 2, var)
  sample.var <- apply(sample.Data, 2, var)
  D_ratio <- sqrt(QC.var) / sqrt(QC.var + sample.var)
  names(D_ratio) <- colnames(QC.Data)

  D_ratio.0_50 <- length(which(D_ratio < 0.5))/p

  if (print_plot == TRUE){
    # 绘制直方图
    interval <- cut(D_ratio, breaks = c(0, 0.5, Inf), labels = c("< 50%", "\u2265 50%"),
                    right = FALSE) # 左闭右开

    histogram <- ggplot2::ggplot(data.frame(D_ratio), aes(x = D_ratio, fill = interval)) +
      geom_histogram(binwidth = 0.1, boundary = 0, closed = 'left', # 左闭右开
                     color = "black", show.legend = TRUE) +
      labs(title = "Histogram of D-ratio", x = "D-ratio (%)", y = "Count", fill = "D-ratio") +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1),
                         labels = scales::label_number(scale = 100)) +
      scale_fill_manual(values = c("< 50%" = "lightblue2", "\u2265 50%" = "lightpink2")) +
      theme_minimal()

    # 绘制累积分布图
    cdf <- ggplot2::ggplot(data.frame(D_ratio), aes(x = D_ratio)) +
      stat_ecdf(geom = "step", color = "black") +
      labs(title = "D-ratio distribution", x = "D-ratio (%)", y = "Cumulative distribution (%)") +
      geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.4) +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1),
                         labels = scales::label_number(scale = 100)) +
      scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1),
                         labels = scales::label_number(scale = 100)) +
      theme_minimal()

    print(histogram)
    print(cdf)
  }

  if (D_ratio_value == FALSE){
    return(list("D-ratio: 0 ~ 50%" = D_ratio.0_50))
  }else if (D_ratio_value == TRUE){
    return(list("D-ratio: 0 ~ 50%" = D_ratio.0_50,
                D_ratio = D_ratio))
  }
}

#' @title D-ratio distribution of different corrected data for comparison
#'
#' @description Dispersion-ratio (D-ratio). The metric examines both QC samples and subject samples, whose formula is:\cr
#' \deqn{\mathrm{D\text{-}ratio}_{i}=\sqrt{\frac{\mathrm{Var}(y^{\mathrm{(QC)}}_{i})}{\mathrm{Var}(y^{\mathrm{(QC)}}_{i}) + \mathrm{Var}(y^{\mathrm{(ss)}}_{i})}},i=1,2,...,p,}\cr
#' where \eqn{y^{\mathrm{(QC)}}_{i}} denotes the intensity of QC samples' \eqn{i}th metabolite;\cr
#' \eqn{y^{\mathrm{(ss)}}_{i}} denotes the intensity of subject samples' \eqn{i}th metabolite;\cr
#' \eqn{p} denotes the variable (metabolite) number.
#'
#' @param data1 A dataframe. \strong{Use \code{data(Dataset_I)} for formats.}
#' @param ... Alternatives. For example, \code{data2, data3}.
#' @param BEC_method A character vector. The length must be the same as \code{c(data1, ...)}. For example, \code{c('raw', 'RF', 'XGBoost')}.
#' @param line_alpha A numeric scalar. Default is 0.8. Determines the curve line opacity in the cumulative distribution plot.
#'
#' @return A cumulative distribution plot.
#'
#' @export
#' @import stats ggplot2
#' @importFrom scales label_number
#'
#' @note
#' \itemize{
#'  \item{The acceptance criterion is 50\%.}
#' }
#' @author Zhendong Guo (\email{guozhendong19@mails.ucas.ac.cn}).
#' @references
#' \itemize{
#'  \item{Broadhurst, D.; Goodacre, R.; Reinke, S. N.; Kuligowski, J.; Wilson, I. D.; Lewis, M. R.; Dunn, W. B. Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies. \emph{Metabolomics} \strong{2018}, 14 (6), 72. DOI: 10.1007/s11306-018-1367-3.}
#' }
#' @seealso \code{\link{D_ratio}}.
#' @examples
#' data(Dataset_I)
#' data <- Dataset_I
#' data.RF <- RF.correction(data)
#' compare_D_ratio(data, data.RF, BEC_method = c('raw', 'RF'))

compare_D_ratio <- function(data1, ..., BEC_method, line_alpha = 0.8){
  data <- list(data1, ...)

  if (length(BEC_method) != length(data)){
    stop("'BEC_method' should be correted.")
  }

  data.D_ratio <- lapply(1:length(data), function(x){
    D_ratio(data[[x]], print_plot = FALSE, D_ratio_value = TRUE)$D_ratio
  })

  D_ratio <- do.call(rbind, lapply(1:length(data), function(x){
    data.frame(D_ratio = data.D_ratio[[x]],
               Method = rep(BEC_method[x], length(data.D_ratio[[x]])))
  }))

  D_ratio$Method <- factor(D_ratio$Method, levels = BEC_method)

  # 绘制每列数据的cdf曲线
  ggplot2::ggplot(D_ratio, aes(x = D_ratio, color = Method)) +
    stat_ecdf(geom = "step", alpha = line_alpha) +
    labs(title = "D-ratio distribution", x = "D-ratio (%)", y = "Cumulative distribution (%)") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.4) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1),
                       labels = scales::label_number(scale = 100)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1),
                       labels = scales::label_number(scale = 100)) +
    theme_minimal()
}
