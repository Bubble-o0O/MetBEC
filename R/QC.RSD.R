#' @title RSD distribution of data
#'
#' @description Relative Standard Deviation (RSD). The metric only examines QC samples, whose formula is:\cr
#' \deqn{\mathrm{RSD}_{i}=\frac{\sqrt{\mathrm{Var}(y^{\mathrm{(QC)}}_{i})}}{\mathrm{mean}(y^{\mathrm{(QC)}}_{i})},i=1,2,...,p,}\cr
#' where \eqn{y^{\mathrm{(QC)}}_{i}} denotes the intensity of QC samples' \eqn{i}th metabolite;\cr
#' \eqn{p} denotes the variable (metabolite) number.
#'
#' @param data A dataframe. \strong{Use \code{data(Dataset_I)} for formats.}
#' @param print_plot Logical. Default is \code{TRUE}. Determines whether to plot the histogram and the cumulative distribution plot.
#' @param RSD_value Logical. Default is \code{FALSE}. Determines whether to output the RSD values.
#'
#' @return
#' \item{`QC.RSD: 0 ~ 30\%`}{A numeric scalar. The cumulative frequency of RSD within 30\%.}
#' \item{`QC.RSD: 0 ~ 20\%`}{A numeric scalar. The cumulative frequency of RSD within 20\%.}
#' \item{`QC.RSD: 0 ~ 15\%`}{A numeric scalar. The cumulative frequency of RSD within 15\%.}
#' \item{QC.RSD}{A numeric vector. Output only when \code{RSD_value = TRUE}.}
#' \item{histogram}{Plotted only when \code{print_plot = TRUE}.}
#' \item{cumulative distribution plot}{Plotted only when \code{print_plot = TRUE}.}
#'
#' @export
#' @import stats ggplot2
#' @importFrom scales label_number
#'
#' @note
#' \itemize{
#'  \item{The acceptance criterion is 15\%, 20\%, and 30\%.}
#' }
#' @author Zhendong Guo (\email{guozhendong19@mails.ucas.ac.cn}).
#' @references
#' \itemize{
#'  \item{Broadhurst, D.; Goodacre, R.; Reinke, S. N.; Kuligowski, J.; Wilson, I. D.; Lewis, M. R.; Dunn, W. B. Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies. \emph{Metabolomics} \strong{2018}, 14 (6), 72. DOI: 10.1007/s11306-018-1367-3.}
#'  \item{Dudzik, D.; Barbas-Bernardos, C.; Garcia, A.; Barbas, C. Quality assurance procedures for mass spectrometry untargeted metabolomics. a review. \emph{J Pharm Biomed Anal} \strong{2018}, 147, 149-173. DOI: 10.1016/j.jpba.2017.07.044.}
#'  \item{Han, W.; Li, L. Evaluating and minimizing batch effects in metabolomics. \emph{Mass Spectrometry Reviews} \strong{2020}, 41 (3), 421-442. DOI: 10.1002/mas.21672.}
#'  \item{Martens, A.; Holle, J.; Mollenhauer, B.; Wegner, A.; Kirwan, J.; Hiller, K. Instrumental Drift in Untargeted Metabolomics: Optimizing Data Quality with Intrastudy QC Samples. \emph{Metabolites} \strong{2023}, 13 (5). DOI: 10.3390/metabo13050665.}
#' }
#' @seealso \code{\link{compare_QC.RSD}}.
#' @examples
#' data(Dataset_I)
#' data <- Dataset_I
#' QC.RSD(data)

QC.RSD <- function(data, print_plot = TRUE, RSD_value = FALSE){
  QC <- subset(data, data$type == 'qc')

  QC.Data <- QC[ ,-1:-4]
  p <- ncol(QC.Data) # 变量数

  mean <- apply(QC.Data, 2, mean)
  sd <- apply(QC.Data, 2, sd)
  RSD <- sd/mean
  names(RSD) <- colnames(QC.Data)

  if (sum(RSD >= 1) > 0){
    warning(sprintf("QC samples' RSD of %.d features is \u2265 100%%.",
                    sum(RSD >= 1)))
  }

  QC.RSD.0_30 <- length(which(RSD < 0.3))/p
  QC.RSD.0_20 <- length(which(RSD < 0.2))/p
  QC.RSD.0_15 <- length(which(RSD < 0.15))/p

  if (print_plot == TRUE){
    # 绘制直方图
    interval <- cut(RSD, breaks = c(0, 0.2, 0.3, Inf), labels = c("< 20%", "20% ~ 30%", "\u2265 30%"),
                    right = FALSE) # 左闭右开

    histogram <- ggplot2::ggplot(data.frame(RSD), aes(x = RSD, fill = interval)) +
      geom_histogram(binwidth = 0.1, boundary = 0, closed = 'left', # 左闭右开
                     color = "black", show.legend = TRUE) +
      labs(title = "Histogram of QC samples' RSD", x = "RSD (%)", y = "Count", fill = "RSD") +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1),
                         labels = scales::label_number(scale = 100)) +
      scale_fill_manual(values = c("< 20%" = "lightgreen",
                                   "20% ~ 30%" = "lightblue2",
                                   "\u2265 30%" = "lightpink2")) +
      theme_minimal()

    # 绘制累积分布图
    cdf <- ggplot2::ggplot(data.frame(RSD), aes(x = RSD)) +
      stat_ecdf(geom = "step", color = "black") +
      labs(title = "QC samples's RSD distribution", x = "RSD (%)", y = "Cumulative frequency (%)") +
      geom_vline(xintercept = c(0.15, 0.2, 0.3), linetype = "dashed", color = "black", linewidth = 0.4) +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1),
                         labels = scales::label_number(scale = 100)) +
      scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1),
                         labels = scales::label_number(scale = 100)) +
      theme_minimal()

    print(histogram)
    print(cdf)
  }

  if (RSD_value == FALSE){
    return(list("QC.RSD: 0 ~ 30%" = QC.RSD.0_30,
                "QC.RSD: 0 ~ 20%" = QC.RSD.0_20,
                "QC.RSD: 0 ~ 15%" = QC.RSD.0_15))
  }else if (RSD_value == TRUE){
    return(list("QC.RSD: 0 ~ 30%" = QC.RSD.0_30,
                "QC.RSD: 0 ~ 20%" = QC.RSD.0_20,
                "QC.RSD: 0 ~ 15%" = QC.RSD.0_15,
                QC.RSD = RSD))
  }
}

#' @title RSD distribution of different corrected data for comparison
#'
#' @description Relative Standard Deviation (RSD). The metric only examines QC samples, whose formula is:\cr
#' \deqn{\mathrm{RSD}_{i}=\frac{\sqrt{\mathrm{Var}(y^{\mathrm{(QC)}}_{i})}}{\mathrm{mean}(y^{\mathrm{(QC)}}_{i})},i=1,2,...,p,}\cr
#' where \eqn{y^{\mathrm{(QC)}}_{i}} denotes the intensity of QC samples' \eqn{i}th metabolite;\cr
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
#'  \item{The acceptance criterion is 15\%, 20\%, and 30\%.}
#' }
#' @author Zhendong Guo (\email{guozhendong19@mails.ucas.ac.cn}).
#' @references
#' \itemize{
#'  \item{Broadhurst, D.; Goodacre, R.; Reinke, S. N.; Kuligowski, J.; Wilson, I. D.; Lewis, M. R.; Dunn, W. B. Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies. \emph{Metabolomics} \strong{2018}, 14 (6), 72. DOI: 10.1007/s11306-018-1367-3.}
#'  \item{Dudzik, D.; Barbas-Bernardos, C.; Garcia, A.; Barbas, C. Quality assurance procedures for mass spectrometry untargeted metabolomics. a review. \emph{J Pharm Biomed Anal} \strong{2018}, 147, 149-173. DOI: 10.1016/j.jpba.2017.07.044.}
#'  \item{Han, W.; Li, L. Evaluating and minimizing batch effects in metabolomics. \emph{Mass Spectrometry Reviews} \strong{2020}, 41 (3), 421-442. DOI: 10.1002/mas.21672.}
#'  \item{Martens, A.; Holle, J.; Mollenhauer, B.; Wegner, A.; Kirwan, J.; Hiller, K. Instrumental Drift in Untargeted Metabolomics: Optimizing Data Quality with Intrastudy QC Samples. \emph{Metabolites} \strong{2023}, 13 (5). DOI: 10.3390/metabo13050665.}
#' }
#' @seealso \code{\link{QC.RSD}}.
#' @examples
#' data(Dataset_I)
#' data <- Dataset_I
#' data.RF <- RF.correction(data)
#' compare_QC.RSD(data, data.RF, BEC_method = c('raw', 'RF'))

compare_QC.RSD <- function(data1, ..., BEC_method, line_alpha = 0.8){
  data <- list(data1, ...)

  if (length(BEC_method) != length(data)){
    stop("'BEC_method' should be correted.")
  }

  data.RSD <- lapply(1:length(data), function(x){
    QC.RSD(data[[x]], print_plot = FALSE, RSD_value = TRUE)$QC.RSD
  })

  RSD <- do.call(rbind, lapply(1:length(data), function(x){
    data.frame(RSD = data.RSD[[x]],
               Method = rep(BEC_method[x], length(data.RSD[[x]])))
  }))

  RSD$Method <- factor(RSD$Method, levels = BEC_method)

  # 绘制每列数据的cdf曲线
  ggplot2::ggplot(RSD, aes(x = RSD, color = Method)) +
    stat_ecdf(geom = "step", alpha = line_alpha) +
    labs(title = "QC samples's RSD distribution", x = "RSD (%)", y = "Cumulative frequency (%)") +
    geom_vline(xintercept = c(0.15, 0.2, 0.3), linetype = "dashed", color = "black", linewidth = 0.4) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1),
                       labels = scales::label_number(scale = 100)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1),
                       labels = scales::label_number(scale = 100)) +
    theme_minimal()
}
