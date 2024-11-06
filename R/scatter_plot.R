#' @title Scatter plot
#'
#' @description A scatter plot for univariate visualization.
#'
#' @param data A dataframe.\cr
#' Execute \code{data(Dataset_I)} and \code{View(Dataset_I)} for formats.
#' @param variable_index A numeric scalar.
#' @param ylab Default is \code{NULL}, which denotes the variable index. Otherwise, the variable name can be input.
#' @param y.log10_trans Logical. Default is \code{FALSE}.
#' @param y.lim Logical. Default is \code{FALSE}.
#' @param y_min A numeric scalar. Input only when \code{y.lim = TRUE}.
#' @param y_max A numeric scalar. Input only when \code{y.lim = TRUE}.
#' @param point_size A numeric scalar. Default is 1.5. Determines the point size in the scatter plot.
#' @param point_alpha A numeric scalar. Default is 0.8. Determines the point opacity in the scatter plot.
#'
#' @return A scatter plot.
#'
#' @export
#' @import stats ggplot2
#'
#' @author Zhendong Guo (\email{guozhendong19@mails.ucas.ac.cn}).
#'
#' @seealso \code{\link{compare_scatter_plot}}.
#' @examples
#' data(Dataset_I)
#' data <- Dataset_I
#' scatter_plot(data, variable_index = 107, y.log10_trans = TRUE,
#'              y.lim = TRUE, y_min = 10^3, y_max = 10^5)

scatter_plot <- function(data, variable_index,
                         ylab = NULL, y.log10_trans = FALSE,
                         y.lim = FALSE, y_min, y_max,
                         point_size = 1.5, point_alpha = 0.8){
  batch.level <- unique(data$batch)
  B <- length(batch.level) # 批次数

  vline <- sapply(1:(B-1), function(b){
    data_b <- subset(data, data$batch == batch.level[b])
    return(max(data_b$order))
  })

  data$type[data$type == 'qc'] <- 'QC'
  data$type[data$type == 'sample'] <- 'subject'
  colnames(data)[colnames(data) == 'type'] <- 'Type'

  if (is.null(ylab) == TRUE){
    ylab <- sprintf("Intensity of Variable %.d", variable_index)
  }else {
    ylab <- paste0("Intensity of ", ylab)
  }

  scatter_plot <- ggplot2::ggplot(data, aes(x = order, y = data[, 4+variable_index],
                                            color = Type, shape = Type)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_shape_manual(values = c(19, 1)) +
    geom_vline(xintercept = vline + 0.5, linetype = "dashed", color = "black", linewidth = 0.4) +
    labs(x = "Injection order", y = ylab) +
    theme_classic()

  if (y.lim == TRUE){
    if (y_min >= y_max){
      stop("y_max must be larger than y_min.")
    }

    y.lim <- c(y_min, y_max)

    if (y.log10_trans == TRUE){
      scatter_plot <- scatter_plot +
        scale_y_log10(limits = y.lim)
    }else {
      scatter_plot <- scatter_plot +
        scale_y_continuous(limits = y.lim)
    }
  }else {
    if (y.log10_trans == TRUE){
      scatter_plot <- scatter_plot +
        scale_y_log10()
    }
  }

  print(scatter_plot)
}

#' @title Scatter plot of different data for comparison
#'
#' @description A scatter plot for univariate visualization.
#'
#' @param data1 A dataframe.\cr
#' Execute \code{data(Dataset_I)} and \code{View(Dataset_I)} for formats.
#' @param ... Alternatives. More data (e.g., \code{data2}, \code{data3}) can be input.
#' @param BEC_method A character vector. The length must be the same as \code{c(data1, ...)}. For example, \code{c("raw", "RF", "XGBoost")}.
#' @param variable_index A numeric scalar.
#' @param sample_type \code{qc} denotes QC samples; \code{sample} denotes subject samples.
#' @param ylab Default is \code{NULL}, which denotes the variable index. Otherwise, the variable name can be input.
#' @param y.log10_trans Logical. Default is \code{FALSE}.
#' @param y.lim Logical. Default is \code{FALSE}.
#' @param y_min A numeric scalar. Input only when \code{y.lim = TRUE}.
#' @param y_max A numeric scalar. Input only when \code{y.lim = TRUE}.
#' @param point_size A numeric scalar. Default is 1.5. Determines the point size in the scatter plot.
#' @param point_alpha A numeric scalar. Default is 0.8. Determines the point opacity in the scatter plot.
#'
#' @return A scatter plot.
#'
#' @export
#' @import stats ggplot2
#'
#' @author Zhendong Guo (\email{guozhendong19@mails.ucas.ac.cn}).
#'
#' @seealso \code{\link{scatter_plot}}.
#' @examples
#' data(Dataset_I)
#' data <- Dataset_I
#' data.RF <- RF.correction(data)
#' compare_scatter_plot(data, data.RF, BEC_method = c("raw", "RF"),
#'                      variable_index = 107, sample_type = "qc",
#'                      y.log10_trans = TRUE, y.lim = TRUE, y_min = 10^3, y_max = 10^5)

compare_scatter_plot <- function(data1, ...,
                                 BEC_method, variable_index,
                                 sample_type = c('qc', 'sample'),
                                 ylab = NULL, y.log10_trans = FALSE,
                                 y.lim = FALSE, y_min, y_max,
                                 point_size = 1.5, point_alpha = 0.8){
  data <- list(data1, ...)

  if (length(BEC_method) != length(data)){
    stop("'BEC_method' should be correted.")
  }

  batch.level <- unique(data1$batch)
  B <- length(batch.level) # 批次数

  vline <- sapply(1:(B-1), function(b){
    data_b <- subset(data1, data1$batch == batch.level[b])
    return(max(data_b$order))
  })

  if (sample_type == 'qc'){
    title <- "BEC for QC samples"
    Data <- lapply(1:length(data), function(x){
      subset(data[[x]], data[[x]]$type == 'qc')
    })
  }else {
    title <- "BEC for subject samples"
    Data <- lapply(1:length(data), function(x){
      subset(data[[x]], data[[x]]$type == 'sample')
    })
  }

  DATA <- as.data.frame(do.call(rbind, Data))
  Method <- rep(BEC_method, each = nrow(Data[[1]]))
  DATA <- cbind(DATA, Method = Method)
  DATA$Method <- factor(Method, levels = BEC_method)

  if (is.null(ylab) == TRUE){
    ylab <- sprintf("Intensity of Variable %.d", variable_index)
  }else {
    ylab <- paste0("Intensity of ", ylab)
  }

  scatter_plot <- ggplot2::ggplot(DATA, aes(x = order, y = DATA[, 4+variable_index],
                                  color = Method)) +
    geom_point(size = point_size, alpha = point_alpha) +
    geom_vline(xintercept = vline + 0.5, linetype = "dashed", color = "black", linewidth = 0.4) +
    labs(title = title, x = "Injection order", y = ylab) +
    theme_classic()

  if (y.lim == TRUE){
    if (y_min >= y_max){
      stop("y_max must be larger than y_min.")
    }

    y.lim <- c(y_min, y_max)

    if (y.log10_trans == TRUE){
      scatter_plot <- scatter_plot +
        scale_y_log10(limits = y.lim)
    }else {
      scatter_plot <- scatter_plot +
        scale_y_continuous(limits = y.lim)
    }
  }else {
    if (y.log10_trans == TRUE){
      scatter_plot <- scatter_plot +
        scale_y_log10()
    }
  }

  print(scatter_plot)
}
