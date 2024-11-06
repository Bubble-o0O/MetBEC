PC_num.TW <- function(X, PCA_scale = TRUE,
                          TW_type = c('2001', '2006'), alpha = 0.01,
                          p_value.adjust = c('Bonferroni-Holm', 'none')){
  ##### Tracy-Widom test #####
  # 原始文献：Johnstone (2001)；Johnstone (2006)
  # 相关文献：Saccenti (2015b)；Saccenti (2017)；Fockman(2019)
  # 备注：Johnstone (2001)、Saccenti (2015b)与Saccenti (2017)均提及TW检验在处理标准化数据时存在缺陷
  # 此外，Patterson (2006)将TW检验应用于基因组学

  # 经测试，Johnstone (2001)与Johnstone (2006)的结果无明显差异，Fockman (2019)也得到了类似结论
  # Johnstone (2001)和Saccenti (2015b)采用alpha = 0.01，
  # Saccenti (2017)和Fockman (2019)采用alpha = 0.05，
  # 但如果进行p-value校正，测试结果将更好

  ##### 默认参数：数据标准化；2001；0.01；Bonferroni-Holm #####
  n <- nrow(X)
  p <- ncol(X)
  r <- min(n-1, p) # 主成分的最大数量（即：协方差矩阵的秩）

  if (length(TW_type) > 1){ # 默认选择Johnstone (2001)
    TW_type <- '2001'
  }

  if (length(p_value.adjust) > 1){ # 默认选择Bonferroni-Holm校正法
    p_value.adjust <- 'Bonferroni-Holm'
  }

  if (PCA_scale == TRUE){
    X_autoscale <- scale(X) # 标准化数据
    l <- (svd(X_autoscale)$d^2)[1:r]/(n - 1)
    # 协方差矩阵的非零特征值（从大到小排列）；对于标准化数据，sum(l) = p
    pca <- prcomp(X, scale. = TRUE)
  }else {
    X_center <- scale(X, scale = FALSE) # 中心化数据
    l <- (svd(X_center)$d^2)[1:r]/(n - 1)
    pca <- prcomp(X, scale. = FALSE)
  }

  if ('2001' %in% TW_type){ # Johnstone (2001) #
    mu_np <- function(n, p){
      if (n > p){
        mu_np <- (sqrt(n-1) + sqrt(p))^2

      }else {
        mu_np <- (sqrt(p-1) + sqrt(n))^2
      }

      return(mu_np)
    }

    sigma_np <- function(n, p){
      if (n > p){
        sigma_np <- (sqrt(n-1) + sqrt(p)) *
          (1/sqrt(n-1) + 1/sqrt(p))^(1/3)

      }else {
        sigma_np <- (sqrt(p-1) + sqrt(n)) *
          (1/sqrt(p-1) + 1/sqrt(n))^(1/3)
      }

      return(sigma_np)
    }

  }else if ('2006' %in% TW_type){ # Johnstone (2006) #
    mu_np <- function(n, p){
      mu_np <- (sqrt(n-0.5) + sqrt(p-0.5))^2

      return(mu_np)
    }

    sigma_np <- function(n, p){
      sigma_np <- (sqrt(n-0.5) + sqrt(p-0.5)) *
        (1/sqrt(n-0.5) + 1/sqrt(p-0.5))^(1/3)

      return(sigma_np)
    }
  }

  sigma_sq_hat_k <- function(n, p, k){
    if (n > p){
      sigma_sq_hat_k <- 0
      for (i in (k+1):r){ # 注意：这里的k+1必须加括号
        sigma_sq_hat_k <- sigma_sq_hat_k + l[i] / (n*(p-k))
      }

    }else {
      sigma_sq_hat_k <- 0
      for (i in (k+1):r){
        sigma_sq_hat_k <- sigma_sq_hat_k + l[i] / (p*(n-k))
      }
    }

    return(sigma_sq_hat_k)
  }

  # 进行Tracy-Widom检验（右尾检验）
  if (PCA_scale == TRUE){ # 与开头的if语句连通
    set.seed(123)
    r.sq <- rchisq(p, n) # 生成p个服从自由度为n的卡方随机数
    X_tilde <- X_autoscale %*% diag(sqrt(r.sq)) # 将r.sq开方后，组成对角阵
    l <- (svd(X_tilde)$d^2)[1:r]/(n - 1)
    # 更新autoscale条件时的l；PCA_scale == FALSE时跳过该if语句
  }

  for (k in 1:(r-1)){
    if ('Bonferroni-Holm' %in% p_value.adjust){
      alpha_k <- alpha / k # Bonferroni-Holm法校正（参考Dray (2008)）
      # 校正后，越往后的主成分，s_alpha越大，检验越严格

    }else if ('none' %in% p_value.adjust){
      alpha_k <- alpha
    }

    s_alpha <- RMTstat::qtw(alpha_k, lower.tail = FALSE) # alpha_k表示上分位数

    if (n > p){
      if (l[k] <= sigma_sq_hat_k(n, p, k) *
          (mu_np(n, p-k) + s_alpha * sigma_np(n, p-k))){
        PC_num.TW <- k - 1
        break
      }else {
        PC_num.TW <- r - 1
      }

    }else {
      if (l[k] <= sigma_sq_hat_k(n, p, k) *
          (mu_np(n-k, p) + s_alpha * sigma_np(n-k, p))){
        PC_num.TW <- k - 1
        break
      }else {
        PC_num.TW <- r - 1
      }
    }
  }

  if (PC_num.TW > 0){
    R2 <- pca$sdev^2/sum(pca$sdev^2)
    R2.TW <- sum(R2[1:PC_num.TW])
  }else {
    R2.TW <- 0
  }

  #### 输出结果 ####
  return(list(PC_num.TW = PC_num.TW,
              R2.TW = R2.TW))
}

# 基于PCA结合Hotelling's T2统计量（D统计量）与平方预测误差（SPE；DModX；Q统计量）检测离群点
# 由于当n <= p时，样本协方差矩阵不可逆，因此需要PCA降维，使得n > k
# 根据模拟数据的试验，在此推荐Tracy-Widom分布选择最佳主成分数的算法（参考Saccenti, 2015)

#' @title Data-preprocessing: outlier detection by PCA for QC samples
#'
#' @description For outlier detection, multivariate statistical methods, such as Principal Component Analysis (PCA), are more powerful than univariate statistical methods due to the consideration for the correlations between variables.\cr
#' Here, we use the Tracy-Widom distribution to determine the principal component number of PCA.\cr
#' Then, we use the 95\% or 99\% upper control limit (UCL) of Hotelling's \eqn{T^2} statistic (also termed \eqn{D}-statistic) and Squared Prediction Errors (SPEs; also termed \eqn{Q}-statistic or DModX) to detect outliers.
#'
#' @param data A dataframe.\cr
#' Execute \code{data(Dataset_I)} and \code{View(Dataset_I)} for formats.
#' @param title Default is \code{NULL}. Otherwise, the batch name can be input.
#' @param PCA_scale Logical. Default is \code{TRUE}, which determines whether to autoscale (standardize) \code{data} and is unnecessary to be tuned for most of the time.
#' @param T2.distribution Default is \code{chisq}. Hotelling's \eqn{T^2} statistic follows four kinds of distribution in different fields, whose strictness ranking (contrary to the size of the confidence ellipse) is:\cr
#' \code{Beta} > \code{chisq} > \code{T2.1} > \code{T2.2}.
#' @param print_plot Logical. Default is \code{TRUE}. Determines whether to plot the Hotelling's \eqn{T^2} chart and the SPE chart.
#'
#' @return
#' \item{PC_num}{A numeric scalar. Determined by the Tracy-Widom distribution.}
#' \item{R2}{A numeric scalar. Cumulative \eqn{R^2} of the top \code{PC_num} principal components.}
#' \item{Hotelling_T2}{A list. Includes "95\%UCL ~ 99\%UCL" and "> 99\%UCL".}
#' \item{SPE}{A list. Includes "95\%UCL ~ 99\%UCL" and "> 99\%UCL".}
#' \item{intersect}{A list. The intersection of \code{Hotelling_T2} and \code{SPE}. Includes "95\%UCL ~ 99\%UCL" and "> 99\%UCL".}
#' \item{union}{A list. The union of \code{Hotelling_T2} and \code{SPE}. Includes "95\%UCL ~ 99\%UCL" and "> 99\%UCL".}
#'
#' @export
#' @import stats ggplot2
#' @importFrom RMTstat qtw
#'
#' @note
#' \itemize{
#'  \item{The principle is based on the premise of the multivariate normal distribution. Hence, the function is suitable for QC samples instead of subject samples. Additionally, it should be used for QC samples batchwise instead of all the QC samples.}
#'  \item{It is an open problem to determine the principal component number of PCA. One of the feasible methods is the Tracy-Widom distribution, which is based on the Random Matrix Theory (RMT).}
#'  \item{Although Hotelling's \eqn{T^2} statistic strictly follows the Beta distribution, the Chi-squared distribution and the Hotelling's \eqn{T^2} distribution are still widely applied.}
#'  \item{By default, the samples whose \code{Hotelling_T2} and \code{SPE} are both larger than the 99\% UCL, namely \code{union}, are regarded as outliers.}
#' }
#' @author Zhendong Guo (\email{guozhendong19@mails.ucas.ac.cn}).
#' @references
#' \itemize{
#'  \item{Wikstrom, C.; Albano, C.; Eriksson, L.; Friden, H.; Johansson, E.; Nordahl, A.; Rannar, S.; Sandberg, M.; Kettaneh-Wold, N.; Wold, S. Multivariate process and quality monitoring applied to an electrolysis process Part I. Process supervision with multivariate control charts. \emph{Chemometrics and Intelligent Laboratory Systems} \strong{1998}, 42, 221-231.}
#'  \item{Saccenti, E.; Camacho, J. Determining the number of components in principal components analysis: A comparison of statistical, crossvalidation and approximated methods. \emph{Chemometrics and Intelligent Laboratory Systems} \strong{2015}, 149, 99-116. DOI: 10.1016/j.chemolab.2015.10.006.}
#'  \item{Blaise, B. J.; Correia, G. D. S.; Haggart, G. A.; Surowiec, I.; Sands, C.; Lewis, M. R.; Pearce, J. T. M.; Trygg, J.; Nicholson, J. K.; Holmes, E.; et al. Statistical analysis in metabolic phenotyping. \emph{Nature Protocols} \strong{2021}, 16 (9), 4299-4326. DOI: 10.1038/s41596-021-00579-1.}
#' }
#' @seealso \code{\link{PCA_plot}}, \code{stats::\link[stats]{prcomp}}.
#' @examples
#' data(Dataset_I.with_outliers)
#' data <- Dataset_I.with_outliers
#'
#' if (length(data$order) != length(unique(data$order))){
#'   warning("The repeated values in 'order' should be corrected.")
#' }
#'
#' # Obtain the batch names and the batch number.
#' batch.level <- unique(data$batch)
#' B <- length(batch.level)
#'
#' # The variance of each variable should not be 0. Otherwise, eliminate them.
#' zero_var.index <- unique(unlist(lapply(1:B, function(b){
#'   QC_b <- subset(data, data$batch == batch.level[b] & data$type == "qc")
#'   return(which(apply(QC_b[, -1:-4], 2, var) == 0))
#' })))
#'
#' if (length(zero_var.index) > 1){
#'   data <- data[, -(4 + zero_var.index)]
#' }
#'
#' # Outlier detection for QC samples batchwise.
#' outlier_result <- lapply(1:B, function(b){
#'   QC_b <- subset(data, data$batch == batch.level[b] & data$type == "qc")
#'   return(outlier_detection_by_PCA(QC_b, title = batch.level[b]))
#' })
#' names(outlier_result) <- batch.level
#'
#' outliers <- unlist(sapply(1:B, function(b){
#'   which(data$order == outlier_result[[b]]$union$`> 99%UCL`$order)
#' }))
#'
#' # Obtain the data without outliers.
#' data <- data[-outliers, ]

outlier_detection_by_PCA <- function(data, title = NULL,
                                     PCA_scale = TRUE,
                                     T2.distribution = c("Beta", "chisq", "T2.1", "T2.2"),
                                     print_plot = TRUE){
  order <- data$order
  name <- data$name

  X <- data[, -1:-4]

  if (any(apply(X, 2, var) == 0) == TRUE){
    stop("Some variances are 0.")
  }

  n <- nrow(X) # 样本量
  p <- ncol(X) # 变量数

  X_scale <- scale(X, scale = PCA_scale)

  result.TW <- PC_num.TW(X, PCA_scale = PCA_scale)

  k <- result.TW$PC_num.TW
  R2 <- result.TW$R2.TW

  if (PCA_scale == TRUE){ # 标准化
    pca <- prcomp(X, scale. = TRUE)
  }else { # 中心化
    pca <- prcomp(X, scale. = FALSE)
  }

  scores <- as.matrix(pca$x[, 1:k]) # 如果没有as.matrix，当k=1时后续将报错

  if (k == 0){
    stop("The PC number is 0, which suggests that
        the correlations between variables are too weak so that it fails to detect outliers.")
  }else if (k == 1){
    S_inv <- 1 / pca$sdev[1]^2
  }else {
    S_inv <- solve(diag(pca$sdev[1:k]^2))
  }

  #### Hotelling's T2（D统计量；马氏距离的平方） ####
  if (length(T2.distribution) > 1){
    T2.distribution <- "chisq"
    # 备注1：其余三种分布的极限分布（n→∞）均为卡方分布
    # 备注2：四种分布的严格程度从大到小依次为：Beta > chisq > T2.1 > T2.2
    # 主要参考文献：Wise (1990)、Tracy (1992)
  }

  T2 <- sapply(1:n, function(i)
    t(scores[i, 1:k]) %*% S_inv %*% scores[i, 1:k]
  )

  if (T2.distribution == 'chisq'){
    T2_UCL.99 <- qchisq(0.99, k)
    T2_UCL.95 <- qchisq(0.95, k)
  }else if (T2.distribution == "Beta"){ # 严格成立；其余三种分布均为近似成立
    T2_UCL.99 <- (n-1)^2/n * qbeta(0.99, k/2, (n-k-1)/2)
    T2_UCL.95 <- (n-1)^2/n * qbeta(0.95, k/2, (n-k-1)/2)
  }else if (T2.distribution == "T2.1"){
    T2_UCL.99 <- k*(n-1)/(n-k) * qf(0.99, k, n-k)
    T2_UCL.95 <- k*(n-1)/(n-k) * qf(0.95, k, n-k)
  }else if (T2.distribution == "T2.2"){ # 对于MSPC的Phase II严格成立
    T2_UCL.99 <- k*(n+1)*(n-1)/(n*(n-k)) * qf(0.99, k, n-k)
    T2_UCL.95 <- k*(n+1)*(n-1)/(n*(n-k)) * qf(0.95, k, n-k)
  }

  D_statistic.Data <- data.frame(order, T2)

  T2_fill <- rep('lightblue2', n)
  for (i in 1:n){
    if (T2[i] > T2_UCL.95 & T2[i] <= T2_UCL.99){
      T2_fill[i] <- 'blue'
    }else if (T2[i] > T2_UCL.99){
      T2_fill[i] <- 'red'
    }
  }

  if (is.null(title) == TRUE){
    title.T2 <- "Hotelling's T\u00B2 control chart"
  }else {
    title.T2 <- paste0("Hotelling's T\u00B2 control chart for Batch ", title)
  }

  plot.T2 <- ggplot2::ggplot(D_statistic.Data, aes(x = factor(order, levels = order), y = T2)) + # 按顺序排列
    geom_bar(stat = "identity", fill = T2_fill) +
    labs(x = "Injection order", y = "T\u00B2 (D-statistic)", title = title.T2) +
    geom_hline(yintercept = T2_UCL.95, linetype = "dashed", color = "blue") +
    geom_hline(yintercept = T2_UCL.99, linetype = "dashed", color = "red") +
    annotate("text", x = 0, y = T2_UCL.95, label = "95%UCL", color = "blue", size = 3, hjust = -0.2, vjust = 1.2) +
    annotate("text", x = 0, y = T2_UCL.99, label = "99%UCL", color = "red", size = 3, hjust = -0.2, vjust = 1.2) +
    theme_minimal()
  # hjust与vjust是根据比例调整

  T2_index.95_99 <- which(T2 > T2_UCL.95 & T2 <= T2_UCL.99)
  T2_count.95_99 <- length(T2_index.95_99)
  T2_order.95_99 <- order[T2_index.95_99]
  T2_name.95_99 <- name[T2_index.95_99]
  D_statistic.95_99 <- list(count = T2_count.95_99,
                            order = T2_order.95_99,
                            name = T2_name.95_99)

  T2_index.99 <- which(T2 > T2_UCL.99)
  T2_count.99 <- length(T2_index.99)
  T2_order.99 <- order[T2_index.99]
  T2_name.99 <- name[T2_index.99]
  D_statistic.99 <- list(count = T2_count.99,
                         order = T2_order.99,
                         name = T2_name.99)

  D_statistic <- list("95%UCL ~ 99%UCL" = D_statistic.95_99,
                      "> 99%UCL" = D_statistic.99)

  #### 平方预测误差（Q统计量；DModX） ####
  # Squared prediction error (SPE)
  # 主要参考文献：MacGregor (1995)、Ferrer (2008)
  # 在实际应用中，D统计量与Q统计量均近似服从于卡方分布；参考：Kerkhof (2013)与Braga (2018)

  # 经济型SVD
  SVD <- svd(X_scale)

  U <- SVD$u # n×k矩阵；U'U = I_k（k阶单位阵）
  S <- diag(SVD$d) # k×k对角阵
  P <- SVD$v # k×k正交矩阵（载荷矩阵）

  # 截断型SVD
  U_nk <- matrix(c(U[, 1:k]), ncol = k)
  S_kk <- matrix(c(S[1:k, 1:k]), ncol = k)
  P_pk <- matrix(c(P[, 1:k]), ncol = k)

  X_scale_hat <- U_nk %*% S_kk %*% t(P_pk)

  # 残差矩阵
  E <- X_scale - X_scale_hat

  SPE <- sapply(1:n, function(i){
    e_i <-  E[i,]
    return(sum(e_i^2))
  })

  S2 <- var(SPE)
  mean <- mean(SPE)

  SPE_UCL.99 <- S2/(2*mean) * qchisq(0.99, 2*mean^2/S2)
  SPE_UCL.95 <- S2/(2*mean) * qchisq(0.95, 2*mean^2/S2)

  Q_statistic.Data <- data.frame(order, SPE)

  SPE_fill <- rep('lightblue2', n)
  for (i in 1:n){
    if (SPE[i] > SPE_UCL.95 & SPE[i] <= SPE_UCL.99){
      SPE_fill[i] <- 'blue'
    }else if (SPE[i] > SPE_UCL.99){
      SPE_fill[i] <- 'red'
    }
  }

  if (is.null(title) == TRUE){
    title.SPE <- "SPE control chart"
  }else {
    title.SPE <- paste0("SPE control chart for Batch ", title)
  }

  plot.SPE <- ggplot2::ggplot(Q_statistic.Data, aes(x = factor(order, levels = order), y = SPE)) + # 按顺序排列
    geom_bar(stat = "identity", fill = SPE_fill) +
    labs(x = "Injection order", y = "SPE (Q-statistic)", title = title.SPE) +
    geom_hline(yintercept = SPE_UCL.95, linetype = "dashed", color = "blue") +
    geom_hline(yintercept = SPE_UCL.99, linetype = "dashed", color = "red") +
    annotate("text", x = 0, y = SPE_UCL.95, label = "95%UCL", color = "blue", size = 3, hjust = -0.2, vjust = 1.2) +
    annotate("text", x = 0, y = SPE_UCL.99, label = "99%UCL", color = "red", size = 3, hjust = -0.2, vjust = 1.2) +
    theme_minimal()
  # hjust与vjust是根据比例调整

  SPE_index.95_99 <- which(SPE > SPE_UCL.95 & SPE <= SPE_UCL.99)
  SPE_count.95_99 <- length(SPE_index.95_99)
  SPE_order.95_99 <- order[SPE_index.95_99]
  SPE_name.95_99 <- name[SPE_index.95_99]
  Q_statistic.95_99 <- list(count = SPE_count.95_99,
                            order = SPE_order.95_99,
                            name = SPE_name.95_99)

  SPE_index.99 <- which(SPE > SPE_UCL.99)
  SPE_count.99 <- length(SPE_index.99)
  SPE_order.99 <- order[SPE_index.99]
  SPE_name.99 <- name[SPE_index.99]
  Q_statistic.99 <- list(count = SPE_count.99,
                         order = SPE_order.99,
                         name = SPE_name.99)

  Q_statistic <- list("95%UCL ~ 99%UCL" = Q_statistic.95_99,
                      "> 99%UCL" = Q_statistic.99)

  # 交集
  intersect_order.95_99 <- intersect(T2_order.95_99, SPE_order.95_99)
  intersect_order.99 <- intersect(T2_order.99, SPE_order.99)

  intersect_name.95_99 <- intersect(T2_name.95_99, SPE_name.95_99)
  intersect_name.99 <- intersect(T2_name.99, SPE_name.99)

  intersect.95_99 = list(order = intersect_order.95_99,
                         name = intersect_name.95_99)
  intersect.99 = list(order = intersect_order.99,
                      name = intersect_name.99)

  intersect <- list("95%UCL ~ 99%UCL" = intersect.95_99,
                    "> 99%UCL" = intersect.99)

  # 并集
  union_order.95_99 <- union(T2_order.95_99, SPE_order.95_99)
  union_order.99 <- union(T2_order.99, SPE_order.99)

  union_name.95_99 <- union(T2_name.95_99, SPE_name.95_99)
  union_name.99 <- union(T2_name.99, SPE_name.99)

  union.95_99 = list(order = union_order.95_99,
                     name = union_name.95_99)
  union.99 = list(order = union_order.99,
                  name = union_name.99)

  union <- list("95%UCL ~ 99%UCL" = union.95_99,
                "> 99%UCL" = union.99)

  if (print_plot == TRUE){
    print(plot.T2)
    print(plot.SPE)
  }

  return(list(PC_num = k,
              R2 = R2,
              Hotelling_T2 = D_statistic,
              SPE = Q_statistic,
              intersect = intersect,
              union = union))
}
