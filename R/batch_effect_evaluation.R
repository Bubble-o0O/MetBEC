Sigma_F.sq <- function(X){
  n <- nrow(X) # 样本量

  sum <- sum(sapply(1:n, function(i)
    sum((X[i,] - colMeans(X))^2)^2))

  K <- sum/(n-1)

  Sigma_F.sq <- (n-1)/(n*(n-2)*(n-3)) *
    ((n-1)*(n-2) * sum(diag(cov(X) %*% cov(X))) +
       sum(diag(cov(X)))^2 - n * K)

  return(Sigma_F.sq)
}

HN_2018 <- function(X, Y){
  n <- nrow(X)
  m <- nrow(Y)
  p <- ncol(X)

  TrS2X <- Sigma_F.sq(X)
  TrS2Y <- Sigma_F.sq(Y)

  diff_mean2 <- sum((colMeans(X)-colMeans(Y))^2) - sum(diag(cov(X)))/n - sum(diag(cov(Y)))/m

  diff_cov2 <- TrS2X + TrS2Y - 2*sum(diag(cov(X) %*% cov(Y)))

  sigma102 <- 2*TrS2X/(n*(n-1)) + 2*TrS2Y/(m*(m-1)) + 4*sum(diag(cov(X)%*%cov(Y)))/(n*m)
  # sigma102 <- 2*TrS2X/(n^2)+2*TrS2Y/(m^2)+4*sum(diag(cov(X)%*%cov(Y)))/(n*m)
  # 推荐使用n*(n-1)的版本，参考Chen & Qin (2010)与Yu (2023)

  sigma202 <- 4*(TrS2X^2)/(n^2) + 4*(TrS2Y^2)/(m^2) + 8*(sum(diag(cov(X)%*%cov(Y)))^2)/(n*m)
  # 该估计量与Li & Chen (2012)与Yu (2023)不一致，但均为一致估计量

  T1 <- diff_mean2/sqrt(sigma102)
  T2 <- diff_cov2/sqrt(sigma202)
  Tn <- (T1+T2)/sqrt(2)

  return(list(T1 = T1,
              mean_pval = 1 - pnorm(T1,0,1),
              T1 = T2,
              cov_pval = 1 - pnorm(T2,0,1),
              StatTn = Tn,
              pTn = 1 - pnorm(Tn,0,1)))
}



simul_test <- function(X, Y,
                      simul_method = c(NULL, 'cauchy', 'fisher', 'HN')){
  if (nrow(X) < 4 || nrow(Y) < 4){
    stop("The sample size of each group must be larger than 3.")
  }

  if (length(simul_method) > 1){ # fisher的功效略高于cauchy
    simul_method <- NULL # 当各批次的QC样本量的平均值不大于10时，HN的功效更高；其余情形时fisher的功效更高
  }

  if (is.null(simul_method) == TRUE){
    if (mean(nrow(X), nrow(Y)) < 10){
      HN <- HN_2018(X, Y)

      mean_pval <- HN$mean_pval
      cov_pval <- HN$cov_pval
      simul_pval <- HN$pTn
    }else {
      mean_pval <- PEtests::meantest.cq(X, Y)$pval
      cov_pval <- PEtests::covtest.lc(X, Y)$pval
      simul_pval <- PEtests::simultest(X, Y, method = 'fisher')$pval
    }
  }else if (simul_method == 'HN'){
    HN <- HN_2018(X, Y)

    mean_pval <- HN$mean_pval
    cov_pval <- HN$cov_pval
    simul_pval <- HN$pTn
  }else {
    mean_pval <- PEtests::meantest.cq(X, Y)$pval
    cov_pval <- PEtests::covtest.lc(X, Y)$pval
    simul_pval <- PEtests::simultest(X, Y, method = simul_method)$pval
  }

  return(list(mean_p_value = mean_pval,
              cov_p_value = cov_pval,
              simul_p_value = simul_pval))
}



#' @title QC-based simultaneous tests (QC-ST)
#'
#' @description The metric only examines QC samples.\cr
#' Considering the homogeneity of all the QC samples, we assume that if without batch effects, QC samples across different batches should follow the same multivariate normal distribution.\cr
#' QC-ST employs the simultaneous test of high-dimensional mean vectors and covariance matrices, which can be applied to the QC samples across different batches.\cr
#' We have verified that QC-ST is competent for batch effect evaluation and correction assessment in metabolomics.
#'
#' @param data A dataframe. \strong{Use \code{data(Dataset_I)} for formats.}
#' @param print_plot Logical. Default is \code{TRUE}. Determines whether to plot the heatmap and the undirected graph.
#' @param sig_level.alpha A numeric scalar. Default is 0.05, which is unnecessary to be tuned for most of the time.
#' @param simul_method Default is \code{NULL}, which selects \code{fisher} when \eqn{\mathrm{mean}\{n_{1},n_{2}\} \ge 10} or else selects \code{HN}.\cr
#' \code{cauchy}, namely "Yu-Cauchy"; \code{fisher}, namely "Yu-Fisher".
#' @param p_value.adjust Default is \code{fdr}, which is equivalent to \code{BH} and is unnecessary to be tuned for most of the time.
#'
#' @return
#' \item{adjacency_matrix}{A matrix. 0 denotes no significant batch effects between the corresponding two batches. 1 denotes significant batch effects between the corresponding two batches.}
#' \item{`the pairs without significant batch effects`}{A character vector. Corresponds to the 0 elements of \code{adjacency_matrix}.}
#' \item{mean_test}{A dataframe. From the homogeneity tests of mean vectors.}
#' \item{cov_test}{A dataframe. From the homogeneity tests of covariance matrices.}
#' \item{simul_test}{A dataframe. From the simultaneous tests.}
#' \item{heatmap}{Plotted only when \code{print_plot = TRUE}. The red block means significant batch effects between the two batches; The blue block means no significant batch effects between the two batches.}
#' \item{undirected graph}{Plotted only when \code{print_plot = TRUE}. Existing no edge means significant batch effects between the two batches; Existing an edge means no significant batch effects between the two batches.}
#'
#' @export
#' @import stats PEtests
#' @importFrom igraph graph_from_adjacency_matrix
#'
#' @note
#' \itemize{
#'  \item{We recommend that the QC sample size of each batch should not be less than 10 to ensure high statistical powers.}
#'  \item{If \code{simul_test$q_value} \eqn{\ge \alpha_{\mathrm{sig}}}, take it as the standard, and ignore \code{mean_test} and \code{cov_test}.\cr
#'  Actually, the similar situation occurs in ANalysis Of Variance (ANOVA), where the \eqn{F}-test should be considered prior to the multiple tests.}
#'  \item{If \code{simul_test$q_value} \eqn{< \alpha_{\mathrm{sig}}}, \code{mean_test$q_value} and \code{cov_test$q_value} can be used to further determine which parameter does.}
#'  \item{If QC-ST suggests that batch effects between some two batches are significant after \code{\link{batch_ratio.correction}} (\code{method = ratio-A} or \code{mean}), the covariance matrices must have statistical significance.}
#'  \item{The principles of \code{mean_test} and \code{cov_test} corresponding to \code{cauchy} and \code{fisher} are the same, but slightly different from those corresponding to \code{HN}.}
#' }
#' @author Zhendong Guo (\email{guozhendong19@mails.ucas.ac.cn}).
#' @references
#' \itemize{
#'  \item{Chen, S. X.; Qin, Y.-L. A two-sample test for high-dimensional data with applications to gene-set testing. \emph{The Annals of Statistics} \strong{2010}, 38 (2). DOI: 10.1214/09-aos716.}
#'  \item{Li, J.; Chen, S. X. Two sample tests for high-dimensional covariance matrices. \emph{The Annals of Statistics} \strong{2012}, 40 (2). DOI: 10.1214/12-aos993.}
#'  \item{Hyodo, M.; Nishiyama, T. A simultaneous testing of the mean vector and the covariance matrix among two populations for high-dimensional data. \emph{Test} \strong{2017}, 27 (3), 680-699. DOI: 10.1007/s11749-017-0567-x.}
#'  \item{Miao, R.; Xu, K. Joint test for homogeneity of high-dimensional means and covariance matrices using maximum-type statistics. \emph{Communications in Statistics - Simulation and Computation} \strong{2022}, 53 (2), 972-992. DOI: 10.1080/03610918.2022.2037641.}
#'  \item{Yu, X.; Li, D.; Xue, L.; Li, R. Power-Enhanced Simultaneous Test of High-Dimensional Mean Vectors and Covariance Matrices with Application to Gene-Set Testing. \emph{Journal of the American Statistical Association} \strong{2022}, 118 (544), 2548-2561. DOI: 10.1080/01621459.2022.2061354.}
#' }
#' @seealso \code{PEtests::\link[PEtests]{meantest.cq}}, \code{PEtests::\link[PEtests]{covtest.lc}}, \code{PEtests::\link[PEtests]{simultest}}.
#' @examples
#' data(Dataset_I)
#' data <- Dataset_I
#' data.QC_ST <- QC_ST(data)

QC_ST <- function(data, print_plot = TRUE,
                  sig_level.alpha = 0.05,
                  simul_method = c(NULL, 'cauchy', 'fisher', 'HN'),
                  p_value.adjust = c("bonferroni", "holm", "hochberg", "hommel",
                                     "fdr", "BH", "BY", "none")){
  if (length(simul_method) > 1){
    simul_method <- NULL
  }

  if (length(p_value.adjust) > 1){ # 默认选择BH法（FDR）校正
    p_value.adjust <- 'fdr' # 与'BH'一致
  }

  batch.level <- unique(data$batch)
  B <- length(batch.level) # 批次数

  QC <- subset(data, data$type == 'qc')
  QC[, -1:-4] <- scale(QC[, -1:-4]) # 标准化

  QC_list <- lapply(1:B, function(b) subset(QC, batch == batch.level[b])[, -1:-4])

  mean_pval <- c()
  cov_pval <- c()
  simul_pval <- c()

  pairs <- c()
  count <- 1

  for (i in 1:(B-1)){
    for (j in (i+1):B){
      X <- as.matrix(QC_list[[i]])
      Y <- as.matrix(QC_list[[j]])

      I <- batch.level[i]
      J <- batch.level[j]

      pairs[count] <- paste0(I, ' vs ', J)

      result <- simul_test(X, Y, simul_method = simul_method)

      mean_pval <- c(mean_pval, result$mean_p_value)
      cov_pval <- c(cov_pval, result$cov_p_value)
      simul_pval <- c(simul_pval, result$simul_p_value)

      count <- count + 1
    }
  }

  # p值校正
  mean_pval.adjust <- p.adjust(mean_pval, method = p_value.adjust)
  cov_pval.adjust <- p.adjust(cov_pval, method = p_value.adjust)
  simul_pval.adjust <- p.adjust(simul_pval, method = p_value.adjust)

  simul_count <- which(simul_pval.adjust >= sig_level.alpha)

  mean_test <- data.frame(pairs = pairs, q_value = mean_pval.adjust)
  cov_test <- data.frame(pairs = pairs, q_value = cov_pval.adjust)
  simul_test <- data.frame(pairs = pairs, q_value = simul_pval.adjust)



  ##### 绘制批次关系图（热图） #####
  entries <- rep(NA, B^2)
  mat <- matrix(entries, nrow = B,
                ncol = B)

  rownames(mat) <- batch.level
  colnames(mat) <- batch.level

  COUNT <- 1
  for (i in 1:(B-1)){
    for (j in (i+1):B){
      ifelse(simul_pval.adjust[COUNT] < sig_level.alpha,
             mat[i,j] <- 0, mat[i,j] <- 1)
      # 若批次效应显著，赋值为0；这样设置就与图论的邻接矩阵保持一致，即1为有边，0为无边
      COUNT <- COUNT + 1
    }
  }

  if (print_plot == TRUE){
    # 定义颜色映射
    if (any(mat == 0, na.rm = TRUE)){ # 检查矩阵中是否存在0元素，并忽略NA值
      my_colors <- c("red", "blue") # 存在0元素
    }else {
      my_colors <- "blue" # 不存在0元素
    }

    Mat <- mat[-B, -1]
    colnames(Mat) <- batch.level[-B]
    rownames(Mat) <- batch.level[-1]

    heatmap(t(Mat), # 只能是矩阵
            Rowv = NA, Colv = NA, # 不显示行和列的聚类
            col = my_colors, # 颜色映射
            symm = TRUE,
            cexRow = 1, # 设置行标签的大小
            cexCol = 1) # 设置列标签的大小

    # 添加颜色注释
    # legend("topleft", legend = c("Without significant batch effects", "With significant batch effects"),
    #        fill = c("blue", "red"), cex = 0.7)
  }

  ##### 绘制批次关系图（无向图） #####
  # 邻接矩阵（用于绘制无向图）
  mat[is.na(mat) == TRUE] <- 0
  adj_mat <- mat + t(mat)

  if (print_plot == TRUE){
    # 将邻接矩阵转换为无向图
    graph <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected")

    vertex_labels <- batch.level

    set.seed(123)
    # 可视化图
    plot(graph,
         vertex.label = vertex_labels,  # 设置顶点标签
         vertex.label.cex = 2,          # 设置顶点标签的大小
         vertex.label.color = 'black',  # 设置顶点标签的颜色
         vertex.color = "lightblue1",   # 设置顶点颜色
         vertex.size = 50,              # 设置顶点大小
         vertex.frame.color = "black",  # 设置顶点边框颜色
         vertex.shape = "circle",       # 设置顶点形状
         edge.width = 2,                # 设置边的宽度
         edge.color = 'black')          # 设置边的颜色
  }

  return(list(adjacency_matrix = adj_mat,
              "the pairs without significant batch effects" = pairs[simul_count],
              mean_test = mean_test,
              cov_test = cov_test,
              simul_test = simul_test))
}
