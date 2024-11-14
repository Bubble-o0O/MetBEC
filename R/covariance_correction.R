#' @title QC-based BEC algorithm: Covariance Correction (CoCo)
#'
#' @description CoCo is used for inter-BEC, which employs the Graphical Elastic Net (GELNET) algorithm of the Gaussian Graphical Model (GGM).\cr
#' For the multivariate normal distribution \eqn{\mathrm{N}_{p}(\boldsymbol{\mu},\boldsymbol{\Sigma})}, the penalized Maximum Likelihood Estimator (MLE) of the precision matrix (inverse covariance matrix, namely \eqn{\boldsymbol{\Theta}=\boldsymbol{\Sigma}^{-1}}) by GELNET is:\cr
#' \deqn{\hat{\boldsymbol{\Theta}} = \arg\min_{\boldsymbol{\Theta}}\{-\mathrm{log}|\boldsymbol{\Theta}| + \mathrm{tr}(\boldsymbol{S}\boldsymbol{\Theta}) + \lambda(\alpha\|\boldsymbol{\Theta}-\boldsymbol{T}\|_{1} + \frac{1-\alpha}{2}\|\boldsymbol{\Theta}-\boldsymbol{T}\|^2_{2})\},}
#' where \eqn{\boldsymbol{S}} denotes the MLE of \eqn{\boldsymbol{\Sigma}};\cr
#' \eqn{\|.\|_{1}} denotes the matrix \eqn{L_{1}}-norm;\cr
#' \eqn{\|.\|_{2}} denotes the matrix \eqn{L_{2}}-norm (Frobenius norm);\cr
#' the hyperparameter \eqn{\alpha \in [0,1]};\cr
#' the hyperparameter \eqn{\lambda \in (0,+\infty)};\cr
#' \eqn{\boldsymbol{T}} denotes the target matrix, which is a positive semi-definite matrix. Here, we select \eqn{\boldsymbol{T}=\boldsymbol{I}_{p}} (see \strong{Note} for details).
#'
#' @param data A dataframe.\cr
#' Execute \code{data(Dataset_I)} and \code{View(Dataset_I)} for formats.
#' @param alpha_range A two-length numeric vector. Default is \code{c(0, 1)}. If the two values are the same, \eqn{\alpha} will be a fixed value.
#' @param lambda_range A two-length numeric vector. Default is \code{c(0, 10)}. If the two values are the same, \eqn{\lambda} will be a fixed value.
#' @param search_iterations A numeric scalar. Default is 500. Applied to random search for hyperparameter optimization.
#' @param random_seed A numeric scalar. Default is 123. Applied to random search for hyperparameter optimization.
#' @param continue Logical. Default is \code{FALSE}. Determines whether to implement CoCo when \code{QC_ST(data)} has indicated no significant batch effects across different batches.
#' @param sig_level.alpha A numeric scalar. Default is 0.05, which is unnecessary to be tuned for most of the time.
#' @param simul_method Default is \code{NULL}, which selects \code{fisher} when \eqn{\frac{n_{1}+n_{2}}{2} \ge 10} or else selects \code{HN}.\cr
#' \code{cauchy}, namely "Yu-Cauchy"; \code{fisher}, namely "Yu-Fisher".
#' @param p_value.adjust Default is \code{fdr}, which is equivalent to \code{BH} and is unnecessary to be tuned for most of the time.
#' @param cl Default is \code{NULL}, which uses all the CPU cores for parallel computing. Otherwise, it should be a numeric scalar.
#'
#' @return
#' \item{hyperparameters}{A dataframe. The optimal hyperparameters.}
#' \item{VarFC_matrix}{A dataframe. The subject samples' variance fold changes before and after CoCo. Denoted as \eqn{\boldsymbol{V}}.}
#' \item{VarFC.all_mean}{A numeric scalar. The mean value of \code{VarFC_matrix}, namely \eqn{\mathrm{mean}(\boldsymbol{V})}, which determines whether the subject samples are overcorrected. We provide an acceptance range as reference: \eqn{(0.25,4)}.}
#' \item{corrected_data}{A dataframe. The format is the same as \code{data}.}
#'
#' @export
#' @import stats parallel GLassoElnetFast
#' @importFrom expm sqrtm
#'
#' @details See our paper for details.
#'
#' @note
#' \itemize{
#'  \item{GELNET is a sophisticated algorithm, which combines the Graphical Lasso (equivalent to \eqn{\alpha=1}) and the Graphical Ridge (equivalent to \eqn{\alpha=0}). When \eqn{\lambda \to 0}, \eqn{\hat{\boldsymbol{\Theta}} \to \boldsymbol{S}^{-1}}; When \eqn{\lambda \to +\infty}, \eqn{\hat{\boldsymbol{\Theta}} \to \boldsymbol{T}}.}
#'  \item{Several references have indicated that \eqn{\boldsymbol{T}=\boldsymbol{I}_{p}} performs well for most of the time. Additionally, \eqn{\boldsymbol{T}=\boldsymbol{I}_{p}} has an another advantage in CoCo (see the next note) comparing to other alternatives.}
#'  \item{If \code{QC_ST(data)} has indicated no significant batch effects across different batches, CoCo is no longer necessary. This is because redundant CoCo might be time-consuming, considering that the algorithm complexity is \eqn{O(Bp^3)}. Further, hyperparameter optimization of CoCo will obtain larger \eqn{\lambda}, which means that \code{data} and \code{CoCo(data)$corrected_data} will be almost the same, namely \eqn{\mathrm{mean}(\boldsymbol{V}) \approx 1}.}
#' }
#' @author Zhendong Guo (\email{guozhendong19@mails.ucas.ac.cn}).
#' @references
#' \itemize{
#'  \item{Kovács, S.; Ruckstuhl, T.; Obrist, H.; Bühlmann, P. Graphical Elastic Net and Target Matrices: Fast Algorithms and Software for Sparse Precision Matrix Estimation. \emph{arXiv} \strong{2021}. DOI: arXiv:2101.02148.}
#'  \item{Kheyri, A.; Bekker, A.; Arashi, M. High-Dimensional Precision Matrix Estimation through GSOS with Application in the Foreign Exchange Market. \emph{Mathematics} \strong{2022}, 10 (22). DOI: 10.3390/math10224232.}
#'  \item{Bekker, A.; Kheyri, A.; Arashi, M. A Computational Note on the Graphical Ridge in High-dimension. \emph{arXiv} \strong{2023}. DOI: arXiv:2312.15781.}
#' }
#' @seealso \code{\link{QC_ST}}, \code{GLassoElnetFast::\link[GLassoElnetFast]{gelnet}}.
#' @examples
#' data(Dataset_I)
#' data <- Dataset_I
#' data.RF <- RF.correction(data)
#' data.RF.QC_ST <- QC_ST(data.RF)
#'
#' RF_CoCo <- CoCo(data.RF)
#' data.RF_CoCo <- RF_CoCo$corrected_data
#' data.RF_CoCo.QC_ST <- QC_ST(data.RF_CoCo)

CoCo <- function(data,
                 alpha_range = c(0, 1), # Kovacs (2021)与Bernardini (2022)选取alpha = 0.5
                 lambda_range = c(0, 10), # 参考Bekker (2023)
                 # 其它优化范围：[0.01,10]参考Kheyri (2022)；[0.9^40,1]参考Kovacs (2021)；[exp(-4),2]参考Bernardini (2022)
                 search_iterations = 500, random_seed = 123,
                 continue = FALSE,
                 sig_level.alpha = 0.05, simul_method = c(NULL, 'cauchy', 'fisher', 'HN'),
                 p_value.adjust = c("bonferroni", "holm", "hochberg", "hommel",
                                    "fdr", "BH", "BY", "none"),
                 cl = NULL){
  batch.level <- unique(data$batch)
  B <- length(batch.level) # 批次数
  if (B == 1){
    stop("There is only one batch.")
  }

  if (any(alpha_range > 1) == TRUE || any(alpha_range < 0) == TRUE){
    stop("'alpha' must be in [0,1].") # 若alpha = 0，即为Graphical Lasso；若alpha = 1，即为Graphical Ridge
  }
  alpha_range <- sort(alpha_range) # 上下限可以相等

  if (any(lambda_range < 0) == TRUE || all(lambda_range == 0) == TRUE){
    stop("'lambda' must be in (0,+\u221E).")
  }
  lambda_range <- sort(lambda_range) # 上下限可以相等，但不能同时为0

  if (length(simul_method) > 1){
    simul_method <- NULL
  }

  if (length(p_value.adjust) > 1){ # 默认选择BH法（FDR）校正
    p_value.adjust <- 'BH' # 与'fdr'一致
  }

  if (is.null(cl) == TRUE){
    cl <- parallel::makeCluster(parallel::detectCores())
  }else {
    cl <- parallel::makeCluster(cl)
  }

  data.QC_ST <- QC_ST(data, print_plot = FALSE,
                      sig_level.alpha = sig_level.alpha,
                      simul_method = simul_method,
                      p_value.adjust = p_value.adjust)

  if (continue == FALSE){
    if (all(data.QC_ST$cov_test$q_value >= sig_level.alpha) == TRUE ||
        all(data.QC_ST$simul_test$q_value >= sig_level.alpha) == TRUE){
      stop("It is no necessary to implement CoCo.
           If you want to continue, please let 'continue = TRUE'.")
    }
  }

  QC <- subset(data, data$type == 'qc')
  sample <- subset(data, data$type == 'sample')

  QC.sd <- apply(QC[, -1:-4], 2, sd)

  n <- nrow(QC) # QC的总样本量
  p <- ncol(QC[, -1:-4]) # 变量数

  QC.Data <- lapply(1:B, function(b){
    QC_b.Data <- subset(QC, batch == batch.level[b])[, -1:-4]

    QC_b.var <- apply(QC_b.Data, 2, var)
    if (any(QC_b.var < 1e-6) == TRUE){
      var_zero <- colnames(QC_b.Data)[which(QC_b.var < 1e-6)]
      stop(sprintf("For QC samples in Batch %s, %d variables' variances are 0 (or close to 0), including: %s.",
                   batch.level[b], length(var_zero), paste(var_zero, collapse = ", ")))
    }else {
      return(QC_b.Data)
    }
  })

  QC.info <- QC[, 1:4]
  QC.scale <- lapply(1:B, function(b){
    subset(cbind(QC.info, scale(QC[, -1:-4])), batch == batch.level[b])[, -1:-4]
  })

  QC.n <- sapply(1:B, function(b) nrow(QC.Data[[b]]))

  QC.centriod <- lapply(QC.Data, colMeans)



  sample.Data <- lapply(1:B, function(b)
    subset(sample, batch == batch.level[b])[, -1:-4]
  )

  sample.centriod <- lapply(sample.Data, colMeans)



  # 超参数的随机搜索
  set.seed(random_seed)
  alpha_mat <- matrix(runif(search_iterations * B, alpha_range[1], alpha_range[2]),
                      nrow = search_iterations, ncol = B)
  lambda_mat <- matrix(runif(search_iterations * B, lambda_range[1], lambda_range[2]),
                       nrow = search_iterations, ncol = B)

  search.result <- parallel::parLapply(cl, 1:search_iterations, function(k){
    alpha_vec <- alpha_mat[k,]
    lambda_vec <- lambda_mat[k,]

    cov.MLE <- function(X){(nrow(X) - 1) / nrow(X) * cov(X)} # 协方差矩阵的极大似然估计

    gelnet.result <- lapply(1:B, function(b){
      alpha <- alpha_vec[b]
      lambda <- lambda_vec[b]

      QC_b.scale <- QC.scale[[b]]

      gelnet <- GLassoElnetFast::gelnet(S = cov.MLE(QC_b.scale),
                                        lambda = lambda, alpha = alpha,
                                        Target = GLassoElnetFast::target(Y = QC_b.scale, type = 'Identity'))

      QC.Sigma_b.scale <- gelnet$W # 标准化数据估计的协方差矩阵


      QC.Sigma_b <- diag(QC.sd) %*% QC.Sigma_b.scale %*% diag(QC.sd)

      return(list(QC.Sigma = QC.Sigma_b,
                  QC.Theta = solve(QC.Sigma_b)))
    })

    QC.Sigma_list <- lapply(gelnet.result, function(x) x[[1]])
    QC.Theta_list <- lapply(gelnet.result, function(x) x[[2]])

    Sigma_tilde <- Reduce("+", lapply(1:B, function(b)
      QC.n[b] * QC.Sigma_list[[b]])) / n # 加权平均

    sample.var <- lapply(1:B, function(b){
      sample_b.Data <- sample.Data[[b]]
      return(sapply(1:p, function(j) var(sample_b.Data[, j])))
    })

    trans_mat_list <- lapply(QC.Theta_list, function(x)
      expm::sqrtm(x) %*% expm::sqrtm(Sigma_tilde))

    # 校正
    correction.result <- lapply(1:B, function(b){
      QC_b.Data <- QC.Data[[b]]
      QC_b.centriod <- QC.centriod[[b]]

      sample_b.Data <- sample.Data[[b]]
      sample_b.centriod <- sample.centriod[[b]]

      QC_b.1 <- as.matrix(QC_b.Data) %*% trans_mat_list[[b]] # 协方差矩阵的校正
      sample_b.1 <- as.matrix(sample_b.Data) %*% trans_mat_list[[b]] # 协方差矩阵的校正

      QC_b.2 <- t(t(scale(QC_b.1, scale = FALSE)) + QC_b.centriod) # 还原初始的均值向量；广播原则
      # 等价形式：QC_b.2 <- t(t(QC_b.1) - colMeans(QC_b.1) + QC_b.centriod)
      sample_b.2 <- t(t(scale(sample_b.1, scale = FALSE)) + sample_b.centriod) # 还原初始的均值向量；广播原则



      # 替换非正值
      all_b <- rbind(QC_b.2, sample_b.2)
      if (any(all_b <= 0) == TRUE){
        # warning(sprintf("%.d samples' intensity in Batch %s is non-positive and has been replaced.",
        #                 sum(all_b <= 0), batch.level[b]))
        QC_b.2[QC_b.2 <= 0] <- QC_b.Data[QC_b.2 <= 0]
        sample_b.2[sample_b.2 <= 0] <- sample_b.Data[sample_b.2 <= 0]
      }

      QC[QC$batch == batch.level[b], -1:-4] <- QC_b.2
      sample[sample$batch == batch.level[b], -1:-4] <- sample_b.2

      sample_b.var.correction <- sapply(1:p, function(j) var(sample_b.2[, j]))

      return(list(rbind(QC[QC$batch == batch.level[b],], sample[sample$batch == batch.level[b],]),
                  sample_b.var.correction))
    })

    all <- do.call(rbind, lapply(correction.result, function(x) x[[1]]))
    all <- all[order(all$order),] # 根据order从小到大排序

    sample.var.correction <- lapply(correction.result, function(x) x[[2]])



    sample.var_fold_change <- lapply(1:B, function(b) sample.var.correction[[b]]/sample.var[[b]]) # 受试样本经协方差校正前后各批次中各变量的方差变化倍数

    VarFC_matrix <- do.call(rbind, sample.var_fold_change)

    VarFC.all_mean <- mean(VarFC_matrix) # 矩阵格式

    hyperparameters.result <- data.frame(alpha = alpha_vec,
                                         lambda = lambda_vec)
    rownames(hyperparameters.result) <- batch.level



    QC_ST <- QC_ST(all, print_plot = FALSE,
                   sig_level.alpha = sig_level.alpha,
                   simul_method = simul_method,
                   p_value.adjust = p_value.adjust)

    if (all(QC_ST$cov_test$q_value >= sig_level.alpha) == TRUE ||
        all(QC_ST$simul_test$q_value >= sig_level.alpha) == TRUE){ # 协方差矩阵的齐性检验比同时检验更严格
      return(list(hyperparameters.result,
                  VarFC_matrix,
                  VarFC.all_mean,
                  all))
    }else {
      return(NULL)
    }
  })

  # 关闭集群
  parallel::stopCluster(cl)

  # 剔除NULL
  search.result <- search.result[sapply(search.result, Negate(is.null))]
  count <- length(search.result)
  if (count == 0){
    if (search_iterations > 5000){
      stop("The covariance matrices between some two batches still have statistical significance.
           When 'search_iterations' > 5000, it suggests that CoCo might fail for the data.
           We provide several alternatives:
           (1) another prepositive BEC algorithm should be used;
           (2) some hyperparameters of the prepositive BEC algorithm should be tuned;
           (3) outliers from the QC samples after prepositive BEC should be excluded.")
    }else {
      stop("The covariance matrices between some two batches still have statistical significance.
         Some arguments must be tuned to ensure effective correction.
         For example, 'search_iterations' should be larger or 'random_seed' should be tuned.")
    }
  }

  hyperparameters <- lapply(search.result, function(x) x[[1]])
  VarFC_matrix <- lapply(search.result, function(x) x[[2]])
  VarFC.all_mean <- sapply(search.result, function(x) x[[3]])
  all.correction <- lapply(search.result, function(x) x[[4]])

  if (min(VarFC.all_mean) >= 4 || min(VarFC.all_mean) <= 0.25){ # 防止受试样本在校正前后的数据变化较大
    if (search_iterations > 5000){
      warning("The mean value of the variance fold change matrix \u2265 4 or \u2264 0.25.
              If 'search_iterations' > 5000 and other arguments are default,
              we recommend that 'p_value.adjust' should be tuned (such as 'bonferroni') or
              'sig_level.alpha' should be tuned (such as 0.01) to prevent overcorrection for the subject samples.
              After that, if the situation does not improve yet, it suggests that CoCo might
              fail for the data, and the prepositive BEC algorithm should be tuned.")
    }else {
      warning("The mean value of variance fold change matrix \u2265 4 or \u2264 0.25, which suggests that
                some arguments should be tuned to prevent overcorrection for the subject samples.
                For example, 'search_iterations' should be larger or 'random_seed' should be tuned.")
    }
  }

  best <- which.min(VarFC.all_mean)

  VarFC_matrix.best <- as.data.frame(VarFC_matrix[[best]])
  rownames(VarFC_matrix.best) <- batch.level
  colnames(VarFC_matrix.best) <- colnames(data[, -1:-4])

  all.best <- all.correction[[best]]

  return(list(hyperparameters = hyperparameters[[best]],
              VarFC_matrix = VarFC_matrix.best,
              VarFC.all_mean = min(VarFC.all_mean),
              corrected_data = all.best))
}
