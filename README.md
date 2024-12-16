# MetBEC
***Quality Control (QC)-based Batch Effect Correction (BEC) in Metabolomics***
# Contact
E-mail to <guozhendong19@mails.ucas.ac.cn> for any questions (**preferably in Chinese**).
## Description
Batch effects are inevitable in large-scale metabolomics. Prior to formal data analysis, BEC is applied to prevent from obscuring biological variations. Here, we provide six QC-based BEC algorithms and several relevant functions.
- Intra-BEC methods:
  - Support Vector Regression (SVR)
  - Random Forest (RF) regression
  - Technical variation elImination with ensemble learninG architEcturR (TIGER)
  - **eXtreme Gradient Boost (XGBoost) regression**
- Inter-BEC methods:
  - batch-ratio
  - **Covariance Correction (CoCo)**
- Metrics for correction assessment:
  - Relative Standard Deviation (RSD)
  - Dispersion-ratio (D-ratio)
  - **QC-based simultaneous tests (QC-ST)**
  - scatter plot for univariate visualization
  - Principal Component Analysis (PCA) score plot for multivariate visualization
- Data pre-processing:
  - outlier detection by PCA with Hotelling's $$T^2$$ statistic and Squared Prediction Errors (SPEs) 
## Install
- The R package depends on [GLassoElnetFast](https://github.com/TobiasRuckstuhl/GLassoElnetFast) or [GLassoElnetFast](https://github.com/Bubble-o0O/GLassoElnetFast) (my version). If necessary, execute
```R
if (require("remotes", quietly = TRUE) == FALSE){
  install.packages("remotes")
}

remotes::install_github("TobiasRuckstuhl/GLassoElnetFast")
## Alternative: My version is also allowed, which has corrected some minor errors.
# remotes::install_github("Bubble-o0O/GLassoElnetFast")
```
where `1: All` is recommended to select. 
- Execute
```R
remotes::install_github("Bubble-o0O/MetBEC")
```
to install ***MetBEC***, where `1: All` is recommended to select.
## Details
- Please read **Help** of the required function in R Studio.
  - Search for ***MetBEC*** in **Packages** and enter.
  - Select the required function in **Help Pages**.
- Please read my preprint: 
  - Guo, Z. High-dimensional Statistics Applications to Batch Effects in Metabolomics. *arXiv* **2024**. DOI: [arXiv:2412.10196](https://arxiv.org/pdf/2412.10196).
## Supplements (in Chinese)
执行**Install**时如报错，请参考以下解决方案：
- 如有必要，安装RBuildTools ([pkgbuild](https://cran.r-project.org/web/packages/pkgbuild/index.html))
```R
if (require("pkgbuild", quietly = TRUE) == FALSE){
  install.packages("pkgbuild")
}
```
- **任务栏 - 搜索 - 编辑系统环境变量**
- **高级 - 环境变量**
- **系统变量 - Path - 编辑**
- **浏览 - C:\RBuildTools\4.3\mingw64\bin**
  - RBuildTools的地址请自行更改
  - R与RBuildTools的版本应相同
  - 地址中不建议出现中文字符
- **确定 - 确定 - 确定**
- 重启R Studio，再次执行**Install**
