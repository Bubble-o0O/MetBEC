# MetBEC
Batch effects are inevitable in large-scale metabolomics. Prior to data analysis, **batch effect correction (BEC)** is imperative to prevent from obscuring biological variations. Here, we provide six quality control (QC)-based BEC algorithms with several metrics for correction assessment.
- Intra-BEC methods:
  - Support Vector Regression (SVR)
  - Random Forest (RF) regression
  - Technical variation elImination with ensemble learninG architEcturR (TIGER)
  - **eXtreme Gradient Boost (XGBoost) regression**
- Inter-BEC methods:
  - batch-ratio
  - **Covariance Correction (CoCo)**
- Metrics for correction assessment:
  - Principal Component Analysis (PCA)
  - Relative Standard Deviation (RSD)
  - Dispersion-ratio (D-ratio)
  - **QC-based simultaneous tests (QC-ST)**
- Data pre-processing
  - outlier detection by PCA with Hotelling's $$T^2$$ statistic and Squared Prediction Errors (SPEs) 
## Install
- The R package depends on [GLassoElnetFast](https://github.com/TobiasRuckstuhl/GLassoElnetFast). If necessary, execute
```R
remotes::install_github("TobiasRuckstuhl/GLassoElnetFast")

# Alternative: My version is also allowed, which has corrected some minor errors.
# remotes::install_github("Bubble-o0O/GLassoElnetFast")
```
([GLassoElnetFast](https://github.com/Bubble-o0O/GLassoElnetFast)), and `1: All` is recommended to select. 
- Execute
```R
remotes::install_github("Bubble-o0O/MetBEC")
```
to install ***MetBEC***, and `1: All` is recommended to select.
## Details
Please read **Help** of the required function in R studio.
- Search for ***MetBEC*** in **Packages** and enter.
- Select the required function in **Help Pages**.
## Contact
E-mail to <guozhendong19@mails.ucas.ac.cn> for any questions (**preferably in Chinese**).
## References
- A new statistical method for evaluating batch effects based on quality control samples with the matching batch effect correction strategy in metabolomics. *Submitting*.
## Supplements
执行**Install**时如报错，请参考以下解决方案：
- 安装RBuildTools（若已安装，即跳过）
```R
if (require("pkgbuild", quietly = TRUE) == FALSE){
  install.packages("pkgbuild")
}
```
- **任务栏 - 搜索 - 编辑系统环境变量**
- **高级 - 环境变量**
- **系统变量 - Path - 编辑**
- **浏览 - C:\RBuildTools\4.3\mingw64\bin**
  - **RBuildTools**文件夹的地址请自行更改
  - 地址中不建议出现中文字符
  - R与RBuildTools的版本应相同
- **确定 - 确定 - 确定**
- 重启R Studio，再执行**Install**的相关代码
