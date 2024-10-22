# MetBEC
Batch effects are inevitable in large-scale metabolomics. Prior to data analysis, batch effect correction (BEC) is imperative to prevent from obscuring biological variations. Here, we provide six quality control (QC)-based BEC algorithms with several metrics for correction assessment.
- Intra-BEC methods:
  - support vector regression (SVR)
  - random forest (RF) regression
  - technical variation elimination with ensemble learning architecture (TIGER)
  - **extreme gradient boost (XGBoost) regression**
- Inter-BEC methods:
  - batch-ratio
  - **covariance correction (CoCo)**
- Metrics for correction assessment:
  - principal component analysis (PCA)
  - relative standard deviation (RSD)
  - dispersion-ratio (D-ratio)
  - **QC-based simultaneous tests (QC-ST)**
- Data pre-processing
  - Outlier detection by PCA with Hotelling's $$T^2$$ statistic and Squared Prediction Errors (SPEs) 
## Install
The R package depends on [GLassoElnetFast](https://github.com/TobiasRuckstuhl/GLassoElnetFast). If necessary, execute
```R
remotes::install_github("TobiasRuckstuhl/GLassoElnetFast")
```
or 
```R
remotes::install_github("Bubble-o0O/GLassoElnetFast")
```
([GLassoElnetFast](https://github.com/Bubble-o0O/GLassoElnetFast)) firstly. Then, execute
```R
remotes::install_github("Bubble-o0O/MetBEC")
```
to install ***MetBEC***.
## Contact
E-mail to <guozhendong19@mails.ucas.ac.cn> for any questions (**preferably in Chinese**).
## References
- A new statistical method for evaluating batch effects based on quality control samples with the matching batch effect correction strategy in metabolomics. *Submitting*.
## Supplement
安装***GLassoElnetFast***包时，可能出现报错。如报错，请参考以下解决方案：
```R
if (require("pkgbuild", quietly = TRUE) == FALSE){
  install.packages("pkgbuild")
}
```
- **任务栏 - 搜索 - 编辑系统环境变量**
- **高级 - 环境变量**
- **系统变量 - Path - 编辑**
- **浏览 - C:\RBuildTools\4.3\mingw64\bin**（备注：**RBuildTools**文件夹的地址请自行更改）
- **确定 - 确定 - 确定**

重启R Studio，再执行**Install**的指令。
