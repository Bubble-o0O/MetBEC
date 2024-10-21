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
## Install
The R package depends on [GLassoElnetFast](https://github.com/TobiasRuckstuhl/GLassoElnetFast). If necessary, execute
```R
remotes::install_github("TobiasRuckstuhl/GLassoElnetFast")
```
or
```R
remotes::install_github("Bubble-o0O/GLassoElnetFast")
```
firstly. Then, execute
```R
remotes::install_github("Bubble-o0O/MetBEC")
```
to install ***MetBEC***.
## Contact
<guozhendong19@mails.ucas.ac.cn> for any questions (**preferably in Chinese**).
## References
A new statistical method for evaluating batch effects based on quality control samples with the matching batch effect correction strategy in metabolomics (*Submitting*).
