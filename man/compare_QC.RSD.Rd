% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QC.RSD.R
\name{compare_QC.RSD}
\alias{compare_QC.RSD}
\title{Relative Standard Deviation (RSD) distribution of different data for comparison}
\usage{
compare_QC.RSD(data1, ..., BEC_method, line_alpha = 0.8)
}
\arguments{
\item{data1}{A dataframe.\cr
Execute \code{data(Dataset_I)} and \code{View(Dataset_I)} for formats.}

\item{...}{Alternatives. More data (e.g., \code{data2}, \code{data3}) can be input.}

\item{BEC_method}{A character vector. The length must be the same as \code{c(data1, ...)}. For example, \code{c("raw", "RF", "XGBoost")}.}

\item{line_alpha}{A numeric scalar. Default is 0.8. Determines the curve line opacity in the cumulative distribution plot.}
}
\value{
A cumulative distribution plot.
}
\description{
The metric only examines QC samples, whose formula is:\cr
\deqn{\mathrm{RSD}_{i}=\frac{\sqrt{\mathrm{Var}(y^{\mathrm{(QC)}}_{i})}}{\mathrm{mean}(y^{\mathrm{(QC)}}_{i})},i=1,2,...,p,}\cr
where \eqn{y^{\mathrm{(QC)}}_{i}} denotes the intensity of QC samples' \eqn{i}th metabolite;\cr
\eqn{p} denotes the variable (metabolite) number.
}
\note{
\itemize{
 \item{The acceptance criterion is 15\%, 20\%, and 30\%.}
}
}
\examples{
data(Dataset_I)
data <- Dataset_I
data.RF <- RF.correction(data)
compare_QC.RSD(data, data.RF, BEC_method = c("raw", "RF"))
}
\references{
\itemize{
 \item{Broadhurst, D.; Goodacre, R.; Reinke, S. N.; Kuligowski, J.; Wilson, I. D.; Lewis, M. R.; Dunn, W. B. Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies. \emph{Metabolomics} \strong{2018}, 14 (6), 72. DOI: 10.1007/s11306-018-1367-3.}
 \item{Dudzik, D.; Barbas-Bernardos, C.; Garcia, A.; Barbas, C. Quality assurance procedures for mass spectrometry untargeted metabolomics. a review. \emph{J Pharm Biomed Anal} \strong{2018}, 147, 149-173. DOI: 10.1016/j.jpba.2017.07.044.}
 \item{Han, W.; Li, L. Evaluating and minimizing batch effects in metabolomics. \emph{Mass Spectrometry Reviews} \strong{2020}, 41 (3), 421-442. DOI: 10.1002/mas.21672.}
 \item{Martens, A.; Holle, J.; Mollenhauer, B.; Wegner, A.; Kirwan, J.; Hiller, K. Instrumental Drift in Untargeted Metabolomics: Optimizing Data Quality with Intrastudy QC Samples. \emph{Metabolites} \strong{2023}, 13 (5). DOI: 10.3390/metabo13050665.}
}
}
\seealso{
\code{\link{QC.RSD}}.
}
\author{
Zhendong Guo (\email{guozhendong19@mails.ucas.ac.cn}).
}
