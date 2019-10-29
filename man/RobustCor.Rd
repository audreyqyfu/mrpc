\name{RobustCor}
\alias{RobustCor}
\title{
Robust correlation matrix
}
\description{
Robust correlation matrix based on beta value. The value of beta plays a key role in the performance of the robust method, which controls the tradeoff between the robustness and efficiency of the estimators.
}
\usage{
RobustCor(xx, Beta, plot = FALSE)
}

\arguments{
  \item{xx}{
Data matrix
}
\item{Beta}{
Tuning parameter, between 0 and 1, if 0 then equal to nonrobust, classical method. We suggest using, Beta = 0.005 in both without and with outliers in simulation study. This value should reflect the amount of outliers in the data. Whereas a large value increases robustness, it reduces sensitivity of identifying an edge. We need a more principled way to determine this value.
}
\item{plot}{
To set no plotting as the default for weight vs gene index.
}
}
\details{
We take a robust approach and calculate the robust correlation matrix (Badsha et al., 2013) on which the series of hypothesis testing is performed. The performance of the robust correlation method depends on the values of the tuning parameter beta. It controls the tradeoff between robustness and efficiency of estimators. This method shows high performance for a wide range of beta. The values of beta lies between 0 and 1, such that a large value of beta decreases the efficiency, while it increases the robustness of an estimator, and vice-versa for a small value of beta. Thus, we need to select an optimal beta to obtain both high robustness and efficiency, while it depends on the initialization of model parameters, data contamination rates, types of data contamination, types of datasets, and so on. We used the beta value from Badsha et al., 2013. The robust method reduces to the classical method (Biased estimator) with the tuning parameter beta -->0. When the data matrix contains missing values, we perform imputation using the R package mice (Buuren and Groothuis-Oudshoorn, 2011).

}
\value{
\link{list} of objects as follows:

\itemize{
\item \code{RR}: Robust correlation matrix.
\item \code{M}: Robust mean vector.
\item \code{V}: Robust covariance matrix.
\item \code{Wt}: Weight for each observation.
}
}

\references{
1. Badsha, M. B. Mollah, M. N. Jahan, N. and Kurata, H. (2013). Robust complementary hierarchical clustering for gene expression data analysis by beta-divergence. J Biosci Bioeng 116(3): 397-407.

2. Van Buuren, S., Groothuis-Oudshoorn, K. (2011). mice: Multivariate Imputation by Chained Equations in R. Journal of Statistical Software, 45(3), 1-67. http://www.jstatsoft.org/v45/i03/
}
\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}

\examples{
RobustCor_objects <- RobustCor(simu_data_M0,
               Beta = 0.005,
               plot = FALSE)
               
Rcorr <- RobustCor_objects $RR  # Correlation matrix
}
