\name{log_mvn_lkh}
\alias{log_mvn_lkh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood for measurements in active intervals.
}
\description{
Calculate the log likelihood of multivariate normal distribution (MVN).
}
\usage{
log_mvn_lkh(y_input, t_input, mu_input, sigma2_input, l_input, g_input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y_innput}{Scalar or vector of y of interest. }
  \item{t_innput}{Length of y_input.}
  \item{mu_innput}{Mean parameter of MVN.}
  \item{sigma2_innput}{Scale parameter of MVN.}
  \item{l_innput}{Range parameter of MVN.}
  \item{g_innput}{Nugget term of MVN.}

}

\value{
\item{value}{Calculated log likelihood}
}
