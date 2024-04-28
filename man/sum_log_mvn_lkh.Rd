\name{sum_log_mvn_lkh}
\alias{sum_log_mvn_lkh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sum of the log likelihood of MVN
}
\description{
Return the sum of the log likelihood of MVN from every active interval.
}
\usage{
sum_log_mvn_lkh(y_input, t_input, mu_input, sigma2_input, l_input, g_input, active_input, rest_input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y_input}{T-dimensional y vector.}
\item{active_input}{B-dimensional vector.}
\item{rest_input}{B-dimensional vector. 'active_input' and 'rest_input' construct active intervals to filter the elements of 'y_input' belonging to the active states.}
}

\value{
\item{sum_value}{Calculated sum of log likelihood.}
}
