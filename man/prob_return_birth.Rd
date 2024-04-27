\name{prob_return_birth}
\alias{prob_return_birth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probabilities for choosing an active (resting) interval
}
\description{
Calculate the ]probabilities for choosing an active (resting) interval in which a new resting (active) interval will be proposed.
}
\usage{
prob_return_birth(active_input, rest_input, n_active_input, n_rest_input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{active_input}{Same as the one in \code{sum_log_mvn_lkh}}
  \item{rest_input}{Same as the one in \code{sum_log_mvn_lkh}}
  \item{n_active_input}{Total number of active intervals. It is the same as the length of 'active_input' but the separate parameter 'n_active_input' is defined to accommodate more general cases
}
  \item{n_rest_input}{total number of resting intervals. It is the same as the length of 'rest_input' - 1.}
}

\value{
  \item{active_prob}{'n_active_input'-dimensional vector with each element representing probability proportional to the length of the corresponding active interval.}
  \item{rest_prob}{'n_rest_input'-dimensional vector with each element representing probability proportional to the length of the corresponding resting interval.}
}
