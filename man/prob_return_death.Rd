\name{prob_return_death}
\alias{prob_return_death}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probabilities for choosing an active (resting) interval
}
\description{
Calculate the probabilities for choosing an active (resting) interval to be removed.}
\usage{
prob_return_death(active_input, rest_input, n_active_input, n_rest_input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{active_input}{All the arguments are the same as \code{prob_return_birth}.}
  \item{rest_input}{}
  \item{n_active_input}{}
  \item{n_rest_input}{}
}

\value{
\item{active_prob}{('n_active_input'-2)-dimensional vector with each element representing a probability inversely proportional to the length of the corresponding active interval. Since we fix the state of the first and last interval as active, they cannot be removed.}
\item{rest_prob}{'n_rest_input'-dimensional vector with each element representing a probability inversely proportional to the length of the corresponding resting interval.}

}
