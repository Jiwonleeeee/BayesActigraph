\name{acc_birth_ar}
\alias{acc_birth_ar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Metropolis-Hastings for the birth of an active interval (a,r)
}
\description{
Return the updated parameters (if the birth of (a,r) is accepted), return the same parameters otherwise.
}
\usage{
acc_birth_ar(y_input, t_input, mvn_inputs, Pi_input, Lambda_input, p_birth_input, p_death_input, active_input, rest_input, n_active_input, n_rest_input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y_input}{
T-dimensional y vector.
}
  \item{t_input}{
T-dimensional time index vector \eqn{(1,\ldots,T)}.
}
  \item{mvn_inputs}{
4-dimensional vector of MVN parameters \eqn{(\mu, \sigma^2, l, g)}.

}
  \item{Pi_input}{
K-dimensional vector of mixing proportions \eqn{\pi_1,\ldots,\pi_K}.
}
  \item{Lambda_input}{
K-dimensional vector of parameters of Poisson distributions.
}
  \item{p_birth_input}{
Probability of choosing birth move.
}
  \item{p_death_input}{
Probability of choosing death move.
}
  \item{active_input}{
B-dimensional vector.
}
  \item{rest_input}{
B-dimensional vector.
}
  \item{n_active_input}{
Total number of active intervals.
}
  \item{n_rest_input}{
Total number of resting intervals.
}
}

\value{
\item{active}{(B+1)-dimensional vector if the birth of (a,r) is accepted, the same B-dimensional vector otherwise.}
\item{rest}{(B+1)-dimensional vector if the birth of (a,r) is accepted, the same B-dimensional vector otherwise.}
\item{accept}{Indicator if accept = 1 then the proposed is accepted, and accept = 0 otherwise. (for tracking the acceptance probability).}
\item{n_active}{The number of active intervals after updating.}
\item{n_rest}{The number of resting intervals after updating.}
}
