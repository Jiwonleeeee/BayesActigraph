\name{acc_death_ar}
\alias{acc_death_ar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Metropolis-Hastings for the death of an active interval (a,r).
}
\description{
Return the updated parameters (if the death of (a,r) is accepted), return the same parameters otherwise. Its structure is the same as \code{acc_birth.ar} except for the part calculating the acceptance probability.
}
\usage{
acc_death_ar(y_input, t_input, mvn_inputs, Pi_input, Lambda_input, p_birth_input, p_death_input, active_input, rest_input, n_active_input, n_rest_input)
}
%- maybe also 'usage' for other objects documented here.

\details{
Used arguments are the same as \code{acc_birth_ar}
}
\value{
\item{active}{(B-1)-dimensional vector if the death of (a,r) is accepted, the same B-dimensional vector otherwise.}
\item{rest}{(B-1)-dimensional vector if the death of (a,r) is accepted, the same B-dimensional vector otherwise.}
\item{accept}{Indicator if accept = 1 then the proposed is accepted, and accept = 0 otherwise. (for tracking the acceptance probability).}
\item{n_active}{The number of active intervals after updating}
\item{n_rest}{The number of resting intervals after updating}
}
