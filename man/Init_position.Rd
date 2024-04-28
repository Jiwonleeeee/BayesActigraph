\name{initial_position}
\alias{initial_position}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Setting initial positions}
\description{
 Set initial positions of the change points \eqn{a_1,\ldots,a_B} and \eqn{r_1,\ldots,r_B} used to construct active and resting intervals.
}
\usage{
initial_position(
y_input, cutoff, end_input
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y_input}{T-dimensional y vector.}
  \item{cutoff}{User-specific cutoff (threshold) value. Measurements less than the cutoff are classified to the resting states.}
  \item{end_input}{The last element of the vector 'rest'. We fix it as \eqn{T + 0.5} throughout the entire algorithm and it varies by subjects.}
}

\value{
\item{active}{B-dimensional vector, denoted by \eqn{a_1,\ldots,a_B} in the paper.}
\item{rest}{B-dimensional vector, denoted by \eqn{r_1,\ldots,r_B}. Each element in 'active' and 'rest' constructs the active interval by (active[b], rest[b]) and the resting interval by (rest[b], active[b+1]) for \eqn{b=1,\ldots,B}.}
%%  ~Describe the value returned
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
y <- runif(10, 0, 10)
cutoff <- 3
end <- length(y) + 0.5
Init <- initial_position(y, cutoff, end)
active <- Init$active
rest <- Init$rest
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
