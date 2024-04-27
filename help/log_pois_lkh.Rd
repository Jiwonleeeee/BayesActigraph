\name{log_pois_lkh}
\alias{log_pois_lkh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood for measurements in resting intervals.
}
\description{
Calculate the sum of the log likelihood of the mixture of Poisson distributions.
}
\usage{
log_pois_lkh(y_input, Pi_input, Lambda_input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y_input}{Scalar or vector of y of interest.} 
  \item{Pi_input}{K-dimensional vector of mixing proportions pi.} 
  \item{Lambda_input}{K-dimensional vector of parameters of Poisson distributions. } 
%%     ~~Describe \code{x} here~~
}
}

\value{
\item{sum_value}{Calculated sum of the log likelihood.} 
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
