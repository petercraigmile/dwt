\name{dwt}
\alias{dwt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the discrete wavelet transform (DWT)
}
\description{
Perform a DWT decomposition of the data 'x' to
          level 'nlevels' using a given 'wavelet'.

Assumes : Can compose to 'nlevels' levels;
          i.e., length of 'x' is divisible by 2^'nlevels'. 
}
\usage{
dwt(x, wavelet = "LA8", nlevels = floor(log2(length(x))), use.C = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{wavelet}{
%%     ~~Describe \code{wavelet} here~~
}
  \item{nlevels}{
%%     ~~Describe \code{nlevels} here~~
}
  \item{use.C}{
%%     ~~Describe \code{use.C} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ts }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
