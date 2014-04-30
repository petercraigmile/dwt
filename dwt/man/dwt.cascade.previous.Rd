\name{dwt.cascade.previous}
\alias{dwt.cascade.previous}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the previous level of DWT scaling coefficients given the current scaling and wavelet coefficients.}
\description{
Given the scaling coefficients 'V' and wavelet coefficients 'W'
on scale j+1, calculate the level j scaling coefficients using
the dwt filter object 'dwt.filter.obj'.
}
\usage{
dwt.cascade.previous(V, W, dwt.filter.obj, use.C = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{V}{
%%     ~~Describe \code{V} here~~
}
  \item{W}{
%%     ~~Describe \code{W} here~~
}
  \item{dwt.filter.obj}{
%%     ~~Describe \code{dwt.filter.obj} here~~
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
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
