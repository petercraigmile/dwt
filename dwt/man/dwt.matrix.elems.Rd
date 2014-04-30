\name{dwt.matrix.elems}
\alias{dwt.matrix.elems}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the unique DWT matrix rows (modulo circular shifts) -- need to explain this better}
\description{
Calculate a subset of wavelet transform matrix for the dwt
          object 'dx'.  Only include the rows corresponding
          to the first wavelet coefficient per level, and the first
          level dx$nlevels scaling coefficient.
}
\usage{
dwt.matrix.elems(dx, use.C = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dx}{
%%     ~~Describe \code{dx} here~~
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
