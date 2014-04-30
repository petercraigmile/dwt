\name{dwt.next.filters}
\alias{dwt.next.filters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the next set of DWT wavelet and scaling filters, given the previous set}
\description{
Given the wavelet and scaling filters at level 'j-1'
           'prev.filters' calculate the level 'j' filters
           using the dwt.filter object 'filter'.

Assumes : 'j' is an integer > 1.
}
\usage{
dwt.next.filters(prev.filters, filter, j)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prev.filters}{
%%     ~~Describe \code{prev.filters} here~~
}
  \item{filter}{
%%     ~~Describe \code{filter} here~~
}
  \item{j}{
%%     ~~Describe \code{j} here~~
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
