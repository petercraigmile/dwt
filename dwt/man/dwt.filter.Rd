\name{dwt.filter}
\alias{dwt.filter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a DWT filter object given a filter name or set of scaling coefficients}
\description{
   If the argument 'scaling' is missing, creates a 'dwt.filter'
   object containing the scaling and wavelet filters of a given name
   'name', if it exists in the list of 'common.scaling.filters'.
   If the vector defining the scaling filter, 'scaling', is provided
   then define a 'dwt.filter' using that vector instead.
}
\usage{
dwt.filter(name, scaling)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
%%     ~~Describe \code{name} here~~
}
  \item{scaling}{
%%     ~~Describe \code{scaling} here~~
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
