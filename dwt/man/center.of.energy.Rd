\name{center.of.energy}
\alias{center.of.energy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the center of energy of a filter}
\description{
Calculate the center of energy of the filter 'f'

References: Wickerhauser (1994, p.171 and p.341)
             Percival and Walden (2000, p.118)
}
\usage{
center.of.energy(f)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{
%%     ~~Describe \code{f} here~~
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
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.

% ## The function is currently defined as
% function (f)
%   ## ======================================================================
%   ## Calculate the center of energy of the filter 'f'
%   ## References: Wickerhauser (1994, p.171 and p.341)
%   ##             Percival and Walden (2000, p.118)
%   ## ======================================================================
% {
%   f.sq <- f^2
%   sum((0:(length(f)-1)) * f.sq) / sum(f.sq)
%   }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
