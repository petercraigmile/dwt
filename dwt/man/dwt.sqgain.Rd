\name{dwt.sqgain}
\alias{dwt.sqgain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the square gain function for a wavelet or scaling filter at a given level, at the specified frequencies}
\description{
Calculates the square gain function for the level 'level'
          wavelet (if wavelet=T) or scaling (if wavelet=F) filter
          at frequencies 'freqs', using the dwt filter object, 'filter'.

Note    : This is a brute force method -- does not use the FFT
          of the filter since we can potentially evaluate on
          any set of frequencies.
}
\usage{
dwt.sqgain(freqs, level, filter, wavelet = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{freqs}{
%%     ~~Describe \code{freqs} here~~
}
  \item{level}{
%%     ~~Describe \code{level} here~~
}
  \item{filter}{
%%     ~~Describe \code{filter} here~~
}
  \item{wavelet}{
%%     ~~Describe \code{wavelet} here~~
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
