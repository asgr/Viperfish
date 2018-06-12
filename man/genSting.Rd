\name{genSting}
\alias{genSting}
\alias{getSFH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert Stingray Outputs to SEDs
}
\description{
Takes Stingray light cones and converts to apparent and absolute SEDs.
}
\usage{
genSting(mocksurvey = "mocksurvey.hdf5", path_shark = ".", h = 0.678, cores = 4,
snapmax = 199, filters = c("FUV", "NUV", "u_SDSS", "g_SDSS", "r_SDSS", "i_SDSS",
"Z_VISTA", "Y_VISTA", "J_VISTA", "H_VISTA", "K_VISTA", "W1", "W2", "W3", "W4", "P100",
"P160", "S250", "S350", "S500"), SFHlist = NULL)
getSFH(mocksurvey = "mocksurvey.hdf5", path_shark = ".", cores = 4, snapmax = 199)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mocksurvey}{
Character scalar; path to Stingray HDF5 file.
}
  \item{path_shark}{
Character scalar; path to top directory containing Shark outputs.
}
  \item{h}{
Numeric scalar; little-h assumed in Shark (needed to convert SFRs to actual SFR since they use Msun/yr/h).
}
  \item{cores}{
Integer scalar; number of cores to process on.
}
  \item{snapmax}{
Integer scalar; maximum snapshot output by Shark.
}
  \item{filters}{
Character vector; filters to use for the photometry.
}
  \item{SFHlist}{
List; output from \code{getSFH}. If present then these SFH matrices will be used directly (rather than value being stripped out of Shark, which can be slow).
}
}
\details{
This interface is pretty high level, so should work fine if the paths are correct.
}
\value{
\code{getSFH} returns a list of the SFRs and the Z for the bulge and disk (which can be passed to \code{genSting} directly).

\code{genSting} returns a data.frame of named absolute and apparent magnitudes for bulge, disk and total photometry.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{genSED}}, \code{\link{genShark}}
}
\examples{
# Nothing yet!
}
