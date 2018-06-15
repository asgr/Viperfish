\name{genShark}
\alias{genShark}
\title{
Convert Shark Outputs to SEDs
}
\description{
Takes Shark boxes and converts to apparent and absolute SEDs.
}
\usage{
genShark(path_shark = ".", snapshot, subsnapshot, redshift = 0.1, h = 0.678, cores = 4,
select = "all", filters = c("FUV", "NUV", "u_SDSS", "g_SDSS", "r_SDSS", "i_SDSS",
"Z_VISTA", "Y_VISTA", "J_VISTA", "H_VISTA", "K_VISTA", "W1", "W2", "W3", "W4", "P100",
"P160", "S250", "S350", "S500"))
}
\arguments{
  \item{path_shark}{
Character scalar; path to top directory containing Shark outputs. The package will check that the files and path exists and we have read access.
}
  \item{snapshot}{
Integer scalar; snapshot number to use.
}
  \item{subsnapshot}{
Integer scalar; sub-volume snapshot number to use.
}
  \item{redshift}{
Numeric vector; can be a single value for all galaxies, or a vector the same length as \option{select}.
}
  \item{h}{
Numeric scalar; little-h assumed in Shark (needed to convert SFRs to actual SFR since they use Msun/yr/h).
}
  \item{cores}{
Integer scalar; number of cores to process on.
}
  \item{select}{
Integer vector; IDs to extract for photometric processing. The default ('all') will process all galaxies in the sub-volume and snapshot requested.
}
  \item{filters}{
Character vector; filters to use for the photometry.
}
}
\details{
By default outputs all SEDs in the box.
}
\value{
Data.frame of named absolute and apparent magnitudes for bulge, disk and total photometry.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{genSED}}, \code{\link{genSting}}
}
\examples{
# Nothing yet!
}
