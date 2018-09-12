\name{genShark}
\alias{genShark}
\title{
Convert Shark Outputs to SEDs
}
\description{
Takes Shark boxes and converts to apparent and absolute SEDs.
}
\usage{
genShark(path_shark = ".", snapshot = NULL, subvolume = NULL, redshift = "get", h = 0.678,
cores = 4, id_galaxy_sam = "all", filters = c("FUV", "NUV", "u_SDSS", "g_SDSS", "r_SDSS",
"i_SDSS", "Z_VISTA", "Y_VISTA", "J_VISTA", "H_VISTA", "K_VISTA", "W1", "W2", "W3", "W4",
"P100", "P160", "S250", "S350", "S500"), tau_birth = 1.5, tau_screen = 0.5, sparse = 5,
intSFR=TRUE, verbose=TRUE, write.csv=FALSE)
}
\arguments{
  \item{path_shark}{
Character scalar; path to top directory containing Shark outputs. The package will check that the files and path exists and we have read access.
}
  \item{snapshot}{
Integer scalar; snapshot number to use.
}
  \item{subvolume}{
Integer scalar; sub-volume snapshot number to use.
}
  \item{redshift}{
Numeric vector; can be a single value for all galaxies, or a vector the same length as \option{select}. If set to "get" (default) it will inherit the redshift from the target \option{snapshot} and \option{subvolume}, which is probably what you want most of the time.
}
  \item{h}{
Numeric scalar; little-h assumed in Shark (needed to convert SFRs to actual SFR since they use Msun/yr/h).
}
  \item{cores}{
Integer scalar; number of cores to process on.
}
  \item{id_galaxy_sam}{
Integer vector; Shark IDs to extract for photometric processing. The default ('all') will process all galaxies in the sub-volume and snapshot requested.
}
  \item{filters}{
Character vector; filters to use for the photometry.
}
  \item{tau_birth}{
Numeric scalar; dust tau for birth clouds. Associated with sub 10 Myr star formation. See \code{\link{SFH}}.
}
  \item{tau_screen}{
Numeric scalar; dust tau for screen. See \code{\link{SFH}}.
}
  \item{sparse}{
Numeric scalar; amount of sparse sampling of the spectra to make. Higher values mean coarser spectra is processed, which means less accurate photometry, but faster processing. Generally done when computing higher resolution libraries for AB mags.
}
  \item{intSFR}{
Logical; should the \option{massfunc} be intergrated between ages? This might be necessary if the SFH is quite bursty, but it is more expensive to compute and not required if the SFH is quite smooth.
}
  \item{verbose}{
Logical; should helpful verbose outputs be printed regarding the processing?
}
  \item{write.csv}{
Logical; whether to write the output of this function into a .csv file within the \option{path_shark} directory.
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
