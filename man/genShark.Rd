\name{genShark}
\alias{genShark}
\title{
Convert Shark Outputs to SEDs
}
\description{
Takes Shark boxes and converts to apparent and absolute SEDs.
}
\usage{
genShark(path_shark = ".", snapshot = NULL, subvolume = NULL, redshift = "get", h = 'get',
cores = 4, id_galaxy_sam = "all", filters = c('FUV_GALEX', 'NUV_GALEX', 'u_SDSS',
'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA',
'W1_WISE', 'W2_WISE', 'W3_WISE', 'W4_WISE', 'P100_Herschel', 'P160_Herschel',
'S250_Herschel', 'S350_Herschel', 'S500_Herschel'), tau_birth = 1.5, tau_screen = 0.5,
sparse = 5, final_file_output = 'Shark-SED.csv', intSFR = TRUE, verbose = TRUE,
write_final_file = FALSE)
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
Numeric vector; can be a single value for all galaxies, or a vector the same length as \option{select}. If set to "get" (default) it will inherit the redshift from the target Shark \option{snapshot} and \option{subvolume}, which is probably what you want most of the time.
}
  \item{h}{
Numeric scalar; little-h assumed in Shark (needed to convert SFRs to actual SFR since they use Msun/yr/h). If set to "get" (default) it will inherit the redshift from the target Shark \option{snapshot} and \option{subvolume}, which is probably what you want most of the time.
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
  \item{final_file_output}{
Character scalar;  name for the final output file to write to. This should not be an existing file and end ".csv".
}
  \item{intSFR}{
Logical; should the \option{massfunc} be intergrated between ages? This might be necessary if the SFH is quite bursty, but it is more expensive to compute and not required if the SFH is quite smooth.
}
  \item{verbose}{
Logical; should helpful verbose outputs be printed regarding the processing?
}
  \item{write_final_file}{
Logical; whether to write the output of this function into a .csv file within the \option{path_shark} directory.
}
}
\details{
By default outputs all SEDs in the box.

Note the output absolute magnitudes are rest frame, whilst the apparent magnitudes are observed frame.

The apparent and absolute magnitudes are physical assuming a fixed H0=67.8 km/s/Mpc (Planck 2018) cosmology. This cannot be changed by the user (to remove silly H0 mistakes), so you should add 5.log10(0.678) (which is -0.844) to the absolute magnitudes if you want H0=100 km/s/Mpc type magnitudes. You should not adjust the apparent magnitudes since the Universe does not look different just because we change our mind about H0, and nor should our apparent magnitudes!
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
