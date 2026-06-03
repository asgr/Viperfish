\name{genShark}
\alias{genShark}
\title{
Convert Shark Snapshot Outputs to SED Products
}
\description{
Reads Shark \code{star_formation_histories.hdf5} outputs and computes SED-derived
photometry or spectra for selected galaxies.
}
\usage{
genShark(path_shark = ".", path_out = ".", snapshot = NULL,
    subvolume = NULL, redshift = "get", h = "get", cores = 4,
    id_galaxy_sam = "all", filters = c("FUV_GALEX", "NUV_GALEX",
        "u_SDSS", "g_SDSS", "r_SDSS", "i_SDSS", "Z_VISTA", "Y_VISTA",
        "J_VISTA", "H_VISTA", "K_VISTA", "W1_WISE", "W2_WISE",
        "W3_WISE", "W4_WISE", "P100_Herschel", "P160_Herschel",
        "S250_Herschel", "S350_Herschel", "S500_Herschel"), tau_birth = 1,
    tau_screen = 0.3, tau_AGN = 1, pow_birth = -0.7, pow_screen = -0.7,
    pow_AGN = -0.7, alpha_SF_birth = 1, alpha_SF_screen = 3,
    alpha_SF_AGN = 0, emission = FALSE, stellarpop = "BC03lr",
    speclib = NULL, read_extinct = FALSE, sparse = 5,
    final_file_output = "Shark-SED.csv", extinction_file = "extinction.hdf5",
    intSFR = TRUE, addradio_SF = FALSE, waveout = seq(2, 30, by = 0.01),
    ff_frac_SF = 0.1, ff_power_SF = -0.1, sy_power_SF = -0.8,
    verbose = TRUE, write_final_file = FALSE, mode = "photom")
}
\arguments{
  \item{path_shark}{
Character scalar; top-level path containing Shark outputs organised by snapshot/subvolume.
}
  \item{path_out}{
Character scalar; destination directory for written outputs.
}
  \item{snapshot}{
Integer scalar; Shark snapshot number.
}
  \item{subvolume}{
Integer scalar; Shark subvolume number.
}
  \item{redshift}{
Numeric scalar/vector; if \code{"get"}, inherited from Shark run metadata for the selected snapshot/subvolume.
}
  \item{h}{
Numeric scalar; little-h conversion factor used to transform Shark SFR histories from Msun/yr/h to physical units. If \code{"get"}, read from Shark cosmology metadata.
}
  \item{cores}{
Integer scalar; number of parallel workers.
}
  \item{id_galaxy_sam}{
Integer vector or \code{"all"}; subset of Shark galaxy IDs to process.
}
  \item{filters}{
Character vector of photometric filters or a precomputed list of filter response functions.
}
  \item{tau_birth, tau_screen, tau_AGN}{
Numeric attenuation optical depth parameters for birth-cloud, screen, and AGN components.
}
  \item{pow_birth, pow_screen, pow_AGN}{
Numeric attenuation slope modifiers for birth-cloud, screen, and AGN components.
}
  \item{alpha_SF_birth, alpha_SF_screen, alpha_SF_AGN}{
Numeric \eqn{\alpha_{SF}} values controlling Dale dust templates for each component.
}
  \item{emission}{
Logical; include emission-line features in generated spectra.
}
  \item{stellarpop}{
Character scalar; one of \code{BC03lr}, \code{BC03hr}, \code{EMILES}, or \code{BPASS}.
}
  \item{speclib}{
Optional stellar population library object to avoid repeated lazy-loading.
}
  \item{read_extinct}{
Logical; if \code{TRUE}, read per-galaxy attenuation values from \code{extinction_file}.
}
  \item{sparse}{
Integer scalar; spectral sparsity factor used for faster photometric integration.
}
  \item{final_file_output}{
Character scalar; filename for optional written output.
}
  \item{extinction_file}{
Character scalar; HDF5 file containing extinction parameters.
}
  \item{intSFR}{
Logical; integrate mass formation over age bins.
}
  \item{addradio_SF}{
Logical; include star-formation radio continuum.
}
  \item{waveout}{
Numeric vector; log10 wavelength grid (Angstrom) for spectrum calculations.
}
  \item{ff_frac_SF, ff_power_SF, sy_power_SF}{
Numeric radio-continuum parameters passed to ProSpect radio models.
}
  \item{verbose}{
Logical; print progress information.
}
  \item{write_final_file}{
Logical; if \code{TRUE}, write outputs to \code{path_out}.
}
  \item{mode}{
Character scalar; \code{"photom"} for filter magnitudes or \code{"spectrum"} for sampled spectra.
}
}
\details{
By default all galaxies in the selected snapshot/subvolume are processed.

This function directly interfaces with the Shark SFH HDF5 schema and consumes
bulge and disk star-formation/metallicity histories, split by bulge growth channel.

Absolute magnitudes are rest-frame quantities, while apparent magnitudes are observer-frame quantities.
}
\value{
Data table of named photometric outputs (mode \code{"photom"}) or matrix-like spectral outputs
(mode \code{"spectrum"}).
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
