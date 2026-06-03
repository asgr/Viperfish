\name{genSting}
\alias{genSting}
\alias{getSFHfull}
\alias{getSFHsing}
\alias{mockcone}
\alias{mockcone_extract}
\alias{mocksubsets}
\alias{extinction_extract}
\alias{process_spec_dump}
\title{
Convert Stingray Light-Cones to SED Products
}
\description{
Reads Stingray mock light-cone catalogs, matches galaxies back to Shark star-formation
history outputs, and computes SED-derived photometry or spectra.
}
\usage{
genSting(file_sting = NULL, path_shark = ".", path_out = ".",
    h = "get", cores_per_subvolume = 1, cores_per_snapshot = 4,
    children_outfile = "/dev/null", filters = c("FUV_GALEX",
        "NUV_GALEX", "u_SDSS", "g_SDSS", "r_SDSS", "i_SDSS",
        "Z_VISTA", "Y_VISTA", "J_VISTA", "H_VISTA", "K_VISTA",
        "W1_WISE", "W2_WISE", "W3_WISE", "W4_WISE", "P100_Herschel",
        "P160_Herschel", "S250_Herschel", "S350_Herschel", "S500_Herschel"),
    tau_birth = 1.5, tau_screen = 0.5, pow_birth = -0.7, pow_screen = -0.7,
    alpha_SF_birth = 1, alpha_SF_screen = 3, alpha_SF_AGN = 0,
    emission = FALSE, stellarpop = "BC03lr", read_extinct = FALSE,
    sparse = 5, time = NULL, mockcone = NULL, intSFR = TRUE,
    final_file_output = "Stingray-SED.csv", temp_file_output = "temp.csv",
    extinction_file = "extinction.hdf5", addradio_SF = FALSE,
    waveout = seq(2, 30, by = 0.01), ff_frac_SF = 0.1, ff_power_SF = -0.1,
    sy_power_SF = -0.8, reorder = TRUE, restart = FALSE, verbose = TRUE,
    write_final_file = FALSE, mode = "photom", spec_range = c(3600, 9600),
    spec_bin = 1)

getSFHfull(file_sting = "mocksurvey.hdf5", path_shark = ".", snapmax = 199,
  cores_per_subvolume = 4, verbose = TRUE)

getSFHsing(id_galaxy_sam, snapshot = NULL, subvolume = NULL, path_shark = ".")

mockcone_extract(file_sting = "mocksurvey.hdf5", reorder = TRUE)

extinction_extract(extinction_file = "extinction.hdf5", reorder = TRUE)

mocksubsets(mockcone)

process_spec_dump(temp_file_input, final_file_output = "Stingray-spec.csv",
  path_out = ".", wavegrid = 3000:10000)
}
\arguments{
  \item{file_sting}{
Character scalar; full path to a Stingray HDF5 file.
}
  \item{path_shark}{
Character scalar; path to the top-level Shark output directory containing snapshot/subvolume folders.
}
  \item{path_out}{
Character scalar; destination directory for final written outputs.
}
  \item{mockcone}{
Data table; output of \code{mockcone_extract}, containing IDs, snapshot/subvolume mapping, and redshifts.
}
\item{id_galaxy_sam}{
Integer vector; Shark galaxy IDs to extract for a specific snapshot/subvolume.
}
  \item{h}{
Numeric scalar; little-h value used to convert Shark SFR units (Msun/yr/h) into physical units. If set to \code{"get"}, values are read from Shark HDF5 cosmology metadata.
}
  \item{cores_per_subvolume}{
Integer scalar; number of worker processes over sub-snapshot groups.
}
  \item{cores_per_snapshot}{
Integer scalar; number of worker processes within each snapshot group.
}
  \item{children_outfile}{
Character scalar; file path passed to SNOW workers for child process logging.
}
  \item{snapmax}{
Integer scalar; maximum snapshot output by Shark.
}
  \item{filters}{
Character vector; filter names used to produce synthetic photometry.
}
  \item{tau_birth}{
Numeric scalar; dust optical depth for birth-cloud attenuation.
}
  \item{tau_screen}{
Numeric scalar; dust optical depth for diffuse screen attenuation.
}
  \item{pow_birth}{
Numeric scalar; power-law slope modifier for birth-cloud attenuation.
}
  \item{pow_screen}{
Numeric scalar; power-law slope modifier for screen attenuation.
}
  \item{alpha_SF_birth}{
Numeric scalar; Dale template \eqn{\alpha_{SF}} value for birth-cloud dust emission.
}
  \item{alpha_SF_screen}{
Numeric scalar; Dale template \eqn{\alpha_{SF}} value for screen dust emission.
}
  \item{alpha_SF_AGN}{
Numeric scalar; Dale template \eqn{\alpha_{SF}} value for AGN dust emission.
}
  \item{emission}{
Logical; should emission-line contributions be included.
}
  \item{stellarpop}{
Character scalar; one of \code{BC03lr}, \code{BC03hr}, \code{EMILES}, or \code{BPASS}.
}
  \item{read_extinct}{
Logical; should per-galaxy attenuation values be read from \code{extinction_file}.
}
  \item{sparse}{
Numeric scalar; spectral sampling factor used by ProSpect (higher values are faster and less accurate).
}
  \item{time}{
Numeric vector; lookback times (years). If omitted, values are read from Shark SFH outputs.
}
  \item{intSFR}{
Logical; integrate SFR between age grid points when evaluating mass formation.
}
  \item{final_file_output}{
Character scalar; output file name for the final catalog (CSV/HDF5 for photometry, CSV for spectra).
}
  \item{temp_file_output}{
Character scalar; temporary CSV accumulation file used during parallel processing.
}
  \item{temp_file_input}{
Character scalar; temporary CSV used as input to \code{process_spec_dump}.
}
  \item{extinction_file}{
Character scalar; HDF5 file containing per-galaxy extinction parameters.
}
  \item{addradio_SF}{
Logical; include star-formation radio continuum in outputs.
}
  \item{waveout}{
Numeric vector; log10 wavelength grid (Angstrom) used for spectral calculations.
}
  \item{wavegrid}{
Numeric vector; linear wavelength grid (Angstrom) used when re-binning dumped spectra.
}
  \item{ff_frac_SF}{
Numeric scalar; thermal free-free fraction at 1.4 GHz for star-formation radio emission.
}
  \item{ff_power_SF}{
Numeric scalar; free-free power-law slope.
}
  \item{sy_power_SF}{
Numeric scalar; synchrotron power-law slope.
}
  \item{write_final_file}{
Logical; if \code{TRUE}, write outputs to \code{path_out/final_file_output}.
}
  \item{mode}{
Character scalar; \code{"photom"} for filter magnitudes or \code{"spectrum"} for spectral outputs.
}
  \item{spec_range}{
Numeric length-two vector; wavelength range (Angstrom) retained in spectrum mode.
}
  \item{spec_bin}{
Integer scalar; binning factor used when writing spectrum-mode temporary output.
}
  \item{reorder}{
Logical; reorder mock-cone rows by \code{subsnapID} and \code{id_galaxy_sam} for faster grouped processing.
}
  \item{restart}{
Logical; continue from existing \code{temp_file_output} rather than starting from scratch.
}
  \item{verbose}{
Logical; print processing progress information.
}
  \item{snapshot}{
Integer scalar; snapshot index for \code{getSFHsing}.
}
  \item{subvolume}{
Integer scalar; subvolume index for \code{getSFHsing}.
}
}
\details{
\code{genSting} provides the high-level Stingray interface in Viperfish.
It groups Stingray rows by Shark snapshot/subvolume, reads matching Shark SFH datasets,
computes bulge/disk component SEDs through ProSpect, and emits either magnitudes or spectra.

The underlying Shark SFH interface expects the \code{star_formation_histories.hdf5} layout used by Shark,
including disk and bulge (disk-instability and merger) SFR and metallicity histories.

Absolute magnitudes are rest-frame quantities; apparent magnitudes are observer-frame quantities.
}
\value{
\code{genSting} returns either photometric columns (mode \code{"photom"}) or spectral samples (mode \code{"spectrum"}).

\code{getSFHfull} and \code{getSFHsing} return SFR/metallicity history lists.

\code{mockcone_extract} and \code{extinction_extract} return data tables extracted from HDF5 inputs.

\code{mocksubsets} returns grouped per-subsnapshot lookup metadata.

\code{process_spec_dump} runs for side effects: it reads temporary spectral dumps,
re-bins wavelengths, and writes the requested final file.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{genSED}}, \code{\link{genShark}}
}
\examples{
\dontrun{
test_def = genSting('localtest/mocksky_19.hdf5', 'localtest/')

test_spec = genSting('localtest/mocksky_19.hdf5', 'localtest/', emission = TRUE,
  mode = 'spec', stellarpop = 'BC03hr')
}
}
