\name{genSting}
\alias{genSting}
\alias{getSFHfull}
\alias{getSFHsing}
\alias{mockcone}
\alias{mockcone_extract}
\alias{mocksubsets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert Stingray Outputs to SEDs
}
\description{
Takes Stingray light cones and converts to apparent and absolute SEDs.
}
\usage{
genSting(file_sting = NULL, path_shark = ".", h = 'get', cores_per_subvolume=1,
  cores_per_snapshot=4, children_outfile = '/dev/null', snapmax = 199,
  filters = c('FUV_GALEX', 'NUV_GALEX', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS',
  'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1_WISE', 'W2_WISE', 'W3_WISE',
  'W4_WISE', 'P100_Herschel', 'P160_Herschel', 'S250_Herschel', 'S350_Herschel',
  'S500_Herschel'),
  tau_birth = 1.5, tau_screen = 0.5, sparse = 5, time = NULL, mockcone = NULL,
  intSFR = TRUE, final_file_output = 'Stingray-SED.csv', temp_file_output = 'temp.csv',
  reorder = TRUE, restart = FALSE, verbose = TRUE, write_final_file = FALSE,
  mode = 'photom')

getSFHfull(file_sting = "mocksurvey.hdf5", path_shark = ".", snapmax = 199,
  cores_per_subvolume = 4, verbose = TRUE)

getSFHsing(id_galaxy_sam, snapshot = NULL, subvolume = NULL, path_shark = ".")

mockcone_extract(file_sting = "mocksurvey.hdf5", reorder = TRUE)

mocksubsets(mockcone)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_sting}{
Character scalar; full path to Stingray HDF5 file. The package will check that the file exists and we have read access.
}
  \item{mockcone}{
data.table; the output of \code{mockcone}.
}
  \item{path_shark}{
Character scalar; path to top directory containing Shark outputs. The package will check that the files and path exists and we have read access.
}
\item{id_galaxy_sam}{
Integer vector; a list of Shark SAM galaxies to be extracted from a given snapshot and subsnapshot.
}
  \item{h}{
Numeric scalar; little-h assumed in Shark (needed to convert SFRs to actual SFR since they use Msun/yr/h). If set to "get" (default) it will inherit the redshift from the target Stingray HDF5 file, which is probably what you want most of the time.
}
  \item{cores_per_subvolume}{
Integer scalar; number of cores to process a subvolume's snapshot on. This multiplied by \code{cores_per_snapshot} yields the total number of processes that will work in parallel.
}
  \item{cores_per_snapshot}{
Integer scalar; number of cores to process a snapshot's galaxies on. This multiplied by \code{cores_per_subvolume} yields the total number of processes that will work in parallel.
}
  \item{children_output}{
Character scalar; file where parallel worker's outputs will be dumped to. It is given as-is to SNOW's \code{makeCluster}.
}
  \item{snapmax}{
Integer scalar; maximum snapshot output by Shark.
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
  \item{time}{
Numeric vector; if provided this should be the full vector of lookback time corresponding to the time outputs in \option{SFHlist}. It should have the same length as the matrices in \option{SFHlist} have columns.
}
  \item{sparse}{
Numeric scalar; amount of sparse sampling of the spectra to make. Higher values mean coarser spectra is processed, which means less accurate photometry, but faster processing. Generally done when computing higher resolution libraries for AB mags.
}
  \item{intSFR}{
Boolean scalar; should the \option{massfunc} be intergrated between ages? This might be necessary if the SFH is quite bursty, but it is more expensive to compute and not required if the SFH is quite smooth.
}
  \item{final_file_output}{
Character scalar;  name for the final output file to write to. This should not be an existing file and end ".csv" or ".hdf5".
}
  \item{temp_file_output}{
Character scalar; name for the temporary file to write to. This should not be an existing file and end ".csv".
}
  \item{write_final_file}{
Logical; whether to write the output of this function into a csv file (with name specified by \option{final_file_output}) within the same directory as \option{file_sting}.
}
  \item{reorder}{
Logical; Should the output mock cone be re-ordered so to enforce that it is first ascending in subsnapID, followed by id_galaxy_sam within groups of subsnapID? It might improve speed for this to be TRUE (meaning subsets are contiguous), but it means the ordering will be different compared to the target mockcone.
}
  \item{restart}{
Logical; should \code{genSting} read in the current temporary file (as specified by \option{temp_file_output}) and carry on from where the code terminated? This is useful if the process is killed accidentally, or hits a machine wall time for running.
}
  \item{verbose}{
Logical; should helpful verbose outputs be printed regarding the processing?
}
  \item{snapshot}{
Integer scalar; snapshot number to use.
}
  \item{subvolume}{
Integer scalar; sub-volume snapshot number to use.
}
}
\details{
This interface is pretty high level, so should work fine if the paths are correct.

Note the output absolute magnitudes are rest frame, whilst the apparent magnitudes are observed frame.

The apparent and absolute magnitudes are physical assuming a fixed H0=67.8 km/s/Mpc (Planck 2018) cosmology. This cannot be changed by the user (to remove silly H0 mistakes), so you should add 5.log10(0.678) (which is -0.844) to the absolute magnitudes if you want H0=100 km/s/Mpc type magnitudes. You should not adjust the apparent magnitudes since the Universe does not look different just because we change our mind about H0, and nor should our apparent magnitudes!
}
\value{
\code{genSting} returns a list called \option{outSED} containing a data.frame of named absolute and apparent magnitudes for bulge, disk and total photometry.

\code{getSFHfull} returns a full list of the SFRs and the Z for the bulge (disk instability and merger driven) and disk for a given Stingray mock cone. This functions use is fairly deprecated now since we instead process on the fly in \code{genSting} (maintains processing speed for large datasets), but it might be useful for other purposes so it is being left available.

\code{getSFHsing} returns a list of the SFRs and the Z for the bulge (disk instability and merger driven) and disk for a subset of galaxies for a given snapshot and subsnapshot. Also returns the id_galaxy_sam actually present, and the location of the present galaxies relative to the input \option{id_galaxy_sam} (keep).

\code{mockcone_extract} outputs the contents of a target Stingray HDF5 into an R data.table. Mostly present for use inside \code{genSting}, but currently deprecated.

\code{mocksubsets} outputs the minimally required sub-volumes and SAM galaxy ID list from a target mockcone. Mostly present for use inside an older version of \code{genSting}, but currently deprecated.
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
