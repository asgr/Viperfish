\name{genSting}
\alias{genSting}
\alias{getSFH}
\alias{getSFHlist}
\alias{mockcone}
\alias{mocksubsets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert Stingray Outputs to SEDs
}
\description{
Takes Stingray light cones and converts to apparent and absolute SEDs.
}
\usage{
genSting(file_sting = "mocksurvey.hdf5", path_shark = ".", h = 0.678, cores = 4,
snapmax = 199, filters = c("FUV", "NUV", "u_SDSS", "g_SDSS", "r_SDSS", "i_SDSS",
"Z_VISTA", "Y_VISTA", "J_VISTA", "H_VISTA", "K_VISTA", "W1", "W2", "W3", "W4", "P100",
"P160", "S250", "S350", "S500"), tau_birth = 1.5, tau_screen = 0.5, sparse = 5,
SFHlist = NULL, time = NULL, mockcone = NULL, intSFR = TRUE, doSFHbatch = TRUE,
verbose = TRUE)
getSFH(file_sting = "mocksurvey.hdf5", path_shark = ".", snapmax = 199, cores = 4,
verbose = TRUE)
getSFHlist(id_galaxy_sam, snapshot = 199, subsnapshot = 0, path_shark = ".")
mockcone(file_sting = "mocksurvey.hdf5", galsname='galaxies')
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
  \item{tau_birth}{
Numeric scalar; dust tau for birth clouds. Associated with sub 10 Myr star formation. See \code{\link{SFH}}.
}
  \item{tau_screen}{
Numeric scalar; dust tau for screen. See \code{\link{SFH}}.
}
  \item{SFHlist}{
List; output from \code{getSFH}. If present then these SFH matrices will be used directly (rather than value being stripped out of Shark, which can be slow). Only relevant if \option{doSFHbatch}=FALSE.
}
  \item{time}{
Numeric vector; if provided this should be the full vector of lookback time corresponding to the time outputs in \option{SFHlist}. It should have the same length as the matrices in \option{SFHlist} have columns.
}
  \item{galsname}{
Character scalar; the name of the group containing galaxies in the Stingray HDF5 output (to cope with name changes through versions).
}
  \item{sparse}{
Numeric scalar; amount of sparse sampling of the spectra to make. Higher values mean coarser spectra is processed, which means less accurate photometry, but faster processing. Generally done when computing higher resolution libraries for AB mags.
}
  \item{intSFR}{
Boolean scalar; should the \option{massfunc} be intergrated between ages? This might be necessary if the SFH is quite bursty, but it is more expensive to compute and not required if the SFH is quite smooth.
}
  \item{doSFHbatch}{
Boolean scalar; should the star formation histories be computed in batches as \code{genSting} is run? It might be faster to set this to FALSE since then all SFHs are extracted in one go, but for large cones (millions of galaxies) this can be too memory itensive, so it is better to set this to TRUE. Setting to TRUE means the SFHlist output component of \code{genSting} is NULL. \option{doSFHbatch} overrides the provision of a matrix of SFHs to \option{SFHlist}.
}
  \item{snapshot}{
Integer scalar; snapshot number to use.
}
  \item{subsnapshot}{
Integer scalar; sub-volume snapshot number to use.
}
  \item{verbose}{
Logical; should helpful verbose outputs be printed regarding the processing?
}
}
\details{
This interface is pretty high level, so should work fine if the paths are correct.
}
\value{
\code{getSFH} returns a full list of the SFRs and the Z for the bulge and disk (which can be passed to \code{genSting} directly) for a given Stingray mock cone.

\code{getSFH} returns a list of the SFRs and the Z for the bulge and disk for a subset of galaxies for a given snapshot and subsnapshot.

\code{genSting} returns a list called \option{outSED} containing a data.frame of named absolute and apparent magnitudes for bulge, disk and total photometry, and a list called \option{SFHlist} containing the output of \code{getSFH} (or the input \option{SFHlist} if provided.

\code{mockcone} outputs the contents of a target Stingray HDF5 into an R data.table. Mostly present for use inside \code{genSting}.

\code{mocksubsets} outputs the minimally required sub-volumes and SAM galaxy ID list from a target mockcone. Mostly present for use inside \code{genSting}.
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
