\name{Viperfish-package}
\alias{Viperfish-package}
\alias{Viperfish}
\alias{viperfish}
\docType{package}
\title{
SED Generation for Shark and Stingray Outputs
}
\description{
Makes SEDs from ICRAR Shark SAM snapshots and Stingray light-cone products using ProSpect.
}
\details{
Named "Viperfish" due to its ability to make dark things (SAMs and dark-cones in this case) shine, a bit like the bioluminescent viperfish.

The usual pipeline is: SURFS N-body volumes are processed by Velociraptor and TreeFrog to produce merger trees, Shark evolves galaxies and writes snapshot products (including star-formation histories), Stingray projects those outputs into mock light-cones, and Viperfish turns those mock catalogs into photometric and/or spectral products.

For Shark interfaces, Viperfish reads the \code{star_formation_histories.hdf5} layout expected by Shark (e.g. \code{bulges_diskins/star_formation_rate_histories}, \code{bulges_mergers/star_formation_rate_histories}, \code{disks/star_formation_rate_histories}, and their metallicity counterparts), converts SFRs from native Msun/yr/h units to physical units via \code{h}, and forwards component histories to ProSpect.

\tabular{ll}{
Package: \tab Viperfish\cr
Type: \tab Package\cr
Version: \tab 0.6.0\cr
Date: \tab 2026-01-13\cr
License: \tab LGPL-3\cr
Depends: \tab ProSpect (>= 0.8.2), data.table, hdf5r\cr
Imports: \tab celestial (>= 1.4.5), bit64, foreach, doSNOW, checkmate, snow, ini, bigmemory\cr
}
}
\author{
Aaron Robotham: aaron.robotham@uwa.edu.au
Claudia Lagos
Rodrigo Tobar
}
\keyword{SED}
