\name{Viperfish-package}
\alias{Viperfish-package}
\alias{Viperfish}
\alias{viperfish}
\docType{package}
\title{
SEDs with Shark or Stingray
}
\description{
Makes SEDs with outputs from Shark and/or Stingray using ProSpect. Named "Viperfish"" due to its ability to make dark things (SAMs and dark-cones in this case) shine, a bit like the bioluminescent viperfish.

The basic flow is usually that SURFS makes an N-body suite, this is processed by Velociraptor and TreeFrog to make trees, Shark uses this to make a SAM output, Stingray turns this into a light-cone (but really a dark-cone at this stage), and Viperfish turns this into a true light-cone with SEDs attached to all of the galaxies and photometric outputs.
}
\details{
\tabular{ll}{
Package: \tab Viperfish\cr
Type: \tab Package\cr
Version: \tab 0.4.2\cr
Date: \tab 2019-06-14\cr
License: \tab LGPL-3\cr
Depends: \tab ProSpect (>= 0.6.0), data.table, hdf5r\cr
Imports: \tab celestial, bit64, foreach, doSNOW, checkmate, snow, ini\cr
}
}
\author{
Aaron Robotham: aaron.robotham@uwa.edu.au
Claudia Lagos
Rodrigo Tobar
}
\keyword{SED}
