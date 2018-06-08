\name{write.SED}
\alias{write.SED}
\title{
Write SED to HDF5 file
}
\description{
Write SED output of \code{\link{genShark}} and/or \code{\link{genSting}}.
}
\usage{
write.SED(SED, filename = "temp.hdf5", overwrite = FALSE, filters = c("FUV", "NUV", "u_SDSS", "g_SDSS", "r_SDSS", "i_SDSS", "Z_VISTA", "Y_VISTA", "J_VISTA", "H_VISTA", "K_VISTA", "W1", "W2", "W3", "W4", "P100", "P160", "S250", "S350", "S500"))
}
\arguments{
  \item{SED}{
Numeric data.frame; output of \code{\link{genShark}} or \code{\link{genSting}}.
}
  \item{filename}{
Character scalar; the path name of the HDF5 file requested.
}
  \item{overwrite}{
Logical; should dataset be overwritten if already present?
}
  \item{filters}{
Character vector; filters to use for the photometry.
}
}
\details{
Writes out a sensible HDF5 file for SEDs.
}
\value{
Run for the side effect of making the desired HDF5 file.

Listing:
          name    obj_type dataset.dims dataset.type_class
           SED   H5I_GROUP         <NA>               <NA>
       ab_dust   H5I_GROUP         <NA>               <NA>
     ab_nodust   H5I_GROUP         <NA>               <NA>
       ap_dust   H5I_GROUP         <NA>               <NA>
     ap_nodust   H5I_GROUP         <NA>               <NA>
       filters H5I_DATASET           20         H5T_STRING
 id_galaxy_sky H5I_DATASET        32558       H5T_COMPOUND
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{write.custom.dataset}}, \code{\link{write.group.safe}}
}
\examples{
# Nothing yet!
}
