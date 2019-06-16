\name{write.SED}
\alias{write.SED}
\alias{write.SFH}
\title{
Write SED to HDF5 file
}
\description{
Write SED output of \code{\link{genShark}} and/or \code{\link{genSting}}.
}
\usage{
write.SED(SED, filters, outdir, filename, verbose = FALSE)

write.SFH(SFHlist, filename = 'temp.hdf5', overwrite = FALSE)
}
\arguments{
  \item{SED}{
Numeric data.frame; output of \code{\link{genShark}} or \code{\link{genSting}}.
}
  \item{filters}{
Character vector; filters to use for the photometry.
}
 \item{outdir}{
Character scalar; output directory.
}
  \item{filename}{
Character scalar; the path name of the HDF5 file requested.
}
  \item{verbose}{
Logical; should verbose output be provided?
}
  \item{SFHlist}{
List; output of \code{\link{getSFHfull}}.
}
  \item{overwrite}{
Logical; should dataset/s be overwritten if already present?
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
