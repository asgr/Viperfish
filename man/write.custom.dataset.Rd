\name{write.custom.dataset}
\alias{write.custom.dataset}
\alias{write.group.safe}
\title{
Write HDF5 Groups and Datasets
}
\description{
Writes hierarchical groups and datasets safely into an HDF5 file.
}
\usage{
write.custom.dataset(filename = "temp.hdf5", group = "group", object, mode = "a",
dataset.name = "matrix", dataset.type = h5types$H5T_IEEE_F32LE, compression.level = 0,
overwrite = FALSE)
write.group.safe(filename = "temp.hdf5", group = "group")
}
\arguments{
  \item{filename}{
Character scalar; the path name of the HDF5 file requested.
}
  \item{group}{
Character scalar; the full path of the group. This will be created if not present.
}
  \item{object}{
Numeric matrix; the object to be written out to HDF5.
}
  \item{mode}{
Character scalar; mode to use when writing the HDF5 (see\code{\link{H5File}}).
}
  \item{dataset.name}{
Character scalar; the name of dataset within the HDF5 group.
}
  \item{dataset.type}{
HDF5 data type as per \code{\link{h5types}}.
}
  \item{compression.level}{
Integer scalar; 0-9 gzip compression level, where 0 is not compressed and 9 is highly compressed.
}
  \item{overwrite}{
Logical; should dataset be overwritten if already present?
}
}
\details{
\code{write.group.safe} will safely write a hierarchical group to a HDF5 file.

\code{write.custom.dataset} is mostly used to save compact 32 bit single numeric matrices to HDF5 datasets, rather than the native \code{R} doubles. Can save other types too by changing \option{dataset}.
}
\value{
Run for the side effect of making the desired HDF5 file.
}
\author{
Danail Obreschkow & Aaron Robotham
}

\seealso{
\code{\link{write.SED}}
}
\examples{
# Nothing yet!
}
