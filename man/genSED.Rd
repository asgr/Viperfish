\name{genSED}
\alias{genSED}
\title{
Generate SED
}
\description{
Convert target bulge and disk SFHs to SEDs.
}
\usage{
genSED(SFRbulge_d, SFRbulge_m, SFRdisk, redshift = 0.1, time = NULL, tau_birth = 1,
tau_screen = 0.3, pow_birth = -0.7, pow_screen = -0.7, Zbulge_d = 5, Zbulge_m = 5,
Zdisk = 5, alpha_SF = 1, AGNfrac = 0, ab_nodust = TRUE, ap_nodust = TRUE, ab_dust = TRUE,
ap_dust = TRUE, emitdust = TRUE, unimax = 1.38e+10, speclib = NULL, Dale = NULL,
filtout = NULL, H0 = 67.8, sparse = 1, intSFR = TRUE)
}
\arguments{
  \item{SFRbulge_d}{
Numeric vector; star formation rate of the bulge (formed via disk instabilities) in units of physical Msun/yr (i.e. scaled by the correct value of H0, so a value near 70 km/s/Mpc). This is key since most simulations will natively output Msun/yr/h, so be careful here!
}
  \item{SFRbulge_m}{
Numeric vector; star formation rate of the bulge (formed via mergers) in units of physical Msun/yr (i.e. scaled by the correct value of H0, so a value near 70 km/s/Mpc). This is key since most simulations will natively output Msun/yr/h, so be careful here!
}
  \item{SFRdisk}{
Numeric vector; star formation rate of the disk in units of physical Msun/yr (i.e. scaled by the correct value of H0, so a value near 70 km/s/Mpc). This is key since most simulations will natively output Msun/yr/h, so be careful here!
}
  \item{redshift}{
Numeric scalar; redshift of the galaxy for computing apparent magnitudes. Must be above 0!
}
  \item{time}{
Numeric vector; lookback time in years. Must be same length as \option{SFRbulge} and \option{SFRdisk}.
}
  \item{tau_birth}{
Numeric scalar; dust tau for birth clouds. Associated with sub 10 Myr star formation. See \code{\link{SFH}}.
}
  \item{tau_screen}{
Numeric scalar; dust tau for screen. See \code{\link{SFH}}.
}
  \item{pow_birth}{
Numeric scalar; power to further raise the dust attenuation of birth clouds. See \code{\link{CF_birth}}.
}
  \item{pow_screen}{
Numeric scalar; power to further raise the dust attenuation of the dust screen. See \code{\link{CF_screen}}.
}
  \item{Zbulge_d}{
Numeric vector; metallicity of star formation in the bulge (formed via disk instabilities). Can be a single value or a vector that matches the history, in which case it must be same length as \option{SFRbulge}.
}
  \item{Zbulge_m}{
Numeric vector; metallicity of star formation in the bulge (formed via mergers). Can be a single value or a vector that matches the history, in which case it must be same length as \option{SFRbulge}.
}
  \item{Zdisk}{
Numeric vector; metallicity of star formation in the disk. Can be a single value or a vector that matches the history, in which case it must be same length as \option{SFRdisk}.
}
  \item{alpha_SF}{
Numeric scalar; desired interpolated alpha slope of the star forming population. See \code{\link{Dale_interp}}.
}
  \item{AGNfrac}{
Numeric scalar; desired interpolated AGN fraction. See \code{\link{Dale_interp}}.
}
  \item{ab_nodust}{
Logical; should absolute magnitudes for intrinsic SEDs be conmputed?
}
  \item{ap_nodust}{
Logical; should apparent magnitudes for intrinsic SEDs be conmputed?
}
  \item{ab_dust}{
Logical; should absolute magnitudes for dust attenuated SEDs be conmputed?
}
  \item{ap_dust}{
Logical; should apparent magnitudes for dust attenuated SEDs be conmputed?
}
  \item{emitdust}{
Logical; should attenuated light be re-emitted via Dale dust templates?
}
  \item{unimax}{
Numeric scalar; maximum allowed age of any stellar population relative to z=0 (i.e. today). Any star formation implied by massfunc that is older than this will be set to 0.
}
  \item{speclib}{
Pass in the spectral library directly. Must be one of \code{\link{BC03lr}}, \code{\link{BC03hr}}, \code{\link{EMILES}}. Doing this speeds up the compute time, since there is no need to lazy load from the package.
}
  \item{Dale}{
Pass in the Dale dust library directly. Must be one of \code{\link{Dale_Orig}}, \code{\link{Dale_Msol}}, \code{\link{Dale_NormTot}}, \code{\link{Dale_NormAGN}}, \code{\link{Dale_NormSFR}}. Doing this speeds up the compute time, since there is no need to lazy load from the package.
}
  \item{filtout}{
Pass in the photometric filters directly. If \code{filters} is a vector of filter names, this can be achieved efficiently with a command like: filtout=foreach(i = filters)\%do\%{getfilt(i)} (see \code{\link{getfilt}}).
}
  \item{H0}{
Numeric scalar; Hubble constant as defined at z=0 (default is H0=67.8 (km/s)/Mpc, i.e. Planck 2018). This only adjusts the output absolute magnitudes, i.e. it acts like it would do in observations where apparent magniatudes of a galaxy do not change just because we update our knowledge of H0, but the absolute values computed will change.
}
  \item{sparse}{
Numeric scalar; amount of sparse sampling of the spectra to make. Higher values mean coarser spectra is processed, which means less accurate photometry, but faster processing. Generally done when computing higher resolution libraries for AB mags.
}
  \item{intSFR}{
Logical; should the \option{massfunc} be intergrated between ages? This might be necessary if the SFH is quite bursty, but it is more expensive to compute and not required if the SFH is quite smooth.
}
}
\details{
Lower level function for computing bulge and disk SEDs. Most users will not use this function directly.

Note the output absolute magnitudes are rest frame, whilst the apparent magnitudes are observed frame.
}
\value{
List of named absolute and apparent magnitudes for bulge, disk and total photometry.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{genShark}}, \code{\link{genSting}}
}
\examples{
# Nothing yet!
}
