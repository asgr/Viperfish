\name{genSED}
\alias{genSED}
\title{
Generate Bulge/Disk SED Components
}
\description{
Converts bulge and disk star-formation and metallicity histories into photometric
or spectral outputs via ProSpect.
}
\usage{
genSED(SFRbulge_d, SFRbulge_m, SFRdisk, redshift = 0.1, time = NULL,
    tau_birth = 1, tau_screen = 0.3, tau_AGN = 1, pow_birth = -0.7,
    pow_screen = -0.7, pow_AGN = -0.7, alpha_SF_birth = 1, alpha_SF_screen = 3,
    alpha_SF_AGN = 0, emission = FALSE, IGMabsorb = FALSE, Zbulge_d = 5,
    Zbulge_m = 5, Zdisk = 5, AGNlum = 0, ab_nodust = TRUE, ap_nodust = TRUE,
    ab_dust = TRUE, ap_dust = TRUE, emitdust = TRUE, unimax = 1.38e+10,
    speclib = NULL, Dale = NULL, AGN = NULL, LKL10 = NULL, filtout = NULL,
    H0 = 67.8, sparse = 5, intSFR = TRUE, addradio_SF = FALSE,
    waveout = seq(2, 10, by = 0.01), ff_frac_SF = 0.1, ff_power_SF = -0.1,
    sy_power_SF = -0.8, mode = "photom", spec_range = c(3600, 9600))
}
\arguments{
  \item{SFRbulge_d, SFRbulge_m, SFRdisk}{
Numeric vectors of star-formation histories in physical Msun/yr units for bulge (disk-instability),
bulge (merger), and disk components.
}
  \item{redshift}{
Numeric scalar; redshift used for observer-frame calculations.
}
  \item{time}{
Numeric vector; lookback times in years, matched to SFR history lengths.
}
  \item{tau_birth, tau_screen, tau_AGN}{
Numeric optical depth values for birth-cloud, screen, and AGN attenuation components.
}
  \item{pow_birth, pow_screen, pow_AGN}{
Numeric power-law attenuation slopes for birth-cloud, screen, and AGN components.
}
  \item{alpha_SF_birth, alpha_SF_screen, alpha_SF_AGN}{
Numeric \eqn{\alpha_{SF}} controls for Dale dust templates in each component.
}
  \item{emission}{
Logical; include emission-line features.
}
  \item{IGMabsorb}{
Logical; include intergalactic medium attenuation approximation.
}
  \item{Zbulge_d, Zbulge_m, Zdisk}{
Numeric scalar or vectors for stellar metallicity histories of each component.
}
  \item{AGNlum}{
Numeric scalar/vector; AGN luminosity scaling per component.
}
  \item{ab_nodust, ap_nodust, ab_dust, ap_dust}{
Logical flags controlling which magnitude products are returned.
}
  \item{emitdust}{
Logical; re-emit absorbed stellar energy using dust templates.
}
  \item{unimax}{
Numeric scalar; maximum stellar age allowed in years.
}
  \item{speclib}{
Stellar library object (e.g. \code{BC03lr}, \code{BC03hr}, \code{EMILES}, \code{BPASS}).
}
  \item{Dale}{
Dust template object (typically \code{Dale_NormTot}).
}
  \item{AGN}{
Optional AGN template object.
}
  \item{LKL10}{
Optional emission line template object.
}
  \item{filtout}{
List of filter response functions.
}
  \item{H0}{
Numeric scalar Hubble constant used to report absolute magnitudes and luminosities.
}
  \item{sparse}{
Integer scalar controlling spectral sparsity/performance trade-off.
}
  \item{intSFR}{
Logical; integrate star-formation history between age bins.
}
  \item{addradio_SF}{
Logical; include star-formation radio continuum component.
}
  \item{waveout}{
Numeric vector; log10 wavelength grid (Angstrom).
}
  \item{ff_frac_SF, ff_power_SF, sy_power_SF}{
Numeric parameters for star-formation radio continuum modelling.
}
  \item{mode}{
Character scalar; \code{"photom"} returns filter products, \code{"spectrum"} returns spectra.
}
  \item{spec_range}{
Numeric length-two vector; wavelength range retained for spectral outputs.
}
}
\details{
Lower-level function used by \code{genShark} and \code{genSting}.
Most users should prefer those interfaces.

Absolute magnitudes are rest-frame quantities, while apparent magnitudes are observer-frame quantities.
}
\value{
For \code{mode = "photom"}: data frame of component magnitudes and IR luminosity diagnostics.

For \code{mode = "spectrum"}: wavelength/flux table restricted to \code{spec_range}.
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
