genSED=function(SFRbulge_d, SFRbulge_m, SFRdisk, redshift=0.1, time=NULL, tau_birth=1, 
                tau_screen=0.3, tau_AGN=1, pow_birth=-0.7, pow_screen=-0.7, pow_AGN=-0.7, 
                alpha_SF_birth=1, alpha_SF_screen=3, alpha_SF_AGN=0, emission=FALSE, IGMabsorb=FALSE, Zbulge_d=5, Zbulge_m=5, 
                Zdisk=5, AGNlum=0, ab_nodust=TRUE, ap_nodust=TRUE, ab_dust=TRUE, ap_dust=TRUE, 
                emitdust=TRUE, unimax=13.8e9, speclib=NULL, Dale=NULL, AGN=NULL, LKL10=NULL, filtout=NULL, 
                H0=67.8, sparse=5, intSFR=TRUE, addradio_SF=FALSE, waveout=seq(2,10,by=0.01), 
                ff_frac_SF=0.1, ff_power_SF=-0.1, sy_power_SF=-0.8, mode='photom', spec_range=c(3600,9600)){

  if(is.null(time)){stop('Need time input!')}
  if(is.null(speclib)){stop('Need speclib (e.g. BC03lr)')}
  if(is.null(filtout)){stop('Need filtout input!')}

  if(emitdust & is.null(Dale)){stop('Need Dale input (e.g. Dale_NormTot)')}

  if(missing(SFRbulge_d)){SFRbulge_d=rep(0,length(time))}
  if(missing(SFRbulge_m)){SFRbulge_m=rep(0,length(time))}
  if(missing(SFRdisk)){SFRdisk=rep(0,length(time))}

  if(length(SFRbulge_d)!=length(time)){stop('SFRbulge_d does not have the same length as time!')}
  if(length(SFRbulge_m)!=length(time)){stop('SFRbulge_m does not have the same length as time!')}
  if(length(SFRdisk)!=length(time)){stop('SFRdisk does not have the same length as time!')}

  if(length(tau_birth)==1){tau_birth=rep(tau_birth,3)}
  if(length(tau_screen)==1){tau_screen=rep(tau_screen,3)}
  if(length(pow_birth)==1){pow_birth=rep(pow_birth,3)}
  if(length(pow_screen)==1){pow_screen=rep(pow_screen,3)}

  if(length(alpha_SF_birth)==1){alpha_SF_birth=rep(alpha_SF_birth,3)}
  if(length(alpha_SF_screen)==1){alpha_SF_screen=rep(alpha_SF_screen,3)}
  if(length(alpha_SF_AGN)==1){alpha_SF_AGN=rep(alpha_SF_AGN,3)}

  if(length(AGNlum)==1){AGNlum=rep(AGNlum,3)}

  assertNumeric(SFRbulge_d)
  assertNumeric(SFRbulge_m)
  assertNumeric(SFRdisk)
  assertScalar(redshift)
  assertNumeric(time)
  assertNumeric(tau_birth, len=3)
  assertNumeric(tau_screen, len=3)
  assertNumeric(pow_birth, len=3)
  assertNumeric(pow_screen, len=3)
  assertNumeric(Zbulge_d)
  assertNumeric(Zbulge_m)
  assertNumeric(Zdisk)
  assertNumeric(alpha_SF_birth, len=3)
  assertNumeric(alpha_SF_screen, len=3)
  assertNumeric(alpha_SF_AGN, len=3)
  assertNumeric(AGNlum, len=3)
  assertFlag(ab_nodust)
  assertFlag(ap_nodust)
  assertFlag(ab_dust)
  assertFlag(ap_dust)
  assertFlag(emitdust)
  assertScalar(unimax)
  assertScalar(H0)
  assertInt(sparse)
  assertFlag(intSFR)
  assertFlag(addradio_SF)
  assertNumeric(ff_frac_SF)
  assertNumeric(ff_power_SF)
  assertNumeric(sy_power_SF)


  SFRbulge_d[!is.finite(SFRbulge_d)] = 0
  SFRbulge_m[!is.finite(SFRbulge_m)] = 0
  SFRdisk[!is.finite(SFRdisk)] = 0
  
  SFRbulge_d[which(SFRbulge_d < 0)] = 0
  SFRbulge_m[which(SFRbulge_m < 0)] = 0
  SFRdisk[which(SFRdisk < 0)] = 0

  SFRbulge_dfunc=approxfun(time, SFRbulge_d, rule=2, yleft=SFRbulge_d[which.min(time)], yright=SFRbulge_d[which.max(time)])
  SFRbulge_mfunc=approxfun(time, SFRbulge_m, rule=2, yleft=SFRbulge_m[which.min(time)], yright=SFRbulge_m[which.max(time)])
  SFRdiskfunc=approxfun(time, SFRdisk, rule=2, yleft=SFRdisk[which.min(time)], yright=SFRdisk[which.max(time)])

  if(length(Zbulge_d)>1){
    if(length(Zbulge_d)!=length(time)){stop('Zbulge_d does not have the same length as time!')}
    Zbulge_d[!is.finite(Zbulge_d)] = 0
    Zbulge_d[which(Zbulge_d < 0)] = 0
    
    Zbulge_d_approx=approxfun(time, Zbulge_d, rule=2, yleft=Zbulge_d[which.min(time)], yright=Zbulge_d[which.max(time)])
    Zbulge_d=function(age, ...){Zbulge_d_approx(age)}
  }
  if(length(Zbulge_m)>1){
    if(length(Zbulge_m)!=length(time)){stop('Zbulge_m does not have the same length as time!')}
    Zbulge_m[!is.finite(Zbulge_m)] = 0
    Zbulge_m[which(Zbulge_m < 0)] = 0
    Zbulge_m_approx=approxfun(time, Zbulge_m, rule=2, yleft=Zbulge_m[which.min(time)], yright=Zbulge_m[which.max(time)])
    Zbulge_m=function(age, ...){Zbulge_m_approx(age)}
  }
  if(length(Zdisk)>1){
    if(length(Zdisk)!=length(time)){stop('Zdisk does not have the same length as time!')}
    Zdisk[!is.finite(Zdisk)] = 0
    Zdisk[which(Zdisk < 0)] = 0
    Zdisk_approx=approxfun(time, Zdisk, rule=2, yleft=Zdisk[which.min(time)], yright=Zdisk[which.max(time)])
    Zdisk=function(age, ...){Zdisk_approx(age)}
  }

  Z=c(Zbulge_d, Zbulge_m, Zdisk)

  if((ap_nodust==FALSE & ap_dust==FALSE) | redshift<=0){
    redshift=0
  }
  
  if(IGMabsorb){
    IGMabsorb = pnorm(redshift, mean=3.8, sd=1.2)
  }else{
    IGMabsorb = 0
  }

  if(mode == 'spectrum' | mode == 'spectra' | mode == 'spec' | mode == 'spectral'){
    waveout = NULL
  }
  
  bulge_d=ProSpectSED(massfunc=SFRbulge_dfunc, tau_birth=tau_birth[1], tau_screen=tau_screen[1], 
                      pow_birth=pow_birth[1], pow_screen=pow_screen[1], pow_AGN=pow_AGN[1], 
                      alpha_SF_birth=alpha_SF_birth[1], alpha_SF_screen=alpha_SF_screen[1], 
                      alpha_SF_AGN=alpha_SF_AGN[1], AGNlum=AGNlum[1], speclib=speclib, 
                      Dale=Dale, AGN=AGN, filters=NULL, filtout=NULL, z=redshift, Z=Z[[1]], outtype=NULL, 
                      unimax=unimax, intSFR=intSFR, sparse=sparse, emission=emission, LKL10=LKL10, IGMabsorb=IGMabsorb,
                      addradio_SF=addradio_SF, waveout=waveout,
                      ff_frac_SF=ff_frac_SF, ff_power_SF=ff_power_SF, sy_power_SF=sy_power_SF)

  bulge_m=ProSpectSED(massfunc=SFRbulge_mfunc, tau_birth=tau_birth[2], tau_screen=tau_screen[2], 
                      pow_birth=pow_birth[2], pow_screen=pow_screen[2], pow_AGN=pow_AGN[2], 
                      alpha_SF_birth=alpha_SF_birth[2], alpha_SF_screen=alpha_SF_screen[2], 
                      alpha_SF_AGN=alpha_SF_AGN[2], AGNlum=AGNlum[2], speclib=speclib, 
                      Dale=Dale, AGN=AGN, filters=NULL, filtout=NULL, z=redshift, Z=Z[[2]], outtype=NULL, 
                      unimax=unimax, intSFR=intSFR, sparse=sparse, emission=emission, LKL10=LKL10, IGMabsorb=IGMabsorb,
                      addradio_SF=addradio_SF, waveout=waveout,
                      ff_frac_SF=ff_frac_SF, ff_power_SF=ff_power_SF, sy_power_SF=sy_power_SF)

  disk=ProSpectSED(massfunc=SFRdiskfunc, tau_birth=tau_birth[3], tau_screen=tau_screen[3], 
                   pow_birth=pow_birth[3], pow_screen=pow_screen[3], pow_AGN=pow_AGN[3], 
                   alpha_SF_birth=alpha_SF_birth[3], alpha_SF_screen=alpha_SF_screen[3], 
                   alpha_SF_AGN=alpha_SF_AGN[3], AGNlum=AGNlum[3], speclib=speclib, 
                   Dale=Dale, AGN=AGN, filters=NULL, filtout=NULL, z=redshift, Z=Z[[3]], outtype=NULL, 
                   unimax=unimax, intSFR=intSFR, sparse=sparse, emission=emission, LKL10=LKL10, IGMabsorb=IGMabsorb,
                   addradio_SF=addradio_SF, waveout=waveout,
                   ff_frac_SF=ff_frac_SF, ff_power_SF=ff_power_SF, sy_power_SF=sy_power_SF)

  if(mode == 'spectrum' | mode == 'spectra' | mode == 'spec' | mode == 'spectral'){
    total = bulge_d$FinalFlux
    total$flux = total$flux + bulge_m$FinalFlux$flux
    total$flux = total$flux + disk$FinalFlux$flux
    
    total = total[total$wave >= spec_range[1] & total$wave <= spec_range[2],]
    return(total)
  }
  
  lir_dust_b_d=bulge_d$Stars$lumtot_atten
  lir_dust_b_m=bulge_m$Stars$lumtot_atten
  lir_dust_b=lir_dust_b_d + lir_dust_b_m
  lir_dust_d=disk$Stars$lumtot_atten
  lir_dust_t=lir_dust_b + lir_dust_d
  flux_ratio_ir_b_d=bulge_d$Stars$lumtot_birth / lir_dust_b_d
  flux_ratio_ir_b_m=bulge_m$Stars$lumtot_birth / lir_dust_b_m
  flux_ratio_ir_b_b=(bulge_d$Stars$lumtot_birth + bulge_m$Stars$lumtot_birth) / lir_dust_b
  flux_ratio_ir_d=disk$Stars$lumtot_birth / lir_dust_d
  flux_ratio_ir_t=(bulge_d$Stars$lumtot_birth + bulge_m$Stars$lumtot_birth + disk$Stars$lumtot_birth) / lir_dust_t

  if(ab_nodust){
    ab_mag_nodust_b_d=photom_flux(bulge_d$StarsUnAtten$wave, bulge_d$StarsUnAtten$lum*.lsol_to_absolute, filters = filtout)
    ab_mag_nodust_b_m=photom_flux(bulge_m$StarsUnAtten$wave, bulge_m$StarsUnAtten$lum*.lsol_to_absolute, filters = filtout)
    ab_mag_nodust_b=-2.5*log10(10^(-0.4*ab_mag_nodust_b_d)+10^(-0.4*ab_mag_nodust_b_m))
    ab_mag_nodust_d=photom_flux(disk$StarsUnAtten$wave, disk$StarsUnAtten$lum*.lsol_to_absolute, filters = filtout)
    ab_mag_nodust_t=-2.5*log10(10^(-0.4*ab_mag_nodust_b)+10^(-0.4*ab_mag_nodust_d))
  }else{
    ab_mag_nodust_b_d=NA
    ab_mag_nodust_b_m=NA
    ab_mag_nodust_b=NA
    ab_mag_nodust_d=NA
    ab_mag_nodust_t=NA
  }

  if(ab_dust){
    ab_mag_dust_b_d=photom_flux(bulge_d$FinalLum$wave, bulge_d$FinalLum$lum*.lsol_to_absolute, filters = filtout)
    ab_mag_dust_b_m=photom_flux(bulge_m$FinalLum$wave, bulge_m$FinalLum$lum*.lsol_to_absolute, filters = filtout)
    ab_mag_dust_b=-2.5*log10(10^(-0.4*ab_mag_dust_b_d)+10^(-0.4*ab_mag_dust_b_m))
    ab_mag_dust_d=photom_flux(disk$FinalLum$wave, disk$FinalLum$lum*.lsol_to_absolute, filters = filtout)
    ab_mag_dust_t=-2.5*log10(10^(-0.4*ab_mag_dust_b)+10^(-0.4*ab_mag_dust_d))
  }else{
    ab_mag_dust_b_d=NA
    ab_mag_dust_b_m=NA
    ab_mag_dust_b=NA
    ab_mag_dust_d=NA
    ab_mag_dust_t=NA
  }

  if(ap_nodust){
    ap_mag_nodust_b_d=photom_lum(bulge_d$StarsUnAtten$wave, bulge_d$StarsUnAtten$lum, filters = filtout, z = redshift)
    ap_mag_nodust_b_m=photom_lum(bulge_m$StarsUnAtten$wave, bulge_m$StarsUnAtten$lum, filters = filtout, z = redshift)
    ap_mag_nodust_b=-2.5*log10(10^(-0.4*ap_mag_nodust_b_d)+10^(-0.4*ap_mag_nodust_b_m))
    ap_mag_nodust_d=photom_lum(disk$StarsUnAtten$wave, disk$StarsUnAtten$lum, filters = filtout, z = redshift)
    ap_mag_nodust_t=-2.5*log10(10^(-0.4*ap_mag_nodust_b)+10^(-0.4*ap_mag_nodust_d))
  }else{
    ap_mag_nodust_b_d=NA
    ap_mag_nodust_b_m=NA
    ap_mag_nodust_b=NA
    ap_mag_nodust_d=NA
    ap_mag_nodust_t=NA
  }

  if(ap_dust){
    ap_mag_dust_b_d=photom_lum(bulge_d$FinalLum$wave, bulge_d$FinalLum$lum, filters = filtout, z = redshift)
    ap_mag_dust_b_m=photom_lum(bulge_m$FinalLum$wave, bulge_m$FinalLum$lum, filters = filtout, z = redshift)
    ap_mag_dust_b=-2.5*log10(10^(-0.4*ap_mag_dust_b_d)+10^(-0.4*ap_mag_dust_b_m))
    ap_mag_dust_d=photom_lum(disk$FinalLum$wave, disk$FinalLum$lum, filters = filtout, z = redshift)
    ap_mag_dust_t=-2.5*log10(10^(-0.4*ap_mag_dust_b)+10^(-0.4*ap_mag_dust_d))
  }else{
    ap_mag_dust_b_d=NA
    ap_mag_dust_b_m=NA
    ap_mag_dust_b=NA
    ap_mag_dust_d=NA
    ap_mag_dust_t=NA
  }

  output=data.frame(
    ab_mag_nodust_b_d=ab_mag_nodust_b_d+5*log10(H0/67.8),
    ab_mag_nodust_b_m=ab_mag_nodust_b_m+5*log10(H0/67.8),
    ab_mag_nodust_b=ab_mag_nodust_b+5*log10(H0/67.8),
    ab_mag_nodust_d=ab_mag_nodust_d+5*log10(H0/67.8),
    ab_mag_nodust_t=ab_mag_nodust_t+5*log10(H0/67.8),

    ap_mag_nodust_b_d=ap_mag_nodust_b_d,
    ap_mag_nodust_b_m=ap_mag_nodust_b_m,
    ap_mag_nodust_b=ap_mag_nodust_b,
    ap_mag_nodust_d=ap_mag_nodust_d,
    ap_mag_nodust_t=ap_mag_nodust_t,

    ab_mag_dust_b_d=ab_mag_dust_b_d+5*log10(H0/67.8),
    ab_mag_dust_b_m=ab_mag_dust_b_m+5*log10(H0/67.8),
    ab_mag_dust_b=ab_mag_dust_b+5*log10(H0/67.8),
    ab_mag_dust_d=ab_mag_dust_d+5*log10(H0/67.8),
    ab_mag_dust_t=ab_mag_dust_t+5*log10(H0/67.8),

    ap_mag_dust_b_d=ap_mag_dust_b_d,
    ap_mag_dust_b_m=ap_mag_dust_b_m,
    ap_mag_dust_b=ap_mag_dust_b,
    ap_mag_dust_d=ap_mag_dust_d,
    ap_mag_dust_t=ap_mag_dust_t,

    lir_dust_b_d=lir_dust_b_d * (H0/67.8)**2.0,
    lir_dust_b_m=lir_dust_b_m * (H0/67.8)**2.0,
    lir_dust_b=lir_dust_b * (H0/67.8)**2.0,
    lir_dust_d=lir_dust_d * (H0/67.8)**2.0,
    lir_dust_t=lir_dust_t * (H0/67.8)**2.0,

    flux_ratio_ir_b_d=flux_ratio_ir_b_d,
    flux_ratio_ir_b_m=flux_ratio_ir_b_m,
    flux_ratio_ir_b_b=flux_ratio_ir_b_b,
    flux_ratio_ir_d=flux_ratio_ir_d,
    flux_ratio_ir_t=flux_ratio_ir_t)

  #cbind(output, lir_dust_b_d, lir_dust_b_m, lir_dust_b, lir_dust_d)
  rownames(output)=names(filtout)
  invisible(output)
}
