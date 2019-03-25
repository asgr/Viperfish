genSED=function(SFRbulge_d, SFRbulge_m, SFRdisk, redshift=0.1, time=NULL, tau_birth=1, tau_screen=0.3, Zbulge_d=5, Zbulge_m=5, Zdisk=5, alpha_SF=1, AGNfrac=0, ab_nodust=TRUE, ap_nodust=TRUE, ab_dust=TRUE, ap_dust=TRUE, emitdust=TRUE, unimax=13.8e9, speclib=NULL, Dale=NULL, filtout=NULL, H0=67.8, sparse=1, intSFR=TRUE){

  if(is.null(time)){stop('Need time input!')}
  if(is.null(speclib)){stop('Need speclib (e.g. BC03lr)')}
  if(is.null(filtout)){stop('Need filtout input!')}

  if(emitdust & is.null(Dale)){stop('Need Dale input (e.g. Dale_Msol)')}

  if(missing(SFRbulge_d)){SFRbulge_d=rep(0,length(time))}
  if(missing(SFRbulge_m)){SFRbulge_m=rep(0,length(time))}
  if(missing(SFRdisk)){SFRdisk=rep(0,length(time))}

  if(length(SFRbulge_d)!=length(time)){stop('SFRbulge_d does not have the same length as time!')}
  if(length(SFRbulge_m)!=length(time)){stop('SFRbulge_m does not have the same length as time!')}
  if(length(SFRdisk)!=length(time)){stop('SFRdisk does not have the same length as time!')}

  if(length(tau_birth)==1){tau_birth=rep(tau_birth,3)}
  if(length(tau_screen)==1){tau_screen=rep(tau_screen,3)}
  if(length(alpha_SF)==1){alpha_SF=rep(alpha_SF,3)}
  if(length(AGNfrac)==1){AGNfrac=rep(AGNfrac,3)}

  assertNumeric(SFRbulge_d)
  assertNumeric(SFRbulge_m)
  assertNumeric(SFRdisk)
  assertScalar(redshift)
  assertNumeric(time)
  assertNumeric(tau_birth, len=3)
  assertNumeric(tau_screen, len=3)
  assertNumeric(Zbulge_d)
  assertNumeric(Zbulge_m)
  assertNumeric(Zdisk)
  assertNumeric(alpha_SF, len=3)
  assertNumeric(AGNfrac, len=3)
  assertFlag(ab_nodust)
  assertFlag(ap_nodust)
  assertFlag(ab_dust)
  assertFlag(ap_dust)
  assertFlag(emitdust)
  assertScalar(unimax)
  assertScalar(H0)
  assertInt(sparse)
  assertFlag(intSFR)

  SFRbulge_d[SFRbulge_d<0]=0
  SFRbulge_m[SFRbulge_m<0]=0
  SFRdisk[SFRdisk<0]=0

  SFRbulge_dfunc=approxfun(time, SFRbulge_d, rule=2, yleft=SFRbulge_d[which.min(time)], yright=SFRbulge_d[which.max(time)])
  SFRbulge_mfunc=approxfun(time, SFRbulge_m, rule=2, yleft=SFRbulge_m[which.min(time)], yright=SFRbulge_m[which.max(time)])
  SFRdiskfunc=approxfun(time, SFRdisk, rule=2, yleft=SFRdisk[which.min(time)], yright=SFRdisk[which.max(time)])

  if(length(Zbulge_d)>1){
    if(length(Zbulge_d)!=length(time)){stop('Zbulge_d does not have the same length as time!')}
    Zbulge_d[Zbulge_d<0]=0
    Zbulge_d=approxfun(time, Zbulge_d, rule=2, yleft=Zbulge_d[which.min(time)], yright=Zbulge_d[which.max(time)])
  }
  if(length(Zbulge_m)>1){
    if(length(Zbulge_m)!=length(time)){stop('Zbulge_m does not have the same length as time!')}
    Zbulge_m[Zbulge_m<0]=0
    Zbulge_m=approxfun(time, Zbulge_m, rule=2, yleft=Zbulge_m[which.min(time)], yright=Zbulge_m[which.max(time)])
  }
  if(length(Zdisk)>1){
    if(length(Zdisk)!=length(time)){stop('Zdisk does not have the same length as time!')}
    Zdisk[Zdisk<0]=0
    Zdisk=approxfun(time, Zdisk, rule=2, yleft=Zdisk[which.min(time)], yright=Zdisk[which.max(time)])
  }

  Z=c(Zbulge_d, Zbulge_m, Zdisk)

  if(sparse>1){
    sparse=seq(1,dim(speclib$Zspec[[1]])[2],by=sparse)
    for(i in 1:length(speclib$Z)){
      speclib$Zspec[[i]]=speclib$Zspec[[i]][,sparse]
    }
    speclib$Wave=speclib$Wave[sparse]
  }

  lumwave=speclib$Wave

  if(ap_nodust){
    bulge_d_nodust=SFHfunc(SFRbulge_dfunc, tau_birth=0, tau_screen=0, speclib=speclib, z=redshift, Z=Z[[1]], outtype=NULL, unimax=unimax, intSFR=intSFR, sparse=sparse)
    bulge_m_nodust=SFHfunc(SFRbulge_mfunc, tau_birth=0, tau_screen=0, speclib=speclib, z=redshift, Z=Z[[2]], outtype=NULL, unimax=unimax, intSFR=intSFR, sparse=sparse)
    disk_nodust=SFHfunc(SFRdiskfunc, tau_birth=0, tau_screen=0, speclib=speclib, z=redshift, Z=Z[[3]], outtype=NULL, unimax=unimax, intSFR=intSFR, sparse=sparse)
  }else if(ab_nodust | emitdust){
    bulge_d_nodust=SFHfunc(SFRbulge_dfunc, tau_birth=0, tau_screen=0, speclib=speclib, z=0, Z=Z[[1]], outtype=NULL, intSFR=intSFR, sparse=sparse)
    bulge_m_nodust=SFHfunc(SFRbulge_mfunc, tau_birth=0, tau_screen=0, speclib=speclib, z=0, Z=Z[[2]], outtype=NULL, intSFR=intSFR, sparse=sparse)
    disk_nodust=SFHfunc(SFRdiskfunc, tau_birth=0, tau_screen=0, speclib=speclib, z=0, Z=Z[[3]], outtype=NULL, intSFR=intSFR, sparse=sparse)
  }

  if(ap_dust){
    bulge_d_dust=SFHfunc(SFRbulge_dfunc, tau_birth=tau_birth[1], tau_screen=tau_screen[1], speclib=speclib, z=redshift, Z=Z[[1]], outtype=NULL, unimax=unimax, intSFR=intSFR, sparse=sparse)
    bulge_m_dust=SFHfunc(SFRbulge_mfunc, tau_birth=tau_birth[2], tau_screen=tau_screen[2], speclib=speclib, z=redshift, Z=Z[[2]], outtype=NULL, unimax=unimax, intSFR=intSFR, sparse=sparse)
    disk_dust=SFHfunc(SFRdiskfunc, tau_birth=tau_birth[3], tau_screen=tau_screen[3], speclib=speclib, z=redshift, Z=Z[[3]], outtype=NULL, unimax=unimax, intSFR=intSFR, sparse=sparse)
  }else if(ab_dust){
    bulge_d_dust=SFHfunc(SFRbulge_dfunc, tau_birth=tau_birth[1], tau_screen=tau_screen[1], speclib=speclib, z=0, Z=Z[[1]], outtype=NULL, intSFR=intSFR, sparse=sparse)
    bulge_m_dust=SFHfunc(SFRbulge_mfunc, tau_birth=tau_birth[2], tau_screen=tau_screen[2], speclib=speclib, z=0, Z=Z[[2]], outtype=NULL, intSFR=intSFR, sparse=sparse)
    disk_dust=SFHfunc(SFRdiskfunc, tau_birth=tau_birth[3], tau_screen=tau_screen[3], speclib=speclib, z=0, Z=Z[[3]], outtype=NULL, intSFR=intSFR, sparse=sparse)
  }

  if((ab_dust | ap_dust) & emitdust & !is.null(Dale)){
    Dale_interp_b_d=Dale_interp(alpha_SF=alpha_SF[1], AGNfrac=AGNfrac[1], Dale=Dale)
    Dale_interp_b_m=Dale_interp(alpha_SF=alpha_SF[2], AGNfrac=AGNfrac[2], Dale=Dale)
    Dale_interp_d=Dale_interp(alpha_SF=alpha_SF[3], AGNfrac=AGNfrac[3], Dale=Dale)

    dustout_b_d=dustmass(speclib$Wave, bulge_d_nodust$lum_atten, bulge_d_dust$lum_atten, Dale$Wave, Dale_interp_b_d[,2])
    dustout_b_m=dustmass(speclib$Wave, bulge_m_nodust$lum_atten, bulge_m_dust$lum_atten, Dale$Wave, Dale_interp_b_m[,2])
    dustout_d=dustmass(speclib$Wave, disk_nodust$lum_atten, disk_dust$lum_atten, Dale$Wave, Dale_interp_d[,2])
    #dustout_t=dustout_b_d+dustout_b_m+dustout_d
    #print(dustout_b)
    #print(dustout_d)

    if(ab_dust){
      #magplot(lumwave,bulge_dust$lum_atten,log='xy',xlim=c(1e2,1e7),ylim=c(1e2,2e6),type='l')
      bulge_d_dust$lum=addspec(Dale$Wave, Dale_interp_b_d[,2]*dustout_b_d[1], lumwave, bulge_d_dust$lum_atten)[,2]
      bulge_m_dust$lum=addspec(Dale$Wave, Dale_interp_b_m[,2]*dustout_b_m[1], lumwave, bulge_m_dust$lum_atten)[,2]
      disk_dust$lum=addspec(Dale$Wave, Dale_interp_d[,2]*dustout_d[1], lumwave, disk_dust$lum_atten)[,2]
      #lines(sort(c(Dale$Wave, lumwave)), bulge_dust$lum,col='red')

      lumwave=sort(c(Dale$Wave, lumwave))
    }

    if(ap_dust){
      dustflux_b_d=Lum2Flux(Dale$Wave, Dale_interp_b_d[,2]*dustout_b_d[1], z=redshift)
      dustflux_b_m=Lum2Flux(Dale$Wave, Dale_interp_b_m[,2]*dustout_b_m[1], z=redshift)
      dustflux_d=Lum2Flux(Dale$Wave, Dale_interp_d[,2]*dustout_d[1], z=redshift)

      bulge_d_dust$flux=addspec(dustflux_b_d[,1], dustflux_b_d[,2], bulge_d_dust$flux[,1], bulge_d_dust$flux[,2])
      bulge_m_dust$flux=addspec(dustflux_b_m[,1], dustflux_b_m[,2], bulge_m_dust$flux[,1], bulge_m_dust$flux[,2])
      disk_dust$flux=addspec(dustflux_d[,1], dustflux_d[,2], disk_dust$flux[,1], disk_dust$flux[,2])
    }
  }else{
    dustout_b_d=NA
    dustout_b_m=NA
    dustout_d=NA
    #dustout_t=NA
  }

  if(ab_nodust){
    bulge_d_nodust_lum=bulge_d_nodust$lum
    bulge_m_nodust_lum=bulge_m_nodust$lum
    disk_nodust_lum=disk_nodust$lum
  }else{
    bulge_d_nodust_lum=NA
    bulge_d_nodust_lum=NA
    disk_nodust_lum=NA
  }

  if(ap_nodust & redshift>0){
    bulge_d_nodust_flux=bulge_d_nodust$flux
    bulge_m_nodust_flux=bulge_m_nodust$flux
    disk_nodust_flux=disk_nodust$flux
  }else{
    bulge_d_nodust_flux=NA
    bulge_m_nodust_flux=NA
    disk_nodust_flux=NA
  }

  if(ab_dust){
    bulge_d_dust_lum=bulge_d_dust$lum
    bulge_m_dust_lum=bulge_m_dust$lum
    disk_dust_lum=disk_dust$lum
  }else{
    bulge_d_dust_lum=NA
    bulge_d_dust_lum=NA
    disk_dust_lum=NA
  }

  if(ap_dust & redshift>0){
    bulge_d_dust_flux=bulge_d_dust$flux
    bulge_m_dust_flux=bulge_m_dust$flux
    disk_dust_flux=disk_dust$flux
  }else{
    bulge_d_dust_flux=NA
    bulge_m_dust_flux=NA
    disk_dust_flux=NA
  }

  #To get to absolute magnitude

  if(ab_nodust){
    ab_fluxnu_nodust_b_d=convert_wave2freq(bulge_d_nodust_lum*3e-07, speclib$Wave) #3e-07 factor to get to absolute magnitude
    ab_fluxnu_nodust_b_m=convert_wave2freq(bulge_m_nodust_lum*3e-07, speclib$Wave)
    ab_fluxnu_nodust_d=convert_wave2freq(disk_nodust_lum*3e-07, speclib$Wave)
  }else{
    ab_fluxnu_nodust_b_d=NA
    ab_fluxnu_nodust_b_m=NA
    ab_fluxnu_nodust_d=NA
  }

  if(ap_nodust & redshift>0){
    ap_fluxnu_nodust_b_d=convert_wave2freq(bulge_d_nodust_flux[,2], bulge_d_nodust_flux[,1])
    ap_fluxnu_nodust_b_m=convert_wave2freq(bulge_m_nodust_flux[,2], bulge_m_nodust_flux[,1])
    ap_fluxnu_nodust_d=convert_wave2freq(disk_nodust_flux[,2], disk_nodust_flux[,1])
  }else{
    ap_fluxnu_nodust_b_d=NA
    ap_fluxnu_nodust_b_m=NA
    ap_fluxnu_nodust_d=NA
  }

  if(ab_dust){
    ab_fluxnu_dust_b_d=convert_wave2freq(bulge_d_dust_lum*3e-07, lumwave)
    ab_fluxnu_dust_b_m=convert_wave2freq(bulge_m_dust_lum*3e-07, lumwave)
    ab_fluxnu_dust_d=convert_wave2freq(disk_dust_lum*3e-07, lumwave)
  }else{
    ab_fluxnu_dust_b_d=NA
    ab_fluxnu_dust_b_m=NA
    ab_fluxnu_dust_d=NA
  }

  if(ap_dust & redshift>0){
    ap_fluxnu_dust_b_d=convert_wave2freq(bulge_d_dust_flux[,2], bulge_d_dust_flux[,1])
    ap_fluxnu_dust_b_m=convert_wave2freq(bulge_m_dust_flux[,2], bulge_m_dust_flux[,1])
    ap_fluxnu_dust_d=convert_wave2freq(disk_dust_flux[,2], disk_dust_flux[,1])
  }else{
    ap_fluxnu_dust_b_d=NA
    ap_fluxnu_dust_b_m=NA
    ap_fluxnu_dust_d=convert_wave2freq(disk_dust_flux[,2], disk_dust_flux[,1])
  }

  if(ab_nodust){
    ab_mag_nodust_b_d={}
    ab_mag_nodust_b_m={}
    ab_mag_nodust_d={}
  }

  if(ap_nodust & redshift>0){
    ap_mag_nodust_b_d={}
    ap_mag_nodust_b_m={}
    ap_mag_nodust_d={}
  }

  if(ab_dust){
    ab_mag_dust_b_d={}
    ab_mag_dust_b_m={}
    ab_mag_dust_d={}
  }

  if(ap_dust & redshift>0){
    ap_mag_dust_b_d={}
    ap_mag_dust_b_m={}
    ap_mag_dust_d={}
  }

  for(i in 1:length(filtout)){
    if(ab_nodust){
      ab_mag_nodust_b_d=c(ab_mag_nodust_b_d, bandpass(flux=ab_fluxnu_nodust_b_d, wave=speclib$Wave, filter=filtout[[i]], lum=TRUE)*1e23)
      ab_mag_nodust_b_m=c(ab_mag_nodust_b_m, bandpass(flux=ab_fluxnu_nodust_b_m, wave=speclib$Wave, filter=filtout[[i]], lum=TRUE)*1e23)
      ab_mag_nodust_d=c(ab_mag_nodust_d, bandpass(flux=ab_fluxnu_nodust_d, wave=speclib$Wave, filter=filtout[[i]], lum=TRUE)*1e23)
    }

    if(ap_nodust & redshift>0){
      ap_mag_nodust_b_d=c(ap_mag_nodust_b_d, bandpass(flux=ap_fluxnu_nodust_b_d, wave=bulge_d_nodust_flux[,1], filter=filtout[[i]], lum=TRUE)*1e23)
      ap_mag_nodust_b_m=c(ap_mag_nodust_b_m, bandpass(flux=ap_fluxnu_nodust_b_m, wave=bulge_m_nodust_flux[,1], filter=filtout[[i]], lum=TRUE)*1e23)
      ap_mag_nodust_d=c(ap_mag_nodust_d, bandpass(flux=ap_fluxnu_nodust_d, wave=disk_nodust_flux[,1], filter=filtout[[i]], lum=TRUE)*1e23)
    }

    if(ab_dust){
      ab_mag_dust_b_d=c(ab_mag_dust_b_d, bandpass(flux=ab_fluxnu_dust_b_d, wave=lumwave, filter=filtout[[i]], lum=TRUE)*1e23)
      ab_mag_dust_b_m=c(ab_mag_dust_b_m, bandpass(flux=ab_fluxnu_dust_b_m, wave=lumwave, filter=filtout[[i]], lum=TRUE)*1e23)
      ab_mag_dust_d=c(ab_mag_dust_d, bandpass(flux=ab_fluxnu_dust_d, wave=lumwave, filter=filtout[[i]], lum=TRUE)*1e23)
    }

    if(ap_dust & redshift>0){
      ap_mag_dust_b_d=c(ap_mag_dust_b_d, bandpass(flux=ap_fluxnu_dust_b_d, wave=bulge_d_dust_flux[,1], filter=filtout[[i]], lum=TRUE)*1e23)
      ap_mag_dust_b_m=c(ap_mag_dust_b_m, bandpass(flux=ap_fluxnu_dust_b_m, wave=bulge_m_dust_flux[,1], filter=filtout[[i]], lum=TRUE)*1e23)
      ap_mag_dust_d=c(ap_mag_dust_d, bandpass(flux=ap_fluxnu_dust_d, wave=disk_dust_flux[,1], filter=filtout[[i]], lum=TRUE)*1e23)
    }
  }

  if(ab_nodust){
    ab_mag_nodust_b_d=Jansky2magAB(ab_mag_nodust_b_d)
    ab_mag_nodust_b_m=Jansky2magAB(ab_mag_nodust_b_m)
    ab_mag_nodust_d=Jansky2magAB(ab_mag_nodust_d)
    ab_mag_nodust_b=-2.5*log10(10^(-0.4*ab_mag_nodust_b_d)+10^(-0.4*ab_mag_nodust_b_m))
    ab_mag_nodust_t=-2.5*log10(10^(-0.4*ab_mag_nodust_b_d)+10^(-0.4*ab_mag_nodust_b_m)+10^(-0.4*ab_mag_nodust_d))
  }else{
    ab_mag_nodust_b_d=NA
    ab_mag_nodust_b_m=NA
    ab_mag_nodust_d=NA
    ab_mag_nodust_b=NA
    ab_mag_nodust_t=NA
  }

  if(ap_nodust & redshift>0){
    ap_mag_nodust_b_d=Jansky2magAB(ap_mag_nodust_b_d)
    ap_mag_nodust_b_m=Jansky2magAB(ap_mag_nodust_b_m)
    ap_mag_nodust_d=Jansky2magAB(ap_mag_nodust_d)
    ap_mag_nodust_b=-2.5*log10(10^(-0.4*ap_mag_nodust_b_d)+10^(-0.4*ap_mag_nodust_b_m))
    ap_mag_nodust_t=-2.5*log10(10^(-0.4*ap_mag_nodust_b_d)+10^(-0.4*ap_mag_nodust_b_m)+10^(-0.4*ap_mag_nodust_d))
  }else{
    ap_mag_nodust_b_d=NA
    ap_mag_nodust_b_m=NA
    ap_mag_nodust_d=NA
    ap_mag_nodust_b=NA
    ap_mag_nodust_t=NA
  }

  if(ab_dust){
    ab_mag_dust_b_d=Jansky2magAB(ab_mag_dust_b_d)
    ab_mag_dust_b_m=Jansky2magAB(ab_mag_dust_b_m)
    ab_mag_dust_d=Jansky2magAB(ab_mag_dust_d)
    ab_mag_dust_b=-2.5*log10(10^(-0.4*ab_mag_dust_b_d)+10^(-0.4*ab_mag_dust_b_m))
    ab_mag_dust_t=-2.5*log10(10^(-0.4*ab_mag_dust_b_d)+10^(-0.4*ab_mag_dust_b_m)+10^(-0.4*ab_mag_dust_d))
  }else{
    ab_mag_dust_b_d=NA
    ab_mag_dust_b_m=NA
    ab_mag_dust_d=NA
    ab_mag_dust_b=NA
    ab_mag_dust_t=NA
  }

  if(ap_dust & redshift>0){
    ap_mag_dust_b_d=Jansky2magAB(ap_mag_dust_b_d)
    ap_mag_dust_b_m=Jansky2magAB(ap_mag_dust_b_m)
    ap_mag_dust_d=Jansky2magAB(ap_mag_dust_d)
    ap_mag_dust_b=-2.5*log10(10^(-0.4*ap_mag_dust_b_d)+10^(-0.4*ap_mag_dust_b_m))
    ap_mag_dust_t=-2.5*log10(10^(-0.4*ap_mag_dust_b_d)+10^(-0.4*ap_mag_dust_b_m)+10^(-0.4*ap_mag_dust_d))
  }else{
    ap_mag_dust_b_d=NA
    ap_mag_dust_b_m=NA
    ap_mag_dust_d=NA
    ap_mag_dust_b=NA
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
    ap_mag_dust_t=ap_mag_dust_t)
  rownames(output)=names(filtout)
  invisible(output)
}
