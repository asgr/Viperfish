genSED=function(SFRbulge, SFRdisk, redshift=0.1, time=NULL, tau_birth=c(1,1), tau_screen=c(0.3,0.3), Zbulge=5, Zdisk=5, alpha_SF=c(1,1), AGNfrac=c(0,0), ab_nodust=TRUE, ap_nodust=TRUE, ab_dust=TRUE, ap_dust=TRUE, emitdust=TRUE, unimax=13.8e9, speclib=NULL, Dale=NULL, filtout=NULL, H0=67.8, sparse=1){

  if(is.null(time)){stop('Need time input!')}
  if(is.null(speclib)){stop('Need speclib (e.g. BC03lr)')}
  if(is.null(filtout)){stop('Need filtout input!')}

  if(emitdust & is.null(Dale)){stop('Need Dale input (e.g. Dale_Msol)')}

  if(missing(SFRbulge)){SFRbulge=rep(0,length(time))}
  if(missing(SFRdisk)){SFRdisk=rep(0,length(time))}

  if(length(SFRbulge)!=length(time)){stop('SFRbulge does not have the same length as time!')}
  if(length(SFRdisk)!=length(time)){stop('SFRdisk does not have the same length as time!')}

  if(length(tau_birth)==1){tau_birth=rep(tau_birth,2)}
  if(length(tau_screen)==1){tau_screen=rep(tau_screen,2)}
  if(length(alpha_SF)==1){alpha_SF=rep(alpha_SF,2)}
  if(length(AGNfrac)==1){AGNfrac=rep(AGNfrac,2)}

  assertNumeric(SFRbulge)
  assertNumeric(SFRdisk)
  assertScalar(redshift)
  assertNumeric(time)
  assertNumeric(tau_birth, len=2)
  assertNumeric(tau_screen, len=2)
  assertNumeric(Zbulge)
  assertNumeric(Zdisk)
  assertNumeric(alpha_SF, len=2)
  assertNumeric(AGNfrac, len=2)
  assertLogical(ab_nodust, len=1)
  assertLogical(ap_nodust, len=1)
  assertLogical(ab_dust, len=1)
  assertLogical(ap_dust, len=1)
  assertLogical(emitdust, len=1)

  SFRbulge[SFRbulge<0]=0
  SFRdisk[SFRdisk<0]=0

  SFRbulgefunc=approxfun(time, SFRbulge, rule=2, yleft=0, yright=0)
  SFRdiskfunc=approxfun(time, SFRdisk, rule=2, yleft=0, yright=0)

  if(length(Zbulge)>1){
    if(length(Zbulge)!=length(time)){stop('Zbulge does not have the same length as time!')}
    Zbulge[Zbulge<0]=0
    Zbulge=approxfun(time, Zbulge, rule=2, yleft=0, yright=0)
  }
  if(length(Zdisk)>1){
    if(length(Zdisk)!=length(time)){stop('Zdisk does not have the same length as time!')}
    Zdisk[Zdisk<0]=0
    Zdisk=approxfun(time, Zdisk, rule=2, yleft=0, yright=0)
  }

  Z=c(Zbulge, Zdisk)

  if(sparse>1){
    sparse=seq(1,dim(speclib$Zspec[[1]])[2],by=sparse)
    for(i in 1:length(speclib$Z)){
      speclib$Zspec[[i]]=speclib$Zspec[[i]][,sparse]
    }
    speclib$Wave=speclib$Wave[sparse]
  }

  lumwave=speclib$Wave

  if(ap_nodust){
    bulge_nodust=SFHfunc(SFRbulgefunc, tau_birth=0, tau_screen=0, speclib=speclib, z=redshift, Z=Z[[1]], outtype=NULL, unimax=unimax)
    disk_nodust=SFHfunc(SFRdiskfunc, tau_birth=0, tau_screen=0, speclib=speclib, z=redshift, Z=Z[[2]], outtype=NULL, unimax=unimax)
  }else if(ab_nodust | emitdust){
    bulge_nodust=SFHfunc(SFRbulgefunc, tau_birth=0, tau_screen=0, speclib=speclib, z=0, Z=Z[[1]], outtype=NULL)
    disk_nodust=SFHfunc(SFRdiskfunc, tau_birth=0, tau_screen=0, speclib=speclib, z=0, Z=Z[[2]], outtype=NULL)
  }

  if(ap_dust){
    bulge_dust=SFHfunc(SFRbulgefunc, tau_birth=tau_birth[1], tau_screen=tau_screen[1], speclib=speclib, z=redshift, Z=Z[[1]], outtype=NULL, unimax=unimax)
    disk_dust=SFHfunc(SFRdiskfunc, tau_birth=tau_birth[2], tau_screen=tau_screen[2], speclib=speclib, z=redshift, Z=Z[[2]], outtype=NULL, unimax=unimax)
  }else if(ab_dust){
    bulge_dust=SFHfunc(SFRbulgefunc, tau_birth=tau_birth[1], tau_screen=tau_screen[1], speclib=speclib, z=0, Z=Z[[1]], outtype=NULL)
    disk_dust=SFHfunc(SFRdiskfunc, tau_birth=tau_birth[2], tau_screen=tau_screen[2], speclib=speclib, z=0, Z=Z[[2]], outtype=NULL)
  }

  if((ab_dust | ap_dust) & emitdust & !is.null(Dale)){
    Dale_interp_b=Dale_interp(alpha_SF=alpha_SF[1], AGNfrac=AGNfrac[1], Dale=Dale)
    Dale_interp_d=Dale_interp(alpha_SF=alpha_SF[2], AGNfrac=AGNfrac[2], Dale=Dale)

    dustout_b=dustmass(speclib$Wave, bulge_nodust$lum, bulge_dust$lum, Dale$Wave, Dale_interp_b)
    dustout_d=dustmass(speclib$Wave, disk_nodust$lum, disk_dust$lum, Dale$Wave, Dale_interp_d)
    dustout_t=dustout_b+dustout_d
    #print(dustout_b)
    #print(dustout_d)

    if(ab_dust){
      #magplot(lumwave,bulge_dust$lum,log='xy',xlim=c(1e2,1e7),ylim=c(1e2,2e6),type='l')

      bulge_dust$lum=addspec(Dale$Wave, Dale_interp_b*dustout_b[1], lumwave, bulge_dust$lum)[,2]
      disk_dust$lum=addspec(Dale$Wave, Dale_interp_d*dustout_d[1], lumwave, disk_dust$lum)[,2]

      #lines(sort(c(Dale$Wave, lumwave)), bulge_dust$lum,col='red')

      lumwave=sort(c(Dale$Wave, lumwave))
    }

    if(ap_dust){
      dustflux_b=Lum2Flux(Dale$Wave, Dale_interp_b*dustout_b[1], z=redshift)
      dustflux_d=Lum2Flux(Dale$Wave, Dale_interp_d*dustout_d[1], z=redshift)

      bulge_dust$flux=addspec(dustflux_b[,1], dustflux_b[,2], bulge_dust$flux[,1], bulge_dust$flux[,2])
      disk_dust$flux=addspec(dustflux_d[,1], dustflux_d[,2], disk_dust$flux[,1], disk_dust$flux[,2])
    }
  }else{
    dustout_b=NA
    dustout_d=NA
    dustout_t=NA
  }

  if(ab_nodust){
    bulgelum_nodust=bulge_nodust$lum
    disklum_nodust=disk_nodust$lum
  }else{
    bulgelum_nodust=NA
    disklum_nodust=NA
  }

  if(ap_nodust & redshift>0){
    bulgeflux_nodust=bulge_nodust$flux
    diskflux_nodust=disk_nodust$flux
  }else{
    bulgeflux_nodust=NA
    diskflux_nodust=NA
  }

  if(ab_dust){
    bulgelum_dust=bulge_dust$lum
    disklum_dust=disk_dust$lum
  }else{
    bulgelum_dust=NA
    disklum_dust=NA
  }

  if(ap_dust & redshift>0){
    bulgeflux_dust=bulge_dust$flux
    diskflux_dust=disk_dust$flux
  }else{
    bulgeflux_dust=NA
    diskflux_dust=NA
  }

  #To get to absolute magnitude

  if(ab_nodust){
    ab_fluxnu_nodust_b=convert_wave2freq(bulgelum_nodust*3e-07, speclib$Wave) #3e-07 factor to get to absolute magnitude
    ab_fluxnu_nodust_d=convert_wave2freq(disklum_nodust*3e-07, speclib$Wave)
  }else{
    ab_fluxnu_nodust_b=NA
    ab_fluxnu_nodust_d=NA
  }

  if(ap_nodust & redshift>0){
    ap_fluxnu_nodust_b=convert_wave2freq(bulgeflux_nodust[,2], bulgeflux_nodust[,1])
    ap_fluxnu_nodust_d=convert_wave2freq(diskflux_nodust[,2], diskflux_nodust[,1])
  }else{
    ap_fluxnu_nodust_b=NA
    ap_fluxnu_nodust_d=NA
  }

  if(ab_dust){
    ab_fluxnu_dust_b=convert_wave2freq(bulgelum_dust*3e-07, lumwave)
    ab_fluxnu_dust_d=convert_wave2freq(disklum_dust*3e-07, lumwave)
  }else{
    ab_fluxnu_dust_b=NA
    ab_fluxnu_dust_d=NA
  }

  if(ap_dust & redshift>0){
    ap_fluxnu_dust_b=convert_wave2freq(bulgeflux_dust[,2], bulgeflux_dust[,1])
    ap_fluxnu_dust_d=convert_wave2freq(diskflux_dust[,2], diskflux_dust[,1])
  }else{
    ap_fluxnu_dust_b=NA
    ap_fluxnu_dust_d=NA
  }

  if(ab_nodust){
    ab_mag_nodust_b={}
    ab_mag_nodust_d={}
  }

  if(ap_nodust & redshift>0){
    ap_mag_nodust_b={}
    ap_mag_nodust_d={}
  }

  if(ab_dust){
    ab_mag_dust_b={}
    ab_mag_dust_d={}
  }

  if(ap_dust & redshift>0){
    ap_mag_dust_b={}
    ap_mag_dust_d={}
  }

  for(i in 1:length(filtout)){
    if(ab_nodust){
      ab_mag_nodust_b=c(ab_mag_nodust_b, bandpass(flux=ab_fluxnu_nodust_b, wave=speclib$Wave, filter=filtout[[i]], lum=TRUE)*1e23)
      ab_mag_nodust_d=c(ab_mag_nodust_d, bandpass(flux=ab_fluxnu_nodust_d, wave=speclib$Wave, filter=filtout[[i]], lum=TRUE)*1e23)
    }

    if(ap_nodust & redshift>0){
      ap_mag_nodust_b=c(ap_mag_nodust_b, bandpass(flux=ap_fluxnu_nodust_b, wave=bulgeflux_nodust[,1], filter=filtout[[i]], lum=TRUE)*1e23)
      ap_mag_nodust_d=c(ap_mag_nodust_d, bandpass(flux=ap_fluxnu_nodust_d, wave=diskflux_nodust[,1], filter=filtout[[i]], lum=TRUE)*1e23)
    }

    if(ab_dust){
      ab_mag_dust_b=c(ab_mag_dust_b, bandpass(flux=ab_fluxnu_dust_b, wave=lumwave, filter=filtout[[i]], lum=TRUE)*1e23)
      ab_mag_dust_d=c(ab_mag_dust_d, bandpass(flux=ab_fluxnu_dust_d, wave=lumwave, filter=filtout[[i]], lum=TRUE)*1e23)
    }

    if(ap_dust & redshift>0){
      ap_mag_dust_b=c(ap_mag_dust_b, bandpass(flux=ap_fluxnu_dust_b, wave=bulgeflux_dust[,1], filter=filtout[[i]], lum=TRUE)*1e23)
      ap_mag_dust_d=c(ap_mag_dust_d, bandpass(flux=ap_fluxnu_dust_d, wave=diskflux_dust[,1], filter=filtout[[i]], lum=TRUE)*1e23)
    }
  }

  if(ab_nodust){
    ab_mag_nodust_b=Jansky2magAB(ab_mag_nodust_b)
    ab_mag_nodust_d=Jansky2magAB(ab_mag_nodust_d)
    ab_mag_nodust_t=-2.5*log10(10^(-0.4*ab_mag_nodust_b)+10^(-0.4*ab_mag_nodust_d))
  }else{
    ab_mag_nodust_b=NA
    ab_mag_nodust_d=NA
    ab_mag_nodust_t=NA
  }

  if(ap_nodust & redshift>0){
    ap_mag_nodust_b=Jansky2magAB(ap_mag_nodust_b)
    ap_mag_nodust_d=Jansky2magAB(ap_mag_nodust_d)
    ap_mag_nodust_t=-2.5*log10(10^(-0.4*ap_mag_nodust_b)+10^(-0.4*ap_mag_nodust_d))
  }else{
    ap_mag_nodust_b=NA
    ap_mag_nodust_d=NA
    ap_mag_nodust_t=NA
  }

  if(ab_dust){
    ab_mag_dust_b=Jansky2magAB(ab_mag_dust_b)
    ab_mag_dust_d=Jansky2magAB(ab_mag_dust_d)
    ab_mag_dust_t=-2.5*log10(10^(-0.4*ab_mag_dust_b)+10^(-0.4*ab_mag_dust_d))
  }else{
    ab_mag_dust_b=NA
    ab_mag_dust_d=NA
    ab_mag_dust_t=NA
  }

  if(ap_dust & redshift>0){
    ap_mag_dust_b=Jansky2magAB(ap_mag_dust_b)
    ap_mag_dust_d=Jansky2magAB(ap_mag_dust_d)
    ap_mag_dust_t=-2.5*log10(10^(-0.4*ap_mag_dust_b)+10^(-0.4*ap_mag_dust_d))
  }else{
    ap_mag_dust_b=NA
    ap_mag_dust_d=NA
    ap_mag_dust_t=NA
  }

  output=data.frame(
    ab_mag_nodust_b=ab_mag_nodust_b+5*log10(H0/67.8),
    ab_mag_nodust_d=ab_mag_nodust_d+5*log10(H0/67.8),
    ab_mag_nodust_t=ab_mag_nodust_t+5*log10(H0/67.8),

    ap_mag_nodust_b=ap_mag_nodust_b,
    ap_mag_nodust_d=ap_mag_nodust_d,
    ap_mag_nodust_t=ap_mag_nodust_t,

    ab_mag_dust_b=ab_mag_dust_b+5*log10(H0/67.8),
    ab_mag_dust_d=ab_mag_dust_d+5*log10(H0/67.8),
    ab_mag_dust_t=ab_mag_dust_t+5*log10(H0/67.8),

    ap_mag_dust_b=ap_mag_dust_b,
    ap_mag_dust_d=ap_mag_dust_d,
    ap_mag_dust_t=ap_mag_dust_t)
  rownames(output)=names(filtout)
  invisible(output)
}
