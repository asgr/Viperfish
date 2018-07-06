genShark=function(path_shark='.', snapshot=199, subsnapshot=0, redshift=0.1, h=0.678, cores=4, select='all', filters=c('FUV', 'NUV', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1', 'W2', 'W3', 'W4', 'P100', 'P160', 'S250', 'S350', 'S500'), intSFR=TRUE, verbose=TRUE){

  timestart=proc.time()[3]

  if(verbose){
    message('Running Viperfish on Shark')
  }

  assertCharacter(path_shark, max.len=1)
  assertAccess(path_shark, access='r')
  assertInt(snapshot)
  assertInt(subsnapshot)
  assertScalar(h)
  assertInt(cores)
  assertCharacter(filters)

  BC03lr=Dale_Msol=SFH=i=subsnapID=snapshot=id_galaxy_sam=Nid=idlist=subsnapID=NULL

  data("BC03lr", envir = environment())
  data("Dale_Msol", envir = environment())

  filtout=foreach(i = filters)%do%{getfilt(i)}
  names(filtout)=filters

  if(!missing(snapshot)){path_shark=paste(path_shark,snapshot,sep='/')}
  if(!missing(subsnapshot)){path_shark=paste(path_shark,subsnapshot,sep='/')}

  assertAccess(paste(path_shark,'star_formation_histories.hdf5',sep='/'), access='r')
  Shark_SFH=h5file(paste(path_shark,'star_formation_histories.hdf5',sep='/'), mode='r')
  time=Shark_SFH[['LBT_mean']][]*1e9

  if(select[1]=='all'){
    select=1:SFH[['Galaxies/id_galaxy']]$dims
  }

  SFRbulge=Shark_SFH[['Bulges/StarFormationRateHistories']][,select,drop=FALSE]
  SFRdisk=Shark_SFH[['Disks/StarFormationRateHistories']][,select,drop=FALSE]
  Zbulge=Shark_SFH[['Bulges/MetallicityHistories']][,select,drop=FALSE]
  Zdisk=Shark_SFH[['Disks/MetallicityHistories']][,select,drop=FALSE]

  if(length(redshift)==1){
    redshift=rep(redshift, length(select))
  }

  if(length(select)!=length(redshift)){
    stop("Length of select does not equal length of redshift!")
  }

  cl=makeCluster(cores)
  registerDoSNOW(cl)

  iterations=length(select)

  if(verbose){
    pb = txtProgressBar(max = iterations, style = 3)
    progress = function(n) setTxtProgressBar(pb, n)
    opts = list(progress=progress)
  }

  outSED=foreach(i=1:iterations, .combine='rbind', .options.snow = if(verbose){opts})%dopar%{
  unlist(genSED(SFRbulge=SFRbulge[,i]/h, SFRdisk=SFRdisk[,i]/h, redshift=redshift[i], time=time, speclib=BC03lr, Zbulge=Zbulge[,i], Zdisk=Zdisk[,i], filtout=filtout, Dale=Dale_Msol, sparse=5, tau_birth = 1.5, tau_screen = 0.5, intSFR = intSFR))
  }

  stopCluster(cl)

  if(verbose){
    close(pb)
  }

  outSED=as.data.table(rbind(outSED))
  colnamesSED=c(
    paste0('ab_mag_nodust_b_',filters),
    paste0('ab_mag_nodust_d_',filters),
    paste0('ab_mag_nodust_t_',filters),
    paste0('ap_mag_nodust_b_',filters),
    paste0('ap_mag_nodust_d_',filters),
    paste0('ap_mag_nodust_t_',filters),
    paste0('ab_mag_dust_b_',filters),
    paste0('ab_mag_dust_d_',filters),
    paste0('ab_mag_dust_t_',filters),
    paste0('ap_mag_dust_b_',filters),
    paste0('ap_mag_dust_d_',filters),
    paste0('ap_mag_dust_t_',filters)
    )
  colnames(outSED)=colnamesSED
  outSED=cbind(id_galaxy=Shark_SFH[['Galaxies/id_galaxy']][select], outSED)

  Shark_SFH$close()

  if(verbose){
    message(paste('Finished Viperfish on Shark -',round(proc.time()[3]-timestart,3),'sec'))
  }

  class(outSED)='Viperfish-Shark'
  invisible(outSED)
}
