genSting=function(file_sting='mocksurvey.hdf5', path_shark='.', h=0.678, cores=4, snapmax=199, filters=c('FUV', 'NUV', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1', 'W2', 'W3', 'W4', 'P100', 'P160', 'S250', 'S350', 'S500'), SFHlist=NULL, time=NULL, verbose=TRUE){

  timestart=proc.time()[3]

  if(verbose){
    message('Running Viperfish on Stingray')
  }

  assertCharacter(file_sting, max.len=1)
  assertAccess(file_sting, access='r')
  assertCharacter(path_shark, max.len=1)
  assertAccess(path_shark, access='r')
  assertScalar(h)
  assertInt(cores)
  assertInt(snapmax)
  assertCharacter(filters)
  assertList(SFHlist, null.ok=TRUE)
  assertFlag(verbose)

  BC03lr=Dale_Msol=Nid=id_galaxy_sam=idlist=snapshot=subsnapID=subsnapshot=z=i=mocksubsets=mockcone=Ntime=zobs=NULL

  data("BC03lr", envir = environment())
  data("Dale_Msol", envir = environment())

  filtout=foreach(i = filters)%do%{getfilt(i)}
  names(filtout)=filters

  if(is.null(SFHlist)){
    SFHlist=getSFH(file_sting=file_sting, path_shark=path_shark, snapmax=snapmax, cores=cores, verbose=verbose)
  }

  SFRbulge=SFHlist$SFRbulge
  SFRdisk=SFHlist$SFRdisk
  Zbulge=SFHlist$Zbulge
  Zdisk=SFHlist$Zdisk

  #Make mock subsets:

  mockcone=mockcone(file_sting=file_sting)
  mocksubsets=mocksubsets(mockcone=mockcone)

  if(is.null(time)){
    assertAccess(paste(path_shark,snapmax,'0/star_formation_histories.hdf5', sep='/'), access='r')
    SFH=h5file(paste(path_shark,snapmax,'0/star_formation_histories.hdf5', sep='/'), mode='r')
    time=SFH[['LBT_mean']][]*1e9
    SFH$close()
  }
  Ntime=length(time)

  SEDlookup=data.table(id=unlist(mocksubsets$idlist), subsnapID=rep(mocksubsets$subsnapID, mocksubsets$Nid))

  cl=makeCluster(cores)
  registerDoSNOW(cl)

  iterations=dim(mockcone)[1]

  if(verbose){
    pb = txtProgressBar(max = iterations, style = 3)
    progress = function(n) setTxtProgressBar(pb, n)
    opts = list(progress=progress)
  }

  outSED=foreach(i=1:iterations, .combine='rbind', .options.snow = if(verbose){opts})%dopar%{
    rowuse=which(SEDlookup$id==mockcone[i,id_galaxy_sam] & SEDlookup$subsnapID==mockcone[i,subsnapID])
    offset=(snapmax-mockcone[i,snapshot])
    unlist(genSED(SFRbulge=SFRbulge[rowuse,]/h, SFRdisk=SFRdisk[rowuse,]/h, redshift=mockcone[i,zobs], time=time-cosdistTravelTime(mockcone[i,zobs], ref='planck')*1e9, speclib=BC03lr, Zbulge=Zbulge[rowuse,], Zdisk=Zdisk[rowuse,], filtout=filtout, Dale=Dale_Msol, sparse=5, tau_birth=1.5, tau_screen=0.5))
  }

  stopCluster(cl)

  if(verbose){
    close(pb)
  }

  outSED=as.data.frame(outSED)
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
  outSED=cbind(id_galaxy_sky=mockcone$id_galaxy_sky, outSED)

  if(verbose){
    message(paste('Finished Viperfish on Stingray -',round(proc.time()[3]-timestart,3),'sec'))
  }

  output=list(outSED=outSED, SFHlist=SFHlist)
  class(output)='Viperfish-Sting'
  invisible(output)
}

mockcone=function(file_sting="mocksurvey.hdf5", galsname='galaxies'){
  subsnapID=snapshot=subsnapshot=NULL
  assertCharacter(file_sting, max.len=1)
  assertAccess(file_sting, access='r')
  mocksurvey=h5file(file_sting, mode='r')[[galsname]]
  extract_col=list.datasets(mocksurvey, recursive = TRUE)
  mockcone=as.data.table(lapply(extract_col, function(x) mocksurvey[[x]][]))
  colnames(mockcone)=extract_col
  mocksurvey$close()
  mockcone[,subsnapID:=snapshot*100+subsnapshot]
  invisible(mockcone)
}

mocksubsets=function(mockcone){
  id_galaxy_sam=subsnapID=Nid=idlist=snapshot=subsnapshot=NULL
  mocksubsets=mockcone[,list(idlist=list(unique(id_galaxy_sam))),by=subsnapID]
  mocksubsets[,Nid:=length(unlist(idlist)),by=subsnapID]
  mocksubsets[,snapshot:=floor(subsnapID/100)]
  mocksubsets[,subsnapshot:=subsnapID%%100]
  invisible(mocksubsets)
}
