genSting=function(file_sting='mocksurvey.hdf5', path_shark='.', h=0.678, cores=4, snapmax=199, filters=c('FUV', 'NUV', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1', 'W2', 'W3', 'W4', 'P100', 'P160', 'S250', 'S350', 'S500'), tau_birth=1.5, tau_screen=0.5, sparse=5, time=NULL, mockcone=NULL, intSFR=TRUE, file_output='temp.csv', verbose=TRUE){

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
  assertScalar(tau_birth)
  assertScalar(tau_screen)
  assertInt(sparse)
  assertNumeric(time, null.ok=TRUE)
  assertDataTable(mockcone, null.ok=TRUE)
  assertFlag(intSFR)
  assertFlag(verbose)

  BC03lr=Dale_Msol=Nid=id_galaxy_sam=idlist=snapshot=subsnapID=subvolume=z=i=j=mocksubsets=Ntime=zobs=NULL

  data("BC03lr", envir = environment())
  data("Dale_Msol", envir = environment())

  filtout=foreach(i = filters)%do%{getfilt(i)}
  names(filtout)=filters

  # if(!is.null(SFHfull) & doSFHbatch==TRUE){
  #   stop('You should not provide an input to SFHfull and have doSFHbatch=TRUE set, since the requested behaviour is ambiguous!')
  # }
  #
  # if(is.null(SFHfull) & doSFHbatch==FALSE){
  #   SFHfull=getSFHfull(file_sting=file_sting, path_shark=path_shark, snapmax=snapmax, cores=cores, verbose=verbose)
  # }
  #
  # if(!is.null(SFHfull)){
  #   SFRbulge=SFHfull$SFRbulge
  #   SFRdisk=SFHfull$SFRdisk
  #   Zbulge=SFHfull$Zbulge
  #   Zdisk=SFHfull$Zdisk
  # }

  #Make mock subsets:

  if(is.null(mockcone)){
    mockcone=mockcone(file_sting=file_sting)
  }else{
    assertDataTable(mockcone)
  }
  mocksubsets=mocksubsets(mockcone=mockcone)

  if(is.null(time)){
    assertAccess(paste(path_shark,snapmax,'0/star_formation_histories.hdf5', sep='/'), access='r')
    Shark_SFH=h5file(paste(path_shark,snapmax,'0/star_formation_histories.hdf5', sep='/'), mode='r')
    time=Shark_SFH[['lbt_mean']][]*1e9
    Shark_SFH$close()
  }
  Ntime=length(time)

  SEDlookup=data.table(id=unlist(mocksubsets$idlist), subsnapID=rep(mocksubsets$subsnapID, mocksubsets$Nid))

  cl=makeCluster(cores)
  registerDoSNOW(cl)

  #iterations=dim(mockcone)[1]
  subsnapIDs=unique(mockcone$subsnapID)

  if(verbose){
    pb = txtProgressBar(max = length(subsnapIDs), style = 3)
    progress = function(n) setTxtProgressBar(pb, n)
    opts = list(progress=progress)
  }

  file(file_output)
  assertAccess(file_output, access='w')

  outSED=foreach(i=1:length(subsnapIDs), .combine=.dumpout, .init=file_output, .final=.dumpin, .inorder=FALSE, .options.snow = if(verbose){opts})%dopar%{
    use=subsnapIDs[i]
    select=which(mockcone$subsnapID==use)
    snapshot=mockcone[select[1],snapshot]
    subvolume=mockcone[select[1],subvolume]
    id_galaxy=mockcone[select,id_galaxy]
    id_galaxy_sam=mockcone[select,id_galaxy_sam]
    zcos=mockcone[select,zcos]
    zobs=mockcone[select,zobs]
    SFHsing_subsnap=getSFHsing(id_galaxy_sam=id_galaxy_sam, snapshot=snapshot, subvolume=subvolume, path_shark=path_shark)

    SFRbulge_d_subsnap=SFHsing_subsnap$SFRbulge_d/h
    SFRbulge_m_subsnap=SFHsing_subsnap$SFRbulge_m/h
    SFRdisk_subsnap=SFHsing_subsnap$SFRdisk/h
    Zbulge_d_subsnap=SFHsing_subsnap$Zbulge_d/h
    Zbulge_m_subsnap=SFHsing_subsnap$Zbulge_m/h
    Zdisk_subsnap=SFHsing_subsnap$Zdisk/h

    tempout=foreach(j=1:length(select), .combine='rbind')%do%{
      unlist(genSED(SFRbulge_d=SFRbulge_d_subsnap[j,], SFRbulge_m=SFRbulge_m_subsnap[j,], SFRdisk=SFRdisk_subsnap[j,], redshift=zobs[j], time=time[1:dim(SFRdisk_subsnap)[2]]-cosdistTravelTime(zcos[j], ref='planck')*1e9, speclib=BC03lr, Zbulge_d=Zbulge_d_subsnap[j,], Zbulge_m=Zbulge_m_subsnap[j,], Zdisk=Zdisk_subsnap[j,], filtout=filtout, Dale=Dale_Msol, sparse=sparse, tau_birth=tau_birth, tau_screen=tau_screen, intSFR = intSFR))
    }
    as.data.table(cbind(id_galaxy,tempout))
  }

  stopCluster(cl)

  if(verbose){
    close(pb)
  }

  outSED=as.data.frame(outSED)
  colnamesSED=c(
    'id_galaxy',
    paste0('ab_mag_nodust_b_d_',filters),
    paste0('ab_mag_nodust_b_m_',filters),
    paste0('ab_mag_nodust_b_',filters),
    paste0('ab_mag_nodust_d_',filters),
    paste0('ab_mag_nodust_t_',filters),
    paste0('ap_mag_nodust_b_d_',filters),
    paste0('ap_mag_nodust_b_m_',filters),
    paste0('ap_mag_nodust_b_',filters),
    paste0('ap_mag_nodust_d_',filters),
    paste0('ap_mag_nodust_t_',filters),
    paste0('ab_mag_dust_b_d_',filters),
    paste0('ab_mag_dust_b_m_',filters),
    paste0('ab_mag_dust_b_',filters),
    paste0('ab_mag_dust_d_',filters),
    paste0('ab_mag_dust_t_',filters),
    paste0('ap_mag_dust_b_d_',filters),
    paste0('ap_mag_dust_b_m_',filters),
    paste0('ap_mag_dust_b_',filters),
    paste0('ap_mag_dust_d_',filters),
    paste0('ap_mag_dust_t_',filters)
    )
  colnames(outSED)=colnamesSED
  #outSED=cbind(id_galaxy=mockcone$id_galaxy, outSED)

  if(verbose){
    message(paste('Finished Viperfish on Stingray -',round(proc.time()[3]-timestart,3),'sec'))
  }

  #output=list(outSED=outSED, SFHfull=SFHfull)

  class(outSED)=c(class(outSED),'Viperfish-Shark')
  invisible(outSED)
}

mockcone=function(file_sting="mocksurvey.hdf5", galsname='galaxies', reorder=TRUE){
  subsnapID=snapshot=subvolume=id_galaxy_sam=NULL
  assertCharacter(file_sting, max.len=1)
  assertAccess(file_sting, access='r')
  mocksurvey=h5file(file_sting, mode='r')[[galsname]]
  extract_col=list.datasets(mocksurvey, recursive = TRUE)
  mockcone=as.data.table(lapply(extract_col, function(x) mocksurvey[[x]][]))
  colnames(mockcone)=extract_col
  mocksurvey$close()
  mockcone[,subsnapID:=snapshot*100+subvolume]
  if(reorder){
    mockcone=mockcone[order(subsnapID,id_galaxy_sam),]
  }
  invisible(mockcone)
}

mocksubsets=function(mockcone){
  id_galaxy_sam=subsnapID=Nid=idlist=snapshot=subvolume=NULL
  mocksubsets=mockcone[,list(idlist=list(unique(id_galaxy_sam))),by=subsnapID]
  mocksubsets[,Nid:=length(unlist(idlist)),by=subsnapID]
  mocksubsets[,snapshot:=floor(subsnapID/100)]
  mocksubsets[,subvolume:=subsnapID%%100]
  invisible(mocksubsets)
}

.dumpout = function(file_output='temp.csv', ...) {
  for(r in list(...))
    fwrite(x=r, file=file_output, append=TRUE)
  file_output
}

.dumpin = function(file_output='temp.csv') {
  fread(file_output)
}


# .filedump=function(file_output='temp.csv', data){
#   fwrite(x=data, file='temp.csv', append=TRUE)
# }
