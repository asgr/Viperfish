getSFHfull=function(file_sting='mocksurvey.hdf5', path_shark='.', snapmax=199, cores=4, verbose=TRUE){

  timestart=proc.time()[3]

  if(verbose){
    message('Running getSFHfull on Stingray')
  }

  assertCharacter(file_sting, max.len=1)
  assertAccess(file_sting, access='r')
  assertCharacter(path_shark, max.len=1)
  assertAccess(path_shark, access='r')
  assertInt(snapmax)
  assertInt(cores)
  assertFlag(verbose)

  BC03lr=Dale_Msol=Nid=id_galaxy_sam=idlist=snapshot=subsnapID=subvolume=z=i=mocksubsets=mockcone=Ntime=time=NULL

  timestart=proc.time()[3]

  sfh_fname = paste(path_shark,snapmax,'0/star_formation_histories.hdf5', sep='/')
  assertAccess(sfh_fname, access='r')
  Shark_SFH=h5file(sfh_fname, mode='r')
  Ntime=Shark_SFH[['lbt_mean']]$dims
  Shark_SFH$close()

  mockcone=mockcone(file_sting=file_sting)
  mocksubsets=mocksubsets(mockcone=mockcone)

  Nunique=length(unlist(mocksubsets$idlist))

  SFRbulge_d=matrix(0,Nunique,Ntime)
  SFRbulge_m=matrix(0,Nunique,Ntime)
  SFRdisk=matrix(0,Nunique,Ntime)
  Zbulge_d=matrix(0,Nunique,Ntime)
  Zbulge_m=matrix(0,Nunique,Ntime)
  Zdisk=matrix(0,Nunique,Ntime)

  cl=makeCluster(cores)
  registerDoSNOW(cl)

  Nstart=1
  iterations=dim(mocksubsets)[1]

  if(verbose){
    pb = txtProgressBar(max = iterations, style = 3)
    progress = function(n) setTxtProgressBar(pb, n)
    opts = list(progress=progress)
  }

  output=foreach(i=1:iterations, .combine='rbind', .options.snow = if(verbose){opts})%dopar%{
    if(verbose){progress(i)}
    Nid=mocksubsets[i,Nid]
    Nend=Nstart+Nid-1
    sfh_fname = paste(path_shark,mocksubsets[i,snapshot],mocksubsets[i,subvolume],'star_formation_histories.hdf5', sep='/')
    assertAccess(sfh_fname, access='r')
    Shark_SFH=h5file(sfh_fname, mode='r')
    Ndim=Shark_SFH[['disks/star_formation_rate_histories']]$dims[1]
    select=match(mocksubsets[i,unlist(idlist)], Shark_SFH[['Galaxies/id_galaxy']][])

    extract=matrix(0,Nid,Ntime*6)

    extract[,1:Ndim+Ntime*0]=t(Shark_SFH[['bulges_diskins/star_formation_rate_histories']][,select,drop=FALSE])
    extract[,1:Ndim+Ntime*1]=t(Shark_SFH[['bulges_mergers/star_formation_rate_histories']][,select,drop=FALSE])
    extract[,1:Ndim+Ntime*2]=t(Shark_SFH[['disks/star_formation_rate_histories']][,select,drop=FALSE])
    extract[,1:Ndim+Ntime*3]=t(Shark_SFH[['bulges_diskins/metallicity_histories']][,select,drop=FALSE])
    extract[,1:Ndim+Ntime*4]=t(Shark_SFH[['bulges_mergers/metallicity_histories']][,select,drop=FALSE])
    extract[,1:Ndim+Ntime*5]=t(Shark_SFH[['disks/metallicity_histories']][,select,drop=FALSE])

    Shark_SFH$close()

    extract
  }

  stopCluster(cl)

  if(verbose){
    close(pb)
  }

  if(verbose){
    message(paste('Finished getSFHfull on Stingray -',round(proc.time()[3]-timestart,3),'sec'))
  }

  outSFH=list(SFRbulge_d=output[,1:Ntime+Ntime*0], SFRbulge_m=output[,1:Ntime+Ntime*1], SFRdisk=output[,1:Ntime+Ntime*2], Zbulge_d=output[,1:Ntime+Ntime*3], Zbulge_m=output[,1:Ntime+Ntime*4], Zdisk=output[,1:Ntime+Ntime*5])
  class(outSFH)='Viperfish-SFH'
  invisible(outSFH)
}

getSFHsing=function(id_galaxy_sam, snapshot=NULL, subvolume=NULL, path_shark='.'){

  assertNumeric(id_galaxy_sam)
  assertInt(snapshot, null.ok=TRUE)
  assertInt(subvolume, null.ok=TRUE)
  assertCharacter(path_shark, max.len=1)
  assertAccess(path_shark, access='r')

  sfh_fname = paste(path_shark, snapshot, subvolume, 'star_formation_histories.hdf5',sep='/')
  assertAccess(sfh_fname, access='r')
  Shark_SFH=h5file(sfh_fname, mode='r')

  select=match(id_galaxy_sam, Shark_SFH[['galaxies/id_galaxy']][])
  keep=which(!is.na(select))
  id_galaxy_sam=id_galaxy_sam[keep]
  select=select[!is.na(select)]

  SFRbulge_d=t(Shark_SFH[['bulges_diskins/star_formation_rate_histories']][,select,drop=FALSE])
  SFRbulge_m=t(Shark_SFH[['bulges_mergers/star_formation_rate_histories']][,select,drop=FALSE])
  SFRdisk=t(Shark_SFH[['disks/star_formation_rate_histories']][,select,drop=FALSE])
  Zbulge_d=t(Shark_SFH[['bulges_diskins/metallicity_histories']][,select,drop=FALSE])
  Zbulge_m=t(Shark_SFH[['bulges_mergers/metallicity_histories']][,select,drop=FALSE])
  Zdisk=t(Shark_SFH[['disks/metallicity_histories']][,select,drop=FALSE])

  #SFRbulge=t(Shark_SFH[['Bulges/StarFormationRateHistories']][,select])
  #SFRdisk=t(Shark_SFH[['Disks/StarFormationRateHistories']][,select])
  #Zbulge=t(Shark_SFH[['Bulges/MetallicityHistories']][,select])
  #Zdisk=t(Shark_SFH[['Disks/MetallicityHistories']][,select])

  outSFH=list(id_galaxy_sam=id_galaxy_sam, keep=keep, SFRbulge_d=SFRbulge_d, SFRbulge_m=SFRbulge_m, SFRdisk=SFRdisk, Zbulge_d=Zbulge_d, Zbulge_m=Zbulge_m, Zdisk=Zdisk)
  class(outSFH)='Viperfish-SFH'
  invisible(outSFH)
}
