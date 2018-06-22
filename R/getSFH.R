getSFH=function(file_sting='mocksurvey.hdf5', path_shark='.', snapmax=199, cores=4, verbose=TRUE){

  timestart=proc.time()[3]

  if(verbose){
    message('Running getSFH on Stingray')
  }

  assertCharacter(file_sting, max.len=1)
  assertAccess(file_sting, access='r')
  assertCharacter(path_shark, max.len=1)
  assertAccess(path_shark, access='r')
  assertInt(snapmax)

  BC03lr=Dale_Msol=Nid=id_galaxy_sam=idlist=snapshot=subsnapID=subsnapshot=z=i=mocksubsets=mockcone=Ntime=time=NULL

  timestart=proc.time()[3]

  assertAccess(paste(path_shark,snapmax,'0/star_formation_histories.hdf5', sep='/'), access='r')
  SFH=h5file(paste(path_shark,snapmax,'0/star_formation_histories.hdf5', sep='/'), mode='r')
  Ntime=SFH[['LBT_mean']]$dims
  SFH$close()

  mockcone=.mockcone(file_sting=file_sting)
  mocksubsets=.mocksubsets(mockcone=mockcone)

  Nunique=length(unlist(mocksubsets$idlist))

  SFRbulge=matrix(0,Nunique,Ntime)
  SFRdisk=matrix(0,Nunique,Ntime)
  Zbulge=matrix(0,Nunique,Ntime)
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
    assertAccess(paste(path_shark,mocksubsets[i,snapshot],mocksubsets[i,subsnapshot],'star_formation_histories.hdf5', sep='/'), access='r')
    SFH=h5file(paste(path_shark,mocksubsets[i,snapshot],mocksubsets[i,subsnapshot],'star_formation_histories.hdf5', sep='/'), mode='r')
    Ndim=SFH[['Bulges/StarFormationRateHistories']]$dims[1]
    select=match(mocksubsets[i,unlist(idlist)], SFH[['Galaxies/id_galaxy']][])

    #SFRbulge[Nstart:Nend,1:Ndim]=t(SFH[['Bulges/StarFormationRateHistories']][,select])
    #SFRdisk[Nstart:Nend,1:Ndim]=t(SFH[['Disks/StarFormationRateHistories']][,select])
    #Zbulge[Nstart:Nend,1:Ndim]=t(SFH[['Bulges/MetallicityHistories']][,select])
    #Zdisk[Nstart:Nend,1:Ndim]=t(SFH[['Disks/MetallicityHistories']][,select])

    extract=matrix(0,Nid,Ntime*4)

    extract[,1:Ndim]=t(SFH[['Bulges/StarFormationRateHistories']][,select])
    extract[,1:Ndim+Ntime]=t(SFH[['Disks/StarFormationRateHistories']][,select])
    extract[,1:Ndim+Ntime*2]=t(SFH[['Bulges/MetallicityHistories']][,select])
    extract[,1:Ndim+Ntime*3]=t(SFH[['Disks/MetallicityHistories']][,select])

    SFH$close()

    extract
  }

  stopCluster(cl)

  if(verbose){
    close(pb)
  }

  if(verbose){
    message(paste('Finished getSFH on Stingray -',round(proc.time()[3]-timestart,3),'sec'))
  }

  outSFH=list(SFRbulge=output[,1:Ntime], SFRdisk=output[,1:Ntime+Ntime], Zbulge=output[,1:Ntime+Ntime*2], Zdisk=output[,1:Ntime+Ntime*3])
  class(outSFH)='Viperfish-StingSFH'
  invisible(outSFH)
}
