genSting=function(file_sting='mocksurvey.hdf5', path_shark='.', h=0.678, cores=4, snapmax=199, filters=c('FUV', 'NUV', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1', 'W2', 'W3', 'W4', 'P100', 'P160', 'S250', 'S350', 'S500'), SFHlist=NULL){

  assertCharacter(file_sting, max.len=1)
  assertCharacter(path_shark, max.len=1)
  assertScalar(h)
  assertInt(cores)
  assertInt(snapmax)
  assertCharacter(filters)
  assertList(SFHlist, null.ok=TRUE)

  timestart=proc.time()[3]

  BC03lr=Dale_Msol=Nid=id_galaxy_sam=idlist=snapshot=subsnapID=subsnapshot=z=i=mocksubsets=mockcone=Ntime=time=NULL

  data("BC03lr", envir = environment())
  data("Dale_Msol", envir = environment())

  filtout=foreach(i = filters)%do%{getfilt(i)}
  names(filtout)=filters

  if(is.null(SFHlist)){
    SFHlist=getSFH(file_sting = file_sting, path_shark = path_shark, snapmax = snapmax)
  }

  SFRbulge=SFHlist$SFRbulge
  SFRdisk=SFHlist$SFRdisk
  Zbulge=SFHlist$Zbulge
  Zdisk=SFHlist$Zdisk

  #Make mock subsets:

  mocksubsets=.mocksubsets(file_sting=file_sting)

  message(paste('Running ProSpect -',round(proc.time()[3]-timestart,3),'sec'))

  SEDlookup=data.table(id=unlist(mocksubsets$idlist), subsnapID=rep(mocksubsets$subsnapID, mocksubsets$Nid))

  registerDoParallel(cores=cores)

  outSED=foreach(i=1:dim(mockcone)[1], .combine='rbind')%dopar%{
    coluse=which(SEDlookup$id==mockcone[i,id_galaxy_sam] & SEDlookup$subsnapID==mockcone[i,subsnapID])
    offset=(snapmax-mockcone[i,snapshot])
    unlist(genSED(SFRbulge=SFRbulge[coluse,1:(Ntime-offset)]/h, SFRdisk=SFRdisk[coluse,1:(Ntime-offset)]/h, redshift=mockcone[i,z], time=time[1:(Ntime-offset)]-cosdistTravelTime(mockcone[i,z], ref='planck')*1e9, speclib=BC03lr, Zbulge=Zbulge[coluse,1:(Ntime-offset)], Zdisk=Zdisk[coluse,1:(Ntime-offset)], filtout=filtout, Dale=Dale_Msol, sparse=5, tau_birth=1.5, tau_screen=0.5))
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

  message(paste('Finished ProSpect -',round(proc.time()[3]-timestart,3),'sec'))

  return=list(outSED=outSED, SFHlist=SFHlist)
}

.mocksubsets=function(file_sting){
  assertCharacter(file_sting, max.len=1)
  mocksurvey=h5file(file_sting, mode='r')[['Galaxies']]

  extract_col=list.datasets(mocksurvey, recursive = TRUE)
  mockcone=as.data.table(lapply(extract_col, function(x) mocksurvey[[x]][]))
  colnames(mockcone)=extract_col
  mocksurvey$close()

  mockcone[,subsnapID:=snapshot*100+subsnapshot]
  mocksubsets=mockcone[,list(idlist=list(unique(id_galaxy_sam))),by=subsnapID]
  mocksubsets[,Nid:=length(unlist(idlist)),by=subsnapID]
  mocksubsets[,snapshot:=floor(subsnapID/100)]
  mocksubsets[,subsnapshot:=subsnapID%%100]
  return=mocksubsets
}
