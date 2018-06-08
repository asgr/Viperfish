genSting=function(mocksurvey='mocksurvey.hdf5', path_shark='.', h=0.678, cores=4, select='all', snapmax=199, filters=c('FUV', 'NUV', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1', 'W2', 'W3', 'W4', 'P100', 'P160', 'S250', 'S350', 'S500')){

  timestart=proc.time()[3]

  data("BC03lr")
  data("Dale_Msol")

  filtout=foreach(i = filters)%do%{getfilt(i)}
  names(filtout)=filters

  SFH=h5file(paste(path_shark,snapmax,'0/star_formation_histories.hdf5', sep='/'), mode='r')
  time=SFH[['age_mean']][]*1e9
  Ntime=SFH[['age_mean']]$dims
  SFH$close()

  mocksurvey=h5file(mocksurvey, mode='r')[['Galaxies']]

  extract_col=list.datasets(mocksurvey, recursive = TRUE)
  mockcone=as.data.table(lapply(extract_col, function(x) mocksurvey[[x]][]))
  colnames(mockcone)=extract_col
  mocksurvey$close()

  mockcone[,subsnapID:=snapshot*100+subsnapshot]
  mocksubsets=mockcone[,list(idlist=list(unique(id_galaxy_sam))),by=subsnapID]
  mocksubsets[,Nid:=length(unlist(idlist)),by=subsnapID]
  mocksubsets[,snapshot:=floor(subsnapID/100)]
  mocksubsets[,subsnapshot:=subsnapID%%100]

  Nunique=length(unlist(mocksubsets$idlist))

  SFRbulge=matrix(0,Nunique,Ntime)
  SFRdisk=matrix(0,Nunique,Ntime)
  Zbulge=matrix(0,Nunique,Ntime)
  Zdisk=matrix(0,Nunique,Ntime)

  message(paste('Extracting Light Cone SFH -',round(proc.time()[3]-timestart,3),'sec'))

  Nstart=1
  for(i in 1:dim(mocksubsets)[1]){
    if(i%%100==0){message(i,' of ',dim(mocksubsets)[1])}
    Nend=Nstart+mocksubsets[i,Nid]-1
    SFH=h5file(paste0(path_shark,'/',mocksubsets[i,snapshot],'/',mocksubsets[i,subsnapshot],'/star_formation_histories.hdf5'), mode='r')
    Ndim=SFH[['Bulges/StarFormationRateHistories']]$dims[1]
    select=which(SFH[['Galaxies/id_galaxy']][] %in% mocksubsets[i,unlist(idlist)])
    SFRbulge[Nstart:Nend,1:Ndim]=SFH[['Bulges/StarFormationRateHistories']][,select]
    SFRdisk[Nstart:Nend,1:Ndim]=SFH[['Disks/StarFormationRateHistories']][,select]
    Zbulge[Nstart:Nend,1:Ndim]=SFH[['Bulges/MetallicityHistories']][,select]
    Zdisk[Nstart:Nend,1:Ndim]=SFH[['Disks/MetallicityHistories']][,select]
    SFH$close()
    Nstart=Nend+1
  }

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

  return=outSED
}
