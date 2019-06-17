.with_64bit_ints = function(x)
{
    # When reading 64 bit integers let's make sure hdf5r reads them as such;
    # otherwise the default behavior is convert to int32 or double if they don't loose precision
    old_option = getOption('hdf5r.h5tor_default')
    options(hdf5r.h5tor_default = h5const$H5TOR_CONV_NONE)
    result = x
    options(hdf5r.h5tor_default = old_option)
    invisible(result)
}

genSting=function(file_sting=NULL, path_shark='.', h='get', cores=4, snapmax=199, filters=c('FUV_GALEX', 'NUV_GALEX', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1_WISE', 'W2_WISE', 'W3_WISE', 'W4_WISE', 'P100_Herschel', 'P160_Herschel', 'S250_Herschel', 'S350_Herschel', 'S500_Herschel'), tau_birth=1.5, tau_screen=0.5, sparse=5, time=NULL, mockcone=NULL, intSFR=TRUE, final_file_output='Stingray-SED.csv', temp_file_output='temp.csv', reorder=TRUE, restart=FALSE, verbose=TRUE, write_final_file=FALSE){

  timestart=proc.time()[3]

  if(verbose){
    message('Running Viperfish on Stingray')
  }

  assertCharacter(path_shark, max.len=1)
  assertAccess(path_shark, access='r')
  assertInt(cores)
  assertInt(snapmax)
  if(is.list(filters)){
    filterlist=TRUE
  }else{
    filterlist=FALSE
    assertCharacter(filters)
  }
  assertScalar(tau_birth)
  assertScalar(tau_screen)
  assertInt(sparse)
  assertNumeric(time, null.ok=TRUE)
  assertDataTable(mockcone, null.ok=TRUE)
  assertFlag(intSFR)
  assertFlag(verbose)

  assertAccess(paste(path_shark,snapmax,'0/star_formation_histories.hdf5', sep='/'), access='r')

  if(h=='get'){
    h=h5file(paste(path_shark,snapmax,'0/star_formation_histories.hdf5', sep='/'), mode='r')[['cosmology/h']][]
  }

  if(! is.null(file_sting)){
    assertCharacter(file_sting, max.len=1, null.ok=TRUE)
    assertAccess(file_sting, access='r')
    Sting_date=h5file(file_sting, mode='r')[['run_info/shark_timestamp']][]
    Shark_date=h5file(paste(path_shark,snapmax,'0/star_formation_histories.hdf5', sep='/'), mode='r')[['run_info/timestamp']][]
    check = (Shark_date==Sting_date)

    if(check==FALSE){
      stop(paste('Date stamps do not match! Shark',Shark_date,'compared to Sting',Sting_date))
    }

    rm(Sting_date)
    rm(Shark_date)

    #Extract the original id_galaxy_sky for final re-ordering.
    Sting_id_galaxy_sky = .with_64bit_ints(h5file(file_sting, mode='r')[['galaxies/id_galaxy_sky']][])
  }else{
    if(! is.null(mockcone)){
      Sting_id_galaxy_sky=mockcone$id_galaxy_sky
    }else{
      stop('Need either file_sting or mockcone to be input!')
    }
  }

  message(paste0('#Galaxy IDs from mock catalogue: ', length(Sting_id_galaxy_sky),
                 ', type/class: ', typeof(Sting_id_galaxy_sky), '/', class(Sting_id_galaxy_sky)))
  assertScalar(h)

  BC03lr=Dale_Msol=Nid=id_galaxy_sky=id_galaxy_sam=idlist=snapshot=subsnapID=subvolume=z=i=j=Ntime=zobs=NULL

  data("BC03lr", envir = environment())
  data("Dale_Msol", envir = environment())

  if(filterlist==FALSE){
    filtout=foreach(i = filters)%do%{getfilt(i)}
    names(filtout)=filters
  }else{
    filtout=filters
  }

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

  if(file.exists(temp_file_output) & restart==FALSE){
    stop(paste(temp_file_output,'already exists! Use another temporary temp_file_output file name, or set restart=TRUE to continue from where genSting finished.'))
    #assertAccess(temp_file_output, access='w')
    #file.remove(temp_file_output)
  }else if(restart==FALSE){
    file.create(temp_file_output)
  }
  assertAccess(temp_file_output, access='w')

  #Make mock subsets:

  if(is.null(mockcone)){
    mockcone = .with_64bit_ints(mockcone_extract(file_sting=file_sting, reorder=reorder))
  }else{
    assertDataTable(mockcone)
  }
  #mocksubsets=mocksubsets(mockcone=mockcone)

  if(is.null(time)){
    assertAccess(paste(path_shark,snapmax,'0/star_formation_histories.hdf5', sep='/'), access='r')
    Shark_SFH=h5file(paste(path_shark,snapmax,'0/star_formation_histories.hdf5', sep='/'), mode='r')
    time=Shark_SFH[['lbt_mean']][]*1e9
    Shark_SFH$close()
  }
  Ntime=length(time)

  #SEDlookup=data.table(id=unlist(mocksubsets$idlist), subsnapID=rep(mocksubsets$subsnapID, mocksubsets$Nid))

  #iterations=dim(mockcone)[1]
  if(restart){
    outSED=fread(temp_file_output)
    subsnapIDs=base::unique(mockcone[!id_galaxy_sky %in% outSED$V1,subsnapID])
    run_foreach = length(subsnapIDs) > 0
  }else{
    subsnapIDs=base::unique(mockcone$subsnapID)
    run_foreach = TRUE
  }

  mockcone=as.big.matrix(mockcone)
  mockpoint=describe(mockcone)

  if (run_foreach) {

    cl=makeCluster(cores)
    registerDoSNOW(cl)

    if(verbose){
      pb = txtProgressBar(max = length(subsnapIDs), style = 3)
      progress = function(n) setTxtProgressBar(pb, n)
      opts = list(progress=progress)
    }

    outSED=foreach(i=1:length(subsnapIDs), .combine=.dumpout, .init=temp_file_output, .final=.dumpin, .inorder=FALSE, .options.snow = if(verbose){opts}, .packages=c('Viperfish','bigmemory'))%dopar%{
      use=subsnapIDs[i]
      mockloop=attach.big.matrix(mockpoint)
      select=which(mockloop[,'subsnapID']==use)
      snapshot=mockloop[select[1],'snapshot']
      subvolume=mockloop[select[1],'subvolume']
      id_galaxy_sky=mockloop[select,'id_galaxy_sky']
      id_galaxy_sam=mockloop[select,'id_galaxy_sam']
      zcos=mockloop[select,'zcos']
      zobs=mockloop[select,'zobs']
      SFHsing_subsnap=getSFHsing(id_galaxy_sam=id_galaxy_sam, snapshot=snapshot, subvolume=subvolume, path_shark=path_shark)

      SFRbulge_d_subsnap=SFHsing_subsnap$SFRbulge_d
      SFRbulge_m_subsnap=SFHsing_subsnap$SFRbulge_m
      SFRdisk_subsnap=SFHsing_subsnap$SFRdisk
      Zbulge_d_subsnap=SFHsing_subsnap$Zbulge_d
      Zbulge_m_subsnap=SFHsing_subsnap$Zbulge_m
      Zdisk_subsnap=SFHsing_subsnap$Zdisk

      # Here we divide by h since the simulations output SFR in their native Msun/yr/h units.
      tempout=foreach(j=1:length(select), .combine='rbind')%do%{
        tempSED=tryCatch(c(id_galaxy_sky[SFHsing_subsnap$keep[j]], unlist(genSED(SFRbulge_d=SFRbulge_d_subsnap[j,]/h, SFRbulge_m=SFRbulge_m_subsnap[j,]/h, SFRdisk=SFRdisk_subsnap[j,]/h, redshift=zobs[j], time=time[1:dim(SFRdisk_subsnap)[2]]-cosdistTravelTime(zcos[j], ref='planck')*1e9, speclib=BC03lr, Zbulge_d=Zbulge_d_subsnap[j,], Zbulge_m=Zbulge_m_subsnap[j,], Zdisk=Zdisk_subsnap[j,], filtout=filtout, Dale=Dale_Msol, tau_birth=tau_birth, tau_screen=tau_screen, sparse=sparse, intSFR = intSFR))), error = function(e) NULL)
        tempSED
      }
      as.data.table(rbind(tempout))
    }

    if(verbose){
      close(pb)
    }

    stopCluster(cl)
  }


  #if(file.exists(temp_file_output)){
  #  assertAccess(temp_file_output, access='w')
  #  file.remove(temp_file_output)
  #}

  outSED=unique(outSED, by=1)
  outSED=as.data.frame(outSED)
  outSED=outSED[match(Sting_id_galaxy_sky, outSED[,1]),]

  if (write_final_file) {
    write.SED(outSED, filters, dirname(file_sting), final_file_output)
  }

  if(verbose){
    message(paste('Finished Viperfish on Stingray -',round(proc.time()[3]-timestart,3),'sec'))
  }

  class(outSED)=c(class(outSED),'Viperfish-Shark')
  invisible(outSED)
}

mockcone_extract=function(file_sting="mocksurvey.hdf5", reorder=TRUE){
  subsnapID=snapshot=subvolume=id_galaxy_sam=NULL
  assertCharacter(file_sting, max.len=1)
  assertAccess(file_sting, access='r')
  mocksurvey=h5file(file_sting, mode='r')[['galaxies']]
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

.dumpout = function(temp_file_output='temp.csv', ...){
  # Writing like this prevents race conditions between threads. All results from a batch of N-cores are written to the file sequentially, so no chance of over-write etc.
  for(r in list(...))
    fwrite(x=r, file=temp_file_output, append=TRUE)
  temp_file_output
}

.dumpin = function(temp_file_output='temp.csv'){
  fread(temp_file_output)
}


# .filedump=function(temp_file_output='temp.csv', data){
#   fwrite(x=data, file='temp.csv', append=TRUE)
# }
