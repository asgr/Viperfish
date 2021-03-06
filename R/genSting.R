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

genSting=function(file_sting=NULL, path_shark='.', h='get', cores_per_subvolume=1, cores_per_snapshot=4, children_outfile="/dev/null", snapmax=199, filters=c('FUV_GALEX', 'NUV_GALEX', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1_WISE', 'W2_WISE', 'W3_WISE', 'W4_WISE', 'P100_Herschel', 'P160_Herschel', 'S250_Herschel', 'S350_Herschel', 'S500_Herschel'), tau_birth=1.5, tau_screen=0.5, pow_birth=-0.7, pow_screen=-0.7,  alpha_SF_birth=1, alpha_SF_screen=3, alpha_SF_AGN=0, read_extinct=FALSE, sparse=5, time=NULL, mockcone=NULL, intSFR=TRUE, final_file_output='Stingray-SED.csv', temp_file_output='temp.csv',  extinction_file='extinction.hdf5', reorder=TRUE, restart=FALSE, verbose=TRUE, write_final_file=FALSE){

  timestart=proc.time()[3]

  if(verbose){
    message('Running Viperfish on Stingray')
  }

  assertCharacter(path_shark, max.len=1)
  assertAccess(path_shark, access='r')
  assertInt(cores_per_subvolume)
  assertInt(cores_per_snapshot)
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

  BC03lr=Dale_NormTot=Nid=id_galaxy_sky=id_galaxy_sam=idlist=snapshot=subsnapID=subvolume=z=i=j=Ntime=zobs=NULL

  data("BC03lr", envir = environment())
  data("Dale_NormTot", envir = environment())

  if(filterlist==FALSE){
    filtout=foreach(i = filters)%do%{approxfun(getfilt(i))}
    names(filtout)=filters
  }else{
    filtout=filters
  }

  Band9_ALMA=approxfun(cbind(wave=c(4000000.0,5000000.0), response=c(1.0,1.0))) #This is what Anne Klitsch sent me- very easy to define any tophat like this
  Band8_ALMA=approxfun(cbind(wave=c(6000000.0,8000000.0), response=c(1.0,1.0))) #This is what Anne Klitsch sent me- very easy to define any tophat like this
  Band7_ALMA=approxfun(cbind(wave=c(8000000.0,11000000.0), response=c(1.0,1.0))) #This is what Anne Klitsch sent me- very easy to define any tophat like this
  Band6_ALMA=approxfun(cbind(wave=c(11000000.0,14000000.0), response=c(1.0,1.0))) #This is what Anne Klitsch sent me- very easy to define any tophat like this
  Band5_ALMA=approxfun(cbind(wave=c(14000000.0,18000000.0), response=c(1.0,1.0))) #This is what Anne Klitsch sent me- very easy to define any tophat like this
  Band4_ALMA=approxfun(cbind(wave=c(18000000.0,24000000.0), response=c(1.0,1.0))) #This is what Anne Klitsch sent me- very easy to define any tophat like this

  filtout=c(filtout, Band9_ALMA=Band9_ALMA)
  filtout=c(filtout, Band8_ALMA=Band8_ALMA)
  filtout=c(filtout, Band7_ALMA=Band7_ALMA)
  filtout=c(filtout, Band6_ALMA=Band6_ALMA)
  filtout=c(filtout, Band5_ALMA=Band5_ALMA)
  filtout=c(filtout, Band4_ALMA=Band4_ALMA)
  filters = names(filtout)

  # if(!is.null(SFHfull) & doSFHbatch==TRUE){
  #   stop('You should not provide an input to SFHfull and have doSFHbatch=TRUE set, since the requested behaviour is ambiguous!')
  # }
  #
  # if(is.null(SFHfull) & doSFHbatch==FALSE){
  #   SFHfull=getSFHfull(file_sting=file_sting, path_shark=path_shark, snapmax=snapmax, cores_per_subvolume=cores_per_subvolume, verbose=verbose)
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
  print ("will make mock cones")
  if(is.null(mockcone)){
    mockcone = .with_64bit_ints(mockcone_extract(file_sting=file_sting, reorder=reorder))
    print ("will make extinction cones")
    if(read_extinct){
      extinctioncone = .with_64bit_ints(extinction_extract(extinction_file=extinction_file, reorder=reorder))
    }
  }else{
    assertDataTable(mockcone)
    if(read_extinct){
      assertDataTable(extinctioncone)
    }
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

  print ("will prepare mockpoint and extinctionpoint")
  mockcone=as.big.matrix(mockcone)
  mockpoint=describe(mockcone)
  if(read_extinct){
    extinctioncone=as.big.matrix(extinctioncone)
    extinctionpoint=describe(extinctioncone)
  }

  if (run_foreach) {

    cl=makeCluster(cores_per_subvolume, outfile=children_outfile)
    registerDoSNOW(cl)

    if(verbose){
      pb = txtProgressBar(max = length(subsnapIDs), style = 3)
      progress = function(n) setTxtProgressBar(pb, n)
      opts = list(progress=progress)
    }

    outSED=foreach(i=1:length(subsnapIDs), .combine=.dumpout, .init=temp_file_output, .final=.dumpin, .inorder=FALSE, .options.snow = if(verbose){opts}, .packages=c('Viperfish','bigmemory','doSNOW'))%dopar%{
      cat(format(Sys.time(), "%X"), "Processing snapshot", i, "of", length(subsnapIDs), "\n")
      use=subsnapIDs[i]
      mockloop=attach.big.matrix(mockpoint)
      select=which(mockloop[,'subsnapID']==use)
      snapshot=mockloop[select[1],'snapshot']
      subvolume=mockloop[select[1],'subvolume']
      id_galaxy_sky=mockloop[select,'id_galaxy_sky']
      id_galaxy_sam=mockloop[select,'id_galaxy_sam']
      zcos=mockloop[select,'zcos']
      zobs=mockloop[select,'zobs']
      read_start = Sys.time()
      cat(format(read_start, "%X"), "Reading SFH for snapshot", i, "\n")
      SFHsing_subsnap=getSFHsing(id_galaxy_sam=id_galaxy_sam, snapshot=snapshot, subvolume=subvolume, path_shark=path_shark)
      read_end = Sys.time()
      cat(format(read_end, "%X"), "Read SFH for snapshot", i, "in", as.numeric(read_end - read_start, units="secs"),"\n")

      SFRbulge_d_subsnap=SFHsing_subsnap$SFRbulge_d
      SFRbulge_m_subsnap=SFHsing_subsnap$SFRbulge_m
      SFRdisk_subsnap=SFHsing_subsnap$SFRdisk
      Zbulge_d_subsnap=SFHsing_subsnap$Zbulge_d
      Zbulge_m_subsnap=SFHsing_subsnap$Zbulge_m
      Zdisk_subsnap=SFHsing_subsnap$Zdisk

      #define tau in extinction laws
      tau_screen_galaxies = matrix(ncol = 3, nrow = length(select)) #this is ordered as bulge (disk-ins), bulge (mergers), disks
      tau_birth_galaxies = matrix(ncol = 3, nrow = length(select)) #this is ordered as bulge (disk-ins), bulge (mergers), disks
    
      pow_screen_galaxies = matrix(ncol = 3, nrow = length(select)) #this is ordered as bulge (disk-ins), bulge (mergers), disks
      pow_birth_galaxies = matrix(ncol = 3, nrow = length(select)) #this is ordered as bulge (disk-ins), bulge (mergers), disks
    
      if(read_extinct){
         #read in disks
         extloop=attach.big.matrix(extinctionpoint)
         select_ext=which(extloop[,'subsnapID']==use)
         tau_screen_galaxies[,3] = extloop[select_ext,'tau_screen_disk']
         tau_birth_galaxies[,3] = extloop[select_ext,'tau_birth_disk']
         pow_screen_galaxies[,3] = extloop[select_ext,'pow_screen_disk']
         #read in bulges
         tau_screen_galaxies[,1] = extloop[select_ext,'tau_screen_bulge']
         tau_birth_galaxies[,1] = extloop[select_ext,'tau_birth_bulge']
         pow_screen_galaxies[,1] = extloop[select_ext,'pow_screen_bulge']
         #assume the same for bulges regardless of origin of star formation
         tau_screen_galaxies[,2] = tau_screen_galaxies[,1]
         tau_birth_galaxies[,2] = tau_birth_galaxies[,1]
         pow_screen_galaxies[,2] = pow_screen_galaxies[,1]
         #all clumps have the same power law index of the Charlot & Fall model
         pow_birth_galaxies[,] = pow_birth
      }
      else{
        tau_screen_galaxies[,] = tau_screen
        tau_birth_galaxies[,] = tau_birth
        pow_screen_galaxies[,] = pow_screen
        pow_birth_galaxies[,] = pow_birth
      }

      cl = makeCluster(cores_per_snapshot, outfile=children_outfile)
      registerDoSNOW(cl)
      cat(format(Sys.time(), "%X"), "Going into inner loop with", length(select), "elements\n")
      # Here we divide by h since the simulations output SFR in their native Msun/yr/h units.
      tempout=foreach(j=1:length(select), .combine='rbind') %dopar% {
         cat(format(Sys.time(), "%X"), "Calculating SED for galaxy", j, "of", length(select), "in snapshot", i, "\n")
         tempSED=tryCatch(c(id_galaxy_sky[SFHsing_subsnap$keep[j]], 
            unlist(genSED(
	      SFRbulge_d=SFRbulge_d_subsnap[j,]/h, 
              SFRbulge_m=SFRbulge_m_subsnap[j,]/h, 
              SFRdisk=SFRdisk_subsnap[j,]/h, 
              Zbulge_d=Zbulge_d_subsnap[j,], 
              Zbulge_m=Zbulge_m_subsnap[j,], 
              Zdisk=Zdisk_subsnap[j,], 

              redshift=zobs[j],
              time=time[1:dim(SFRdisk_subsnap)[2]]-cosdistTravelTime(zcos[j], ref='planck')*1e9, 

              tau_birth=tau_birth_galaxies[j,], 
              tau_screen=tau_screen_galaxies[j,], 
              pow_birth=pow_birth_galaxies[j,], 
              pow_screen=pow_screen_galaxies[j,], 
              alpha_SF_birth=alpha_SF_birth, 
              alpha_SF_screen=alpha_SF_screen, 
              alpha_SF_AGN=alpha_SF_AGN, 

              speclib=BC03lr, 
              filtout=filtout, 
              Dale=Dale_NormTot, 
              sparse=sparse, 
              intSFR = intSFR))), error = function(e) NULL)
        #if(class(tempSED)=="try-error"){tempSED=NA}
        tempSED
      }
      stopCluster(cl)
      as.data.table(rbind(tempout))
    }

    warnings()
    if(verbose){
      close(pb)
    }

    stopCluster(cl)
  }


  #if(file.exists(temp_file_output)){
  #  assertAccess(temp_file_output, access='w')
  #  file.remove(temp_file_output)
  #}

  print("will sort ids")
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

extinction_extract=function(extinction_file="extinction.hdf5", reorder=TRUE){
  subsnapID=snapshot=subvolume=id_galaxy_sam=NULL
  assertCharacter(extinction_file, max.len=1)
  assertAccess(extinction_file, access='r')
  extinction=h5file(extinction_file, mode='r')[['galaxies']]
  extract_col=list.datasets(extinction, recursive = TRUE)
  extinctioncone=as.data.table(lapply(extract_col, function(x) extinction[[x]][]))
  colnames(extinctioncone)=extract_col
  extinction$close()
  extinctioncone[,subsnapID:=snapshot*100+subvolume]
  if(reorder){
    extinctioncone=extinctioncone[order(subsnapID,id_galaxy_sam),]
  }
  invisible(extinctioncone)
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
  fread(temp_file_output, integer64 = "integer64")
}


# .filedump=function(temp_file_output='temp.csv', data){
#   fwrite(x=data, file='temp.csv', append=TRUE)
# }
