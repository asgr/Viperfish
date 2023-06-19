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

genSting = function(file_sting=NULL, path_shark='.', path_out='.', h='get', cores_per_subvolume=1,
                  cores_per_snapshot=4, children_outfile="/dev/null",
                  filters=c('FUV_GALEX', 'NUV_GALEX', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS',
                  'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1_WISE', 'W2_WISE',
                  'W3_WISE', 'W4_WISE', 'P100_Herschel', 'P160_Herschel', 'S250_Herschel',
                  'S350_Herschel', 'S500_Herschel'), tau_birth=1.5, tau_screen=0.5, pow_birth=-0.7,
                  pow_screen=-0.7,  alpha_SF_birth=1, alpha_SF_screen=3, alpha_SF_AGN=0, emission=FALSE,
                  stellarpop='BC03lr', read_extinct=FALSE, sparse=5, time=NULL, mockcone=NULL,
                  intSFR=TRUE, final_file_output='Stingray-SED.csv', temp_file_output='temp.csv',
                  extinction_file='extinction.hdf5', addradio_SF=FALSE, waveout=seq(2,30,by=0.01),
                  ff_frac_SF=0.1, ff_power_SF=-0.1, sy_power_SF=-0.8, reorder=TRUE, restart=FALSE,
                  verbose=TRUE, write_final_file=FALSE, mode='photom', spec_range=c(3600,9600), spec_bin=1){

  timestart=proc.time()[3]

  if(verbose){
    message('Running Viperfish on Stingray')
  }

  assertCharacter(path_shark, max.len=1)
  assertAccess(path_shark, access='r')
  assertInt(cores_per_subvolume)
  assertInt(cores_per_snapshot)
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

  suppressWarnings({
    snapmax = max(as.integer(list.files(path_shark)), na.rm=TRUE)
    assertAccess(paste(path_shark,snapmax, sep='/'), access='r')
    submin = min(as.integer(list.files(paste(path_shark,snapmax, sep='/'))), na.rm=TRUE)
    assertAccess(paste(path_shark,snapmax,submin,'star_formation_histories.hdf5', sep='/'), access='r')
  })

  if(h=='get'){
    h = h5file(paste(path_shark,snapmax,submin,'star_formation_histories.hdf5', sep='/'), mode='r')[['cosmology/h']][]
  }

  if(! is.null(file_sting)){
    assertCharacter(file_sting, max.len=1, null.ok=TRUE)
    assertAccess(file_sting, access='r')
    Sting_date=h5file(file_sting, mode='r')[['run_info/shark_timestamp']][]
    Shark_date=h5file(paste(path_shark,snapmax,submin,'star_formation_histories.hdf5', sep='/'), mode='r')[['run_info/timestamp']][]
    check = (Shark_date==Sting_date)

    if(check==FALSE){
      message('Date stamps do not match! Shark: ',Shark_date,' compared to Sting: ',Sting_date)
    }

    rm(Sting_date)
    rm(Shark_date)

    #Extract the original id_galaxy_sky for final re-ordering.
    Sting_id_galaxy_sky = .with_64bit_ints(h5file(file_sting, mode='r')[['galaxies/id_galaxy_sky']][])
  }else{
    if(! is.null(mockcone)){
      Sting_id_galaxy_sky = mockcone$id_galaxy_sky
    }else{
      stop('Need either file_sting or mockcone to be input!')
    }
  }

  message(paste0('#Galaxy IDs from mock catalogue: ', length(Sting_id_galaxy_sky),
                 ', type/class: ', typeof(Sting_id_galaxy_sky), '/', class(Sting_id_galaxy_sky)))
  assertScalar(h)

  speclib=Dale_NormTot=Nid=id_galaxy_sky=id_galaxy_sam=idlist=snapshot=subsnapID=subvolume=z=i=j=Ntime=zobs=NULL

  if (stellarpop == 'BC03lr') {
    if (is.null(speclib)) {
      BC03lr = NULL
      data('BC03lr', envir = environment())
      speclib = BC03lr
    }
  }else if (stellarpop == 'BC03hr') {
    if (is.null(speclib)) {
      BC03hr = NULL
      data('BC03hr', envir = environment())
      speclib = BC03hr
    }
  }else if (stellarpop == 'EMILES') {
    if (is.null(speclib)) {
      EMILES = NULL
      data('EMILES', envir = environment())
      speclib = EMILES
    }
  }else if (stellarpop == 'BPASS') {
    if (is.null(speclib)) {
      BPASS = NULL
      data('BPASS', envir = environment())
      speclib = BPASS
    }
  }

  data("Dale_NormTot", envir = environment())

  if(filterlist==FALSE){
    filtout=foreach(i = filters)%do%{approxfun(getfilt(i))}
    names(filtout)=filters
  }else{
    filtout=filters
  }


  Band_ionising_photons=approxfun(cbind(wave=c(0,912), response=c(0.01,0.01))) #Band to compute the luminosity left to the 912Angs wavelength.
  Band9_ALMA=approxfun(cbind(wave=c(4000000.0,5000000.0), response=c(1.0,1.0))) #This is what Anne Klitsch sent me- very easy to define any tophat like this
  Band8_ALMA=approxfun(cbind(wave=c(6000000.0,8000000.0), response=c(1.0,1.0))) #This is what Anne Klitsch sent me- very easy to define any tophat like this
  Band7_ALMA=approxfun(cbind(wave=c(8000000.0,11000000.0), response=c(1.0,1.0))) #This is what Anne Klitsch sent me- very easy to define any tophat like this
  Band6_ALMA=approxfun(cbind(wave=c(11000000.0,14000000.0), response=c(1.0,1.0))) #This is what Anne Klitsch sent me- very easy to define any tophat like this
  Band5_ALMA=approxfun(cbind(wave=c(14000000.0,18000000.0), response=c(1.0,1.0))) #This is what Anne Klitsch sent me- very easy to define any tophat like this
  Band4_ALMA=approxfun(cbind(wave=c(18000000.0,24000000.0), response=c(1.0,1.0))) #This is what Anne Klitsch sent me- very easy to define any tophat like this
  Band3_ALMA=approxfun(cbind(wave=c(26000000.0,36000000.0), response=c(1.0,1.0))) #This is what Anne Klitsch sent me- very easy to define any tophat like this

  BandX_VLA=approxfun(cbind(wave=c(249827050.0,374740570.0), response=c(1000.0,1000.0)))
  BandC_VLA=approxfun(cbind(wave=c(374740570.0,749481150.0), response=c(1000.0,1000.0)))
  BandS_VLA=approxfun(cbind(wave=c(749481150.0,1498962290.0), response=c(1000.0,1000.0)))
  BandL_VLA=approxfun(cbind(wave=c(1498962290.0,2997924580.0), response=c(1000.0,1000.0)))
  Band_610MHz=approxfun(cbind(wave=c(4542309970.0,5353436750.0), response=c(10000.0,10000.0)))
  Band_325MHz=approxfun(cbind(wave=c(7994465550.0,10901543930.0), response=c(10000.0,10000.0)))
  Band_150MHz=approxfun(cbind(wave=c(14989622900.0,29979245800.0), response=c(10000.0,10000.0)))

  filtout=c(filtout, Band_ionising_photons=Band_ionising_photons)
  filtout=c(filtout, Band9_ALMA=Band9_ALMA)
  filtout=c(filtout, Band8_ALMA=Band8_ALMA)
  filtout=c(filtout, Band7_ALMA=Band7_ALMA)
  filtout=c(filtout, Band6_ALMA=Band6_ALMA)
  filtout=c(filtout, Band5_ALMA=Band5_ALMA)
  filtout=c(filtout, Band4_ALMA=Band4_ALMA)
  filtout=c(filtout, Band3_ALMA=Band3_ALMA)
  filtout=c(filtout, BandX_VLA=BandX_VLA)
  filtout=c(filtout, BandC_VLA=BandC_VLA)
  filtout=c(filtout, BandS_VLA=BandS_VLA)
  filtout=c(filtout, BandL_VLA=BandL_VLA)
  filtout=c(filtout, Band_610MHz=Band_610MHz)
  filtout=c(filtout, Band_325MHz=Band_325MHz)
  filtout=c(filtout, Band_150MHz=Band_150MHz)


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
    assertAccess(paste(path_shark,snapmax,submin,'star_formation_histories.hdf5', sep='/'), access='r')
    Shark_SFH=h5file(paste(path_shark,snapmax,submin,'star_formation_histories.hdf5', sep='/'), mode='r')
    time=Shark_SFH[['lbt_mean']][]*1e9
    Shark_SFH$close()
  }
  Ntime=length(time)

  #SEDlookup=data.table(id=unlist(mocksubsets$idlist), subsnapID=rep(mocksubsets$subsnapID, mocksubsets$Nid))

  #iterations=dim(mockcone)[1]
  if(restart){
    outSED=fread(temp_file_output)
    subsnapIDs=base::unique(mockcone[!id_galaxy_sky %in% unique(outSED$V1),subsnapID])
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

  if(run_foreach){

    cl_sub = makeCluster(cores_per_subvolume, outfile=children_outfile)
    registerDoSNOW(cl_sub)

    if(verbose){
      pb = txtProgressBar(max = length(subsnapIDs), style = 3)
      progress = function(n) setTxtProgressBar(pb, n)
      opts = list(progress=progress)
    }

    outSED = foreach(i=1:length(subsnapIDs), .combine=.dumpout, .init=temp_file_output, .final=.dumpin,
                   .inorder=FALSE, .options.snow = if(verbose){opts}, .packages=c('Viperfish','bigmemory','doSNOW'))%dopar%{
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

      cl_snap = makeCluster(cores_per_snapshot, outfile=children_outfile)
      registerDoSNOW(cl_snap)
      cat(format(Sys.time(), "%X"), "Going into inner loop with", length(select), "elements\n")
      # Here we divide by h since the simulations output SFR in their native Msun/yr/h units.
      tempout = foreach(j=1:length(select), .combine='rbind') %dopar% {
         cat(format(Sys.time(), "%X"), "Calculating SED for galaxy", j, "of", length(select), "in snapshot", i, "\n")
        if(mode == 'photom'){
           tempSED = tryCatch(c(id_galaxy_sky[SFHsing_subsnap$keep[j]], 
              unlist(genSED(
  	            SFRbulge_d = SFRbulge_d_subsnap[j,]/h, 
                SFRbulge_m = SFRbulge_m_subsnap[j,]/h, 
                SFRdisk = SFRdisk_subsnap[j,]/h, 
                Zbulge_d = Zbulge_d_subsnap[j,], 
                Zbulge_m = Zbulge_m_subsnap[j,], 
                Zdisk = Zdisk_subsnap[j,], 
  
                redshift = zobs[j],
                time = time[1:dim(SFRdisk_subsnap)[2]]-cosdistTravelTime(zcos[j], ref='planck')*1e9, 
  
                tau_birth = tau_birth_galaxies[j,], 
                tau_screen = tau_screen_galaxies[j,], 
                pow_birth = pow_birth_galaxies[j,], 
                pow_screen = pow_screen_galaxies[j,], 
                alpha_SF_birth = alpha_SF_birth, 
                alpha_SF_screen = alpha_SF_screen, 
                alpha_SF_AGN = alpha_SF_AGN, 
  
  	            emission = emission,
  	            
                speclib = speclib, 
                filtout = filtout, 
                Dale = Dale_NormTot, 
                sparse = sparse, 
                intSFR = intSFR,
  
                addradio_SF = addradio_SF,
                waveout = waveout,
                ff_frac_SF = ff_frac_SF,
                ff_power_SF = ff_power_SF,
                sy_power_SF = sy_power_SF,
  	      
  	            mode = mode)
  	      )), error = function(e) NULL)
          #if(class(tempSED)=="try-error"){tempSED=NA}
          return(tempSED)
        }else if(mode == 'spectrum' | mode == 'spectra' | mode == 'spec' | mode == 'spectral'){
          tempSED = tryCatch(genSED(
            SFRbulge_d = SFRbulge_d_subsnap[j,]/h, 
            SFRbulge_m = SFRbulge_m_subsnap[j,]/h, 
            SFRdisk = SFRdisk_subsnap[j,]/h, 
            Zbulge_d = Zbulge_d_subsnap[j,], 
            Zbulge_m = Zbulge_m_subsnap[j,], 
            Zdisk = Zdisk_subsnap[j,], 
            
            redshift = zobs[j],
            time = time[1:dim(SFRdisk_subsnap)[2]]-cosdistTravelTime(zcos[j], ref='planck')*1e9, 
            
            tau_birth = tau_birth_galaxies[j,], 
            tau_screen = tau_screen_galaxies[j,], 
            pow_birth = pow_birth_galaxies[j,], 
            pow_screen = pow_screen_galaxies[j,], 
            alpha_SF_birth = alpha_SF_birth, 
            alpha_SF_screen = alpha_SF_screen, 
            alpha_SF_AGN = alpha_SF_AGN, 
            
            emission = emission,
            
            speclib = speclib, 
            filtout = filtout, 
            Dale = Dale_NormTot, 
            sparse = sparse, 
            intSFR = intSFR,
            
            addradio_SF = addradio_SF,
            waveout = NULL,
            ff_frac_SF = ff_frac_SF,
            ff_power_SF = ff_power_SF,
            sy_power_SF = sy_power_SF,
            
            mode = mode,
            spec_range = spec_range + c(-100,100)), error = function(e) NULL)
          return(cbind(id_galaxy_sky[SFHsing_subsnap$keep[j]],tempSED))
        }
      }
      stopCluster(cl_snap)
      if(mode == 'photom'){
        return(as.data.table(rbind(tempout)))
      }else if(mode == 'spectrum' | mode == 'spectra' | mode == 'spec' | mode == 'spectral'){
        return(tempout)
      }
    }

    warnings()
    if(verbose){
      close(pb)
    }

    stopCluster(cl_sub)
  }


  #if(file.exists(temp_file_output)){
  #  assertAccess(temp_file_output, access='w')
  #  file.remove(temp_file_output)
  #}

  if(mode == 'photom'){
    print("will sort ids")
    outSED = unique(outSED, by=1)
    outSED = as.data.frame(outSED)
    outSED = outSED[match(Sting_id_galaxy_sky, outSED[,1]),]
  
    if (write_final_file){
      write.SED(outSED, filters, path_out, final_file_output)
    }
  
    if(verbose){
      message(paste('Finished Viperfish on Stingray -',round(proc.time()[3]-timestart,3),'sec'))
    }
  
    class(outSED)=c(class(outSED),'Viperfish-Sting')
    return(invisible(outSED))
  }else if(mode == 'spectrum' | mode == 'spectra' | mode == 'spec'){
    #waveN = dim(outSED)[1]/length(Sting_id_galaxy_sky)
    wavegrid = seq(spec_range[1], spec_range[2], by=spec_bin)
    outSED_split = foreach(i = unique(outSED$V1), .combine='rbind')%do%{
      temp_spec = outSED[outSED$V1 == i,list(V2,V3)]
      temp_spec = temp_spec[!duplicated(temp_spec$V2),]
      setkey(temp_spec, V2)
      return(specReBin(temp_spec, wavegrid=wavegrid)$flux)
    }
    outSED_split = outSED_split[match(Sting_id_galaxy_sky, unique(outSED$V1)),]
    outSED_split = data.frame(rbind(wavegrid, outSED_split))
    outSED_split = cbind(data.frame(c(as.integer64(0), Sting_id_galaxy_sky)), outSED_split)
    
    fwrite(outSED_split, paste(path_out, final_file_output, sep='/'), row.names = FALSE, col.names = FALSE)
    
    return(invisible(outSED_split))
  }
}

mockcone_extract = function(file_sting="mocksurvey.hdf5", reorder=TRUE){
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

extinction_extract = function(extinction_file="extinction.hdf5", reorder=TRUE){
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


mocksubsets = function(mockcone){
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
