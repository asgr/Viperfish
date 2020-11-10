genShark = function(path_shark='.', snapshot=NULL, subvolume=NULL, redshift="get",
                  h='get', cores=4, id_galaxy_sam='all', filters=c('FUV_GALEX', 'NUV_GALEX',
                  'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA', 'Y_VISTA', 'J_VISTA',
                  'H_VISTA', 'K_VISTA', 'W1_WISE', 'W2_WISE', 'W3_WISE', 'W4_WISE',
                  'P100_Herschel', 'P160_Herschel', 'S250_Herschel', 'S350_Herschel',
                  'S500_Herschel'), tau_birth=1, tau_screen=0.3, tau_AGN=1, pow_birth=-0.7,
                  pow_screen=-0.7, pow_AGN=-0.7, alpha_SF_birth=1, alpha_SF_screen=3,
                  alpha_SF_AGN=0, emission=FALSE, IGMabsorb=FALSE, read_extinct=FALSE, sparse=5, final_file_output='Shark-SED.csv',
                  extinction_file='extinction.hdf5', intSFR=TRUE, verbose=TRUE, write_final_file=FALSE){

  timestart = proc.time()[3]

  if(verbose){
    message('Running Viperfish on Shark')
  }

  assertCharacter(path_shark, max.len=1)
  assertAccess(path_shark, access='r')
  assertInt(snapshot, null.ok=TRUE)
  assertInt(subvolume, null.ok=TRUE)
  assertInt(cores)
  if(is.list(filters)){
    filterlist = TRUE
  }else{
    filterlist = FALSE
    assertCharacter(filters)
  }
  assertScalar(tau_birth)
  assertScalar(tau_screen)

  assertInt(sparse)
  assertFlag(intSFR)
  assertFlag(verbose)

  BC03lr = Dale_NormTot = AGN_UnOb_Sparse = LKL10_NormAll = SFH = i = subsnapID = NULL
  Nid = idlist = subsnapID = NULL

  data("BC03lr", envir = environment())
  data("Dale_NormTot", envir = environment())
  data("AGN_UnOb_Sparse", envir = environment())
  if(emission){
    data("LKL10_NormAll", envir = environment())
  }else{
    LKL10_NormAll = NULL
  }

  if(filterlist==FALSE){
    filtout=foreach(i = filters)%do%{approxfun(getfilt(i))}
    names(filtout) = filters
  }else{
    filtout = filters
  }

  # FUV_Nathan=approxfun(cbind(wave=c(1450,1550), response=c(0.01,0.01))) #This is what Claudia sent me- very easy to define any tophat like this
  # Band9_ALMA=approxfun(cbind(wave=c(4000000.0,5000000.0), response=c(1.0,1.0))) #This is what Anne Klitsch sent me- very easy to define any tophat like this
  # Band8_ALMA=approxfun(cbind(wave=c(6000000.0,8000000.0), response=c(1.0,1.0))) #This is what Anne Klitsch sent me- very easy to define any tophat like this
  # Band7_ALMA=approxfun(cbind(wave=c(8000000.0,11000000.0), response=c(1.0,1.0))) #This is what Anne Klitsch sent me- very easy to define any tophat like this
  # Band6_ALMA=approxfun(cbind(wave=c(11000000.0,14000000.0), response=c(1.0,1.0))) #This is what Anne Klitsch sent me- very easy to define any tophat like this
  # Band4_ALMA=approxfun(cbind(wave=c(18000000.0,24000000.0), response=c(1.0,1.0))) #This is what Anne Klitsch sent me- very easy to define any tophat like this
  # 
  # filtout=c(filtout, FUV_Nathan=FUV_Nathan)
  # filtout=c(filtout, Band9_ALMA=Band9_ALMA)
  # filtout=c(filtout, Band8_ALMA=Band8_ALMA)
  # filtout=c(filtout, Band7_ALMA=Band7_ALMA)
  # filtout=c(filtout, Band6_ALMA=Band6_ALMA)
  # filtout=c(filtout, Band4_ALMA=Band4_ALMA)
  # filters = names(filtout)

  sfh_fname = paste(path_shark, snapshot, subvolume, 'star_formation_histories.hdf5',sep='/')
  assertAccess(sfh_fname, access='r')
  Shark_SFH=h5file(sfh_fname, mode='r')

  if(read_extinct){
      extinct_fname = paste(path_shark, snapshot, subvolume, extinction_file,sep='/')
      assertAccess(extinct_fname, access='r')
      Shark_Extinct = h5file(extinct_fname, mode='r')
  }

  # Read things in from the Shark HDF5 file:
  time = Shark_SFH[['lbt_mean']][]*1e9

  if(redshift[1]=='get'){
    redshift = Shark_SFH[['run_info/redshift']][]
    if(redshift==0){redshift = 1e-10}
  }

  if(h=='get'){
    h = Shark_SFH[['cosmology/h']][]
  }

  assertNumeric(redshift)
  assertScalar(h)

  if(id_galaxy_sam[1]=='all'){
    select = 1:Shark_SFH[['galaxies/id_galaxy']]$dims
  }else{
    select = match(id_galaxy_sam, Shark_SFH[['galaxies/id_galaxy']][])
  }

  SFRbulge_d = Shark_SFH[['bulges_diskins/star_formation_rate_histories']][,select,drop=FALSE]
  SFRbulge_m = Shark_SFH[['bulges_mergers/star_formation_rate_histories']][,select,drop=FALSE]
  SFRdisk = Shark_SFH[['disks/star_formation_rate_histories']][,select,drop=FALSE]
  Zbulge_d = Shark_SFH[['bulges_diskins/metallicity_histories']][,select,drop=FALSE]
  Zbulge_m = Shark_SFH[['bulges_mergers/metallicity_histories']][,select,drop=FALSE]
  Zdisk = Shark_SFH[['disks/metallicity_histories']][,select,drop=FALSE]

  #define tau in extinction laws
  tau_screen_galaxies = matrix(ncol = 3, nrow = length(select)) #this is ordered as bulge (disk-ins), bulge (mergers), disks
  tau_birth_galaxies = matrix(ncol = 3, nrow = length(select)) #this is ordered as bulge (disk-ins), bulge (mergers), disks

  pow_screen_galaxies = matrix(ncol = 3, nrow = length(select)) #this is ordered as bulge (disk-ins), bulge (mergers), disks
  pow_birth_galaxies = matrix(ncol = 3, nrow = length(select)) #this is ordered as bulge (disk-ins), bulge (mergers), disks

  if(read_extinct){
    #read in disks
    tau_screen_galaxies[,3] = Shark_Extinct[['galaxies/tau_screen_disk']][select]
    tau_birth_galaxies[,3] = Shark_Extinct[['galaxies/tau_birth_disk']][select]
    pow_screen_galaxies[,3] = Shark_Extinct[['galaxies/pow_screen_disk']][select]

    tau_screen_galaxies[,1] = Shark_Extinct[['galaxies/tau_screen_bulge']][select]
    tau_birth_galaxies[,1] = Shark_Extinct[['galaxies/tau_birth_bulge']][select]
    pow_screen_galaxies[,1] = Shark_Extinct[['galaxies/pow_screen_bulge']][select]

    #assume the same for bulges regardless of origin of star formation
    tau_screen_galaxies[,2] = tau_screen_galaxies[,1]
    tau_birth_galaxies[,2] = tau_birth_galaxies[,1]
    pow_screen_galaxies[,2] = pow_screen_galaxies[,1]

    #all clumps have the same power law index of the Charlot & Fall model
    pow_birth_galaxies[,] = pow_birth
  }else{
    tau_screen_galaxies[,] = tau_screen
    tau_birth_galaxies[,] = tau_birth
    pow_screen_galaxies[,] = pow_screen
    pow_birth_galaxies[,] = pow_birth
  }

  if(length(redshift)==1){
    redshift = rep(redshift, length(select))
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

  # Here we divide by h since the simulations output SFR in their native Msun/yr/h units.

  outSED=foreach(i=1:iterations, .combine='rbind', .options.snow = if(verbose){opts})%dopar%{
    
    unlist(genSED(
      SFRbulge_d = SFRbulge_d[,i]/h,
      SFRbulge_m = SFRbulge_m[,i]/h,
      SFRdisk = SFRdisk[,i]/h,

      Zbulge_d = Zbulge_d[,i],
      Zbulge_m = Zbulge_m[,i],
      Zdisk = Zdisk[,i],

      redshift = redshift[i],
      time = time-cosdistTravelTime(redshift[i], ref='planck')*1e9,

      tau_birth = tau_birth_galaxies[i,],
      tau_screen = tau_screen_galaxies[i,],
      tau_AGN = 1, #hard coded for now

      pow_birth = pow_birth_galaxies[i,],
      pow_screen = pow_screen_galaxies[i,],
      pow_AGN = -0.7, #hard coded for now

      alpha_SF_birth = alpha_SF_birth, #hard coded for now
      alpha_SF_screen = alpha_SF_screen, #hard coded for now
      alpha_SF_AGN = alpha_SF_AGN, #hard coded for now
      
      emission = emission,
      IGMabsorb = IGMabsorb,
      
      AGNlum = 0, #hard coded for now

      speclib = BC03lr,
      Dale = Dale_NormTot,
      AGN = AGN_UnOb_Sparse,
      LKL10 = LKL10_NormAll,
      
      filtout = filtout,

      sparse = sparse,
      intSFR = intSFR
    ))
  }

  stopCluster(cl)

  if(verbose){
    close(pb)
  }

  outSED = as.data.table(rbind(outSED))
  outSED = cbind(id_galaxy=Shark_SFH[['galaxies/id_galaxy']][select], outSED)

  Shark_SFH$close()

  if(verbose){
    message(paste('Finished Viperfish on Shark -',round(proc.time()[3]-timestart,3),'sec'))
  }

  if (write_final_file) {
    outdir = paste(path_shark, 'Photometry', snapshot, subvolume, sep='/')
    write.SED(outSED, filters, outdir, final_file_output, verbose=TRUE)
  }

  class(outSED) = c(class(outSED),'Viperfish-Shark')
  invisible(outSED)
}
