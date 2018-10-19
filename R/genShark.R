genShark=function(path_shark='.', snapshot=NULL, subvolume=NULL, redshift="get", h='get', cores=4, id_galaxy_sam='all', filters=c('FUV', 'NUV', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1', 'W2', 'W3', 'W4', 'P100', 'P160', 'S250', 'S350', 'S500'), tau_birth=1.5, tau_screen=0.5, sparse=5, final_file_output='Shark-SED.csv', intSFR=TRUE, verbose=TRUE, write_final_file=FALSE){

  timestart=proc.time()[3]

  if(verbose){
    message('Running Viperfish on Shark')
  }

  assertCharacter(path_shark, max.len=1)
  assertAccess(path_shark, access='r')
  assertInt(snapshot, null.ok=TRUE)
  assertInt(subvolume, null.ok=TRUE)
  assertInt(cores)
  if(is.list(filters)){
    filterlist=TRUE
  }else{
    filterlist=FALSE
    assertCharacter(filters)
  }
  assertScalar(tau_birth)
  assertScalar(tau_screen)
  assertInt(sparse)
  assertFlag(intSFR)
  assertFlag(verbose)


  BC03lr=Dale_Msol=SFH=i=subsnapID=Nid=idlist=subsnapID=NULL

  data("BC03lr", envir = environment())
  data("Dale_Msol", envir = environment())

  if(filterlist==FALSE){
    filtout=foreach(i = filters)%do%{getfilt(i)}
    names(filtout)=filters
  }else{
    filtout=filters
  }

  sfh_fname = paste(path_shark, snapshot, subvolume, 'star_formation_histories.hdf5',sep='/')
  assertAccess(sfh_fname, access='r')
  Shark_SFH=h5file(sfh_fname, mode='r')

  # Read things in from the Shark HDF5 file:

  time=Shark_SFH[['lbt_mean']][]*1e9

  if(redshift[1]=='get'){
    redshift=Shark_SFH[['run_info/redshift']][]
    if(redshift==0){redshift=1e-10}
  }

  if(h=='get'){
    h=Shark_SFH[['cosmology/h']][]
  }

  assertNumeric(redshift)
  assertScalar(h)

  if(id_galaxy_sam[1]=='all'){
    select=1:Shark_SFH[['galaxies/id_galaxy']]$dims
  }else{
    select=match(id_galaxy_sam, Shark_SFH[['galaxies/id_galaxy']][])
  }

  SFRbulge_d=Shark_SFH[['bulges_diskins/star_formation_rate_histories']][,select,drop=FALSE]
  SFRbulge_m=Shark_SFH[['bulges_mergers/star_formation_rate_histories']][,select,drop=FALSE]
  SFRdisk=Shark_SFH[['disks/star_formation_rate_histories']][,select,drop=FALSE]
  Zbulge_d=Shark_SFH[['bulges_diskins/metallicity_histories']][,select,drop=FALSE]
  Zbulge_m=Shark_SFH[['bulges_mergers/metallicity_histories']][,select,drop=FALSE]
  Zdisk=Shark_SFH[['disks/metallicity_histories']][,select,drop=FALSE]

  if(length(redshift)==1){
    redshift=rep(redshift, length(select))
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
  unlist(genSED(SFRbulge_d=SFRbulge_d[,i]/h, SFRbulge_m=SFRbulge_m[,i]/h, SFRdisk=SFRdisk[,i]/h, redshift=redshift[i], time=time, speclib=BC03lr, Zbulge_d=Zbulge_d[,i], Zbulge_m=Zbulge_m[,i], Zdisk=Zdisk[,i], filtout=filtout, Dale=Dale_Msol, tau_birth=tau_birth, tau_screen=tau_screen, sparse=sparse, intSFR=intSFR))
  }

  stopCluster(cl)

  if(verbose){
    close(pb)
  }

  outSED=as.data.table(rbind(outSED))
  colnamesSED=c(
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
  outSED=cbind(id_galaxy=Shark_SFH[['galaxies/id_galaxy']][select], outSED)

  Shark_SFH$close()

  if(verbose){
    message(paste('Finished Viperfish on Shark -',round(proc.time()[3]-timestart,3),'sec'))
  }

  if (write_final_file) {
    outdir = paste(path_shark, 'Photometry', snapshot, subvolume, sep='/')
    if (!dir.exists(outdir)) {
      dir.create(outdir, recursive=TRUE)
    }
    outfile = paste(outdir, final_file_output, sep='/')
    if (verbose) {
      message(paste('Writing CSV file on ', outfile))
    }
    fwrite(outSED, file=outfile)
  }

  class(outSED)=c(class(outSED),'Viperfish-Shark')
  invisible(outSED)
}
