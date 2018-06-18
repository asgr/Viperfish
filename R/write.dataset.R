write.group.safe = function(filename='temp.hdf5', group='group'){
  assertCharacter(filename, max.len=1)
  assertCharacter(group, max.len=1)

  assertPathForOutput(filename, overwrite=TRUE)
  file.h5=h5file(filename=filename, mode='a')

  splitgroup=strsplit(group,'/')[[1]]

  for(i in 1:length(splitgroup)){
    current=paste(splitgroup[1:i],collapse='/')
    if(!file.h5$exists(current)){
      file.h5$create_group(current)
    }
  }

  file.h5$close()
}

write.custom.dataset = function(filename='temp.hdf5', # hd5 filename
                         group='group', # name of group
                         object, # array of any dimension
                         mode="a",
                         dataset.name='matrix',
                         dataset.type=h5types$H5T_IEEE_F32LE,
                         compression.level = 0,
                         overwrite=FALSE) {
  assertCharacter(filename, max.len=1)
  assertCharacter(group, max.len=1)
  assertCharacter(mode, max.len=1)
  assertCharacter(dataset.name, max.len=1)
  if(!testEnvironment(dataset.type)){stop('dataset.type must be one of the allowed datatypes provided by h5types accessed with a $')}
  assertInt(compression.level)
  assertFlag(overwrite)

  write.group.safe(filename=filename, group=group)

  if(!missing(object)){

    assertPathForOutput(filename, overwrite=TRUE)
    file.h5=h5file(filename=filename, mode=mode)

    if(file.h5$exists(paste0(group,'/',dataset.name))){
      if(overwrite){
        file.h5$link_delete(paste0(group,'/',dataset.name))
      }else{
        stop('data set is already present with name ',dataset.name,'!')
      }
    }


    dims = dim(object)
    space = H5S$new(dims = dims)
    ds = H5P_DATASET_CREATE$new()
    ds$set_chunk(dims)$set_fill_value(dataset.type, 1)$set_nbit()
    dataset.h5 = file.h5[[group]]$create_dataset(name = dataset.name, space = space,
                                    dtype = dataset.type, dataset_create_pl = ds,
                                    gzip_level = compression.level)
    dataset.h5[,] = object
    file.h5$close()
  }
}

write.SED=function(SED, filename='temp.hdf5', overwrite=FALSE, filters=c('FUV', 'NUV', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1', 'W2', 'W3', 'W4', 'P100', 'P160', 'S250', 'S350', 'S500')){

  assertPathForOutput(filename, overwrite=TRUE)
  file.h5=h5file(filename=filename, mode='a')

  if(file.h5$exists('filters')){
      if(overwrite){
        file.h5$link_delete('filters')
      }else{
        stop('data set is already present with name filters!')
      }
  }

  if(file.h5$exists('id_galaxy_sky')){
      if(overwrite){
        file.h5$link_delete('id_galaxy_sky')
      }else{
        stop('data set is already present with name id_galaxy_sky!')
      }
  }

  file.h5[['filters']]=filters
  file.h5[['id_galaxy_sky']]=SED[,1]
  file.h5$close()

  Ncol=length(filters)
  colrun=1:Ncol

  SED=as.matrix(SED)
  SED[!is.finite(SED)]=-999

  write.custom.dataset(filename=filename, group="SED/ab_nodust", object=SED[,1+colrun], dataset.name='bulge', overwrite=overwrite)
  write.custom.dataset(filename=filename, group="SED/ab_nodust", object=SED[,1+colrun+Ncol], dataset.name='disk', overwrite=overwrite)
  write.custom.dataset(filename=filename, group="SED/ab_nodust", object=SED[,1+colrun+Ncol*2], dataset.name='total', overwrite=overwrite)

  write.custom.dataset(filename=filename, group="SED/ap_nodust", object=SED[,1+colrun+Ncol*3], dataset.name='bulge', overwrite=overwrite)
  write.custom.dataset(filename=filename, group="SED/ap_nodust", object=SED[,1+colrun+Ncol*4], dataset.name='disk', overwrite=overwrite)
  write.custom.dataset(filename=filename, group="SED/ap_nodust", object=SED[,1+colrun+Ncol*5], dataset.name='total', overwrite=overwrite)

  write.custom.dataset(filename=filename, group="SED/ab_dust", object=SED[,1+colrun+Ncol*6], dataset.name='bulge', overwrite=overwrite)
  write.custom.dataset(filename=filename, group="SED/ab_dust", object=SED[,1+colrun+Ncol*7], dataset.name='disk', overwrite=overwrite)
  write.custom.dataset(filename=filename, group="SED/ab_dust", object=SED[,1+colrun+Ncol*8], dataset.name='total', overwrite=overwrite)

  write.custom.dataset(filename=filename, group="SED/ap_dust", object=SED[,1+colrun+Ncol*9], dataset.name='bulge', overwrite=overwrite)
  write.custom.dataset(filename=filename, group="SED/ap_dust", object=SED[,1+colrun+Ncol*10], dataset.name='disk', overwrite=overwrite)
  write.custom.dataset(filename=filename, group="SED/ap_dust", object=SED[,1+colrun+Ncol*11], dataset.name='total', overwrite=overwrite)
}

write.SFH=function(SFHlist, filename='temp.hdf5', overwrite=FALSE){
  assertPathForOutput(filename, overwrite=TRUE)

  write.custom.dataset(filename=filename, group="SFR", object=as.matrix(SFHlist$SFRbulge), dataset.name='bulge', overwrite=overwrite)
  write.custom.dataset(filename=filename, group="SFR", object=as.matrix(SFHlist$SFRdisk), dataset.name='disk', overwrite=overwrite)

  write.custom.dataset(filename=filename, group="Z", object=as.matrix(SFHlist$Zbulge), dataset.name='bulge', overwrite=overwrite)
  write.custom.dataset(filename=filename, group="Z", object=as.matrix(SFHlist$Zdisk), dataset.name='disk', overwrite=overwrite)
}
