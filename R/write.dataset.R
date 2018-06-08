write.group.safe = function(filename='temp.hdf5', group='group'){
  file.h5=h5file(filename, mode='a')

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

  write.group.safe(filename=filename, group=group)

  if(!missing(object)){

    file.h5=h5file(filename, mode=mode)

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

  Ncol=length(filters)
  colrun=1:Ncol

  write.custom.dataset(filename=filename, group="ab_nodust/bulge", object=testSED[,1+colrun], overwrite=overwrite)
  write.custom.dataset(filename=filename, group="ab_nodust/disk", object=testSED[,1+colrun+Ncol], overwrite=overwrite)
  write.custom.dataset(filename=filename, group="ab_nodust/total", object=testSED[,1+colrun+Ncol*2], overwrite=overwrite)

  write.custom.dataset(filename=filename, group="ap_nodust/bulge", object=testSED[,1+colrun+Ncol*3], overwrite=overwrite)
  write.custom.dataset(filename=filename, group="ap_nodust/disk", object=testSED[,1+colrun+Ncol*4], overwrite=overwrite)
  write.custom.dataset(filename=filename, group="ap_nodust/total", object=testSED[,1+colrun+Ncol*5], overwrite=overwrite)

  write.custom.dataset(filename=filename, group="ab_dust/bulge", object=testSED[,1+colrun+Ncol*6], overwrite=overwrite)
  write.custom.dataset(filename=filename, group="ab_dust/disk", object=testSED[,1+colrun+Ncol*7], overwrite=overwrite)
  write.custom.dataset(filename=filename, group="ab_dust/total", object=testSED[,1+colrun+Ncol*8], overwrite=overwrite)

  write.custom.dataset(filename=filename, group="ap_dust/bulge", object=testSED[,1+colrun+Ncol*9], overwrite=overwrite)
  write.custom.dataset(filename=filename, group="ap_dust/disk", object=testSED[,1+colrun+Ncol*10], overwrite=overwrite)
  write.custom.dataset(filename=filename, group="ap_dust/total", object=testSED[,1+colrun+Ncol*11], overwrite=overwrite)
}
