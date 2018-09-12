get_shark_path = function(config_file){
  assertAccess(config_file, access='r')
  cfg = read.ini(config_file)
  path = paste(cfg$execution$output_directory, cfg$simulation$sim_name, cfg$execution$name_model, sep='/')
}
