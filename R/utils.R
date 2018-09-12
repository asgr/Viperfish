get_shark_path = function(path=NULL, config_file=NULL)
{
  if (is.null(path) && is.null(config_file)) {
    path = '.'
  }
  else if (!is.null(path) && !is.null(config_file)) {
    stop('path and config_file cannot be specified together')
  }

  if (is.null(path)) {
    assertAccess(config_file, access='r')
    cfg = configr::read.config(config_file)
    path = paste(cfg$execution$output_directory, cfg$simulation$sim_name, cfg$execution$name_model, sep='/')
  }

  assertAccess(path, access='r')
  path
}
