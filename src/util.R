# Determine repo_dir base path from current working directory
# Must be run inside the repo directory (or a subdirectory)

subdir_list <- strsplit(getwd(), .Platform$file.sep)[[1]]
while (subdir_list[length(subdir_list)] != 'EpiClockInvasiveBRCA') {
  subdir_list <- subdir_list[-length(subdir_list)]
  if (subdir_list[length(subdir_list)] == '') {
    stop('Must be run from inside the local clone of the EpiClockInvasiveBRCA repository...')
  }
}
repo_dir <- paste(subdir_list, collapse=.Platform$file.sep)
