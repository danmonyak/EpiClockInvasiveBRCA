subdir_list <- strsplit(getwd(), .Platform$file.sep)[[1]]
while (subdir_list[length(subdir_list)] != 'EpiClockInvasiveBRCA') {
  subdir_list <- subdir_list[-length(subdir_list)]
}
repo_dir <- paste(subdir_list, collapse=.Platform$file.sep)
