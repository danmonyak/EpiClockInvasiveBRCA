### EpiClockUtil.R
### Utility functions

keepN_fields <- function(x, n, sep='-') {
  # Return first n fields of segmented string
  # x - string
  # n - number of fields
  # sep - separator pattern
  fields <- strsplit(x, sep)[[1]]
  paste(fields[1:n], collapse = sep)
}
getParticipantIDs <- function(fullTumorIDs) {
  # Get TCGA participant IDs - first 3 fields
  # Vectorized function
  sapply(fullTumorIDs, function (x) keepN_fields(x, 3), USE.NAMES=F)
}
getTumorIDs <- function(fullTumorIDs) {
  # Get TCGA tumor IDs - first 4 fields
  # Vectorized function
  sapply(fullTumorIDs, function (x) keepN_fields(x, 4), USE.NAMES=F)
}
########################