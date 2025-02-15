#' Gather degree day outputs from multiple files
#'
#' Collects and combines degree day calculation outputs from multiple CSV files
#' based on a specified file mapping.
#'
#' @param fileMapping Data frame containing columns:
#'   gcm, rcp, start, end
#' @param outDir A string specifying the absolute path to the \code{output} directory.
#'
#' @return Data frame combining all found output files
#'
#' @details
#' Searches for files matching pattern: "hddcdd_[gcm]_[rcp]_[start]-[end].csv"
#' in the "hddcdd" subdirectory of \code{outDir}.
#'
#' @importFrom dplyr select
#' @importFrom utils read.csv
#'
#' @export

gatherData <- function(fileMapping,
                       outDir) {

  # directory of temporary outputs
  hddcddDir <- "hddcdd"

  # extract information to identify output files
  fileMapping <- fileMapping %>%
    select("gcm", "rcp", "start", "end")


  # construct list of file names
  filePaths <- vapply(seq_len(nrow(fileMapping)), function(i) {

    # create pattern
    pattern <- paste0(
      "hddcdd_",
      fileMapping$gcm[i], "_",
      fileMapping$rcp[i], "_",
      fileMapping$start[i], "-",
      fileMapping$end[i],
      "\\.csv$"
    )

    # find files matching the pattern in the directory
    matches <- list.files(path = file.path(outDir, hddcddDir),
                          pattern = pattern,
                          full.names = TRUE)

    # return the matched files (empty string if none found, stop if duplicated)
    if (length(matches) == 0) {
      return("")
    } else if (length(matches) == 1) {
      return(matches)
    } else {
      warning(paste0("The following file seems to be duplicated in the output directory: ", pattern, ".\n",
                     "Only the first match will be used."))
      return(matches[[1]])
    }
  },
  character(1))


  # filter out empty strings (files not found)
  filePaths <- filePaths[filePaths != ""]


  # read data
  data <- do.call("rbind", lapply(filePaths, function(f) {
    read.csv(f)
  }))

  return(data)
}
