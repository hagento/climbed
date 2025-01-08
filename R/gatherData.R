#' Gather Heating and Cooling Degree Day Outputs
#'
#' This function collects and combines degree day calculation outputs (HDDs and CDDs)
#' from multiple CSV files based on a specified file mapping. It is designed to
#' streamline the integration of multiple model outputs into a single data frame for
#' further analysis and post-processing.
#'
#' The function searches for files matching the pattern:
#' \code{"hddcdd_[gcm]_[rcp]_[start]-[end].csv"} in the \code{"hddcdd"} subdirectory
#' of the specified \code{outDir}. Each file corresponds to a specific General Circulation
#' Model (GCM), Representative Carbon Pathway (RCP), and time period.
#'
#' @param fileMapping A data frame containing metadata for locating and processing
#' degree day output files. The following columns are required:
#'   \describe{
#'     \item{\code{"gcm"}}{The General Circulation Model (GCM) name.}
#'     \item{\code{"rcp"}}{The RCP scenario (e.g., RCP2.6, RCP8.5).}
#'     \item{\code{"start"}}{The starting year of the time period.}
#'     \item{\code{"end"}}{The ending year of the time period.}
#'   }
#'
#' @param outDir A string specifying the absolute path to the \code{output} directory
#' where degree day outputs are stored.
#'
#' @return A data frame that combines all gathered degree day outputs, organized by region
#' and scenario.
#'
#' @author Hagen Tockhorn
#'
#' @importFrom dplyr select
#' @importFrom utils read.csv


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
