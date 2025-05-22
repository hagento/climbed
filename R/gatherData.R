#' Gather Heating and Cooling Degree Day Outputs
#'
#' This function collects and combines degree day calculation outputs (HDDs and CDDs)
#' from multiple CSV files based on a specified file mapping. It is designed to
#' streamline the integration of multiple model outputs into a single data frame for
#' further analysis and post-processing.
#'
#' The function searches for files matching the unique hash \code{runTag}.
#'
#' @param runTag A string specifying the unique batch tag for temporary files
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

gatherData <- function(runTag,
                       outDir) {

  # directory of temporary outputs
  hddcddDir <- "hddcdd"


  # construct list of file names
  filePaths <- list.files(path = file.path(outDir, hddcddDir),
                          pattern = runTag,
                          full.names = TRUE)


  # read data
  data <- do.call(rbind, lapply(filePaths, read.csv))

  return(data)
}
