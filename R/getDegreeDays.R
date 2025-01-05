#' Calculate Heating and Cooling Degree Days
#'
#' This function initiates the calculation of annual Heating (HDDs) and Cooling Degree Days (CDDs)
#' for historical and future scenarios. The calculation supports raw near-surface
#' air temperature and the bias-adjusted internal temperature (BAIT) which incorporates
#' additional climate variables.
#' The user can specify Shared Socioeconomic Pathways (SSPs) to define future population
#' development. The corresponding climate scenarios, represented by Representative Carbon Pathways (RCPs),
#' are chosen according to the IPCC scenario matrix.
#' Annual degree days are calculated per data set (typically 10-year period) in an
#' individually defined and submitted SLURM job. After successful completion, the job
#' outputs are gathered, post-processed and saved as a .csv file in the \code{output} folder.
#'
#' @param mappingFile A string specifying the path to the mapping file containing information about the climate data
#' Must include columns like \code{"gcm"}, \code{"rcp"}, \code{"start"}, \code{"end"}, \code{"tas"}, \code{"rsds"},
#' \code{"sfcwind"} and \code{"huss"}.
#' @param bait Logical. If \code{TRUE}, bias-adjusted internal temperature (BAIT)
#' is included in the calculations. Defaults to \code{TRUE}.
#' @param tLim A list defining temperature limits for HDD and CDD calculations.
#' Defaults to \code{list("HDD" = seq(9, 19), "CDD" = seq(15, 25))}.
#' @param std A named vector of standardization parameters for temperature limits and ambient temperatures.
#' Defaults to \code{c("tLim" = 2, "tAmb" = 2)}.
#' @param ssp A character vector specifying the SSP scenarios to include. Defaults to \code{c("historical", "SSP2")}.
#' @param outDir A string specifying the absolute path to the \code{output} directory. This will also
#' contain the \code{logs} and \code{tmp} directories for log and temporary files.
#' @param globalPars \code{logical} indicating whether to use global or gridded BAIT parameters
#' (required if \code{bait} is TRUE).
#'
#' @returns Saves a \code{.csv} file containing the calculated degree days.
#'
#' @author Hagen Tockhorn
#'
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @importFrom utils read.csv
#' @importFrom stats setNames
#' @importFrom piamutils getSystemFile
#'
#' @export

getDegreeDays <- function(mappingFile = NULL,
                          bait = TRUE,
                          tLim = list("HDD" = seq(9, 19), "CDD" = seq(15, 25)),
                          std  = c("tLim" = 2, "tAmb" = 2),
                          ssp  = c("historical", "SSP2"),
                          outDir = "output",
                          globalPars = FALSE) {
  # CHECKS ---------------------------------------------------------------------

  # check for mapping file
  if (is.null(mappingFile)) {
    stop("No mapping file was provided.")
  }

  # absolute path
  mappingFile <- getSystemFile("extdata", "mappings", mappingFile, package = "climbed")

  if (!file.exists(mappingFile)) {
    stop("Provided mapping file does not exist in /extdata/mappings.")
  }



  # PARAMETERS -----------------------------------------------------------------

  # process arguments
  args <- processArgs(tLim = tLim, std = std, ssp = ssp)

  # se standardized arguments
  tLim <- args$tLim
  std  <- args$std
  ssp  <- args$ssp


  # range of pre-calculated HDD/CDD-values, e.g. [173, 348] K, converted to [C]
  tLow <- 173 - 273.15
  tUp  <- 348 - 273.15


  # track submitted jobs
  allJobs <- list()



  # READ-IN DATA ---------------------------------------------------------------

  # file mapping
  fileMapping <- read.csv2(mappingFile)


  # external BAIT weights and parameters
  wBAIT <- read.csv2(getSystemFile("extdata", "mappings", "BAITweights.csv",
                                   package = "climbed"),
                     stringsAsFactors = FALSE)


  # population files
  popMapping <- read.csv2(getSystemFile("extdata", "mappings", "populationMapping.csv",
                                        package = "climbed"),
                          stringsAsFactors = FALSE)


  # scenario matrix
  scenMatrix <- read.csv2(getSystemFile("extdata", "mappings", "scenarioMatrix.csv",
                                        package = "climbed"),
                          stringsAsFactors = FALSE)



  # CHECKS ---------------------------------------------------------------------

  # check if mapping file contains correct columns
  mappingCols <- c("gcm", "rcp", "start", "end", "tas", "rsds", "sfc", "huss")
  missingCols <- setdiff(mappingCols, colnames(fileMapping)) # Identify missing columns
  if (length(missingCols) > 0) { # Check if there are any missing columns
    stop("Please provide file mapping with correct columns.\nMissing columns:\n",
         paste(missingCols, collapse = ", "))
  }

  # create output directory if it doesn't exist
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

  # create output logs directory if it doesn't exist
  dir.create(file.path(outDir, "logs"), showWarnings = FALSE, recursive = TRUE)



  # PROCESS DATA ---------------------------------------------------------------

  # --- Prepare Mappings

  # BAIT weights
  wBAIT <- wBAIT %>%
    pull("value", "variable") %>%
    as.list()

  # population data
  popMapping <- popMapping  %>%
    pull("file", "scenario") %>%
    as.list()

  # scenario matrix
  scenMatrix <- scenMatrix %>%
    pull("rcp", "ssp") %>%
    strsplit(",") %>%
    lapply(trimws)


  # calculate HDD/CDD-factors
  hddcddFactor <- compFactors(tLow = tLow, tUp = tUp, tLim, std[["tAmb"]], std[["tLim"]])


  # get the main job ID if running as a slurm job
  mainJobId <- getMainJobId()
  message("Main job ID: ", mainJobId)



  # --- CALCULATE DEGREE DAYS

  for (s in ssp) {
    message("\nProcessing SSP scenario: ", s)

    # read in population data
    pop <- importData(subtype = popMapping[[s]])

    # filter compatible RCP scenarios
    files <- fileMapping %>%
      filter(.data[["rcp"]] %in% scenMatrix[[s]])

    if (nrow(files) == 0) {
      stop("Provided SSP scenario not in file mapping.")
    }


    message("Submitting ", nrow(files), " jobs...")

    # submit jobs and collect job details
    for (i in seq(1, nrow(files))) {
      message("\nSubmitting job ", i, " of ", nrow(files))
      tryCatch(
        {
          job <- createSlurm(fileRow = files[i, ],
                             pop = pop,
                             ssp = s,
                             bait = bait,
                             tLim = tLim,
                             hddcddFactor = hddcddFactor,
                             wBAIT = wBAIT,
                             outDir = outDir,
                             globalPars = globalPars)

          allJobs <- c(allJobs, job)
          message("Job submitted successfully")
        },
        error = function(e) {
          warning("Failed to submit job ", i, ": ", e$message)
        }
      )
    }
  }

  if (length(allJobs) == 0) {
    stop("No jobs were successfully submitted")
  }

  message("\nSubmitted ", length(allJobs), " jobs successfully")
  message("Waiting for jobs to complete...")

  # extract all job IDs
  jobIds <- lapply(allJobs, function(x) x$jobId)

  # wait for our specific jobs to complete (max. 12hrs)
  waitForSlurm(jobIds, maxWaitTime = 12 * 60 * 60)



  # now gather all jobs and do some further processing ...
}
