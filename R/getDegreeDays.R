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
#' @param mappingFile A string specifying the path to the mapping file containing information about the climate data.
#' The file must include the following columns:
#'   \describe{
#'     \item{\code{"gcm"}}{The General Circulation Model (GCM) name of the dataset.}
#'     \item{\code{"rcp"}}{The RCP scenario.}
#'     \item{\code{"start"}}{The starting year of the time period.}
#'     \item{\code{"end"}}{The ending year of the time period.}
#'     \item{\code{"tas"}}{The filename of the near-surface air temperature data. Needs
#'     to exist in the \code{"ISIMIPbuildings"} directory within your madrat \code{"sourcefolder"}.}
#'     \item{\code{"rsds"}}{The filename of the near-surface air temperature data. Needs
#'     to exist in the \code{"ISIMIPbuildings"} directory within your madrat \code{"sourcefolder"}.}
#'     \item{\code{"sfcwind"}}{The filename of the near-surface wind speed data. Needs
#'     to exist in the \code{"ISIMIPbuildings"} directory within your madrat \code{"sourcefolder"}.}
#'     \item{\code{"huss"}}{The filename of the near-surface relative humidity data. Needs
#'     to exist in the \code{"ISIMIPbuildings"} directory within your madrat \code{"sourcefolder"}.}
#'   }
#'
#' @param bait Logical. If \code{TRUE}, bias-adjusted internal temperature (BAIT)
#' is included in the calculations. Defaults to \code{TRUE}.
#'
#' @param tLim A list defining temperature limits for HDD and CDD calculations.
#' Defaults to:
#'   \code{list("HDD" = seq(9, 19), "CDD" = seq(15, 25))}.
#'
#' @param std A named vector of standard deviations for temperature limits
#' and ambient temperatures in K. Defaults to:
#'   \code{c("tLim" = 2, "tAmb" = 2)}.
#'
#' @param ssp A character vector specifying the SSP scenarios for population data. Defaults to:
#'   \code{c("historical", "SSP2")}.
#'
#' @param outDir A string specifying the absolute or relative path to the \code{output} directory.
#' The output directory will also include \code{logs} and \code{tmp} subdirectories for log
#' and temporary files. If a relative path is provided, the directory is created relative
#' to the current working directory.
#'
#' @param fileRev (Optional) A string specifying the revision number to identify the output file.
#'
#' @param globalPars Logical. Indicates whether to use global or gridded BAIT parameters
#' (required if \code{bait} is TRUE).
#'
#' @param endOfHistory An integer specifying the upper temporal limit for historical data.
#' Defaults to 2025.
#'
#' @param noCC \code{logical} indicating whether to compute a no-climate-change scenario.
#'        If \code{TRUE}, the function will calculate degree days assuming constant climate conditions.
#'        Default is \code{FALSE}.
#'
#' @param packagePath (Optional) A string specifying the path to the package for development mode.
#' If provided, the SLURM jobs will use devtools::load_all() instead of library() to load the package.
#'
#' @author Hagen Tockhorn
#'
#' @importFrom dplyr filter pull anti_join
#' @importFrom magrittr %>%
#' @importFrom utils read.csv
#' @importFrom stats setNames
#' @importFrom piamutils getSystemFile
#' @importFrom utils write.csv
#'
#' @export

getDegreeDays <- function(mappingFile = NULL,
                          bait = TRUE,
                          tLim = list("HDD" = seq(9, 19), "CDD" = seq(15, 25)),
                          std  = c("tLim" = 2, "tAmb" = 2),
                          ssp  = c("historical", "SSP2"),
                          outDir = "output",
                          fileRev = NULL,
                          globalPars = TRUE,
                          endOfHistory = 2025,
                          noCC = FALSE,
                          packagePath = NULL) {

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


  # indicate whether SSP2 was added for historical fill-up
  sspAddedForFillup <- FALSE



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

  if (!is.null(fileMapping)) {
    # check if mapping file contains correct columns
    mappingCols <- c("gcm", "rcp", "start", "end", "tas", "rsds", "sfcwind", "huss")
    missingCols <- setdiff(mappingCols, colnames(fileMapping)) # Identify missing columns

    if (length(missingCols) > 0) { # Check if there are any missing columns
      stop("Please provide file mapping with correct columns.\nMissing columns:\n",
           paste(missingCols, collapse = ", "))
    }
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

  # mutate values to numerics
  wBAIT <- lapply(wBAIT, as.numeric)


  # population data
  popMapping <- popMapping  %>%
    pull("file", "scenario") %>%
    as.list()

  # scenario matrix
  scenMatrix <- scenMatrix %>%
    pull("rcp", "ssp") %>%
    strsplit(",") %>%
    lapply(trimws)


  # if no file mapping is given for noCC-case, use default
  if (is.null(fileMapping) && isTRUE(noCC)) {
    fileMapping <- read.csv2(getSystemFile("extdata", "mappings", "ISIMIPbuildings_fileMapping.csv",
                                           package = "climbed")) %>%
      filter(.data$rcp == "historical")
  }


  # indicate original mapping entries
  fileMapping <- fileMapping %>%
    mutate(originalFile = TRUE)

  # if necessary, add files to fill up missing historical data points
  if (max(fileMapping$end[fileMapping$rcp == "historical"]) < endOfHistory &&
        !(any(fileMapping$rcp == "2.6" &
                fileMapping$start <= max(fileMapping$end[fileMapping$rcp == "historical"]) + 1 &
                fileMapping$end >= endOfHistory) &&
            length(unique(fileMapping$gcm[fileMapping$rcp == "2.6"])) == 5)) {

    fullMappingFile <- getSystemFile("extdata", "mappings", "ISIMIPbuildings_fileMapping.csv",
                                     package = "climbed")
    fileMapping <- fullMappingFile %>%
      read.csv2() %>%
      # ensure the files cover the period from the end of historical data to endOfHistory
      filter(.data$rcp == "2.6",
             (.data$start <= max(fileMapping$end[fileMapping$rcp == "historical"]) + 1 &
                .data$end >= endOfHistory)) %>%
      # filter out duplicated entries
      anti_join(fileMapping, by = c("gcm", "rcp", "start", "end")) %>%
      mutate(originalFile = FALSE) %>%
      rbind(fileMapping)

    # If we added fill-up files and SSP2 is not in the original list, add it now
    if (nrow(fileMapping[fileMapping$originalFile == FALSE, ]) > 0 &&
          !("SSP2" %in% ssp)) {
      ssp <- c(ssp, "ssp2")
      sspAddedForFillup <- TRUE
    }
  }


  # calculate HDD/CDD-factors
  hddcddFactor <- compFactors(tLow = tLow, tUp = tUp, tLim, std[["tAmb"]], std[["tLim"]])


  # get the main job ID if running as a slurm job
  mainJobId <- getMainJobId()
  message("Main job ID: ", mainJobId)



  # --- CALCULATE DEGREE DAYS

  for (s in ssp) {
    message("\nProcessing SSP scenario: ", s)

    # filter compatible RCP scenarios and ensure that added SSPs are only used for fill-up
    if (s == "ssp2" && isTRUE(sspAddedForFillup)) {
      files <- fileMapping %>%
        filter(!.data$originalFile)
    } else {
      files <- fileMapping %>%
        filter(.data[["rcp"]] %in% scenMatrix[[s]])
    }

    if (nrow(files) == 0) {
      stop("Provided SSP scenario not in file mapping.")
    }


    message("Submitting ", nrow(files), " jobs...")

    # submit jobs and collect job details
    for (i in seq(1, nrow(files))) {
      tryCatch(
        {
          job <- createSlurm(fileRow = files[i, ],
                             pop = popMapping[[s]],
                             ssp = s,
                             bait = bait,
                             tLim = tLim,
                             hddcddFactor = hddcddFactor,
                             wBAIT = wBAIT,
                             outDir = outDir,
                             globalPars = globalPars,
                             noCC = noCC,
                             packagePath = packagePath)

          allJobs[[job$jobId]] <- job
          message("Job ", job$jobId, " submitted successfully")
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
  jobIds <- as.numeric(names(allJobs))

  # wait for our specific jobs to complete (max. 6hrs)
  waitForSlurm(jobIds, maxWaitTime = 6 * 60 * 60)



  # --- GATHER AND SMOOTH DEGREE DAYS

  # gather outputs from slurm jobs
  data <- gatherData(fileMapping = fileMapping, outDir = outDir)


  # if specified, calculate degree days for constant climate
  if (isTRUE(noCC)) {
    # extract directory for grid data
    gridDataDir <- lapply(allJobs, function(x) x$gridDataDir) %>%
      unique() %>%
      unlist()

    dataNoCC <- computeConstantClimate(fileMapping = fileMapping,
                                       ssp = ssp,
                                       popMapping = popMapping,
                                       nHistYears = 5,
                                       gridDataDir = gridDataDir)

    data <- rbind(data, dataNoCC)
  }


  # smooth degree days
  dataSmooth <- smoothDegreeDays(data,
                                 fileMapping,
                                 nSmoothIter = 100,
                                 transitionYears = 10,
                                 endOfHistory = endOfHistory,
                                 noCC = noCC,
                                 predictTransition = FALSE)




  # OUTPUT ---------------------------------------------------------------------

  fileName <- "hddcdd"

  if (is.character(fileRev)) {
    fileName <- paste0(fileName, "_", fileRev)
  } else if (!is.null(fileRev)) {
    warning("fileRev must be character")
  }

  outPath <- file.path(outDir, paste0(fileName, ".csv"))

  write.csv(dataSmooth, outPath, row.names = FALSE)

  return(invisible(outPath))
}
