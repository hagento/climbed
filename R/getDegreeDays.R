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
#' and ambient temperatures. Defaults to:
#'   \code{c("tLim" = 2, "tAmb" = 2)}.
#'
#' @param ssp A character vector specifying the SSP scenarios for population data. Defaults to:
#'   \code{c("historical", "SSP2")}.
#'
#' @param outDir A string specifying the absolute path to the \code{output} directory.
#' The output directory will also include \code{logs} and \code{tmp} subdirectories for log
#' and temporary files. If a relative path is provided, the package directory is
#' used as the working directory.
#'
#' @author Hagen Tockhorn
#'
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @importFrom utils read.csv
#' @importFrom stats setNames
#' @importFrom piamutils getSystemFile
#' @importFrom fs is_absolute_path
#'
#' @export

getDegreeDays <- function(mappingFile = NULL,
                          bait = TRUE,
                          tLim = list("HDD" = seq(9, 19), "CDD" = seq(15, 25)),
                          std  = c("tLim" = 2, "tAmb" = 2),
                          ssp  = c("historical", "SSP2"),
                          outDir = "output",
                          fileRev = NULL) {
  # CHECKS ---------------------------------------------------------------------

  # check for mapping file
  if (is.null(mappingFile)) {
    stop("No mapping file was provided.")
  }

  # absolute path
  mappingFile <- getSystemFile("extdata", "sectoral", mappingFile, package = "climbed")

  if (!file.exists(mappingFile)) {
    stop("Provided mapping file does not exist in /extdata/sectoral.")
  }


  # PARAMETERS -----------------------------------------------------------------

  # process arguments
  args <- processArgs(tLim = tLim, std = std, ssp = ssp)

  # se standardized arguments
  tLim <- args$tLim
  std  <- args$std
  ssp  <- args$ssp


  # range of pre-calculated HDD/CDD-values, e.g. [173, 348] K, converted to [C]
  tLow <- 153 - 273.15
  tUp  <- 383 - 273.15


  #--- BAIT parameters

  # The weights (wRSDS, wSFC, wHUSS) for calcBAIT and the smoothing coefficient are assumed to
  # be region-independent and equal to the mean of the values given in Staffell
  # et al. 2023.

  # The weights below give the respective weight of the difference between the
  # climate factor and its counterfactual in calculating BAIT.
  wRSDS <- 0.012
  wSFC <- -0.20
  wHUSS <- 0.05

  # The smoothing coefficient defines the assumed thermal intertia of the building
  # and weighs the influence of the internal temperature of the preceding two days.
  sig <- 0.50

  # The blending parameters for the blending of BAIT and raw temperature are like-wise
  # taken from the paper.
  # At bLower, we consider only BAIT. For higher temperatures, we assume a mix with the ambient
  # temperature reaching the highest contribution of bMax at bUpper following a sigmoid function
  bLower <- 15
  bUpper <- 23
  bMax   <- 0.5

  # concatenate to vector
  wBAIT <- list("wRSDS"  = wRSDS,
                "wSFC"   = wSFC,
                "wHUSS"  = wHUSS,
                "sig"    = sig,
                "bLower" = bLower,
                "bUpper" = bUpper,
                "bMax"   = bMax)


  # scenario matrix
  # TODO: bring this into a mapping file
  sspMapping <- list(
    "historical" = "historical",
    "ssp1"       = c("1.9", "2.6", "3.4", "4.5", "6.0"),
    "ssp2"       = c("1.9", "2.6", "3.4", "4.5", "6.0", "7.0"),
    "ssp3"       = c("3.4", "4.5", "6.0", "7.0"),
    "ssp4"       = c("2.6", "3.4", "4.5", "6.0"),
    "ssp5"       = c("2.6", "3.4", "4.5", "6.0", "7.0", "8.5")
  )


  # population mapping
  # TODO: bring this into a mapping file
  popMapping <- c("historical" = "population_histsoc_30arcmin_annual_1901_2021.nc",
                  "ssp2"       = "population_ssp2_30arcmin_annual_2015_2100.nc")


  # track submitted jobs
  allJobs <- list()



  # READ-IN DATA ---------------------------------------------------------------

  # file mapping
  fileMapping <- read.csv2(mappingFile)



  # CHECKS ---------------------------------------------------------------------

  # check if mapping file contains correct columns
  mappingCols <- c("gcm", "rcp", "start", "end", "tas", "rsds", "sfcwind", "huss")
  if (!all(mappingCols %in% colnames(fileMapping))) {
    missing <- mappingCols[!mappingCols %in% colnames(fileMapping)]
    stop("Please provide file mapping with correct columns.\n Missing columns: ",
         paste(missing, collapse=", "))
  }

  # create output directory if it doesn't exist
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

  # create output logs directory if it doesn't exist
  dir.create(file.path(outDir, "logs"), showWarnings = FALSE, recursive = TRUE)



  # PROCESS DATA ---------------------------------------------------------------

  # calculate HDD/CDD-factors
  hddcddFactor <- compFactors(tLow = tLow, tUp = tUp, tLim, std[["tAmb"]], std[["tLim"]])


  # get the main job ID if running as a slurm job
  mainJobId <- getMainJobId()
  message("Main job ID: ", mainJobId)



  # --- CALCULATE DEGREE DAYS

  for (s in ssp) {
    message("\nProcessing SSP scenario: ", s)

    # read in population data
    message("Reading population data...")
    pop <- importData(subtype = popMapping[[s]])


    # filter compatible RCP scenarios
    files <- fileMapping %>%
      filter(.data[["rcp"]] %in% sspMapping[[s]])

    if (nrow(files) == 0) {
      stop("Provided SSP scenario not in file mapping.")
    }


    message("Submitting ", nrow(files), " jobs...")

    # submit jobs and collect job details
    for (i in seq(1, nrow(files))) {
      tryCatch(
        {
          job <- createSlurm(fileRow = files[i, ],
                             pop = pop,
                             ssp = s,
                             bait = bait,
                             tLim = tLim,
                             hddcddFactor = hddcddFactor,
                             wBAIT = wBAIT,
                             outDir = outDir)

          allJobs[[length(allJobs) + 1]] <- job
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
  jobIds <- sapply(allJobs, function(x) x$jobId) # nolint

  # wait for our specific jobs to complete (max. 6hrs)
  waitForSlurm(jobIds, maxWaitTime = 6 * 60 * 60)



  # --- GATHER AND SMOOTH DEGREE DAYS

  # gather outputs from slurm jobs
  data <- gatherData(fileMapping = fileMapping, outDir = outDir)

  # smooth degree days
  dataSmooth <- smoothDegreeDays(data, nSmoothIter = 50, transitionYears = 5)




  # OUTPUT ---------------------------------------------------------------------

  fileName <- "hddcdd"

  if(is.character(fileRev)) {
    fileName <- paste0(fileName, "_", fileRev)
  } else if(!is.null(fileRev)) {
    warning("fileRev must be character")
  }

  outPath <- file.path(outDir, paste0(fileName, ".csv"))

  write.csv(dataSmooth, outPath)


}
