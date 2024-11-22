start <- function(mappingFile = NULL,
                  bait = TRUE,
                  tLim = list("HDD" = seq(9, 19), "CDD" = seq(15, 25)),
                  std  = c("tLim" = 2, "tAmb" = 2),
                  ssp  = c("SSP2")) {

  # CHECKS ---------------------------------------------------------------------

  # check for mapping file
  if (is.null(mappingFile)) {
    stop("No mapping file was provided.")
  }


  # PARAMETERS -----------------------------------------------------------------

  # process arguments
  args <- processArgs(tLim = tLim, std = std, ssp = ssp)

  # se standardized arguments
  tLim <- args$tLim
  std  <- args$std
  ssp  <- args$ssp

  # add historical
  ssp <- c("historical", ssp) %>%
    tolower()


  # range of pre-calculated HDD/CDD-values, e.g. [173, 348] K, converted to [C]
  tlow <- 173 - 273.15
  tup  <- 348 - 273.15


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

  # BAIT parameter names
  # (note: this is really just a check in case calcBAITpars returns no layer names)
  parNames <- c("aRSDS", "bRSDS", "aSFC", "bSFC", "aHUSS", "bHUSS") # TODO: change this


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
                  "ssp2" = "population_ssp2_30arcmin_annual_2015_2100.nc")


  # Store job details
  jobDetails <- list()



  # READ-IN DATA ---------------------------------------------------------------

  # file mapping
  fileMapping <- read.csv(paste0("inst/extdata/sectoral/", mappingFile))

  # country mask
  countries <- readSource("ISIMIPbuildings",
                          subtype = "countrymasks-fractional_30arcmin.nc",
                          # subtype = "countrymasks.nc",
                          convert = FALSE)


  # bait regression parameters
  if (bait) {
    cacheDir <- "/p/projects/rd3mod/inputdata/sources/BAITpars"

    baitPars <- calcOutput("BAITpars",
                           aggregate = FALSE,
                           model     = unique(fileMapping$gcm),
                           cacheDir  = cacheDir)

    names(baitPars) <- parNames
  }



  # CHECKS ---------------------------------------------------------------------

  # check if mapping file contains correct columns
  mappingCols <- c("gcm", "rcp", "start", "end", "tas", "rsds", "sfc", "huss")
  # if (!mappingCols %in% colnames(mappingFile)) {
  #   stop("Please provide file mapping with correct columns.\n Missing columns: ")
  # }

  # Create output directory if it doesn't exist
  dir.create("output", showWarnings = FALSE, recursive = TRUE)

  # Create output logs directory if it doesn't exist
  dir.create("output/logs", showWarnings = FALSE, recursive = TRUE)



  # PROCESS DATA ---------------------------------------------------------------

  # calculate HDD/CDD-factors
  hddcddFactor <- compHDDCDDFactors(tlow = tlow, tup = tup, tLim, std[["tAmb"]], std[["tLim"]])


  # calculate degree days
  for (s in ssp) {

    # read in population data
    pop <- readSource("ISIMIPbuildings", subtype = popMapping[[s]], convert = FALSE)

    # filter compatible RCP scenarios
    files <- fileMapping %>%
      filter(.data[["rcp"]] %in% sspMapping[[s]])

    # loop over and calculate degree days for each row
    for (i in seq(1, nrow(files))) {
      job <- createSlurm(fileMapping  = files[i],
                         pop          = pop,
                         ssp          = s,
                         bait         = bait,
                         tLim         = tLim,
                         countries    = countries,
                         hddcddFactor = hddcddFactor,
                         wBAIT        = wBAIT,
                         baitPars     = baitPars)

      # Submit the Slurm job
      submitSlurm(job)

      # Store job details
      JobDetails[[length(JobDetails) + 1]] <- job
    }
  }

  # Wait for all jobs to complete
  waitForSlurm()

  print("all jobs done")
}






