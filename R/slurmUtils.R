#### SLURM JOB UTILITIES ####


#' Create a Slurm Job for Running initCalculation
#'
#' @param fileRow A data frame row containing file mapping details
#' @param pop Population data to be passed to initCalculation
#' @param ssp Shared Socioeconomic Pathway scenario
#' @param bait Logical indicating whether to use BAIT calculation
#' @param tLim Temperature limits for degree day calculation
#' @param countries Country mask data
#' @param hddCddFactor Heating/Cooling Degree Day factors
#' @param wBait Weights for BAIT calculation
#' @param baitPars BAIT regression parameters
#' @param jobConfig List of Slurm job configuration parameters
#'
#' @return A list containing job details including job name, script path,
#'         output file path, and Slurm submission command
#'
#' @author Hagen Tockhorn


createSlurm <- function(fileRow,
                        pop,
                        ssp,
                        bait,
                        tLim,
                        countries,
                        hddCddFactor,
                        wBait,
                        baitPars,
                        jobConfig = list()) {

  # default job configuration
  defaultConfig <- list(
    outputDir     = "output",
    logsDir       = "output/logs",
    jobNamePrefix = "hddcdd",
    cpusPerTask   = 16,
    nodes         = 1,
    ntasks        = 1,
    partition     = "standard",
    qos           = "short"
  )

  # merge default config with user-provided config
  config <- modifyList(defaultConfig, jobConfig)

  # create output file path
  outputFile <- file.path(
    config$outputDir,
    paste0("result_", fileRow$gcm, "_", fileRow$rcp,
           "_", fileRow$start, "-", fileRow$end, ".csv")
  )

  # Create unique job script with Slurm directives and R commands
  jobScript <- tempfile(pattern = "initcalc_job_", fileext = ".slurm")

  writeLines(c(
    "#!/bin/bash",
    sprintf("#SBATCH --job-name=%s", jobName),
    sprintf("#SBATCH --output=%s/%s.out", config$logsDir, jobName),
    sprintf("#SBATCH --error=%s/%s.err", config$logsDir, jobName),
    sprintf("#SBATCH --cpus-per-task=%s", config$cpusPerTask),
    sprintf("#SBATCH --nodes=%s", config$nodes),
    sprintf("#SBATCH --ntasks=%s", config$ntasks),
    sprintf("#SBATCH --partition=%s", config$partition),
    sprintf("#SBATCH --qos=%s", config$qos),
    "",
    "module load R", # Load R module, adjust as needed
    "",
    "Rscript <<EOF",
    paste(c(
      "devtools::load_all()",
      "devtools::load_all(\"/p/tmp/hagento/dev/mredgebuildings\")",
      "library(magrittr)",
      "source('R/initCalculation.R')",
      "library(madrat)",
      sprintf("fileMapping <- data.frame(
      gcm = '%s',
      rcp = '%s',
      start = %s,
      end = %s,
      tas = '%s',
      rsds = '%s',
      sfc = '%s',
      huss = '%s',
      stringsAsFactors = FALSE)",
              fileRow$gcm, fileRow$rcp, fileRow$start, fileRow$end,
              fileRow$tas, fileRow$rsds, fileRow$sfc, fileRow$huss
      ),
      "result <- initCalculation(",
      "  fileMapping = fileMapping,",
      "  pop = pop,",
      sprintf("  ssp = '%s',", ssp),
      sprintf("  bait = %s,", bait),
      "  tLim = tLim,",
      "  countries = countries,",
      "  hddCddFactor = hddCddFactor,",
      "  wBAIT = wBait,",
      "  baitPars = baitPars",
      ")",
      sprintf("write.csv(result, '%s', row.names = FALSE)", outputFile)
    ), collapse = "\n"),
    "EOF"
  ), jobScript)

  # Create job name
  jobName <- paste0(
    config$jobNamePrefix,
    "_", fileRow$gcm,
    "_", fileRow$rcp
  )

  # Construct Slurm submission command with flexible configuration
  slurmCommand <- sprintf("sbatch %s", jobScript)

  # Submit job
  submissionResult <- tryCatch({
    system(slurmCommand, intern = TRUE)
  }, error = function(e) {
    stop("Failed to submit Slurm job: ", e$message)
  })

  # Extract job ID
  jobId <- str_extract(submissionResult, "\\d+")

  if (is.na(jobId)) {
    warning("Unable to extract job ID from submission result.")
  }

  list(
    jobName = jobName,
    jobScript = jobScript,
    outputFile = outputFile,
    slurmCommand = slurmCommand,
    jobId = jobId
  )
}



#' Wait for All Slurm Jobs to Complete
#'
#' @param checkInterval Interval in seconds between job status checks
#'
#' @return NULL
waitForSlurm <- function(checkInterval = 60) {
  repeat {
    jobStatus <- system("squeue -u $USER", intern = TRUE)
    if (length(jobStatus) <= 1) break
    Sys.sleep(checkInterval)
  }
}
