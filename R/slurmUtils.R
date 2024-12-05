#' Create a Slurm Job for Running initCalculation
#'
#' @param fileRow \code{data.frame} row containing file mapping details
#' @param pop \code{SpatRaster} Population data to be passed to initCalculation
#' @param ssp \code{character} Shared Socioeconomic Pathway scenario
#' @param bait \code{logical} indicating whether to use BAIT calculation
#' @param tLim \code{numeric} Temperature limits for degree day calculation
#' @param hddcddFactor \code{data.frame} Heating/Cooling Degree Day factors
#' @param wBAIT \code{numeric} Weights for BAIT calculation
#' @param jobConfig \code{list} of Slurm job configuration parameters
#' @param outDir \code{character} Absolute path to the output directory, containing logs/ and tmp/
#' @param globalPars \code{logical} indicating whether to use global or gridded BAIT parameters
#' (required if \code{bait} is TRUE).
#'
#' @returns \code{list} containing job details:
#'   - jobName: Name of the Slurm job
#'   - jobScript: Path to job script
#'   - outputFile: Path to output file
#'   - slurmCommand: Slurm submission command
#'   - jobId: Slurm job ID
#'   - batch_tag: Unique batch identifier
#'
#' @author Hagen Tockhorn
#'
#' @importFrom stringr str_extract
#' @importFrom piamutils getSystemFile
#' @importFrom terra writeCDF
#' @importFrom utils modifyList

createSlurm <- function(fileRow,
                        pop,
                        ssp,
                        bait,
                        tLim,
                        hddcddFactor,
                        wBAIT,
                        jobConfig = list(),
                        outDir = "output",
                        globalPars = FALSE) {
  # PARAMETERS -----------------------------------------------------------------

  # define default slurm job configuration
  defaultConfig <- list(
    logsDir      = file.path(outDir, "logs"),
    jobNamePrefix = "hddcdd",
    cpusPerTask   = 32,
    nodes         = 1,
    ntasks        = 1,
    partition     = "standard",
    qos           = "short"
  )



  # PROCESS DATA ---------------------------------------------------------------

  # match slurm job configs
  config <- modifyList(defaultConfig, jobConfig)

  # create directories output, output/tmp, output/logs and output/hddcdd
  tmpDir <- file.path(outDir, "tmp")
  dir.create(tmpDir, recursive = TRUE, showWarnings = FALSE)
  dir.create(config$logsDir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outDir, "hddcdd"), recursive = TRUE, showWarnings = FALSE)

  # create a unique tag for this batch of files
  batch_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # save files temporarily with a time tag to remove after successful processing
  pop_names <- names(pop)
  saveRDS(pop_names, sprintf("%s/pop_names_%s.rds", tmpDir, batch_tag))
  writeCDF(pop, sprintf("%s/pop_%s.nc", tmpDir, batch_tag), overwrite = TRUE)

  saveRDS(tLim, sprintf("%s/tLim_%s.rds", tmpDir, batch_tag))
  saveRDS(hddcddFactor, sprintf("%s/hddcddFactor_%s.rds", tmpDir, batch_tag))
  saveRDS(wBAIT, sprintf("%s/wBAIT_%s.rds", tmpDir, batch_tag))

  # output file
  outputFile <- file.path(outDir, "hddcdd",
                          paste0("hddcdd_", fileRow$gcm, "_", fileRow$rcp,
                                 "_", fileRow$start, "-", fileRow$end, ".csv"))

  # job name
  jobName <- paste0(config$jobNamePrefix, "_", fileRow$gcm, "_", fileRow$rcp,
                    "_", fileRow$start, "-", fileRow$end)

  # write slurm job script
  jobScript <- tempfile(pattern = "initcalc_job_", fileext = ".slurm")

  writeLines(c(
    "#!/bin/bash",
    sprintf("#SBATCH --job-name=%s", jobName),
    sprintf("#SBATCH --output=%s/log-%%j.out", config$logsDir),
    sprintf("#SBATCH --error=%s/errlog-%%j.out", config$logsDir),
    # sprintf("#SBATCH --output=%s/%s.out", config$logsDir, jobName),
    # sprintf("#SBATCH --error=%s/%s.err", config$logsDir, jobName),
    sprintf("#SBATCH --cpus-per-task=%s", config$cpusPerTask),
    sprintf("#SBATCH --nodes=%s", config$nodes),
    sprintf("#SBATCH --ntasks=%s", config$ntasks),
    sprintf("#SBATCH --partition=%s", config$partition),
    sprintf("#SBATCH --qos=%s", config$qos),
    "",
    "source /p/system/modulefiles/defaults/piam/module_load_piam",
    "",
    "R --no-save <<EOF",
    "library(devtools)",
    sprintf("devtools::load_all(\"/p/tmp/hagento/dev/climbed\")"),   # not sure how to properly manage this yet
    "",
    "# Load and restore raster with names",
    sprintf("pop <- terra::rast('%s/pop_%s.nc')", tmpDir, batch_tag),
    sprintf("pop_names <- readRDS('%s/pop_names_%s.rds')", tmpDir, batch_tag),
    "names(pop) <- pop_names",
    "",
    sprintf("tLim <- readRDS('%s/tLim_%s.rds')", tmpDir, batch_tag),
    sprintf("hddcddFactor <- readRDS('%s/hddcddFactor_%s.rds')", tmpDir, batch_tag),
    sprintf("wBAIT <- readRDS('%s/wBAIT_%s.rds')", tmpDir, batch_tag),
    "",
    sprintf("fileMapping <- data.frame(
      gcm = '%s',
      rcp = '%s',
      start = %s,
      end = %s,
      tas = '%s',
      rsds = '%s',
      sfcwind = '%s',
      huss = '%s',
      stringsAsFactors = FALSE)",
      fileRow$gcm, fileRow$rcp, fileRow$start, fileRow$end,
      fileRow$tas, fileRow$rsds, fileRow$sfcwind, fileRow$huss
    ),
    "",
    "result <- initCalculation(",
    "  fileMapping = fileMapping,",
    sprintf("  ssp = '%s',", ssp),
    sprintf("  bait = '%s',", bait),
    "  tLim = tLim,",
    "  pop = pop,",
    "  hddcddFactor = hddcddFactor,",
    "  wBAIT = wBAIT,",
    sprintf("  globalPars = '%s'", globalPars),
    ")",
    "",
    sprintf("write.csv(result, '%s', row.names = FALSE)", outputFile),
    "",
    "# Clean up all temporary files",
    sprintf("unlink(list.files('%s', pattern='%s', full.names=TRUE))", tmpDir, batch_tag),
    "EOF"
  ), jobScript)


  # submit slurm job
  slurmCommand <- sprintf("sbatch %s", jobScript)
  submissionResult <- system(slurmCommand, intern = TRUE)
  jobId <- str_extract(submissionResult, "\\d+")


  # return relevant job information
  return(invisible(list(jobName = jobName,
                        jobScript = jobScript,
                        outputFile = outputFile,
                        slurmCommand = slurmCommand,
                        jobId = jobId,
                        batch_tag = batch_tag)))
}



#' Wait for SLURM Jobs to Complete
#'
#' Monitors the status of SLURM jobs and waits for their completion. The function
#' periodically checks the status of specified jobs using \code{sacct} and handles
#' different job states including failures, timeouts, and successful completions.
#'
#' @param jobIds A vector of SLURM job IDs to monitor. Can be numeric or character.
#' @param checkInterval Number of seconds to wait between status checks (default: 60).
#' @param maxWaitTime Maximum time in seconds to wait for job completion (default: 3600).
#'
#' @returns Returns TRUE if all jobs complete successfully. Stops with an error if:
#'   - Any job fails (FAILED, CANCELLED, TIMEOUT, OUT_OF_MEMORY, NODE_FAIL)
#'   - Maximum wait time is exceeded
#'
#' @importFrom utils read.table

waitForSlurm <- function(jobIds, checkInterval = 60, maxWaitTime = 3600) {
  startTime <- Sys.time()
  jobIds <- as.character(jobIds)
  jobSet <- unique(jobIds)  # Remove duplicates

  while (difftime(Sys.time(), startTime, units = "secs") < maxWaitTime) {
    # Get detailed job status including job steps
    jobsCommand <- sprintf("sacct -j %s --parsable2 --noheader --format=jobid,state",
                           paste(jobSet, collapse = ","))
    jobsStatus <- system(jobsCommand, intern = TRUE)

    if (length(jobsStatus) == 0) {
      stop("No job information found. Check if jobs exist.")
    }

    # Process status data
    statusMatrix <- do.call(rbind, strsplit(jobsStatus, "|", fixed = TRUE))
    parentJobs <- statusMatrix[!grepl("\\.", statusMatrix[, 1]), , drop = FALSE]

    # Check for failed states
    failedStates <- c("FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY",
                      "NODE_FAIL", "PREEMPTED", "DEADLINE")
    failedJobs <- parentJobs[parentJobs[, 2] %in% failedStates, , drop = FALSE]

    if (nrow(failedJobs) > 0) {
      failMessage <- sprintf("Jobs failed: %s",
                             paste(sprintf("JobID %s: %s",
                                           failedJobs[, 1], failedJobs[, 2]),
                                   collapse = ", "))
      stop(failMessage)
    }

    # Check for active states
    activeStates <- c("PENDING", "RUNNING", "COMPLETING", "CONFIGURING",
                      "SUSPENDED", "REQUEUED", "RESIZING")
    activeJobs <- parentJobs[parentJobs[, 2] %in% activeStates, , drop = FALSE]

    # Verify completion
    if (nrow(activeJobs) == 0) {
      completedJobs <- parentJobs[parentJobs[, 2] == "COMPLETED", , drop = FALSE]
      if (nrow(completedJobs) == length(jobSet)) {
        message(sprintf("All %d jobs completed successfully", length(jobSet)))
        return(TRUE)
      } else {
        # Some jobs are missing or in unexpected states
        unknownJobs <- setdiff(jobSet, parentJobs[, 1])
        if (length(unknownJobs) > 0) {
          stop(sprintf("Jobs with unknown status: %s",
                       paste(unknownJobs, collapse = ", ")))
        }
        stop("Some jobs are in unexpected states. Check sacct manually.")
      }
    }

    Sys.sleep(checkInterval)
  }

  if (difftime(Sys.time(), startTime, units = "secs") > maxWaitTime) {
    stop("Maximum wait time exceeded")
  }
}



#' Get Main Job ID
#' @returns Character string containing the main job ID or NULL if not in a Slurm job
getMainJobId <- function() {
  slurm_job_id <- Sys.getenv("SLURM_JOB_ID")
  if (slurm_job_id == "") return(NULL)
  return(slurm_job_id)
}
