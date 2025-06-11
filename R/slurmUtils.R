#' Create a Slurm Job for Running initCalculation or computeBAITpars
#'
#' @param subtype \code{character} Type of calculation: "initCalculation" or "computeBAITpars"
#' @param fileRow \code{data.frame} row containing file mapping details (required for initCalculation)
#' @param pop \code{SpatRaster} Population data to be passed to initCalculation (required for initCalculation)
#' @param ssp \code{character} Shared Socioeconomic Pathway scenario (required for initCalculation)
#' @param bait \code{logical} indicating whether to use BAIT calculation (required for initCalculation)
#' @param tLim \code{numeric} Temperature limits for degree day calculation (required for initCalculation)
#' @param hddcddFactor \code{data.frame} Heating/Cooling Degree Day factors (required for initCalculation)
#' @param wBAIT \code{numeric} Weights for BAIT calculation (required for initCalculation)
#' @param jobConfig \code{list} of Slurm job configuration parameters
#' @param outDir \code{character} Absolute path to the output directory, containing logs/ and tmp/
#' @param globalPars \code{logical} indicating whether to use global or gridded BAIT parameters
#' (required if \code{bait} is TRUE).
#' @param noCC \code{logical} indicating whether to compute a no-climate-change scenario.
#'        If \code{TRUE}, the function will calculate degree days assuming constant climate conditions.
#'        Default is \code{FALSE}.
#' @param packagePath \code{character} (Optional) Path to the package for development mode.
#' If provided, the SLURM jobs will use devtools::load_all() with this path.
#' If NULL (default), the installed climbed package will be loaded.
#' @param runTag A string specifying the unique batch tag for temporary files
#'
#' @returns \code{list} containing job details:
#'   - jobName: Name of the Slurm job
#'   - jobScript: Path to job script
#'   - outputFile: Path to output file (for initCalculation only)
#'   - slurmCommand: Slurm submission command
#'   - jobId: Slurm job ID
#'   - batchTag: Unique batch identifier
#'
#' @author Hagen Tockhorn
#'
#' @importFrom stringr str_extract
#' @importFrom piamutils getSystemFile
#' @importFrom terra writeCDF
#' @importFrom utils modifyList
#' @importFrom digest digest
#' @importFrom stats runif

createSlurm <- function(subtype,
                        fileRow = NULL,
                        pop = NULL,
                        ssp = NULL,
                        bait = NULL,
                        tLim = NULL,
                        hddcddFactor = NULL,
                        wBAIT = NULL,
                        jobConfig = list(),
                        outDir = "output",
                        globalPars = FALSE,
                        noCC = FALSE,
                        packagePath = NULL,
                        runTag = "") {

  # PARAMETERS -----------------------------------------------------------------

  # determine required memory based on subtype
  mem <- if (subtype == "computeBAITpars") {
    "200G"
  } else {
    ifelse(isTRUE(bait), ifelse(isTRUE(globalPars), "64G", "128G"), "32G")
  }

  # define default slurm job configuration
  defaultConfig <- list(
    logsDir      = file.path(outDir, "logs"),
    jobNamePrefix = if (subtype == "computeBAITpars") "computeBAITpars" else "hddcdd",
    nodes         = 1,
    ntasks        = 1,
    mem           = mem,
    partition     = "standard",
    qos           = "short",
    SLURMtmpDir   = "slurm_tmp"  # Default temporary directory for slurm jobs
  )

  # PROCESS DATA ---------------------------------------------------------------

  # match slurm job configs
  config <- modifyList(defaultConfig, jobConfig)

  # create basic directories
  tmpDir <- file.path(outDir, "tmp")
  dir.create(tmpDir, recursive = TRUE, showWarnings = FALSE)
  dir.create(config$logsDir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(tmpDir, config$SLURMtmpDir), recursive = TRUE, showWarnings = FALSE)

  # Create a unique tag with timestamp + hash
  batchTag <- paste0(
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    "_",
    substr(digest::digest(runif(1)), 1, 6)  # 6-character hash
  )

  # initCalculation-specific setup
  if (subtype == "initCalculation") {
    # create directories
    dir.create(file.path(outDir, "hddcdd"), recursive = TRUE, showWarnings = FALSE)

    # create output directory for grid data in noCC-case
    if (isTRUE(noCC)) {
      gridDataDir <- file.path(outDir, "hddcddGrid")
      dir.create(gridDataDir, recursive = TRUE, showWarnings = FALSE)
    }

    # save files temporarily
    saveRDS(tLim, sprintf("%s/tLim_%s.rds", tmpDir, batchTag))
    saveRDS(hddcddFactor, sprintf("%s/hddcddFactor_%s.rds", tmpDir, batchTag))
    saveRDS(wBAIT, sprintf("%s/wBAIT_%s.rds", tmpDir, batchTag))

    # define output file and job name
    outputFile <- file.path(outDir, "hddcdd",
                            paste0("hddcdd_", fileRow$gcm, "_", ssp, "_", fileRow$rcp,
                                   "_", fileRow$start, "-", fileRow$end, "_", runTag, ".csv"))
    jobName <- paste0(config$jobNamePrefix, "_", fileRow$gcm, "_", fileRow$rcp,
                      "_", fileRow$start, "-", fileRow$end)
  } else {
    # computeBAITpars job name
    jobName <- paste0(config$jobNamePrefix, "_", batchTag)
  }

  # Create job script directory
  scriptDir <- file.path(config$SLURMtmpDir, "job_scripts")
  dir.create(scriptDir, recursive = TRUE, showWarnings = FALSE)

  # Create separate R script file
  rScriptPath <- file.path(scriptDir, paste0("r_script_", batchTag, ".R"))
  jobScript <- file.path(scriptDir, paste0("job_", batchTag, ".slurm"))

  # determine package loading approach
  if (!is.null(packagePath)) {
    # use provided package path with load_all
    packageLoadingCode <- sprintf("devtools::load_all(\"%s\")", packagePath)
  } else {
    # use standard library loading - assuming package is installed
    packageLoadingCode <- "library(climbed)"
  }

  # Write the R script file based on subtype
  if (subtype == "computeBAITpars") {
    rScriptContent <- c(
      "library(devtools)",
      packageLoadingCode,
      "",
      "computeBAITpars()"
    )
  } else {
    # initCalculation R script
    rScriptContent <- c(
      "library(devtools)",
      packageLoadingCode,
      "",
      sprintf("tLim <- readRDS('%s/tLim_%s.rds')", tmpDir, batchTag),
      sprintf("hddcddFactor <- readRDS('%s/hddcddFactor_%s.rds')", tmpDir, batchTag),
      sprintf("wBAIT <- readRDS('%s/wBAIT_%s.rds')", tmpDir, batchTag),
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
      sprintf("  bait = %s,", bait),
      "  tLim = tLim,",
      sprintf("  pop = '%s',", pop),
      "  hddcddFactor = hddcddFactor,",
      "  wBAIT = wBAIT,",
      sprintf("  globalPars = %s,", globalPars),
      sprintf("  noCC = %s,", noCC),
      if (exists("gridDataDir")) sprintf("  gridDataDir = '%s'", gridDataDir) else "  gridDataDir = NULL",
      ")",
      "",
      sprintf("write.csv(result, '%s', row.names = FALSE)", outputFile),
      "",
      "# Clean up all temporary files",
      sprintf("unlink(list.files('%s', pattern='%s', full.names=TRUE))", tmpDir, batchTag)
    )
  }

  writeLines(rScriptContent, rScriptPath)

  # Write the SLURM job script without using here-document
  writeLines(c(
    "#!/bin/bash",
    sprintf("#SBATCH --job-name=%s", jobName),
    sprintf("#SBATCH --output=%s/log-%%j.out", config$logsDir),
    sprintf("#SBATCH --error=%s/errlog-%%j.out", config$logsDir),
    sprintf("#SBATCH --nodes=%s", config$nodes),
    sprintf("#SBATCH --mem=%s", config$mem),
    sprintf("#SBATCH --ntasks=%s", config$ntasks),
    sprintf("#SBATCH --partition=%s", config$partition),
    sprintf("#SBATCH --qos=%s", config$qos),
    "",
    "source /p/system/modulefiles/defaults/piam/module_load_piam",
    "",
    # Execute the R script with Rscript instead of using here-document
    sprintf("Rscript %s", rScriptPath)
  ), jobScript)

  # submit slurm job
  slurmCommand <- sprintf("sbatch %s", jobScript)
  submissionResult <- system(slurmCommand, intern = TRUE)
  jobId <- str_extract(submissionResult, "\\d+")

  # return relevant job information
  result <- list(jobName = jobName,
                 jobScript = jobScript,
                 rScriptPath = rScriptPath,
                 slurmCommand = slurmCommand,
                 jobId = jobId,
                 batchTag = batchTag)

  # Add initCalculation-specific fields
  if (subtype == "initCalculation") {
    result$outputFile <- outputFile
    if (exists("gridDataDir")) result$gridDataDir <- gridDataDir
  }

  return(invisible(result))
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
  message("Monitoring jobs: ", paste(jobSet, collapse = ", "))
  while (difftime(Sys.time(), startTime, units = "secs") < maxWaitTime) {
    # Try squeue first (for running/queued jobs)
    squeueCmd <- sprintf("squeue -j %s -h -o '%%A|%%T'", paste(jobSet, collapse = ","))
    squeueOut <- system(squeueCmd, intern = TRUE, ignore.stderr = TRUE)
    if (length(squeueOut) > 0) {
      Sys.sleep(checkInterval)
      next
    }
    # Check if job exists in accounting
    sacctCmd <- sprintf("sacct -j %s -X -n -o jobid,state", paste(jobSet, collapse = ","))
    sacctOut <- system(sacctCmd, intern = TRUE, ignore.stderr = TRUE)
    # If we got output from sacct
    if (length(sacctOut) > 0) {
      # Check for completed jobs
      allCompleted <- all(vapply(jobSet, function(id) {
        any(grepl(paste0(id, ".*COMPLETED"), sacctOut, ignore.case = TRUE))
      }, logical(1)))
      if (allCompleted) {
        message("All jobs completed successfully")
        return(TRUE)
      }
      # Check for failed jobs
      failedJobs <- character(0)
      failedStates <- c("FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY",
                        "NODE_FAIL", "PREEMPTED", "DEADLINE")
      for (id in jobSet) {
        for (state in failedStates) {
          if (any(grepl(paste0(id, ".*", state), sacctOut, ignore.case = TRUE))) {
            failedJobs <- c(failedJobs, paste0(id, ": ", state))
          }
        }
      }
      if (length(failedJobs) > 0) {
        stop("Jobs failed: ", paste(failedJobs, collapse = ", "))
      }
      # If we get here, jobs are in some other state
      message("Jobs in progress, waiting...")
    } else {
      # Try scontrol as an alternative
      sctrlCmd <- sprintf("scontrol show job %s", paste(jobSet, collapse = " "))
      sctrlOut <- system(sctrlCmd, intern = TRUE, ignore.stderr = TRUE)
      if (length(sctrlOut) > 0 && !any(grepl("Invalid job id", sctrlOut))) {
        if (any(grepl("JobState=COMPLETED", sctrlOut))) {
          message("All jobs completed successfully")
          return(TRUE)
        }
        message("Jobs found but not completed, waiting...")
      } else {
        message("No job information found. Will retry.")
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
