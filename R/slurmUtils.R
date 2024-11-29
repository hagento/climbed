#### SLURM JOB UTILITIES ####


#' Create a Slurm Job for Running initCalculation
#'
#' @param fileRow A data frame row containing file mapping details
#' @param pop Population data to be passed to initCalculation
#' @param ssp Shared Socioeconomic Pathway scenario
#' @param bait Logical indicating whether to use BAIT calculation
#' @param tLim Temperature limits for degree day calculation
#' @param hddCddFactor Heating/Cooling Degree Day factors
#' @param wBAIT Weights for BAIT calculation
#' @param jobConfig List of Slurm job configuration parameters. Currently supports
#' @param outDir A string specifying the absolute path to the \code{output} directory. This will also
#' contain the \code{logs} and \code{tmp} directories for log and temporary files.
#'
#' @return A list containing job details including job name, script path,
#'         output file path, and Slurm submission command
#'
#' @author Hagen Tockhorn
#'
#' @importFrom stringr str_extract
#' @importFrom piamutils getSystemFile

createSlurm <- function(fileRow,
                        pop,
                        ssp,
                        bait,
                        tLim,
                        hddcddFactor,
                        wBAIT,
                        jobConfig = list(),
                        outDir = "output") {
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

  # create temporary file and log file directory
  tmpDir <- file.path(outDir, "tmp")
  dir.create(tmpDir, recursive = TRUE, showWarnings = FALSE)
  dir.create(config$logsDir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outDir, "hddcdd"), recursive = TRUE, showWarnings = FALSE)

  # create a unique tag for this batch of files
  batch_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # save files temporarily with the unique tag
  pop_names <- names(pop)
  saveRDS(pop_names, sprintf("%s/pop_names_%s.rds", tmpDir, batch_tag))
  terra::writeCDF(pop, sprintf("%s/pop_%s.nc", tmpDir, batch_tag), overwrite = TRUE)
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
    sprintf("#SBATCH --output=%s/%s.out", config$logsDir, jobName),
    sprintf("#SBATCH --error=%s/%s.err", config$logsDir, jobName),
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
    "  wBAIT = wBAIT",
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
  return(list(jobName = jobName,
              jobScript = jobScript,
              outputFile = outputFile,
              slurmCommand = slurmCommand,
              jobId = jobId,
              batch_tag = batch_tag))
}



#' Wait for SLURM Jobs to Complete
#'
#' Monitors the status of SLURM jobs and waits for their completion. The function
#' periodically checks the status of specified jobs using `sacct` and handles
#' different job states including failures, timeouts, and successful completions.
#'
#' @param jobIds A vector of SLURM job IDs to monitor. Can be numeric or character.
#' @param checkInterval Number of seconds to wait between status checks (default: 60).
#' @param maxWaitTime Maximum time in seconds to wait for job completion (default: 3600).
#'
#' @return Returns TRUE if all jobs complete successfully. Stops with an error if:
#'   - Any job fails (FAILED, CANCELLED, TIMEOUT, OUT_OF_MEMORY, NODE_FAIL)
#'   - Maximum wait time is exceeded
#'
#' @importFrom utils read.table

waitForSlurm <- function(jobIds, checkInterval = 60, maxWaitTime = 3600) {
  startTime <- Sys.time()
  jobIds <- as.character(jobIds)
  jobSet <- unique(jobIds)  # Remove duplicates

  while (TRUE) {
    # Get detailed job status including job steps
    jobs_cmd <- sprintf("sacct -j %s --parsable2 --noheader --format=jobid,state",
                        paste(jobSet, collapse = ","))
    jobs_status <- system(jobs_cmd, intern = TRUE)

    if (length(jobs_status) == 0) {
      stop("No job information found. Check if jobs exist.")
    }

    # Process status data
    status_matrix <- do.call(rbind, strsplit(jobs_status, "|", fixed = TRUE))
    parent_jobs <- status_matrix[!grepl("\\.", status_matrix[,1]), , drop = FALSE]

    # Check for failed states
    failed_states <- c("FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY",
                       "NODE_FAIL", "PREEMPTED", "DEADLINE")
    failed_jobs <- parent_jobs[parent_jobs[,2] %in% failed_states, , drop = FALSE]

    if (nrow(failed_jobs) > 0) {
      failed_msg <- sprintf("Jobs failed: %s",
                            paste(sprintf("JobID %s: %s",
                                          failed_jobs[,1], failed_jobs[,2]),
                                  collapse = ", "))
      stop(failed_msg)
    }

    # Check for active states
    active_states <- c("PENDING", "RUNNING", "COMPLETING", "CONFIGURING",
                       "SUSPENDED", "REQUEUED", "RESIZING")
    active_jobs <- parent_jobs[parent_jobs[,2] %in% active_states, , drop = FALSE]

    # Verify completion
    if (nrow(active_jobs) == 0) {
      completed_jobs <- parent_jobs[parent_jobs[,2] == "COMPLETED", , drop = FALSE]
      if (nrow(completed_jobs) == length(jobSet)) {
        message(sprintf("All %d jobs completed successfully", length(jobSet)))
        return(TRUE)
      } else {
        # Some jobs are missing or in unexpected states
        unknown_jobs <- setdiff(jobSet, parent_jobs[,1])
        if (length(unknown_jobs) > 0) {
          stop(sprintf("Jobs with unknown status: %s",
                       paste(unknown_jobs, collapse = ", ")))
        }
        stop("Some jobs are in unexpected states. Check sacct manually.")
      }
    }

    # Check timeout
    if (difftime(Sys.time(), startTime, units = "secs") > maxWaitTime) {
      stop("Maximum wait time exceeded")
    }

    Sys.sleep(checkInterval)
  }
}



#' Get Main Job ID
#' @return Character string containing the main job ID or NULL if not in a Slurm job
getMainJobId <- function() {
  slurm_job_id <- Sys.getenv("SLURM_JOB_ID")
  if (slurm_job_id == "") return(NULL)
  return(slurm_job_id)
}
