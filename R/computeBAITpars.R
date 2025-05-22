#' Calculate Global Regression Parameters for BAIT Climate Variables
#'
#' @description Performs linear regression on historical climate data to calculate
#' global regression parameters for key climate variables: surface downwelling shortwave
#' radiation (rsds), near-surface wind speed (sfcwind), and near-surface specific
#' humidity (huss), with respect to near-surface air temperature (tas).
#'
#' The regression utilizes a simple linear model based on historical data from
#' 2000 to 2019-21 across all available models. For rsds and sfcwind, a linear relationship
#' is assumed, while for huss, an exponential relationship is adopted, reflecting
#' the non-linear dependence of water vapor pressure on temperature.
#'
#' @param cacheDir \code{character} string specifying the directory containing pre-calculated
#'   function outputs.
#'
#' @returns \code{terra::SpatRaster} object with layers representing regression parameters for each cell.
#'
#' @author Hagen Tockhorn
#'
#' @importFrom terra regress rast
#' @importFrom madrat toolGetMapping readSource
#' @importFrom magclass as.magpie
#' @importFrom utils read.csv2
#' @importFrom piamutils getSystemFile

computeBAITpars <- function(cacheDir = "intdata/BAITpars") {
  # CHECK CACHE-----------------------------------------------------------------
  # absolut package path
  pkgPath <- getSystemFile(package = "climbed")
  if (!is.null(cacheDir)) {
    fpath <- list.files(file.path(pkgPath, cacheDir), pattern = "globalBaitPars", full.names = TRUE) %>%
      unlist()
    if (length(fpath) > 0 && file.exists(fpath)) {
      print(paste0("Load global parameters from cache: ", basename(fpath)))
      regPars <- rast(fpath)
      return(regPars)
    } else {
      print("Global BAITpars file not in given cache directory - will be re-calculated.")
    }
  }

  # READ-IN DATA----------------------------------------------------------------
  # Read file mapping
  fileMapping <- read.csv2(getSystemFile("extdata", "mappings", "BAITpars_fileMapping.csv", package = "climbed"))

  vars <- list("tas" = "tas", "rsds" = "RSDS", "sfcwind" = "SFC", "huss" = "HUSS")

  # Create empty list to store data
  allData <- setNames(vector("list", length(names(vars))), names(vars))

  # Load all data at once instead of model by model
  for (v in names(vars)) {
    # Get all non-NA file paths for this variable across all models
    varFiles <- fileMapping[[v]]
    varFiles <- varFiles[!is.na(varFiles)]

    if (length(varFiles) > 0) {
      # Import data for each file
      dataList <- list()
      for (f in varFiles) {
        tryCatch({
          importedData <- importData(subtype = f)
          dataList <- c(dataList, list(importedData))
        }, error = function(e) {
          warning(paste("Error importing file:", f, "-", e$message))
        })
      }

      # Only create raster if we have data
      if (length(dataList) > 0) {
        allData[[v]] <- rast(dataList)
      }
    }
  }

  # PROCESS DATA----------------------------------------------------------------
  # convert huss into log scale
  allData$huss <- log(allData$huss)
  # convert tas into [C]
  allData$tas <- allData$tas - 273.15

  regPars <- do.call(c, lapply(names(vars)[names(vars) != "tas"], function(v) {
    print(v)
    x <- allData[["tas"]]
    y <- allData[[v]]
    r <- regress(x = x, y = y, formula = y ~ x)
    names(r) <- c(paste0("a", vars[v]), paste0("b", vars[v]))
    return(r)
  }))

  # Save to cache
  if (!is.null(cacheDir)) {
    dir.create(file.path(pkgPath, cacheDir), showWarnings = FALSE, recursive = TRUE)
    outputPath <- file.path(pkgPath, cacheDir, "globalBaitPars.tif")
    writeRaster(regPars, outputPath, overwrite = TRUE)
    print(paste0("Saved global parameters to cache: ", basename(outputPath)))
  }

  # OUTPUT----------------------------------------------------------------------
  return(regPars)
}
