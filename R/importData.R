#' Read and Pre-process Input Data
#'
#' Reads and processes input data such as region masks, population, and climate data.
#' The file to process is specified via \code{subtype}, which includes the full file name.
#'
#' Files with suffixes in the format \code{_<int>.filetype} are split by year.
#' For example, \code{<filename>_2001_2010_2.nc} returns data for 2002 (the second year in the range).
#'
#' @param subtype Character string specifying the file name.
#'
#' @return Pre-processed dataset based on \code{subtype}.
#'
#' @details
#' The function supports data from ISIMIP3a and ISIMIP3b. It expects the following folder structure
#' under the base directory \code{inputdata/sources/ISIMIPbuildings}:
#' \itemize{
#'   \item \strong{Country masks}: \code{inputdata/sources/ISIMIPbuildings/countrymasks/}
#'   \item \strong{Population data}: \code{inputdata/sources/ISIMIPbuildings/population/<scenario>/}
#'   \item \strong{Other data (e.g., tas, sfcwind)}:
#'   \code{inputdata/sources/ISIMIPbuildings/otherdata/<scenario>/<model>/}
#' }
#'
#' Note that \code{<scenario>} and \code{<model>} are placeholders that need to be replaced
#' with the actual scenario and model names provided as part of the dataset.
#'
#' Ensure that the directory structure matches this expected format to avoid errors.
#'
#' @author Hagen Tockhorn
#'
#' @importFrom stringr str_split str_detect
#' @importFrom terra rast subset aggregate ext res
#' @importFrom ncdf4 nc_open
#' @importFrom madrat getConfig
#' @importFrom utils tail
#'
#' @note The folder structure and correct placeholders for \code{<scenario>} and
#' \code{<model>} are crucial for successful data loading.



importData <- function(subtype) {
  # SOURCE DIRECTORY -----------------------------------------------------------

  # check if source folder exists
  sourceDir <- getConfig("sourcefolder")
  if (!dir.exists(sourceDir)) {
    stop("No madrat source directory found. Please define with `madrat::setConfig(\"sourcefolder\") = ...`")
  }

  sourceDir <- file.path(sourceDir, "ISIMIPbuildings")



  # PARAMETERS -----------------------------------------------------------------

  baitVars <- c("tas", "sfcwind", "rsds", "huss")

  firstHistYear <- 1960



  # FUNCTIONS ------------------------------------------------------------------

  splitSubtype <- function(subtype) {
    vars <- list()

    if (grepl("countrymask", subtype)) {
      vars[["variable"]] <- "countrymask"
    } else if (grepl("population", subtype)) {
      subSplit <- str_split(subtype, "_") %>% unlist()

      vars[["variable"]] <- subSplit[[1]]
      vars[["scenario"]] <- subSplit[[2]]
    } else if (any(vapply(baitVars, grepl, logical(1), x = subtype))) {
      subSplit <- str_split(subtype, "_") %>% unlist()

      # observations have a shorter file name
      if (grepl("obsclim", subtype)) {
        vars[["variable"]] <- subSplit[[3]]
        vars[["scenario"]] <- "historical"
      } else {
        vars[["variable"]] <- subSplit[[5]]
        vars[["scenario"]] <- subSplit[[4]]
      }

      vars[["model"]]    <- subSplit[[1]]

      # raster data will be split into individual years
      if (length(subSplit) > 9 || (length(subSplit) > 7 && grepl("obsclim", subtype))) {
        # split index defines the year
        vars[["idx"]] <- gsub(".nc", "", tail(subSplit, 1)) %>%
          as.numeric()

        # temporal range of data
        if (grepl("obsclim", subtype)) {
          vars[["yStart"]] <- subSplit[[6]]
          vars[["yEnd"]]   <- subSplit[[7]]
        } else {
          vars[["yStart"]] <- subSplit[[8]]
          vars[["yEnd"]]   <- subSplit[[9]]
        }

        # year of interest
        vars[["year"]] <- seq(vars[["yStart"]] %>%
                                as.numeric(),
                              vars[["yEnd"]] %>%
                                as.numeric())[[vars[["idx"]]]] %>%
          as.character()

        vars[["subtype"]] <- sub("_(1[0-9]|\\d)\\.nc$", ".nc", subtype)
      }
    } else {
      stop("Invalid subtype given.")
    }

    return(vars)
  }


  # determine period range of subset
  getRanges <- function(vars) {
    # total temporal range
    dRange <- seq.Date(from = as.Date(paste0(vars[["yStart"]], "-01-01")),
                       to = as.Date(paste0(vars[["yEnd"]],   "-12-31")),
                       by = "day")

    # temporal range of interest
    yRange <- seq.Date(from = as.Date(paste0(vars[["year"]], "-01-01")),
                       to   = as.Date(paste0(vars[["year"]], "-12-31")),
                       by   = "day")

    # indices of range of interest
    idxRange <- match(yRange, dRange)

    return(list("yRange"   = yRange,
                "idxRange" = idxRange))
  }


  # fill dates for unnamed data
  fillDates <- function(r, filename, pop = FALSE) {

    yStart <- stringr::str_sub(filename, -12, -9)
    n <- terra::nlyr(r)

    if (!pop) {
      dStart <- as.Date(paste0(yStart, "-1-1"))
      dates <- seq.Date(dStart, by = "day", length.out = n)
    } else {
      dates <- seq(yStart, by = 1, length.out = n) %>%
        as.character()
    }

    # fill dates
    names(r) <- dates
    return(r)
  }



  # PROCESS DATA ---------------------------------------------------------------

  vars <- splitSubtype(subtype)

  # region mask
  if (vars[["variable"]] == "countrymask") {
    fpath     <- file.path(sourceDir, "countrymasks", subtype)
    varNames  <- names(nc_open(fpath)[["var"]])
    countries <- list()

    for (var in varNames) {
      countries[[var]] <- suppressWarnings(rast(fpath, subds = var))
    }

    r        <- rast(countries)
    names(r) <- gsub("m_", "", varNames)

    x <- r

  } else if (vars[["variable"]] == "population") {
    fpath <- file.path(sourceDir, vars[["variable"]], vars[["scenario"]], subtype)

    if (vars[["scenario"]] == "picontrol") {
      r <- suppressWarnings(rast(fpath))
    } else {
      r <- suppressWarnings(rast(fpath, subds = "total-population"))
    }

    subtype <- gsub(".nc", "", subtype)

    # rename years
    years    <- tail(strsplit(subtype, "_")[[1]], 2)
    names(r) <- years[1]:years[2]

    # filter relevant years
    r <- subset(r, as.numeric(names(r)) >= firstHistYear)

    # aggregate to common resolution of 0.5 deg
    if (any(res(r) != 0.5)) {
      r <- aggregate(r, fun = "sum", fact = round(0.5 / res(r), 3))
    }

    x <- r

  } else if (any(vars[["variable"]] %in% baitVars)) {
    # slice single years
    if (!is.null(vars[["yStart"]])) {
      fpath  <- file.path(sourceDir, vars[["variable"]], vars[["scenario"]], vars[["model"]], vars[["subtype"]])
      ranges <- getRanges(vars)

      r        <- suppressWarnings(rast(fpath, lyrs = ranges[["idxRange"]]))
      names(r) <- ranges[["yRange"]]

    } else {
      fpath <- file.path(sourceDir, vars[["variable"]], vars[["scenario"]], vars[["model"]], subtype)
      r <- suppressWarnings(rast(fpath))
    }

    if (!all(nchar(names(r)) == 10)) {
      # dates have specific length of n = 10
      r <- fillDates(r, subtype)
    }

    if (vars[["variable"]] == "huss") {
      # convert kg/kg -> g/kg
      r <- r * 1e3
    }

    x <- r

  } else {
    stop("Subtype was incorrectly split or invalid subtype given.")
  }

  return(x)
}
