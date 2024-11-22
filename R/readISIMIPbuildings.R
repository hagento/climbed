#' Read relevant ISIMIP data for mredgebuildings
#'
#' Relevant data such as region masks, population and relevant climate data are
#' read in. The relevant file is declared in the subtype with the full file name.
#'
#' If the file name includes a suffix in the form of an integer such as
#' \code{_<int>.filetype}, the file is split into a single year period, e.g.
#' \code{<filename>_2001_2010_2.nc} will return data for the second year of the
#' 2001-2010 period, here 2002.
#'
#' @param subtype filename
#'
#' @author Hagen Tockhorn
#'
#' @importFrom stringr str_split
#' @importFrom terra rast subset aggregate ext res
#' @importFrom ncdf4 nc_open
#'
#' @note
#' folder structure in inputdata/sources/ISIMIPbuildings is expected to be:
#'    country masks : var/
#'    population :    var/scenario
#'    other :         var/scenario/model
#'
#' @note currently, this function only reads data from ISIMIP3b


readISIMIPbuildings <- function(subtype) {

  # PARAMETERS------------------------------------------------------------------

  baitVars <- c("tas", "sfcwind", "rsds", "huss")

  firstHistYear <- 1960


  # FUNCTIONS-------------------------------------------------------------------

  splitSubtype <- function(subtype) {
    vars <- list()

    # nolint start
    if (grepl("countrymask", subtype)) {
      vars[["variable"]] <- "countrymask"
    }

    else if (grepl("population", subtype)) {
      subSplit <- str_split(subtype, "_") %>% unlist()

      vars[["variable"]] <- subSplit[[1]]
      vars[["scenario"]] <- subSplit[[2]]
    }

    else if (any(sapply(baitVars, grepl, x = subtype))) {
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
      if (length(subSplit) > 9) {
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
    # nolint end

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


  # PROCESS DATA----------------------------------------------------------------

  vars <- splitSubtype(subtype)

  # region mask
  if (vars[["variable"]] == "countrymask") {
    fpath     <- file.path("countrymasks", subtype)
    varNames  <- names(nc_open(fpath)[["var"]])
    countries <- list()

    for (var in varNames) {
      countries[[var]] <- suppressWarnings(rast(fpath, subds = var))
    }

    r        <- rast(countries)
    names(r) <- gsub("m_", "", varNames)

    x <- list(x = r, class = "SpatRaster")
  }


  # population
  else if (vars[["variable"]] == "population") { #nolint
    fpath <- file.path(vars[["variable"]], vars[["scenario"]], subtype)

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

    x <- list(x = r, class = "SpatRaster", cache = FALSE)
  }


  # climate data
  else if (any(vars[["variable"]] %in% baitVars)) { #nolint
    # slice single years
    if (!is.null(vars[["yStart"]])) {
      fpath  <- file.path(vars[["variable"]], vars[["scenario"]], vars[["model"]], vars[["subtype"]])
      ranges <- getRanges(vars)

      r        <- suppressWarnings(rast(fpath, lyrs = ranges[["idxRange"]]))
      names(r) <- ranges[["yRange"]]
    }

    # full data set
    else { #nolint
      fpath <- file.path(vars[["variable"]], vars[["scenario"]], vars[["model"]], subtype)
      r <- suppressWarnings(rast(fpath))
    }

    x <- list(x = r, class = "SpatRaster", cache = FALSE)
  }

  else {stop("Subtype was incorrectly split or invalid subtype given.")} #nolint

  return(x)
}
