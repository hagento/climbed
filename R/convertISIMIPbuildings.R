#' convert ISIMIP data for mredgebuildings
#'
#' @param x MAgPIE object on cellular level
#' @param subtype filename
#'
#' @importFrom stringr str_detect
#'
#' @author Hagen Tockhorn
#'
#' @return rasterBrick object

convertISIMIPbuildings <- function(x, subtype) {

  # FUNCTIONS-------------------------------------------------------------------

  # fill dates for unnamed data
  fillDates <- function(r, filename, pop = FALSE) {

    yStart <- stringr::str_sub(filename, -12, -9)
    n <- terra::nlyr(r)

    if (!pop) {
      dStart <- as.Date(paste0(yStart, "-1-1"))
      dates <- seq.Date(dStart, by = "day", length.out = n)
    } else {
      print("pop")
      dates <- seq(yStart, by = 1, length.out = n) %>% as.character()
    }

    # fill dates
    names(r) <- dates
    return(r)
  }


  # PARAMETERS------------------------------------------------------------------

  # variable units
  unitMapping <- list(
    "tas"     = "K",
    "rsds"    = "Wm-2",
    "sfcWind" = "ms-1",
    "huss"    = "gkg-1"
  )

  edgeVars <- c("tas", "rsds", "sfcwind", "huss", "population")


  # PROCESS DATA----------------------------------------------------------------

  if (grepl(paste(edgeVars, collapse = "|"), subtype)) {
    var <- edgeVars[str_detect(subtype, edgeVars)]


    if (var == "population") {
      x <- fillDates(x, subtype, pop = TRUE)
    } else if (!all(nchar(names(x)) == 10)) {
      # dates have specific length of n = 10
      x <- fillDates(x, subtype)
    }

    if (var == "huss") {
      # convert kg/kg -> g/kg
      x <- x * 1e3
    }
  }

  return(list(x = x,
              class = "SpatRaster",
              unit = unitMapping[[var]],
              cache = FALSE))
}
