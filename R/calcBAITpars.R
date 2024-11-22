#' calculate regression parameters for BAIT climate variables
#'
#' @description linear regression on historic data to determine regression
#' parameters for surface downdwelling shortwave radiation (rsds), near-surface
#' wind speed (sfcwind) and near-surface specific humidity (huss) with respect
#' to near-surface air temperature (tas).
#' The regression is done with on a simple linear model, where the historical
#' input data covers the years of 2000-2014. For rsds and sfcwind, a simple linear
#' relationship is assumed where for huss an exponential relation is assumed,
#' buildig upon the non-linear relation between water vapor pressure and temperature.
#'
#' @param model specify GCM responsible for data input
#' @param cacheDir directory containing pre-calculated function output
#'
#' @return terra SpatRaster covering one regression parameter per layer per cell
#'
#' @author Hagen Tockhorn
#'
#' @importFrom terra regress rast
#' @importFrom madrat toolGetMapping
#' @importFrom magclass as.magpie
#' @importFrom madrat readSource


calcBAITpars <- function(model = "GFDL-ESM4",
                         cacheDir = NULL) {

  # CHECK CACHE-----------------------------------------------------------------

  if (!is.null(cacheDir)) {
    fpath <- list.files(cacheDir, pattern = model, full.names = TRUE) %>%
      unlist()

    if (file.exists(fpath)) {
      print(paste0("Load parameters from cache: ", basename(fpath)))
      regPars <- rast(fpath)

      return(list(x = regPars,
                  class = "SpatRaster",
                  unit = "(unit)",
                  description = "Regression parameters for calcHDDCDD"))
    }
    else {
      print("BAITpars file not in given cache directory - will be re-calculated.")
    }
  }


  # READ-IN DATA----------------------------------------------------------------

  files <- toolGetMapping("baitregression-files_test.csv", type = "sectoral") %>%
    filter(.data[["gcm"]] == model)

  vars <- unique(files$variable)

  # nolint start
  data <- sapply(vars, function(v) {
    tmp <- sapply(files[files$variable == v, ]$file,
                  function(f) {
                    return(readSource("ISIMIPbuildings", subtype = f))
                    },
                  USE.NAMES = FALSE) %>%
      rast()

    return(tmp)
  },
  USE.NAMES = TRUE)
  # nolint end

  print("Reading completed")



  # PROCESS DATA----------------------------------------------------------------

  # convert huss into log scale
  data$huss <- log(data$huss)

  # convert tas into [C]
  data$tas <- data$tas - 273.15


  # nolint start
  regPars <- sapply(vars[vars != "tas"], function(v) {
    x <- data[["tas"]]
    y <- data[[v]]

    r <- regress(x = x, y = y, formula = y ~ x)

    names(r) <- c(paste0("a_", v), paste0("b_", v))
    return(r)
  },
  USE.NAMES = FALSE) %>%
    rast()
  # nolint end



  # OUTPUT----------------------------------------------------------------------

  return(list(x = regPars,
              class = "SpatRaster",
              unit = "(unit)",
              description = "Regression parameters for calcHDDCDD"))
}
