#' Calculate Regression Parameters for BAIT Climate Variables
#'
#' @description Performs linear regression on historical climate data to calculate
#' regression parameters for key climate variables: surface downwelling shortwave
#' radiation (rsds), near-surface wind speed (sfcwind), and near-surface specific
#' humidity (huss), with respect to near-surface air temperature (tas).
#'
#' The regression utilizes a simple linear model based on historical data from
#' 2000 to 2019-21 depending on the model. For rsds and sfcwind, a linear relationship
#' is assumed, while for huss, an exponential relationship is adopted, reflecting
#' the non-linear dependence of water vapor pressure on temperature.
#'
#' @param model Character string specifying the GCM responsible for data input.
#' @param cacheDir Character string specifying the directory containing pre-calculated
#' function outputs.
#'
#' @return A `terra::SpatRaster` object with layers representing regression
#' parameters for each cell.
#'
#' @author Hagen Tockhorn
#'
#' @importFrom terra regress rast
#' @importFrom madrat toolGetMapping readSource
#' @importFrom magclass as.magpie
#' @importFrom utils read.csv2



computeBAITpars <- function(model = "20crv3-era5",
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
    } else {
      print("BAITpars file not in given cache directory - will be re-calculated.")
    }
  }


  # READ-IN DATA----------------------------------------------------------------

  files <- read.csv2("inst/extdata/sectoral/BAITpars_fileMapping.csv") %>% # nolint
    filter(.data[["gcm"]] == model)

  vars <- c("tas", "sfcwind", "rsds", "huss")

  # nolint start
  data <- sapply(vars, function(v) {
    tmp <- sapply(files[[v]],
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

  terra::writeCDF(regPars, paste0("/p/tmp/hagento/output/", model, "/baitpars_", model, ".nc")) # nolint


  # OUTPUT----------------------------------------------------------------------

  return(list(x = regPars,
              class = "SpatRaster",
              unit = "(unit)",
              description = "Regression parameters for calcHDDCDD"))
}
