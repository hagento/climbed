#' Read and Process Climate Data for BAIT Calculation or Data Imputation
#'
#' This function reads the necessary climate data to calculate BAIT (Bias-Adjusted Internal Temperature)
#' or computes mean values of the specified climate variables to impute missing data.
#' Imputation is used in cases where there is a temporal mismatch between near-surface atmospheric temperature
#' data and other climate variables.
#'
#' @param fileNames \code{character} A named vector of file paths extracted from \code{fileMapping}. Elements include:
#'   \describe{
#'     \item{tas}{Temperature data file path.}
#'     \item{rsds}{Solar radiation data file path (included if \code{bait} is \code{TRUE}).}
#'     \item{sfc}{Surface wind data file path (included if \code{bait} is \code{TRUE}).}
#'     \item{huss}{Humidity data file path (included if \code{bait} is \code{TRUE}).}
#'   }
#' @param baitInput \code{list} containing \code{terra::SpatRaster} objects for
#' climate data inputs used in BAIT calculation.
#' @param fillWithMean \code{logical}; if \code{TRUE}, the function calculates and
#' returns the mean of the climate data
#'   instead of the raw data.
#'
#' @returns \code{list} containing either the read-in \code{terra::SpatRaster} climate
#' data or the computed mean values of the data.
#'
#' @importFrom terra tapp
#' @importFrom stats setNames

prepBaitInput <- function(fileNames,
                          baitInput = NULL,
                          fillWithMean = FALSE) {

  if (isTRUE(fillWithMean)) {
    # optional: calculate daily means over years to fill missing data
    baitInputMean <- setNames(lapply(names(baitInput), function(var) {
      meanData <- tapp(baitInput[[var]],
                       unique(substr(names(baitInput[[var]]), 6, 11)),
                       fun = "mean")

      names(meanData) <- gsub("\\.", "-", substr(names(meanData), 2, 6))
      return(meanData)
    }),
    names(baitInput))

    return(baitInputMean)
  } else {
    input <- list("rsds" = importData(subtype = fileNames[["rsds"]]),
                  "sfc"  = importData(subtype = fileNames[["sfc"]]),
                  "huss" = importData(subtype = fileNames[["huss"]]))
    return(input)
  }
}



#' Calculate Counterfactual Climate Variables as a Function of Near-Surface Temperature
#'
#' Computes the expected values of solar radiation, wind speed, specific humidity,
#' or near-surface temperature based on raster data of near-surface temperature.
#' The calculation uses regression parameters derived from prior analysis, correlating
#' the respective climate variable with near-surface atmospheric temperature.
#'
#' If cell-specific regression parameters are unavailable, globally averaged
#' parameters from Staffel et al. (2023) are used as a fallback.
#' For further details, see \url{https://doi.org/10.1038/s41560-023-01341-5}.
#'
#' @param t \code{terra::SpatRaster} object representing near-surface atmospheric temperature.
#' @param type \code{character} string specifying the climate variable to calculate. Options include:
#'   \itemize{
#'     \item \code{"rsds"}: Shortwave solar radiation.
#'     \item \code{"sfcwind"}: Surface wind speed.
#'     \item \code{"huss"}: Specific humidity.
#'     \item \code{"tas"}: Near-surface temperature (optional, if applying the function to temperature itself).
#'   }
#' @param params Regression parameters provided as either:
#'   \itemize{
#'     \item A \code{vector} for globally averaged parameters.
#'     \item A \code{terra::SpatRaster} object for cell-specific parameters.
#'   }
#'
#' @returns \code{terra::SpatRaster} object representing the calculated counterfactual
#'   values for the specified climate variable.
#'
#' @importFrom terra rast
#' @importFrom terra subset
#' @importFrom stringr str_sub
#' @importFrom utils read.csv2

cfac <- function(t, type, params) {
  return(switch(type,
                s = params[["aRSDS"]] + params[["bRSDS"]] * t,
                w = params[["aSFC"]] + params[["bSFC"]] * t,
                h = exp(params[["aHUSS"]] + params[["bHUSS"]] * t),
                t = params[["T"]],
                warning("No valid parameter type specified.")))
}



#' Smooth Data Over Preceding Two Days
#'
#' Applies a smoothing operation to raster data using a specified smoothing parameter.
#'
#' @param r \code{terra::SpatRaster} object containing the data to be smoothed.
#' @param weight \code{list} with a named smoothing parameter \code{sig}.
#'
#' @returns \code{terra::SpatRaster} object with smoothed data.
#'
#' @importFrom terra nlyr

smooth <- function(r, weight) {
  # one day indented
  r1D <- r[[c(nlyr(r), 1:(nlyr(r) - 1))]]
  r1D[[1]] <- 0

  # two days indented
  r2D <-  r[[c(nlyr(r) - 1, nlyr(r), 1:(nlyr(r) - 2))]]
  r2D[[1:2]] <- 0

  # smooth
  rSmooth <- (r + weight[["sig"]] * r1D + weight[["sig"]]**2 * r2D) / (1 + weight[["sig"]] + weight[["sig"]]**2)

  return(rSmooth)
}



#' Weighted Blend of BAIT and Near-Surface Air Temperature
#'
#' Blends BAIT and near-surface temperature data to account for the loss of a building's
#' internal memory of previous conditions at high external temperatures (e.g., due to window opening).
#' The blend is applied between a lower and upper temperature threshold, \code{bLower} and \code{bUpper},
#' mapped to a range of -5 to +5 in a sigmoid function, corresponding to 1% and 99% blending.
#' The maximum blending amplitude is controlled by the parameter \code{bMax}.
#'
#' @param bait \code{terra::SpatRaster} object containing BAIT data.
#' @param tas \code{terra::SpatRaster} object containing near-surface air temperature data.
#' @param weight \code{list} with blending parameters \code{bLower}, \code{bUpper}, and \code{bMax}.
#'
#' @returns \code{terra::SpatRaster} object with blended data.

blend <- function(bait, tas, weight) {
  bBar <- (tas - mean(weight[c("bUpper", "bLower")])) * 10 / (weight[["bUpper"]] - weight[["bLower"]])
  b    <- weight[["bMax"]] / (1 + exp(-bBar))

  blend <- bait * (1 - b) + (tas * b)
  return(blend)
}



#' Calculate Bias-Adjusted Internal Temperature (BAIT)
#'
#' Computes BAIT from raster data of near-surface atmospheric temperature (\code{tas}),
#' surface downwelling shortwave radiation (\code{rsds}), near-surface wind speed
#' (\code{sfcwind}), and near-surface specific humidity (\code{huss}).
#' The latter three climate variables are adjusted by calculating the difference
#' between their real values and their expected values relative to near-surface
#' temperature (see \code{\link{cfac}}). These adjusted values are combined in
#' a weighted sum to account for the respective influence of each variable on BAIT.
#'
#' The resulting BAIT data is smoothed to account for thermal inertia (see \code{\link{smooth}})
#' and blended with near-surface temperature (see \code{\link{blend}}) to address the
#' loss of internal memory at high temperatures.
#'
#' @param baitInput \code{list} containing \code{rsds}, \code{sfcwind}, and \code{huss}
#'   climate data as \code{terra::SpatRaster} objects.
#' @param tasData \code{terra::SpatRaster} object representing near-surface atmospheric temperature.
#' @param weight \code{list} containing weights for the climate variables and blending parameters.
#' @param params \code{list} (optional) containing regression parameters derived from \code{calcBAITpars}.
#'
#' @returns \code{terra::SpatRaster} object containing BAIT values.

compBAIT <- function(baitInput, tasData, weight, params = NULL) {
  # extract climate data
  solar <- baitInput$rsds
  wind  <- baitInput$sfc
  hum   <- baitInput$huss

  # calculate respective summands
  s <- solar   - cfac(tasData, type = "s", params = params)
  w <- wind    - cfac(tasData, type = "w", params = params)
  h <- hum     - cfac(tasData, type = "h", params = params)
  t <- tasData - cfac(tasData, type = "t", params = params)

  # calc bait
  bait <- tasData + weight[["wRSDS"]] * s + weight[["wSFC"]] * w + weight[["wHUSS"]] * h * t

  # smooth bait
  bait <- smooth(bait, weight)

  # blend bait
  bait <- blend(bait, tasData, weight)

  return(bait)
}
