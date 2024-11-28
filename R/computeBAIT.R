# --- FUNCTIONS TO CALCULATE BIAS-ADJUSTED INTERNAL TEMPERATURE ----------------



#' Read in necessary climate data to calculate BAIT or calculate mean values of
#' said climate data to fill missing data in case of temporal mismatch between
#' near-surface atmospherical temperature and other considered climate data.
#'
#' @param frsds file path to raster data on surface downdwelling shortwave radiation
#' @param fsfc file path to raster data on near-surface wind speed
#' @param fhuss file path to raster data on near-surface specific humidity
#' @param baitInput named list of climate data
#' @param fillWithMean boolean, only mean is calculated and returned
#'
#' @return named list with read-in or meaned climate data
#'
#' @importFrom terra tapp
#' @importFrom madrat readSource

prepBaitInput <- function(frsds = NULL,
                          fsfc = NULL,
                          fhuss = NULL,
                          baitInput = NULL,
                          fillWithMean = FALSE) {

  if (isTRUE(fillWithMean)) {
    # optional: calculate daily means over years to fill missing data
    baitInputMean <- sapply( #nolint
      names(baitInput),
      function(var) {
        meanData <- tapp(baitInput[[var]],
                         unique(substr(names(baitInput[[var]]), 6, 11)),
                         fun = "mean")
        names(meanData) <- gsub("\\.", "-", substr(names(mean), 2, 6))
        return(meanData)
      }
    )
    return(baitInputMean)
  } else {
    input <- list( #nolint start
      "rsds" = readSource("ISIMIPbuildings", subtype = frsds, convert = TRUE),
      "sfc"  = readSource("ISIMIPbuildings", subtype = fsfc,  convert = TRUE),
      "huss" = readSource("ISIMIPbuildings", subtype = fhuss, convert = TRUE))
    return(input) # nolint end
  }
}



#' Calculate counterfactuals for solar radiation, wind speed, specific humidity
#' and near-surface temperature as function of raster data on near-surface temperature.
#'
#' The expected value of the respective climate variable (except temperature) is
#' calculated from parameters taken from a preceding linear regression done in
#' calcBAITpars where the respective variable is correlated with the near-surface
#' atmospherical temperature.
#' If no cell-resoluted parameters are given, the globally-meaned parameters from
#' Staffel et. all 2023 are taken (see https://doi.org/10.1038/s41560-023-01341-5).
#'
#' @param t raster data on near-surface atmospherical temperature
#' @param type considered climate variable
#' @param params regression parameters as vector or raster object
#'
#' @return counterfactuals for respective climate variable

cfac <- function(t, type, params = NULL) {
  if (is.null(params)) {
    params <- switch(type,
                     s = c(100, 7),
                     w = c(4.5, -0.025),
                     h = c(1.1, 0.06),
                     t = c(16))
  }

  # nolint start
  return(switch(type,
                s = {params[[1]] + params[[2]] * t},
                w = {params[[1]] + params[[2]] * t},
                h = {exp(params[[1]] + params[[2]] * t)},
                t = {params[[1]]},
                warning("No valid parameter type specified.")
  )
  )
  # nolint end
}



#' Smooth data over preceding two days
#'
#' @param r raster data to be smoothed
#' @param weight named list with smoothing parameter sig
#'
#' @return smoothed raster data
#'
#' @importFrom terra nlyr

smooth <- function(r, weight) {
  print("smooth")

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



#' Weighted blend of BAIT and near-surface atmospherical temperature
#'
#' To adress loss of buildings' internal memory of previous conditions at high
#' outside temperatures due to window opening, etc., BAIT and outside temperature
#' are blended. The blend is active between a lower and upper temperature threshold,
#' \code{bLower} and \code{bUpper}, which are mapped to a range of -5 to +5 of a
#' sigmoid function (corresponding to a 1% and 99% blend). The maximum amount of
#' blending (i.e. the amplitude of the sigmoid function) is given by a parameter \code{bMax}.
#'
#' @param bait raster data on BAIT
#' @param tas raster data on near-surface atmospherical temperature
#' @param weight named list with blending parameters bLower, bUpper, bMax
#'
#' @return blended raster data

blend <- function(bait, tas, weight) {
  print("blend")

  bBar <- (tas - 0.5 * (weight[["bUpper"]] + weight[["bLower"]])) * 10 / (weight[["bUpper"]] - weight[["bLower"]])
  b    <- weight[["bMax"]] / (1 + exp(-bBar))

  blend <- bait * (1 - b) + (tas * b)
  return(blend)
}



#' Calculate bias-adjusted internal temperature (BAIT)
#'
#' BAIT is calculated from raster data on near-surface atmospherical temperature
#' (tas), surface downdwelling shortwave radiation (rsds), near-surface wind speed
#' (sfcwind) and near-surface specific humidity (huss). The latter three climate
#' parameters are incorporated in the calculation of BAIT as the difference from
#' their real value to the their expected value w.r.t. the near-surface temperature
#' (see \code{\link{cfac}}). These are then incorporated in a weighted sum to
#' account for the respective influence of each climate parameter on BAIT.
#' The raster data containing BAIT is smoothed to account for the
#' buildings' thermal inertia (see \code{\link{smooth}}) and blended with the
#' near-surface temperature (see \code{\link{blend}}).
#'
#' @param baitInput named list containing rsds, sfcwind, huss climate data
#' @param tasData raster data on near-surface atmospherical temperature
#' @param weight named list with weights
#' @param params optional named list with regression parameters from calcBAITpars
#'
#' @return raster object with BAIT values

compBAIT <- function(baitInput, tasData, weight = NULL, params = NULL) {
  if (is.null(weight)) {
    warning("Please give appropriate weights for the calculation of BAIT.")
    weight <- list("wRSDS"  = 0.012,
                   "wSFC"   = -0.20,
                   "wHUSS"  = 0.05,
                   "sig"    = 0.5,
                   "bLower" = 15,
                   "bUpper" = 23,
                   "bMax"   = 0.5)
  }

  solar <- baitInput$rsds
  wind  <- baitInput$sfc
  hum   <- baitInput$huss

  print("calc s")
  s <- solar -  cfac(tasData, type = "s", params = c(params[["aRSDS"]], params[["bRSDS"]]))
  print("calc w")
  w <- wind  -  cfac(tasData, type = "w", params = c(params[["aSFC"]], params[["bSFC"]]))
  print("calc h")
  h <- hum   -  cfac(tasData, type = "h", params = c(params[["aHUSS"]], params[["bHUSS"]]))
  print("calc t")
  t <- tasData - cfac(tasData, type = "t", params = NULL)

  # calc bait
  print("calc bait")
  bait <- tasData + weight[["wRSDS"]] * s + weight[["wSFC"]] * w + weight[["wHUSS"]] * h * t

  # smooth bait
  bait <- smooth(bait, weight)

  # # blend bait
  bait <- blend(bait, tasData, weight)

  return(bait)
}
