#' Calculate HDD and CDD based on outdoor/indoor temperature difference
#'
#' Heating and cooling degree days based on raw outside temperature
#' or bias-adjusted internal temperature (BAIT), driver for space heating and
#' cooling demand in buildings.
#'
#' @param mappingFile file name of sectoral mapping containing input data file names and directories
#' @param bait boolean, use BAIT instead of raw temperature
#' @param multiscen boolean, does \code{mappingFile} cover more than one scenario?
#' @param rasDir absolute path to directory for saving raster files
#' @param cacheDir absolute path to directory for pre-calculated BAIT regression parameters
#' @param fromSource logical, if TRUE, read the previously calculated data that
#'   is temporarily saved as a source
#'
#' @return magpie object of heating and cooling degree days
#'
#' @author Robin Hasse, Hagen Tockhorn
#'
#' @importFrom madrat toolGetMapping readSource calcOutput toolCountryFill
#' @importFrom tidyr %>%
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @importFrom stringr str_sub
#' @importFrom terra setGDALconfig


calcHDDCDD <- function(mappingFile,
                       bait = FALSE,
                       multiscen = FALSE,
                       rasDir = NULL,
                       cacheDir = NULL,
                       fromSource = FALSE) {

  if (isTRUE(fromSource)) {
    # Read previously calculated data from source
    hddcdd <- readSource("HDDCDDtemp")

    weight <- hddcdd
    weight[] <- 1

    return(list(x = hddcdd,
                min = 0,
                weight = weight,
                unit = "K.d/yr",
                description = "Heating and cooling degree days"))
  }


  # initialize full calculation
  makeCalculations <- function(f, m, n, tLim, countries, pop, hddcddFactor,
                               bait, wBAIT = NULL, params = NULL, suffix = NULL,
                               rasDir = NULL) {

    if (!is.null(suffix)) {
      ftas <- gsub(".nc",
                   paste0("_", suffix, ".nc"),
                   f[f$variable == "tas"       & f$gcm == m, ][[n, "file"]])

      if (bait) {
        frsds <- gsub(".nc",
                      paste0("_", suffix, ".nc"),
                      f[f$variable == "rsds"    & f$gcm == m, ][[n, "file"]])

        fsfc  <- gsub(".nc",
                      paste0("_", suffix, ".nc"),
                      f[f$variable == "sfcwind" & f$gcm == m, ][[n, "file"]])

        fhuss <- gsub(".nc",
                      paste0("_", suffix, ".nc"),
                      f[f$variable == "huss"    & f$gcm == m, ][[n, "file"]])
      }
    } else {
      ftas  <- f[f$variable == "tas"     & f$gcm == m, ][[n, "file"]]
      frsds <- f[f$variable == "rsds"    & f$gcm == m, ][[n, "file"]]
      fsfc  <- f[f$variable == "sfcwind" & f$gcm == m, ][[n, "file"]]
      fhuss <- f[f$variable == "huss"    & f$gcm == m, ][[n, "file"]]
    }

    print(paste("Processing temperature file:", ftas))

    if (bait) {
      hddcddCell <- compStackHDDCDD(ftas,
                                    tLim,
                                    countries,
                                    pop,
                                    hddcddFactor,
                                    bait,
                                    frsds  = frsds,
                                    fsfc   = fsfc,
                                    fhuss  = fhuss,
                                    wBAIT  = wBAIT,
                                    params = params,
                                    rasDir = rasDir)
    } else {
      hddcddCell <- compStackHDDCDD(ftas,
                                    tLim,
                                    countries,
                                    pop,
                                    hddcddFactor,
                                    bait,
                                    rasDir = rasDir)
    }
    return(hddcddCell)
  }



  # PARAMETERS------------------------------------------------------------------

  setGDALconfig(c("BIGTIFF = YES"))


  # threshold temperature for heating and cooling [C]
  # NOTE: Staffel at. al 2023 gives global average of T_heat = 14, T_cool = 20
  tLim <- list("HDD" = seq(12, 22), "CDD" = seq(20, 26))

  # standard deviations for temperature distributions
  tlimStd <- 5   # threshold
  tambStd <- 5   # ambient

  # range of pre-calculated HDD/CDD-values, e.g. [173, 348] K, converted to [C]
  tlow <- 173 - 273.15
  tup  <- 348 - 273.15



  #--- BAIT parameters

  # The weights (wRSDS, wSFC, wHUSS) for calcBAIT and the smoothing coefficient are assumed to
  # be region-independent and equal to the mean of the values given in Staffell
  # et al. 2023.

  # The weights below give the respective weight of the difference between the
  # climate factor and its counterfactual in calculating BAIT.
  wRSDS <- 0.012
  wSFC <- -0.20
  wHUSS <- 0.05

  # The smoothing coefficient defines the assumed thermal intertia of the building
  # and weighs the influence of the internal temperature of the preceding two days.
  sig <- 0.50

  # The blending parameters for the blending of BAIT and raw temperature are like-wise
  # taken from the paper.
  # At bLower, we consider only BAIT. For higher temperatures, we assume a mix with the ambient
  # temperature reaching the highest contribution of bMax at bUpper following a sigmoid function
  bLower <- 15
  bUpper <- 23
  bMax   <- 0.5

  # concatenate to vector
  wBAIT <- list("wRSDS"  = wRSDS,
                "wSFC"   = wSFC,
                "wHUSS"  = wHUSS,
                "sig"    = sig,
                "bLower" = bLower,
                "bUpper" = bUpper,
                "bMax"   = bMax)

  # BAIT parameter names
  # (note: this is really just a check in case calcBAITpars returns no layer names)
  parNames <- c("aRSDS", "bRSDS", "aSFC", "bSFC", "aHUSS", "bHUSS")



  # READ-IN DATA----------------------------------------------------------------

  # list of files that are processed
  files <- toolGetMapping(name  = mappingFile,
                          type  = "sectoral",
                          where = "mredgebuildings") %>%
    filter(.data[["variable"]] != "")

  # cells -> country mask
  fCM <- file.path(files[files$variable == "CountryMask", "file"])
  countries <- readSource("ISIMIPbuildings", subtype = fCM, convert = FALSE)

  # population
  fpop <- files %>% filter(.data[["variable"]] == "pop")
  pop  <- readSource("ISIMIPbuildings", subtype = fpop$file,
                     convert = FALSE)

  # extract run-specific variables
  ssp   <- unique(files$ssp[files$variable == "tas"])
  rcp   <- unique(files$rcp[files$variable == "tas"])
  model <- unique(files$gcm[files$variable == "tas"])


  # bait regression parameters
  if (bait) {
    baitPars <- calcOutput("BAITpars",
                           aggregate = FALSE,
                           model     = model,
                           cacheDir  = cacheDir)

    names(baitPars) <- parNames
  }



  # PROCESS DATA----------------------------------------------------------------

  # calculate HDD/CDD-factors
  hddcddFactor <- compHDDCDDFactors(tlow = tlow, tup = tup, tLim, tambStd, tlimStd)


  # full calculation of degree days
  hddcdd <- do.call( # file iteration
    "rbind",
    lapply(
      seq_len(nrow(filter(files, files$variable == "tas"))),
      function(n) {
        split <- files[files$variable == "tas" & files$gcm == model, ][[n, "split"]]

        # split large raster files to save memory / speed up processing
        if (as.logical(split)) {
          hddcddCell <- do.call(
            "rbind",
            lapply(list("A", "B"), function(suffix) {
              tmp <- makeCalculations(f = files,
                                      m = model,
                                      n = n,
                                      tLim = tLim,
                                      countries = countries,
                                      pop = pop,
                                      hddcddFactor = hddcddFactor,
                                      bait = bait,
                                      wBAIT = wBAIT,
                                      params = baitPars,
                                      suffix = suffix,
                                      rasDir = rasDir)
              return(tmp)
            }
            )
          )
        } else {
          hddcddCell <- makeCalculations(f = files,
                                         m = model,
                                         n = n,
                                         tLim = tLim,
                                         countries = countries,
                                         pop = pop,
                                         hddcddFactor = hddcddFactor,
                                         bait = bait,
                                         wBAIT = wBAIT,
                                         params = baitPars,
                                         rasDir = rasDir)
        }

        hddcddCell <- hddcddCell %>%
          mutate("model" = model,
                 "ssp" = ssp,
                 "rcp" = rcp)

        gc()

        return(hddcddCell)
      }
    )
  )

  rownames(hddcdd) <- c()



  # OUTPUT----------------------------------------------------------------------

  # might be necessary for magpie output
  data <- hddcdd %>%
    mutate(period = gsub("y", "", .data[["period"]])) %>%
    unite("variable", .data[["variable"]], .data[["tlim"]], sep = "_") %>%
    unite("scenario", .data[["ssp"]], .data[["rcp"]], sep = "_")

  # prepare output
  data <- data %>%
    as.quitte() %>%
    as.magpie() %>%
    toolCountryFill()


  return(list(x = data))
}


#--- CALCULATE BAIT (building-adjusted internal temperature)


#' Check if time period of BAIT input data (rsds, sfc, huss) is congruent with
#' near-surface temperature data (tas).
#'
#' @param baitInput list of raster data encompassing different climate variables
#' @param tasData raster data on near-surface atmosphere temperature
#'
#' @return baitInput with congruent time periods w.r.t. tasData
#'
#' @importFrom terra rast

checkDates <- function(baitInput, tasData) {
  datesT <- names(tasData)

  baitInput <- sapply(names(baitInput), function(var) { #nolint
    # fill missing data with means from previous years
    # NOTE: "temp" and "baitInput" have the same global temporal lower
    #       boundary, since "temp" is the constraining dataset, only
    #       "baitInput" needs to be filled up.

    tmp <- baitInput[[var]]

    datesBait <- names(tmp)

    datesFill <- setdiff(datesT, datesBait)        # dates to fill up
    daysFill  <- unique(substr(datesFill, 6, 11))

    datesKeep <- intersect(datesBait, datesT)      # dates to keep
    keep      <- length(datesKeep) > 0

    if (keep) {
      tmp        <- subset(tmp, datesKeep)
      names(tmp) <- datesKeep
    }

    if (length(daysFill) > 0) {
      baitInputMean <- prepBaitInput(fillWithMean = TRUE, baitInput = baitInput)

      # fill up missing dates with yearly-average value for specific day/cell
      baitInputFill <- rast(
        lapply(
          daysFill,
          function(d) {
            idx <- which(grepl(d, stringr::str_sub(datesFill, -5, -1)))
            r   <- rast(replicate(length(idx), baitInputMean[[var]][[d]]))
            names(r) <- datesFill[idx]
            return(r)
          }
        )
      )

      # concatenate data
      if (keep) {
        tmp <- rast(list(tmp, baitInputFill))
      } else {
        tmp <- baitInputFill
      }

      # re-order dates
      tmp <- rast(tmp[[order(names(tmp))]])
    }


    if (!identical(names(tmp), names(tasData))) {
      warning("Dates of Temperature and BAIT Input Data are not aligned.")
    }
    return(tmp)
  },
  USE.NAMES = TRUE
  )
  return(baitInput)
}


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
  print("print h")
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

  rm(baitInput, solar, wind, hum, s, w, h, t, bait)

  return(bait)
}


#--- CALCULATE HDD/CDD

#' Calculate HDD/CDD values for different ambient/limit temperature combinations
#'
#' HDD/CDD values are pre-calculated for an interval \code{tlow}-\code{tup} and
#' for a set of limit temperatures \code{tlim} with a temperature resolution of
#' 0.1C.
#'
#' The respective heating/cooling degree days are calculated as the difference
#' between the assumed ambient and a limit temperature, aggregated to a full day.
#' The latter defines a threshold above/below which cooling/heating is assumed to
#' be initiated.
#'
#' To account for heterogenity in heating/cooling behavior, the ambient and limit
#' temperature, \code{tamb} and \code{tlim}, are assumed to be normally distributed.
#' This changes the calculation of a degree day to a double integration of
#' \code{tLimit - T_ambient_day} with integration boundaries set at 3 standard
#' deviations, \code{tambStd} and \code{tlimStd}, from \code{tamb} and \code{tlim}
#' respectively.
#'
#' As consequence, the ramp function of \code{HDD_day = max(0, tLimit - T_ambient_day)}
#' changes to a curved function that is above zero even if the mean of \code{T_ambient_day}
#' is above the mean of \code{tLimit}.
#'
#' @param tlow lower temperature boundary
#' @param tup upper temperature boundary
#' @param tlim named list of limit temperature sequences for \code{HDD} and \code{CDD}
#' @param tambStd std of ambient temperature
#' @param tlimStd std of limit temperature
#'
#' @return data frame of HDD/CDD values
#'
#' @importFrom stats dnorm
#' @importFrom pracma integral2

compHDDCDDFactors <- function(tlow, tup, tlim, tambStd = 5, tlimStd = 5) {

  # t1 : ambient temperature variable
  # t2 : limit temperature variable

  # HDD
  heatingFactor <- function(t2, t1, tamb, tambStd, tlim, tlimStd) {
    h <- dnorm(t2, mean = tlim, sd = tlimStd) * dnorm(t1, mean = tamb, sd = tambStd) * (t2 - t1)
    return(h)
  }

  # CDD
  coolingFactor <- function(t2, t1, tamb, tambStd, tlim, tlimStd) {
    h <- dnorm(t2, mean = tlim, sd = tlimStd) * dnorm(t1, mean = tamb, sd = tambStd) * (t1 - t2)
    return(h)
  }

  t <- seq(tlow, tup, .1)

  hddcddFactors <- do.call(
    "rbind", lapply(
      c("HDD", "CDD"), function(typeDD) {
        do.call(
          "rbind", lapply(
            t, function(tamb) {
              do.call(
                "rbind", lapply(
                  tlim[[typeDD]], function(.tlim) {

                    # tlim integration boundaries
                    x1 <- .tlim - 3 * tlimStd
                    x2 <- .tlim + 3 * tlimStd

                    # nolint start
                    switch(typeDD,
                           HDD = {
                             fun <- heatingFactor
                             ymin <- tamb - 3 * tambStd
                             ymax <- min(.tlim, tamb + 3 * tambStd)
                           },
                           CDD = {
                             fun  <- coolingFactor
                             ymin <- max(.tlim, tamb - 3 * tambStd)
                             ymax <- tamb + 3 * tambStd
                           }
                    )
                    # nolint end

                    f <- integral2(fun,
                                   xmin = x1,
                                   xmax = x2,
                                   ymin = ymin,
                                   ymax = ymax,
                                   tamb = tamb,
                                   tambStd = tambStd,
                                   tlim = .tlim,
                                   tlimStd = tlimStd,
                                   reltol = 1e-1)

                    tmp <- data.frame("T_amb"        = tamb,
                                      "T_amb_K"      = round(tamb + 273.15, 1),
                                      "tLim"         = .tlim,
                                      "factor"       = f$Q,
                                      "factor_err"   = f$error,
                                      "typeDD"       = typeDD)

                    return(tmp)
                  }
                )
              )
            }
          )
        )
      }
    )
  )
  return(hddcddFactors)
}


#' Assign HDD/CDD values for given ambient/limit temperature
#'
#' @param temp raster data containing temperature/BAIT values
#' @param typeDD type of degree day
#' @param tlim limit temperature
#' @param factors data frame with degree day values for \code{temp/tlim} combination
#'
#' @return raster object with HDD/CDD values
#'
#' @importFrom terra classify tapp

compCellHDDCDD <- function(temp, typeDD, tlim, factors) {
  # extract years
  dates <- names(temp)

  # add tolerance of 0.04K to avoid machine precision errors
  factors <- factors[factors$typeDD == typeDD, ]

  factors <- factors %>%
    filter(.data[["tLim"]] == tlim) %>%
    dplyr::reframe(from = .data[["T_amb_K"]] - 0.049,
                   to = .data[["T_amb_K"]] + 0.049,
                   becomes = .data[["factor"]]) %>%
    data.matrix()

  # swap ambient temperature values with corresponding DD values
  hddcdd <- classify(temp, factors)

  terra::time(hddcdd) <- as.Date(dates)

  # aggregate to yearly HDD/CDD [K.d/a]
  hddcdd <- tapp(hddcdd, "years", fun = sum, na.rm = TRUE)

  names(hddcdd) <- gsub("_", "", names(hddcdd))
  return(hddcdd)
}


#' Aggregate cellular HDD/CDD values to country-wide average (population-weighted)
#'
#' @param data raster object containing HDD/CDD values
#' @param weight raster object containing aggregation weights
#' @param mask raster opbject defining (regional) aggregation boundaries
#'
#' @return data frame containing regionally meaned HDD/CDD values
#'
#' @importFrom terra subset

aggCells <- function(data, weight, mask) {
  yearsData   <- names(data)
  yearsWeight <- names(weight)

  if (!all(yearsData %in% yearsWeight)) {
    stop("Time periods of raster file and aggregation weights do not match.")
  }

  # loop: years in raster file r
  hddcddAgg <- do.call(
    "rbind", lapply(
      yearsData, function(y) {
        # mask data and weights to considered regions
        regData   <- subset(data, y) * subset(weight, y) * mask
        regWeight <- subset(weight, y) * mask

        # aggregate regional data
        regDataAgg   <- terra::global(regData,   "sum", na.rm = TRUE)$sum
        regWeightAgg <- terra::global(regWeight, "sum", na.rm = TRUE)$sum

        # calculate weighted sum
        weightedAgg <- regDataAgg / regWeightAgg

        tmp <- data.frame("region" = names(mask),
                          "period" = y,
                          "value"  = round(weightedAgg, 3))

        rm(regData, regWeight, regDataAgg, regWeightAgg, weightedAgg)

        rownames(tmp) <- c()
        return(tmp)
      }
    )
  )
  return(hddcddAgg)
}


#' Calculate country-wise population-weighted HDD/CDD values
#'
#' This function calculates country-wise population-weighted HDD/CDD values for
#' raw ambient temperature or bias-adjusted internal temperature for a given set
#' of limit temperatures from raster data on (various) climate variables.
#'
#' For further processing, raster objects containing degree day data are written
#' for an interval of ten years.
#'
#' @param ftas file name of data on near-surface atmospherical temperature
#' @param tlim named list of limit temperature sequences for \code{HDD} and \code{CDD}
#' @param countries raster opbject defining (regional) aggregation boundaries
#' @param pop raster object containing population data
#' @param factors data frame with degree day values for \code{temp/tlim} combination
#' @param bait boolean, BAIT is used as ambient temperature
#' @param frsds file name of data on surface downdwelling shortwave radiation (optional)
#' @param fsfc file name of data on near-surface wind speed (optional)
#' @param fhuss file name of data on near-surface specific humidity (optional)
#' @param wBAIT named list containing BAIT weights (optional)
#' @param params raster object containing regression parameters from \code{calcBAITpars} (optional)
#' @param rasDir absolute path to directory for saving raster files
#'
#' @importFrom raster writeRaster
#' @importFrom stringr str_split
#' @importFrom terra writeCDF

compStackHDDCDD <- function(ftas, tlim, countries, pop, factors, bait,
                            frsds = NULL,
                            fsfc  = NULL,
                            fhuss = NULL,
                            wBAIT = NULL,
                            params = NULL,
                            rasDir = NULL) {

  # read cellular temperature
  temp <- readSource("ISIMIPbuildings", subtype = ftas, convert = TRUE)

  dates <- names(temp)

  # optional: transform raw temperature into BAIT
  if (bait) {
    # note: easier to do in [C]
    tasData <- temp - 273.15   # [C]

    # read and prepare bait input data
    baitInput <- prepBaitInput(frsds, fsfc, fhuss) %>%
      checkDates(tasData)

    # calculate bait
    temp <- compBAIT(baitInput, tasData, weight = wBAIT, params = params)

    # convert back to [K]
    temp <- temp + 273.15   # [K]
  }

  # round and assign dates
  temp <- terra::round(temp, digits = 1)
  names(temp) <- dates

  print("Calculating HDD/CDDs per cell.")


  fSplit <- str_split(ftas, "_") %>% unlist()

  # loop: typeDD
  hddcdd <- do.call(
    "rbind", lapply(
      c("HDD", "CDD"), function(typeDD) {
        # loop: threshold temperatures
        do.call(
          "rbind", lapply(
            tlim[[typeDD]], function(t) {
              hddcddAgg <- compCellHDDCDD(temp, typeDD, t, factors)

              # write raster files
              if (!is.null(rasDir)) {
                y <- names(hddcddAgg)

                rname <- paste0(fSplit[[1]], "_", y, "_", fSplit[[4]], "_", typeDD, "_", t)
                rname <- paste0(rname, if (bait) "_bait" else "", ".nc")

                writeCDF(hddcddAgg,
                         file.path(rasDir, fSplit[[1]], rname),
                         overwrite = TRUE)
              }

              hddcddAgg <- hddcddAgg %>%
                aggCells(pop, countries) %>%
                mutate("variable" = typeDD,
                       "tlim"     = t)    # [C]

              return(hddcddAgg)
            }
          )
        )
      }
    )
  )
  rm(tasData, baitInput, temp)

  return(hddcdd)
}
