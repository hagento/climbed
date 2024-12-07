#' Initialize Degree Days Calculation
#'
#' Initiates the calculation of degree days for a single scenario, model, and time period.
#' The calculation is split into yearly periods to improve computational stability.
#'
#' @param fileMapping A named list containing file paths and metadata. Expected elements include:
#'   \describe{
#'     \item{tas}{Path to temperature data file.}
#'     \item{rsds}{Path to solar radiation data file (required if \code{bait} is TRUE).}
#'     \item{sfcwind}{Path to surface wind data file (required if \code{bait} is TRUE).}
#'     \item{huss}{Path to humidity data file (required if \code{bait} is TRUE).}
#'     \item{start}{Start year of the time period.}
#'     \item{end}{End year of the time period.}
#'     \item{rcp}{RCP scenario identifier.}
#'     \item{gcm}{GCM model identifier.}
#'   }
#' @param pop \code{SpatRaster} with annual population data.
#' @param ssp \code{character} SSP scenario.
#' @param bait \code{logical} indicating whether to use raw temperature or BAIT as ambient temperature.
#' @param tLim \code{numeric} Temperature limits for degree day calculations.
#' @param hddcddFactor \code{data.frame} containing pre-computed degree days.
#' @param wBAIT \code{numeric} (Optional) Weights for BAIT adjustments. Default is \code{NULL}.
#'
#' @returns \code{data.frame} containing annual degree days.
#'
#' @author Hagen Tockhorn
#'
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%

initCalculation <- function(fileMapping,
                            pop,
                            ssp,
                            bait,
                            tLim,
                            hddcddFactor,
                            wBAIT = NULL) {
  # extract filenames
  ftas  <- fileMapping[["tas"]]
  frsds <- if (bait) fileMapping[["rsds"]]    else NULL
  fsfc  <- if (bait) fileMapping[["sfcwind"]] else NULL
  fhuss <- if (bait) fileMapping[["huss"]]    else NULL

  # extract temporal interval
  yStart <- fileMapping[["start"]] %>% as.numeric()
  yEnd   <- fileMapping[["end"]] %>% as.numeric()

  # extract RCP scenario + model
  rcp   <- fileMapping[["rcp"]] %>% unique()
  model <- fileMapping[["gcm"]] %>% unique()


  # read country masks
  countries <- importData(subtype = "countrymasks-fractional_30arcmin.nc")


  if (bait) {
    # bait regression parameters
    baitPars <- computeBAITpars(model = unique(fileMapping$gcm))
    names(baitPars) <- c("aRSDS", "bRSDS", "aSFC", "bSFC", "aHUSS", "bHUSS") # TODO: change this
  }


  # loop over single years and compute annual degree days
  hddcdd <- do.call("rbind", lapply(seq(1, yEnd - yStart + 1), function(i) {
    message("Initiating calculating degree days for the year: ", seq(yStart, yEnd)[[i]])

    compStackHDDCDD(ftas  = gsub(".nc", paste0("_", i, ".nc"), ftas),
                    frsds = if (bait) gsub(".nc", paste0("_", i, ".nc"), frsds) else NULL,
                    fsfc  = if (bait) gsub(".nc", paste0("_", i, ".nc"), fsfc)  else NULL,
                    fhuss = if (bait) gsub(".nc", paste0("_", i, ".nc"), fhuss) else NULL,
                    tlim = tLim,
                    pop = pop,
                    countries = countries,
                    factors = hddcddFactor,
                    bait = bait,
                    wBAIT  = wBAIT,
                    baitPars = baitPars)
  }))

  hddcdd <- hddcdd %>%
    mutate("model" = model,
           "ssp" = ssp,
           "rcp" = rcp)

  return(hddcdd)
}




#' Calculate Country-wise Population-weighted HDD/CDD Values
#'
#' This function computes population-weighted Heating Degree Days (HDD) and Cooling Degree Days (CDD)
#' for individual countries, based on raw ambient temperature or bias-adjusted internal temperature.
#' Calculations are performed using raster data for various climate variables and specified limit temperatures.
#'
#' For further processing, raster objects containing degree day data are written for ten-year intervals.
#'
#' @param ftas \code{character} string specifying the file name of data on near-surface atmospheric temperature.
#' @param tlim \code{list} named list of temperature limit sequences for \code{HDD} and \code{CDD}.
#' @param countries \code{terra::SpatRaster} object defining the regional aggregation boundaries.
#' @param pop \code{terra::SpatRaster} object containing population data.
#' @param factors \code{data.frame} containing degree day values for the combinations of \code{temp} and \code{tlim}.
#' @param bait \code{logical}; if \code{TRUE}, BAIT is used as the ambient temperature instead of
#' near-surface temperature.
#' @param frsds (Optional) \code{character} string specifying the file name of data on surface
#' downdwelling shortwave radiation.
#' @param fsfc (Optional) \code{character} string specifying the file name of data on near-surface wind speed.
#' @param fhuss (Optional) \code{character} string specifying the file name of data on near-surface specific humidity.
#' @param wBAIT (Optional) \code{list} named list containing weights for BAIT calculations.
#' @param baitPars (Optional) \code{terra::SpatRaster} object containing regression parameters
#' generated by \code{calcBAITpars}.
#'
#' @returns \code{data.frame} containing regional population-weighted annual degree days (HDD/CDD).
#'
#' @details
#' The function integrates various climate and population datasets to produce a regionalized summary of
#' HDD/CDD values. It supports optional adjustments based on additional climate variables and regression
#' parameters when using BAIT for internal temperature calculations.
#'
#' @author Hagen Tockhorn
#'
#' @importFrom raster writeRaster
#' @importFrom stringr str_split
#' @importFrom terra writeCDF round
#' @importFrom madrat readSource
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%

compStackHDDCDD <- function(ftas, tlim, countries, pop, factors, bait,
                            frsds = NULL,
                            fsfc  = NULL,
                            fhuss = NULL,
                            wBAIT = NULL,
                            baitPars = NULL) {
  # read cellular temperature
  temp <- importData(subtype = ftas)

  dates <- names(temp)

  # optional: transform raw temperature into BAIT
  if (bait) {
    # note: easier to do in [C]
    tempCelsius <- temp - 273.15   # [C]

    # read and prepare bait input data
    baitInput <- prepBaitInput(frsds, fsfc, fhuss) %>%
      checkDates(tempCelsius)

    # calculate bait
    tempBAIT <- compBAIT(baitInput, tempCelsius, weight = wBAIT, params = baitPars)   # [C]

    # convert back to [K]
    temp <- tempBAIT + 273.15   # [K]
  }

  # round and assign dates
  temp <- terra::round(temp, digits = 1)
  names(temp) <- dates

  # loop: type of degree day
  hddcdd <- do.call(rbind, lapply(c("HDD", "CDD"), function(typeDD) {
    # loop: threshold temperatures
    do.call(rbind, lapply(tlim[[typeDD]], function(t) {
      compCellHDDCDD(temp, typeDD, t, factors) %>%
        aggCells(pop, countries) %>%
        mutate("variable" = typeDD,
               "tlim"     = t)    # [C]
    }))
  }))

  return(hddcdd)
}



#' Assign HDD/CDD values for given ambient/limit temperature
#'
#' @param temp \code{terra::SpatRaster} containing temperature or BAIT values.
#' @param typeDD \code{character} string specifying the type of degree day to calculate: \code{"HDD"} or \code{"CDD"}.
#' @param tlim \code{numeric} limit temperature used as a threshold for degree day calculations.
#' @param factors \code{data.frame} containing degree day values for the combinations of \code{temp} and \code{tlim}.
#'
#' @returns \code{terra::SpatRaster} with calculated HDD or CDD values.
#'
#' @author Hagen Tockhorn
#'
#' @importFrom terra classify tapp time
#' @importFrom dplyr filter reframe .data
#' @importFrom magrittr %>%

compCellHDDCDD <- function(temp, typeDD, tlim, factors) {
  # extract years
  dates <- names(temp)

  # add tolerance of 0.049K to avoid machine precision errors
  factors <- factors[factors$typeDD == typeDD, ]

  factors <- factors %>%
    filter(.data[["tLim"]] == tlim) %>%
    dplyr::reframe(from =    .data[["T_amb_K"]] - 0.049,
                   to =      .data[["T_amb_K"]] + 0.049,
                   becomes = .data[["factor"]]) %>%
    data.matrix()

  # swap ambient temperature values with corresponding DD values
  hddcdd <- classify(temp, factors)

  terra::time(hddcdd) <- as.Date(dates)

  # aggregate to yearly HDD/CDD [K.d/a]
  hddcdd <- tapp(hddcdd, "years", fun = sum, na.rm = TRUE)

  names(hddcdd) <- gsub("y_", "", names(hddcdd))
  return(hddcdd)
}



#' Aggregate cellular HDD/CDD values to country-wide average (population-weighted)
#'
#' @param data \code{terra::SpatRaster} object containing HDD or CDD values.
#' @param weight \code{terra::SpatRaster} object containing weights for aggregation (i.e. population).
#' @param mask \code{terra::SpatRaster} object defining the regional aggregation boundaries.
#'
#' @returns \code{data.frame} containing regionally averaged HDD or CDD values.
#'
#' @author Hagen Tockhorn
#'
#' @importFrom terra subset global

aggCells <- function(data, weight, mask) {
  message("Aggregating degree days to regions...")

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

        aggData <- data.frame("region" = names(mask),
                              "period" = y,
                              "value"  = round(weightedAgg, 3))

        rownames(aggData) <- c()
        return(aggData)
      }
    )
  )
  return(hddcddAgg)
}
