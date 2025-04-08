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
#' @param globalPars \code{logical} indicating whether to use global or gridded BAIT parameters.
#' @param noCC \code{logical} indicating whether to compute a no-climate-change scenario.
#'        If \code{TRUE}, the function will calculate degree days assuming constant climate conditions.
#'        Default is \code{FALSE}.
#' @param gridDataDir \code{character} (Optional) path to directory where grid data files will be stored
#'        when \code{noCC} is \code{TRUE}. Required when calculating no-climate-change scenarios.
#'        Default is \code{NULL}.
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
                            wBAIT = NULL,
                            globalPars = TRUE,
                            noCC = FALSE,
                            gridDataDir = NULL) {

  # ensure parameters of type logical are correctly passed
  bait       <- as.logical(bait)
  globalPars <- as.logical(globalPars)
  noCC       <- as.logical(noCC)


  # extract filenames
  fileNames <- c("tas" = fileMapping[["tas"]],
                 "gcm" = fileMapping[["gcm"]])
  if (bait) {
    fileNames <- c(
      fileNames,
      "rsds" = fileMapping[["rsds"]],
      "sfc"  = fileMapping[["sfcwind"]],
      "huss" = fileMapping[["huss"]]
    )
  }

  # extract temporal interval
  yStart <- fileMapping[["start"]] %>% as.numeric()
  yEnd   <- fileMapping[["end"]] %>% as.numeric()

  # extract RCP scenario + model
  rcp   <- fileMapping[["rcp"]] %>% unique()
  model <- fileMapping[["gcm"]] %>% unique()


  # read country masks
  countries <- importData(subtype = "countrymasks-fractional_30arcmin.nc")


  baitPars <- NULL

  if (isTRUE(bait)) {
    if (isFALSE(globalPars)) {
      # gridded bait regression parameters
      baitPars <- computeBAITpars(model = unique(fileMapping$gcm))
      names(baitPars) <- c("aRSDS", "bRSDS", "aSFC", "bSFC", "aHUSS", "bHUSS") # TODO: change this
    } else {
      # load default global parameters
      paramsMap <- read.csv2(getSystemFile("extdata", "mappings", "cfacBAITpars.csv",
                                           package = "climbed"))

      baitPars <- setNames(lapply(paramsMap$value, function(x) eval(parse(text = x))),
                           paramsMap$variable)
    }
  }




  # loop over single years and compute annual degree days
  hddcdd <- do.call("rbind", lapply(seq_len(yEnd - yStart + 1), function(i) {
    message("Initiating calculating degree days for the year: ", seq(yStart, yEnd)[[i]])

    # add year tag
    fileNames <- lapply(fileNames, function(x) sub("\\.nc$", paste0("_", i, ".nc"), x))

    compStackHDDCDD(fileNames = fileNames,
                    tlim = tLim,
                    pop = pop,
                    countries = countries,
                    factors = hddcddFactor,
                    bait = bait,
                    wBAIT  = wBAIT,
                    baitPars = baitPars,
                    noCC = noCC,
                    gridDataDir = gridDataDir)
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
#' @param fileNames \code{character} A named vector of file paths extracted from \code{fileMapping}. Elements include:
#'   \describe{
#'     \item{tas}{Temperature data file path.}
#'     \item{rsds}{(Optional) Solar radiation data file path (included if \code{bait} is \code{TRUE}).}
#'     \item{sfc}{(Optional) Surface wind data file path (included if \code{bait} is \code{TRUE}).}
#'     \item{huss}{(Optional) Humidity data file path (included if \code{bait} is \code{TRUE}).}
#'   }
#' @param tlim \code{list} named list of temperature limit sequences for \code{HDD} and \code{CDD}.
#' @param countries \code{terra::SpatRaster} object defining the regional aggregation boundaries.
#' @param pop \code{terra::SpatRaster} object containing population data.
#' @param factors \code{data.frame} containing degree day values for the combinations of \code{temp} and \code{tlim}.
#' @param bait \code{logical}; if \code{TRUE}, BAIT is used as the ambient temperature instead of
#' near-surface temperature.
#' @param wBAIT (Optional) \code{list} named list containing weights for BAIT calculations.
#' @param baitPars (Optional) \code{terra::SpatRaster} object containing regression parameters
#' generated by \code{calcBAITpars}.
#' @param noCC \code{logical} indicating whether to compute a no-climate-change scenario.
#'        If \code{TRUE}, the function will calculate degree days assuming constant climate conditions.
#'        Default is \code{FALSE}.
#' @param gridDataDir \code{character} (Optional) path to directory where grid data files will be stored
#'        when \code{noCC} is \code{TRUE}. Required when calculating no-climate-change scenarios.
#'        Default is \code{NULL}.
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

compStackHDDCDD <- function(fileNames, tlim, countries, pop, factors, bait,
                            wBAIT = NULL,
                            baitPars = NULL,
                            noCC = FALSE,
                            gridDataDir = NULL) {
  # read cellular temperature
  temp <- importData(subtype = fileNames[["tas"]])

  dates <- names(temp)

  # optional: transform raw temperature into BAIT
  if (isTRUE(bait)) {
    # note: easier to do in [C]
    tempCelsius <- temp - 273.15   # [C]

    # read and prepare bait input data
    baitInput <- prepBaitInput(fileNames) %>%
      checkDates(tempCelsius)

    # calculate bait
    tempBAIT <- compBAIT(baitInput, tempCelsius, weight = wBAIT, params = baitPars)   # [C]

    # convert back to [K]
    temp <- tempBAIT + 273.15   # [K]
  }

  # round and assign dates
  temp <- terra::round(temp, digits = 1)
  names(temp) <- dates

  # Create an empty list to store all layers for the current year if noCC == TRUE
  yearlyRasterLayers <- list()

  # loop: type of degree day
  hddcdd <- do.call(rbind, lapply(c("HDD", "CDD"), function(typeDD) {
    # loop: threshold temperatures
    do.call(rbind, lapply(tlim[[typeDD]], function(t) {
      # Compute annual degree days
      annualDegreeDays <- compCellHDDCDD(temp, typeDD, t, factors)

      if (isTRUE(noCC)) {
        # create and set appropriate layer name
        layerName <- paste0(typeDD, "_", t)
        yearlyRasterLayers[[layerName]] <<- annualDegreeDays
        names(yearlyRasterLayers[[layerName]]) <- layerName
      }

      # Continue with your aggregation as before, using original annualDegreeDays
      annualDegreeDays %>%
        aggCells(pop, countries) %>%
        mutate("variable" = typeDD,
               "tlim"     = t)    # [C]
    }))
  }))


  # Save the combined layers
  if (isTRUE(noCC) && length(yearlyRasterLayers) > 0) {
    currentYear <- hddcdd %>%
      pull("period") %>%
      unique()

    # combine all layers into a single SpatRaster
    combinedYearRaster <- rast(yearlyRasterLayers)

    # Save the combined raster with a filename that includes the year
    outputFileName <- paste0(tolower(fileNames[["gcm"]]), "_", currentYear, ".tif")
    outputFilePath <- file.path(gridDataDir, outputFileName)
    writeRaster(combinedYearRaster, outputFilePath, overwrite = TRUE)
  }

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
    reframe(from    = .data[["T_amb_K"]] - 0.049,
            to      = .data[["T_amb_K"]] + 0.049,
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
#' @param noCC \code{logical} flag indicating whether the data represents a no-climate-change scenario.
#'        When TRUE, the data is treated as a constant value applied to all time periods in weight.
#'
#' @returns \code{data.frame} containing regionally averaged HDD or CDD values.
#'
#' @author Hagen Tockhorn
#'
#' @importFrom terra subset global

aggCells <- function(data, weight, mask, noCC = FALSE) {
  # Determine which years to process
  if (isTRUE(noCC)) {
    # For constant climate, use all years from weight
    yearsToProcess <- names(weight)
  } else {
    # For normal case, use matching years
    yearsData <- names(data)
    yearsWeight <- names(weight)

    if (!all(yearsData %in% yearsWeight)) {
      stop("Time periods of raster file and aggregation weights do not match.")
    }

    yearsToProcess <- yearsData
  }

  # Process each year
  hddcddAgg <- do.call("rbind", lapply(yearsToProcess, function(y) {
    # Get data for this year (either subset or use constant data)
    yearData <- if (noCC) data else subset(data, y)

    # Mask data and weights to considered regions
    regData <- yearData * subset(weight, y) * mask
    regWeight <- subset(weight, y) * mask

    # Aggregate regional data
    regDataAgg <- terra::global(regData, "sum", na.rm = TRUE)$sum
    regWeightAgg <- terra::global(regWeight, "sum", na.rm = TRUE)$sum

    # Calculate weighted sum
    weightedAgg <- regDataAgg / regWeightAgg

    # Create result data frame
    aggData <- data.frame(
      "region" = names(mask),
      "period" = y,
      "value" = round(weightedAgg, 3)
    )

    rownames(aggData) <- c()
    return(aggData)
  }))
  return(hddcddAgg)
}
