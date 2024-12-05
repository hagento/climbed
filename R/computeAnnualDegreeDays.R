# --- FUNCTIONS TO CALCULATE ANNUAL DEGREE DAYS --------------------------------


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
#' @param pop SpatRaster with annual population data.
#' @param ssp Shared Socioeconomic Pathway scenario.
#' @param bait Logical indicating whether to use raw temperature or BAIT as ambient temperature.
#' @param tLim Temperature limits for degree day calculations.
#' @param countries SpatRaster containing the region masks.
#' @param hddcddFactor data.frame containing pre-computed degree days.
#' @param wBAIT (Optional) Weights for BAIT adjustments. Default is \code{NULL}.
#' @param baitPars (Optional) Parameters for BAIT calculations. Default is \code{NULL}.
#'
#' @return A data frame containing annual degree days.
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
  hddcdd <- do.call(
    "rbind",
    lapply(
      seq(1, yEnd - yStart + 1),
      function(i) {
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
      }
    )
  )

  hddcdd <- hddcdd %>%
    mutate("model" = model,
           "ssp" = ssp,
           "rcp" = rcp)

  return(hddcdd)
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
#' @param countries SpatRaster defining (regional) aggregation boundaries
#' @param pop SpatRaster containing population data
#' @param factors data frame with degree day values for \code{temp/tlim} combination
#' @param bait boolean, BAIT is used as ambient temperature
#' @param frsds file name of data on surface downdwelling shortwave radiation (optional)
#' @param fsfc file name of data on near-surface wind speed (optional)
#' @param fhuss file name of data on near-surface specific humidity (optional)
#' @param wBAIT named list containing BAIT weights (optional)
#' @param params SpatRaster containing regression parameters from \code{calcBAITpars} (optional)
#' @param rasDir absolute path to directory for saving raster files
#'
#' @return data frame containing regional population-weighted annual degree days
#'
#' @author Hagen Tockhorn
#'
#' @importFrom raster writeRaster
#' @importFrom stringr str_split
#' @importFrom terra writeCDF

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
  hddcdd <- do.call(
    "rbind", lapply(
      c("HDD", "CDD"), function(typeDD) {
        # loop: threshold temperatures
        do.call(
          "rbind", lapply(
            tlim[[typeDD]], function(t) {
              hddcddAgg <- compCellHDDCDD(temp, typeDD, t, factors)

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

  return(hddcdd)
}



#' Assign HDD/CDD values for given ambient/limit temperature
#'
#' @param temp SpatRaster containing temperature/BAIT values
#' @param typeDD type of degree day
#' @param tlim limit temperature
#' @param factors data frame with degree day values for \code{temp/tlim} combination
#'
#' @return SpatRaster with HDD/CDD values
#'
#' @author Hagen Tockhorn
#'
#' @importFrom terra classify tapp
#' @importFrom dplyr .data

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
#' @param data SpatRaster object containing HDD/CDD values
#' @param weight SpatRaster object containing aggregation weights
#' @param mask SpatRaster object defining (regional) aggregation boundaries
#'
#' @return data frame containing regionally averaged HDD/CDD values
#'
#' @author Hagen Tockhorn
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
