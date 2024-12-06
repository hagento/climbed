#' Align BAIT Input Variables with Near-Surface Temperature Data
#'
#' This function ensures that the time periods of BAIT input data
#' (e.g., shortwave radiation, surface wind, and specific humidity) align
#' with the time period of near-surface temperature data. Missing data in the
#' BAIT input are filled using yearly average values for the same day, derived
#' from the BAIT input data itself.
#'
#' @param baitInput A named list containing \code{terra::SpatRaster} objects
#'   for BAIT variables: \code{rsds} (shortwave radiation), \code{sfcwind} (surface wind),
#'   and \code{huss} (specific humidity).
#' @param tasData A \code{terra::SpatRaster} containing data on near-surface air temperature (\code{tas}).
#'
#' @return A named list of \code{terra::SpatRaster} objects with time periods aligned
#'   to the near-surface temperature data (\code{tasData}).
#'
#' @importFrom terra rast subset

checkDates <- function(baitInput, tasData) {
  # Extract time periods from temperature data
  datesT <- names(tasData)

  # Align each BAIT variable
  baitInput <- lapply(names(baitInput), function(var) {
    # Retrieve the data for the current BAIT variable
    baitLayer <- baitInput[[var]]
    datesBait <- names(baitLayer)

    # Determine dates to keep and fill
    datesFill <- setdiff(datesT, datesBait)  # Dates to fill
    datesKeep <- intersect(datesBait, datesT)  # Dates to keep
    daysFill <- unique(substr(datesFill, 6, 11))  # Unique day strings (MM-DD)

    # Subset to dates that are present
    if (length(datesKeep) > 0) {
      baitLayer <- subset(baitLayer, datesKeep)
      names(baitLayer) <- datesKeep
    }

    # Fill missing dates with yearly-average values
    if (length(daysFill) > 0) {
      baitInputMean <- prepBaitInput(fillWithMean = TRUE, baitInput = baitInput)
      baitFill <- rast(lapply(daysFill, function(d) {
        idx <- which(grepl(d, substr(datesFill, 6, 11)))
        filledRaster <- rast(replicate(length(idx), baitInputMean[[var]][[d]]))
        names(filledRaster) <- datesFill[idx]
        filledRaster
      }))

      # Combine existing and filled data
      baitLayer <- if (length(datesKeep) > 0) rast(list(baitLayer, baitFill)) else baitFill

      # Reorder dates
      baitLayer <- rast(baitLayer[[order(names(baitLayer))]])
    }

    # Alignment check
    if (!identical(names(baitLayer), names(tasData))) {
      warning(sprintf("Dates of Temperature and BAIT Input Data for '%s' are not aligned.", var))
    }

    return(baitLayer)
  })

  return(baitInput)
}
