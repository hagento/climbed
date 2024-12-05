#' Check if time period of BAIT input data (rsds, sfc, huss) is congruent with
#' near-surface temperature data (tas).
#'
#' @param baitInput list of raster data encompassing different climate variables
#' @param tasData raster data on near-surface atmosphere temperature
#'
#' @return baitInput with congruent time periods w.r.t. tasData
#'
#' @importFrom terra rast subset
#' @importFrom stringr str_sub

checkDates <- function(baitInput, tasData) {
  datesT <- names(tasData)

  baitInput <- vapply(
    names(baitInput),
    function(var) {
      # Fill missing data with means from previous years
      baitLayer <- baitInput[[var]]
      datesBait <- names(baitLayer)

      # Determine dates to fill and keep
      datesFill <- setdiff(datesT, datesBait) # Dates to fill up
      datesKeep <- intersect(datesBait, datesT) # Dates to keep
      daysFill  <- unique(substr(datesFill, 6, 11))

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
          r <- rast(replicate(length(idx), baitInputMean[[var]][[d]]))
          names(r) <- datesFill[idx]
          r
        }))

        # Combine existing and filled data
        baitLayer <- if (length(datesKeep) > 0) rast(list(baitLayer, baitFill)) else baitFill

        # Reorder dates
        baitLayer <- rast(baitLayer[[order(names(baitLayer))]])
      }

      # Alignment check
      if (!identical(names(baitLayer), names(tasData))) {
        warning("Dates of Temperature and BAIT Input Data are not aligned.")
      }

      baitLayer
    },
    FUN.VALUE = rast(),
    USE.NAMES = TRUE
  )

  return(baitInput)
}
