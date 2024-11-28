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
  
  baitInput <- sapply(names(baitInput), function(var) { # nolint
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