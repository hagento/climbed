#' Fill up historical climate data with SSP2 RCP2.6 data
#'
#' This function fills gaps in historical climate data (up to endOfHistory) using
#' SSP2 RCP2.6 projections that have been marked with rcp="historical" for this purpose.
#' The fill-up data is relabeled to ssp="historical" and merged into the historical series.
#'
#' Note: Fill-up data comes from different GCMs than historical observations, so smoothing
#' must be done after model averaging in smoothDegreeDays(), not here.
#'
#' @param data data frame with necessary historical and scenario climate data
#' @param endOfHistory upper temporal limit for historical data
#'
#' @importFrom dplyr filter mutate anti_join %>%
#'
#' @return A data frame with filled historical climate data

fillHistory <- function(data, endOfHistory) {

  # check whether SSP2 RCP2.6 fill-up data exists (marked as rcp="historical")
  if (!any(data$ssp == "ssp2" & data$rcp == "historical")) {
    stop("SSP2 RCP2.6 fill-up data not found. Cannot fill missing historical data.")
  }

  # separate existing historical data
  histData <- data %>%
    filter(.data$rcp == "historical",
           .data$ssp == "historical")

  # identify SSP2 RCP2.6 data to fill up missing historical data points
  # (these were added in getDegreeDays with rcp="historical")
  transitionData <- data %>%
    filter(.data$ssp == "ssp2" & .data$rcp == "historical",
           .data$period <= endOfHistory) %>%
    mutate(ssp = "historical",
           rcp = "historical") %>%
    anti_join(histData, by = c("region", "period", "variable", "tlim"))

  # combine: non-historical data + original historical + fill-up data
  filledData <- data %>%
    filter(.data$rcp != "historical") %>%
    rbind(histData,
          transitionData)

  return(filledData)
}
