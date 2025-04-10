#' Fill up historical climate data with RCP2.6
#'
#' @param data data frame with necessary historical and scenario climate data
#' @param endOfHistory upper temporal limit for historical data
#'
#' @importFrom dplyr filter mutate anti_join %>%
#'
#' @return A data frame with filled historical climate data

fillHistory <- function(data, endOfHistory) {

  # check whether the scenario SSP2-2.6 exists
  if (!any(data$ssp == "ssp2" & data$rcp == "2.6")) {
    stop("Please provide necessary climate data on SSP2-2.6 to fill missing historical data.")
  }

  # separate existing historical data
  histData <- data %>%
    filter(.data$rcp == "historical")

  # identify data to fill up missing historical data points
  transitionData <- data %>%
    filter(.data$ssp == "ssp2" & .data$rcp == "2.6",
           .data$period <= endOfHistory) %>%
    mutate(ssp = "historical",
           rcp = "historical") %>%
    anti_join(histData, by = c("region", "period", "variable", "tlim"))

  # fill data
  filledData <- data %>%
    filter(.data$rcp != "historical") %>%
    rbind(histData,
          transitionData)

  return(filledData)
}
