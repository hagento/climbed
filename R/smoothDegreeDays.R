#' Smooth Degree Days Time Series and Create Transition Period
#'
#' This function smooths climate model projections of Heating and Cooling Degree Days (HDDs and CDDs)
#' and establishes a smooth transition from historical observations to future projections.
#' Historical data points remain unchanged, while all scenarios transition seamlessly
#' from the same historical end value to their respective projected trajectories.
#'
#' @param data A data frame containing degree day projections with the following required columns:
#'   \describe{
#'     \item{\code{"period"}}{The time period (e.g., year).}
#'     \item{\code{"value"}}{The degree day value.}
#'     \item{\code{"model"}}{The climate model name.}
#'     \item{\code{"rcp"}}{The Representative Carbon Pathway (RCP) scenario.}
#'     \item{\code{"region"}}{The geographical region.}
#'     \item{\code{"variable"}}{The degree day variable type (e.g., HDD, CDD).}
#'     \item{\code{"ssp"}}{The Shared Socioeconomic Pathway (SSP) scenario.}
#'     \item{\code{"tlim"}}{The temperature limit associated with the degree day calculation.}
#'   }
#'
#' @param nSmoothIter An integer specifying the number of iterations for lowpass smoothing.
#' Defaults to \code{50}.
#'
#' @param transitionYears An integer specifying the number of years for the transition period
#' from historical observations to projections. Defaults to \code{10}.
#'
#' @return A data frame with smoothed degree day values and a seamless transition period
#' between historical and projected data.
#'
#' @details
#' The function applies lowpass filtering to smooth model projections and uses a linear
#' interpolation to create a consistent transition period. The transition ensures that
#' all projections align with the historical data at the transition's start, diverging smoothly
#' into their respective future trajectories.
#'
#' @importFrom dplyr ungroup group_by across all_of mutate rowwise reframe filter
#' left_join select case_when
#' @importFrom magclass lowpass
#'
#' @author Hagen Tockhorn


smoothDegreeDays <- function(data, nSmoothIter = 50, transitionYears = 10) {

  # PARAMETERS -----------------------------------------------------------------

  # upper temporal threshold for historical data points
  endOfHistory <- 2020



  # PROCESS DATA ---------------------------------------------------------------

  dataSmooth <- data %>%

    # smooth model data before averaging to avoid impact of outliers
    group_by(across(-all_of(c("period", "value")))) %>%
    mutate(value = ifelse(.data[["period"]] <= endOfHistory |
                            .data[["rcp"]] == "picontrol",
                          .data[["value"]],
                          lowpass(.data[["value"]], i = nSmoothIter))) %>%
    ungroup() %>%

    # take mean over models
    group_by(across(-all_of(c("model", "value")))) %>%
    reframe(value = mean(.data[["value"]]))



  # get single historical value per region/variable
  lastHistValues <- dataSmooth %>%
    filter(.data[["period"]] == endOfHistory, !is.na(.data[["value"]])) %>%
    group_by(across(all_of(c("region", "variable")))) %>%
    reframe(lastHistValue = mean(.data[["value"]], na.rm = TRUE))

  # filter data w.r.t. periods
  dataSmooth <- rbind(dataSmooth %>%
                        filter(.data[["rcp"]] == "historical",
                               .data[["period"]] <= endOfHistory),
                      dataSmooth %>%
                        filter(.data[["rcp"]] != "historical",
                               .data[["period"]] >= endOfHistory))


  # smooth transition between last historical value and projections to minimize
  # deviations caused by the preceding smoothing
  dataSmooth <- dataSmooth %>%
    left_join(lastHistValues,
              by = c("region", "variable")) %>%
    mutate(
      value = case_when(
        .data[["period"]] >= endOfHistory &
          .data[["period"]] <= (endOfHistory + transitionYears) &
          rcp != "picontrol" & !is.na(.data[["lastHistValue"]]) & !is.na(.data[["value"]]) ~
          .data[["lastHistValue"]] + (.data[["value"]] - .data[["lastHistValue"]]) *
            ((.data[["period"]] - endOfHistory) / transitionYears),
        TRUE ~ .data[["value"]]
      )
    ) %>%
    select(-"lastHistValue")

  return(dataSmooth)
}
