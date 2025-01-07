#' Smooth degree days time series and create transition period
#'
#' Smooths climate model projections of degree days and creates a smooth transition
#' from historical observations to projections. Historical data points remain unchanged.
#' All scenarios transition from the same historical end value to their respective trajectories.
#'
#' @param data Data frame containing degree day projections with columns:
#'   period, value, model, rcp, region, variable, ssp, tlim
#' @param nSmoothIter Number of iterations for lowpass smoothing (default: 50)
#' @param transitionYears Number of years for transition period (default: 10)
#'
#' @return Data frame with smoothed values and transition period
#'
#' @author Hagen Tockhorn
#'
#' @importFrom dplyr ungroup group_by across all_of mutate rowwise reframe filter
#' left_join select case_when
#' @importFrom magclass lowpass

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


  # filter

  return(dataSmooth)
}
