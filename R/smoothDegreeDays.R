#' Smooth Degree Days Time Series and Create Transition Period
#'
#' This function smooths climate model projections of Heating and Cooling Degree Days (HDDs and CDDs)
#' and establishes a smooth transition from historical observations to future projections.
#' Historical data points remain unchanged, while all scenarios transition seamlessly
#' from the same historical end value to their respective projected trajectories.
#'
#' @param data A data frame containing degree day projections.
#' @param fileMapping A data frame containing metadata for locating and processing degree day output files.
#' @param nSmoothIter An integer specifying the number of iterations for lowpass smoothing.
#' @param transitionYears An integer specifying the number of years for the transition period.
#' @param nHistYears An integer specifying the number of years used for the linear regression.
#' @param endOfHistory An integer specifying the upper temporal limit for historical data.
#' @param noCC Logical indicating a no-climate-change scenario.
#' @param predictTransition Logical indicating whether to predict transition values (TRUE)
#'        or use mean-based transition (FALSE).
#'
#' @returns A data frame with smoothed degree day values and a seamless transition period
#' between historical and projected data.
#'
#' @author Hagen Tockhorn
#'
#' @importFrom dplyr filter mutate select anti_join group_by across all_of ungroup reframe left_join group_modify
#' @importFrom magclass lowpass
#' @importFrom purrr map2
#' @importFrom tidyr unnest

smoothDegreeDays <- function(data,
                             fileMapping,
                             nSmoothIter = 50,
                             transitionYears = 10,
                             nHistYears = 20,
                             endOfHistory = 2025,
                             noCC = FALSE,
                             predictTransition = FALSE) {

  # PROCESS DATA ---------------------------------------------------------------

  # extract list of files only used to fill up missing historical data points
  excludeEntries <- fileMapping %>%
    filter(.data$originalFile == FALSE) %>%
    mutate(period = map2(.data$start, .data$end, seq)) %>%
    select("model" = "gcm", "rcp", "period") %>%
    unnest("period")


  dataSmooth <- data %>%

    # fill up historical data
    fillHistory(endOfHistory = endOfHistory) %>%

    # remove entries only added for historical fill-up
    anti_join(excludeEntries, by = c("model", "rcp", "period")) %>%

    # smooth model data before averaging to avoid impact of outliers
    group_by(across(-all_of(c("period", "value")))) %>%
    mutate(value = ifelse(.data[["period"]] <= endOfHistory,
                          .data[["value"]],
                          lowpass(.data[["value"]], i = nSmoothIter))) %>%
    ungroup() %>%

    # take mean over models
    group_by(across(-all_of(c("model", "value")))) %>%
    reframe(value = mean(.data[["value"]])) %>%
    ungroup()


  # filter data w.r.t. periods
  dataSmooth <- rbind(dataSmooth %>%
                        filter(.data[["rcp"]] == "historical",
                               .data[["period"]] <= endOfHistory),
                      dataSmooth %>%
                        filter(.data[["rcp"]] != "historical",
                               .data[["period"]] > endOfHistory))


  # no transition in noCC case
  if (isTRUE(noCC)) {
    return(dataSmooth)
  }

  # Get transition values based on chosen approach
  if (predictTransition) {
    # APPROACH 1: Prediction-based transition
    transitionPreds <- data %>%
      filter(.data[["period"]] >= (endOfHistory - nHistYears),
             .data[["period"]] <= endOfHistory,
             .data[["rcp"]] == "historical",
             !is.na(.data[["value"]])) %>%

      group_by(across(all_of(c("region", "variable", "tlim", "model")))) %>%

      # smooth historical data to better obtain trends
      mutate(value = lowpass(.data[["value"]], i = nSmoothIter)) %>%

      # predict continuation of historical trends
      group_modify(~predTransitionValues(.x, endOfHistory, transitionYears)) %>%
      ungroup() %>%

      # average predictions across models
      group_by(across(-all_of(c("model", "prediction")))) %>%
      reframe(prediction = mean(.data[["prediction"]])) %>%
      ungroup()
  } else {
    # APPROACH 2: Mean-based transition
    transitionPreds <- dataSmooth %>%
      filter(.data$rcp != "historical",
             .data$period > endOfHistory,
             .data$period <= (endOfHistory + transitionYears)) %>%
      group_by(across(all_of(c("region", "variable", "tlim", "period")))) %>%
      reframe(prediction = mean(.data$value, na.rm = TRUE))
  }

  # Apply transition using consistent approach for both methods
  dataSmooth <- dataSmooth %>%
    left_join(transitionPreds, by = c("region", "variable", "tlim", "period")) %>%
    mutate(
      value = ifelse(
        # Condition: In transition period with prediction available
        .data[["period"]] > endOfHistory &
          .data[["period"]] <= (endOfHistory + transitionYears) &
          .data[["rcp"]] != "historical" &
          !is.na(.data[["prediction"]]) &
          !is.na(.data[["value"]]),

        # Apply linear weighted transition
        .data[["prediction"]] + (.data[["value"]] - .data[["prediction"]]) *
          ((.data[["period"]] - endOfHistory) / transitionYears),

        # else: keep original value
        .data[["value"]]
      )
    ) %>%
    select(-"prediction")

  return(dataSmooth)
}



#' Predict temporal trends via linear regression for smoothing transition between
#' historical and projection data point
#'
#' @param data (grouped) data frame with necessary columns \code{period} and \code{value}
#' @param endOfHistory upper temporal limit for historical data
#' @param transitionYears An integer specifying the number of years for the transition period
#' from historical observations to projections.
#'
#' @importFrom stats lm predict

predTransitionValues <- function(data, endOfHistory, transitionYears) {
  # Create prediction data frame
  predPeriods <- data.frame(
    period = seq(endOfHistory + 1, endOfHistory + transitionYears, 1)
  )

  # Linear fit
  fit <- lm(value ~ period, data = data)

  # Predict future values
  predPeriods$prediction <- predict(fit, newdata = predPeriods)

  return(predPeriods)
}
