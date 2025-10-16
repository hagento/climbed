#' Smooth Degree Days Time Series and Create Transition Period
#'
#' This function smooths climate model projections of Heating and Cooling Degree Days (HDDs and CDDs)
#' and establishes a smooth transition from historical observations to future projections.
#' Historical data points remain unchanged, while all scenarios transition seamlessly
#' from the same historical end value to their respective projected trajectories.
#'
#' @param data A data frame containing degree day projections.
#' @param fileMapping A data frame containing metadata for locating and processing degree day output files.
#' @param ssp2Requested Logical indicating whether SSP2 was explicitly requested by the user.
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
#' arrange summarise rename
#' @importFrom magclass lowpass
#' @importFrom purrr map2
#' @importFrom tidyr unnest

smoothDegreeDays <- function(data,
                             fileMapping,
                             ssp2Requested = TRUE,
                             nSmoothIter = 50,
                             transitionYears = 10,
                             nHistYears = 20,
                             endOfHistory = 2025,
                             noCC = FALSE,
                             predictTransition = FALSE) {

  # PROCESS DATA ---------------------------------------------------------------

  # Determine max period of true historical observations (before fill-up)
  maxHistPeriod <- max(data$period[data$rcp == "historical" & data$ssp == "historical"])

  # 1) Fill historical data up to endOfHistory with SSP2 RCP2.6 data
  dataFilled <- data %>%
    fillHistory(endOfHistory = endOfHistory)

  # 2) Smooth future projections per model to reduce impact of outliers
  #    Keep true historical observations unchanged
  dataSmoothedModel <- dataFilled %>%
    group_by(across(all_of(c("region", "variable", "tlim", "model", "rcp", "ssp")))) %>%
    arrange(.data$period) %>%
    mutate(value = ifelse(.data[["period"]] <= endOfHistory,
                          .data[["value"]],
                          lowpass(.data[["value"]], i = nSmoothIter))) %>%
    ungroup()

  # 3) Calculate ensemble mean across climate models
  dataEnsemble <- dataSmoothedModel %>%
    group_by(across(all_of(c("region", "variable", "tlim", "rcp", "ssp", "period")))) %>%
    summarise(value = mean(.data[["value"]], na.rm = TRUE), .groups = "drop")

  # 4) Smooth filled historical timeline
  #    Preserve original values for true historical period (before fill-up)
  dataEnsemble <- dataEnsemble %>%
    group_by(across(all_of(c("region", "variable", "tlim", "rcp", "ssp")))) %>%
    arrange(.data$period) %>%
    mutate(valueSmooth = ifelse(.data[["period"]] <= endOfHistory,
                                lowpass(.data[["value"]], i = nSmoothIter),
                                .data[["value"]]),
           value = ifelse(.data[["period"]] <= maxHistPeriod,
                          .data[["value"]],
                          .data[["valueSmooth"]])) %>%
    ungroup() %>%
    select(-"valueSmooth")

  # 5) Align all scenarios (including noCC) with historical endpoint via delta shift
  #    Calculate historical value at endOfHistory
  histAtEnd <- dataEnsemble %>%
    filter(.data$period == endOfHistory, .data$rcp == "historical") %>%
    select(all_of(c("region", "variable", "tlim", "value"))) %>%
    rename(histValue = "value")

  # Calculate delta for each scenario (including noCC) at endOfHistory
  alignDelta <- dataEnsemble %>%
    filter(.data$period == endOfHistory) %>%
    left_join(histAtEnd, by = c("region", "variable", "tlim")) %>%
    filter(.data$rcp != "historical") %>%
    mutate(delta = .data$histValue - .data$value) %>%
    select(all_of(c("region", "variable", "tlim", "rcp", "ssp", "delta")))

  # Apply delta shift to align all scenarios with historical endpoint
  dataEnsemble <- dataEnsemble %>%
    left_join(alignDelta, by = c("region", "variable", "tlim", "rcp", "ssp")) %>%
    mutate(value = ifelse(!is.na(.data$delta),
                          .data$value + .data$delta,
                          .data$value)) %>%
    select(-"delta")

  # 6) Filter data by period: keep historical <= endOfHistory, scenarios > endOfHistory
  dataSmooth <- rbind(
    dataEnsemble %>%
      filter(.data[["rcp"]] == "historical",
             .data[["period"]] <= endOfHistory),
    dataEnsemble %>%
      filter(.data[["rcp"]] != "historical",
             .data[["period"]] > endOfHistory)
  )

  # 7) Calculate transition predictions by averaging across all scenarios
  #    This creates a convergence point for all scenarios in early transition years
  transitionPreds <- dataSmooth %>%
    filter(.data$rcp != "historical",
           .data$period > endOfHistory,
           .data$period <= (endOfHistory + transitionYears)) %>%
    group_by(across(all_of(c("region", "variable", "tlim", "period")))) %>%
    reframe(prediction = mean(.data$value, na.rm = TRUE))

  # 8) Apply linear weighted transition from mean prediction to scenario-specific values
  #    Weight gradually shifts from prediction (at endOfHistory) to scenario value
  dataSmooth <- dataSmooth %>%
    left_join(transitionPreds, by = c("region", "variable", "tlim", "period")) %>%
    mutate(
      value = ifelse(
        .data[["period"]] > endOfHistory &
          .data[["period"]] <= (endOfHistory + transitionYears) &
          .data[["rcp"]] != "historical" &
          !is.na(.data[["prediction"]]) &
          !is.na(.data[["value"]]),
        # Linear blend: prediction + (scenario - prediction) * weight
        .data[["prediction"]] + (.data[["value"]] - .data[["prediction"]]) *
          ((.data[["period"]] - endOfHistory) / transitionYears),
        .data[["value"]]
      )
    ) %>%
    select(-"prediction")

  return(dataSmooth)


}
