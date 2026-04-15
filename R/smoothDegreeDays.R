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
#' @param smoothMethod Character string specifying the smoothing method: "lowpass" (default)
#'        for magclass lowpass filter, or "spline" for smooth spline interpolation.
#' @param splineSpar Numeric value for the spline smoothing parameter (0-1). Lower values
#'        produce closer fits to the data. Only used when smoothMethod = "spline".
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
#' @importFrom tidyr replace_na
#' @importFrom magclass lowpass
#' @importFrom purrr map2
#' @importFrom tidyr unnest
#' @importFrom stats smooth.spline predict

smoothDegreeDays <- function(data,
                             fileMapping,
                             ssp2Requested = TRUE,
                             nSmoothIter = 50,
                             smoothMethod = "lowpass",
                             splineSpar = 0.5,
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
                          .applySmoothing(x = .data[["value"]],
                                          smoothMethod = smoothMethod,
                                          splineSpar = splineSpar,
                                          nSmoothIter = nSmoothIter,
                                          periods = .data[["period"]]))) %>%
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
                                .applySmoothing(x = .data[["value"]],
                                                smoothMethod = smoothMethod,
                                                splineSpar = splineSpar,
                                                nSmoothIter = nSmoothIter,
                                                periods = .data[["period"]]),
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
    mutate(value = .data$value + replace_na(.data$delta, 0)) %>%
    select(-"delta")

  # 6) Filter data by period: keep historical <= endOfHistory, scenarios > endOfHistory
  #    After delta shift alignment, scenarios naturally diverge from historical endpoint
  dataSmooth <- rbind(
    dataEnsemble %>%
      filter(.data[["rcp"]] == "historical",
             .data[["period"]] <= endOfHistory),
    dataEnsemble %>%
      filter(.data[["rcp"]] != "historical",
             .data[["period"]] > endOfHistory)
  )

  return(dataSmooth)
}


#' Apply smoothing to a time series vector
#'
#' @param x Numeric vector to smooth
#' @param periods Optional period values for spline (same length as x)
#' @param smoothMethod A character specifying the smoothing method (either "lowpass" or "spline")
#' @param nSmoothIter An integer specifying the number of iterations for lowpass smoothing.
#' @param splineSpar Numeric value for the spline smoothing parameter (0-1). Lower values
#'        produce closer fits to the data. Only used when smoothMethod = "spline".
#'
#' @return Smoothed numeric vector

.applySmoothing <- function(x, smoothMethod, splineSpar, nSmoothIter, periods = NULL) {
  if (smoothMethod == "lowpass") {
    # Use magclass lowpass filter
    return(lowpass(x, i = nSmoothIter))
  } else if (smoothMethod == "spline") {
    # Use smooth spline
    if (is.null(periods)) {
      periods <- seq_along(x)
    }
    # Handle edge cases
    if (length(x) < 4) {
      return(x)  # Not enough points for spline
    }
    # Fit smooth spline
    splineFit <- smooth.spline(x = periods, y = x, spar = splineSpar)
    return(predict(splineFit, periods)$y)
  } else {
    stop("smoothMethod must be either 'lowpass' or 'spline'")
  }
}
