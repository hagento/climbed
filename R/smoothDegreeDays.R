#' Smooth Degree Days Time Series and Create Transition Period
#'
#' This function smooths climate model projections of Heating and Cooling Degree Days (HDDs and CDDs)
#' and establishes a smooth transition from historical observations to future projections.
#' Historical data points remain unchanged, while all scenarios transition seamlessly
#' from the same historical end value to their respective projected trajectories.
#'
#' @param data A data frame containing degree day projections.
#' @param splineSpar Numeric value for the spline smoothing parameter (0-1). Lower values
#'        produce closer fits to the data.
#' @param endWeight Numeric weight of the endOfHistory data point in the spline smoothing relative
#'        to all other points (weight 1).
#' @param endOfHistory An integer specifying the upper temporal limit for historical data.
#'
#' @returns A data frame with smoothed degree day values and a seamless transition period
#' between historical and projected data.
#'
#' @author Hagen Tockhorn
#'
#' @importFrom dplyr filter mutate select group_by across all_of ungroup left_join arrange summarise rename
#'  distinct n
#' @importFrom stats smooth.spline predict

smoothDegreeDays <- function(data,
                             splineSpar = 0.5,
                             endWeight = 1000,
                             endOfHistory = 2025) {

  # PROCESS DATA ---------------------------------------------------------------

  # Determine max period of true historical observations (before fill-up)
  maxHistPeriod <- max(data$period[data$rcp == "historical" & data$ssp == "historical"])

  # 1) Fill historical data up to endOfHistory with SSP2 RCP2.6 data
  dataFilled <- data %>%
    fillHistory(endOfHistory = endOfHistory)

  # 2) Calculate ensemble mean across climate models
  dataEnsemble <- dataFilled %>%
    group_by(across(all_of(c("region", "variable", "tlim", "rcp", "ssp", "period")))) %>%
    summarise(value = mean(.data[["value"]], na.rm = TRUE), .groups = "drop")

  # 3) Replace historical fill-up by smoothed SSP2 RCP2.6 projection
  #    Reference combines fill-up data and SSP2 RCP2.6 scenario (if requested), which are identical where overlapping
  #    True historical observations remain unchanged
  histEnsemble <- dataEnsemble %>%
    filter(.data$rcp == "historical", .data$period <= endOfHistory) %>%
    select(-"rcp", -"ssp")

  refEnsemble <- data %>%
    filter(.data$ssp == "ssp2", .data$rcp %in% c("historical", "2.6"), .data$period > endOfHistory) %>%
    group_by(across(all_of(c("region", "variable", "tlim", "period")))) %>%
    summarise(value = mean(.data[["value"]], na.rm = TRUE), .groups = "drop")

  histFillSmooth <- histEnsemble %>%
    rbind(refEnsemble) %>%
    group_by(across(all_of(c("region", "variable", "tlim")))) %>%
    arrange(.data$period) %>%
    mutate(valueSmooth = predict(smooth.spline(.data$period, .data$value, spar = splineSpar),
                                 .data$period)$y) %>%
    ungroup() %>%
    select(all_of(c("region", "variable", "tlim", "period", "valueSmooth")))

  histEnsemble <- histEnsemble %>%
    left_join(histFillSmooth, by = c("region", "variable", "tlim", "period")) %>%
    mutate(value = ifelse(.data$period <= maxHistPeriod, .data$value, .data$valueSmooth)) %>%
    select(-"valueSmooth")

  # 4) Smooth each scenario (excluding noCC) over the combined historical and scenario timeline
  #    Historical data up to endOfHistory is prepended to each scenario to avoid edge effects
  #    The endOfHistory point is weighted strongly to pull the smoothed curve through it
  #    smooth.spline normalizes weights to sum to the number of points, therefore lambda of the
  #    unweighted fit is rescaled to keep the smoothness of the remaining points unchanged
  dataScenSmooth <- dataEnsemble %>%
    filter(!.data$rcp %in% c("historical", "noCC")) %>%
    distinct(across(all_of(c("region", "variable", "tlim", "rcp", "ssp")))) %>%
    left_join(histEnsemble, by = c("region", "variable", "tlim"), relationship = "many-to-many") %>%
    rbind(dataEnsemble %>%
            filter(!.data$rcp %in% c("historical", "noCC"), .data$period > endOfHistory)) %>%
    group_by(across(all_of(c("region", "variable", "tlim", "rcp", "ssp")))) %>%
    arrange(.data$period) %>%
    mutate(w      = ifelse(.data$period == endOfHistory, endWeight, 1),
           lambda = smooth.spline(.data$period, .data$value, spar = splineSpar)$lambda * n() / sum(.data$w),
           value  = predict(smooth.spline(.data$period,
                                          .data$value,
                                          w = .data$w,
                                          lambda = .data$lambda[1]),
                            .data$period)$y) %>%
    ungroup() %>%
    filter(.data$period > endOfHistory) %>%
    select(-"w", -"lambda")

  # 5) Historical value at endOfHistory as common starting point of all scenarios
  histAtEnd <- histEnsemble %>%
    filter(.data$period == endOfHistory) %>%
    select(all_of(c("region", "variable", "tlim", "value"))) %>%
    rename(histValue = "value")

  # 6) Scale noCC scenario to historical endpoint
  #    Relative changes due to population weighting are preserved, no smoothing applied
  dataNoCC <- dataEnsemble %>%
    filter(.data$rcp == "noCC") %>%
    group_by(across(all_of(c("region", "variable", "tlim", "rcp", "ssp")))) %>%
    mutate(ratio = .data$value / .data$value[.data$period == endOfHistory][1]) %>%
    ungroup() %>%
    filter(.data$period > endOfHistory) %>%
    left_join(histAtEnd, by = c("region", "variable", "tlim")) %>%
    mutate(value = ifelse(is.finite(.data$ratio), .data$histValue * .data$ratio, .data$histValue)) %>%
    select(-"ratio", -"histValue")

  # 7) Combine historical data <= endOfHistory with aligned scenarios > endOfHistory
  dataSmooth <- rbind(
    histEnsemble %>%
      mutate(rcp = "historical", ssp = "historical"),
    dataScenSmooth,
    dataNoCC
  )

  return(dataSmooth)
}
