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
#' @param fileMapping A data frame containing metadata for locating and processing
#' degree day output files. The following columns are required:
#'   \describe{
#'     \item{\code{"gcm"}}{The General Circulation Model (GCM) name.}
#'     \item{\code{"rcp"}}{The RCP scenario (e.g., RCP2.6, RCP8.5).}
#'     \item{\code{"start"}}{The starting year of the time period.}
#'     \item{\code{"end"}}{The ending year of the time period.}
#'     \item{\code{"originalFile"}}{Logical indicating subset of original mapping.}
#'   }
#'
#' @param nSmoothIter An integer specifying the number of iterations for lowpass smoothing.
#' Defaults to \code{50}.
#'
#' @param transitionYears An integer specifying the number of years for the transition period
#' from historical observations to projections. Defaults to \code{10}.
#'
#' @param nHistYears An integer specifying the number of years used for the linear regression
#' to blend the continued historical trend into the scenario projections over a specified period
#' in \code{transitionYears}.
#'
#' @param endOfHistory An integer specifying the upper temporal limit for historical data.
#' Defaults to 2025.
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
#' @importFrom dplyr ungroup group_by across all_of mutate reframe filter
#' left_join select case_when group_modify
#' @importFrom magclass lowpass
#' @importFrom stats lm predict
#' @importFrom purrr map2
#' @importFrom tidyr unnest
#'
#' @author Hagen Tockhorn


smoothDegreeDays <- function(data,
                             fileMapping,
                             nSmoothIter = 50,
                             transitionYears = 10,
                             nHistYears = 10,
                             endOfHistory = 2025) {

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


  # fit linear regression to last n (= nHistYears) values per model and create
  # predictions for m (= transitionYears) periods
  transitionPreds <- data %>%
    filter(.data[["period"]] >= (endOfHistory - nHistYears),
           .data[["period"]] <= endOfHistory,
           .data[["rcp"]] == "historical",
           !is.na(.data[["value"]])) %>%

    group_by(across(all_of(c("region", "variable", "tlim", "model")))) %>%
    group_modify(~predTransitionValues(.x, endOfHistory, transitionYears)) %>%
    ungroup() %>%

    # average predictions across models
    group_by(across(-all_of(c("model", "prediction")))) %>%
    reframe(prediction = mean(.data[["prediction"]])) %>%
    ungroup()


  # filter data w.r.t. periods
  dataSmooth <- rbind(dataSmooth %>%
                        filter(.data[["rcp"]] == "historical",
                               .data[["period"]] <= endOfHistory),
                      dataSmooth %>%
                        filter(.data[["rcp"]] != "historical",
                               .data[["period"]] > endOfHistory))


  # merge historical and scenario data
  dataSmooth <- dataSmooth %>%
    left_join(transitionPreds,
              by = c("region", "variable", "period", "tlim")) %>%
    mutate(
      value = ifelse(
        # Condition: In transition period with both prediction and value available
        .data[["period"]] >= endOfHistory &
          .data[["period"]] <= (endOfHistory + transitionYears) &
          !is.na(.data[["prediction"]]) &
          !is.na(.data[["value"]]),

        # apply weighted transition
        {
          # calculate normalized position in transition (0 to 1)
          t <- (.data[["period"]] - endOfHistory) / transitionYears

          # use cosine-based weighting for smoother S-curve transition
          weight <- 0.5 * (1 - cos(t * pi))

          .data[["prediction"]] + (.data[["value"]] - .data[["prediction"]]) * weight
        },

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

predTransitionValues <- function(data, endOfHistory, transitionYears) {
  # Create prediction data frame
  predPeriods <- data.frame(
    period = seq(endOfHistory, endOfHistory + transitionYears, 1)
  )

  # Linear fit
  fit <- lm(value ~ period, data = data)

  # Predict future values
  predPeriods$prediction <- predict(fit, newdata = predPeriods)

  return(predPeriods)
}
