#' Calculate HDD/CDD values for different ambient/limit temperature combinations
#'
#' HDD/CDD values are pre-calculated for an interval \code{tLow}-\code{tUp} and
#' for a set of limit temperatures \code{tLim} with a temperature resolution of
#' 0.1C.
#'
#' The respective heating/cooling degree days are calculated as the difference
#' between the assumed ambient and a limit temperature, aggregated to a full day.
#' The latter defines a threshold above/below which cooling/heating is assumed to
#' be initiated.
#'
#' To account for heterogenity in heating/cooling behavior, the ambient and limit
#' temperatures, \code{tAmb} and \code{tLim}, are assumed to be normally distributed.
#' This changes the calculation of a degree day to a double integration of
#' \code{tLimit - T_ambient_day} with integration boundaries set at 3 standard
#' deviations, \code{tAmbStd} and \code{tLimStd}, from \code{tAmb} and \code{tLim}
#' respectively.
#'
#' As consequence, the ramp function of \code{HDD_day = max(0, tLimit - T_ambient_day)}
#' changes to a curved function that is above zero even if the mean of \code{T_ambient_day}
#' is above the mean of \code{tLimit}.
#'
#' @param tLow lower temperature boundary
#' @param tUp upper temperature boundary
#' @param tLim named list of limit temperature sequences for \code{HDD} and \code{CDD}
#' @param tAmbStd std of ambient temperature
#' @param tLimStd std of limit temperature
#'
#' @returns data frame of HDD/CDD values
#'
#' @author Hagen Tockhorn
#'
#' @importFrom stats dnorm
#' @importFrom pracma integral2

compDegDayFactors <- function(tLow, tUp, tLim, tAmbStd = 2, tLimStd = 2) {
  # t1 : ambient temperature variable
  # t2 : limit temperature variable

  temperatureFactor <- function(t2, t1, tAmb, tAmbStd, tLim, tLimStd, type = "HDD") {
    difference <- if (type == "HDD") (t2 - t1) else (t1 - t2)
    h <- dnorm(t2, mean = tLim, sd = tLimStd) * dnorm(t1, mean = tAmb, sd = tAmbStd) * difference
    return(h)
  }

  t <- seq(tLow, tUp, .1)

  hddcddFactors <-
    do.call("rbind", lapply(c("HDD", "CDD"), function(typeDD) {
      do.call("rbind", lapply(t, function(tAmb) {
        do.call("rbind", lapply(tLim[[typeDD]], function(.tLim) {
          # tLim integration boundaries
          x1 <- .tLim - 3 * tLimStd
          x2 <- .tLim + 3 * tLimStd

          switch(typeDD,
                 HDD = {
                   ymin <- tAmb - 3 * tAmbStd
                   ymax <- min(.tLim, tAmb + 3 * tAmbStd)
                 },
                 CDD = {
                   ymin <- max(.tLim, tAmb - 3 * tAmbStd)
                   ymax <- tAmb + 3 * tAmbStd
                 })

          f <- integral2(temperatureFactor,
                         xmin = x1,
                         xmax = x2,
                         ymin = ymin,
                         ymax = ymax,
                         tAmb = tAmb,
                         tAmbStd = tAmbStd,
                         tLim = .tLim,
                         tLimStd = tLimStd,
                         type = typeDD,
                         reltol = 1e-1)

          tmp <- data.frame("T_amb"        = tAmb,
                            "T_amb_K"      = round(tAmb + 273.15, 1),
                            "tLim"         = .tLim,
                            "factor"       = f$Q,
                            "factor_err"   = f$error,
                            "typeDD"       = typeDD)

          return(tmp)
        }))
      }))
    }))
  return(hddcddFactors)
}
