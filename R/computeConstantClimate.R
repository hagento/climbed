#' Calculate Constant Climate Heating and Cooling Degree Days
#'
#' This function imports pre-calculated heating and cooling degree days from historical periods
#' and uses them to project future degree days under a constant climate scenario (no climate change).
#' It processes grid data from historical years, calculates temporal averages, and then performs
#' population-weighted aggregation for different SSP scenarios. The result represents what
#' future heating and cooling demands would be if climate remained constant at historical levels.
#'
#' @param fileMapping A data frame mapping climate files to their time periods. The data frame
#' must include the following columns:
#'   \describe{
#'     \item{\code{"rcp"}}{The Representative Concentration Pathway scenario.}
#'     \item{\code{"start"}}{The starting year of the time period.}
#'     \item{\code{"end"}}{The ending year of the time period.}
#'   }
#'
#' @param ssp A character vector specifying the Shared Socioeconomic Pathway (SSP) scenarios
#' to process. Should include other scenarios than "historical".
#'
#' @param popMapping A named list mapping SSP scenarios to population data filenames. The names
#' of the list should correspond to the SSP scenarios provided in the \code{ssp} parameter.
#'
#' @param nHistYears An integer specifying the number of historical years to use for
#' calculating the constant climate conditions.
#'
#' @param gridDataDir A string specifying the path to the directory containing grid data files
#' in .tif format. The directory must exist and contain the required degree day data files.
#'
#' @return A data frame containing population-weighted aggregated degree day data for
#' constant climate scenarios, with columns for period, model, variable, temperature limit,
#' SSP scenario, and RCP scenario (set to "noCC" for no climate change).
#'
#' @author Hagen Tockhorn
#'
#' @importFrom dplyr filter select mutate pull
#' @importFrom magrittr %>%
#' @importFrom purrr map2
#' @importFrom tidyr unnest
#' @importFrom stringr str_extract
#' @importFrom terra rast app
#'
#' @export

computeConstantClimate <- function(fileMapping, ssp, popMapping, nHistYears, gridDataDir) {

  # READ IN DATA ---------------------------------------------------------------

  # check directory
  if (!dir.exists(gridDataDir)) {
    stop("Provided directory for stored grid data does not exist.")
  }

  # read country masks
  countries <- importData(subtype = "countrymasks-fractional_30arcmin.nc")

  # get relevant years
  relevantPeriods <- fileMapping %>%
    filter(.data$rcp == "historical") %>%
    select("start", "end") %>%
    mutate(period = map2(.data$start, .data$end, seq)) %>%
    unnest("period") %>%
    pull("period") %>%
    unique() %>%
    tail(nHistYears) %>%
    as.character()

  if (length(relevantPeriods) == 0) {
    stop("Relevant periods for constant climate calculation could not be determined from provided file mapping.")
  }

  # get relevant files
  fileList <- list.files(path = gridDataDir, pattern = "\\.tif$", full.names = TRUE)
  fileList <- fileList[grep(pattern = paste(relevantPeriods, collapse = "|"), fileList)]

  # import and organize data
  resultList <- .importNoCCData(fileList)



  # CALCULATE DEGREE DAYS ------------------------------------------------------

  if (length(ssp[ssp != "historical"]) == 0) {
    stop("Please define future SSP scenarios for constant climate calculation.")
  }


  # loop over SSP (population) scenarios
  hddcddNoCC <- do.call(rbind, lapply(ssp[ssp != "historical"], function(s) {
    # read in population data
    pop <- importData(subtype = popMapping[[s]])

    # process variables (HDD, CDD)
    do.call(rbind, lapply(names(resultList), function(variable) {
      # process temperature limits
      do.call(rbind, lapply(names(resultList[[variable]]), function(tlim) {

        # Combine all rasters into a single SpatRaster (across models and years)
        allRasters <- lapply(names(resultList[[variable]][[tlim]]), function(model) {
          resultList[[variable]][[tlim]][[model]]
        })

        # Stack all rasters and calculate the mean
        dataMean <- do.call(c, allRasters) %>%
          app(fun = mean, na.rm = TRUE)

        # aggregate using population weights
        aggData <- aggCells(dataMean, pop, countries, noCC = TRUE) %>%
          mutate(period = as.numeric(.data$period),
                 model = "model_average",
                 variable = variable,
                 tlim = tlim,
                 ssp = s,
                 rcp = "noCC")

        return(aggData)
      }))
    }))
  }))

  return(hddcddNoCC)
}



#' Import and Organize Degree Day Data for Constant Climate Calculation
#'
#' This internal helper function imports pre-calculated heating and cooling degree day data from
#' a list of .tif files and organizes it into a nested list structure by variable (HDD/CDD),
#' temperature limit, and climate model. The resulting structure allows for easily averaging
#' across all climate models.
#'
#' @param fileList A character vector of full file paths to .tif files containing degree day data.
#' The filenames should follow a specific pattern where the model name is at the beginning
#' of the basename (before the first underscore), and the year is represented as a 4-digit
#' number at the end of the filename (before .tif). Layer names within the raster files
#' should follow the pattern "variable_number", where variable typically represents HDD or CDD
#' and number represents the temperature threshold.
#'
#' @return A nested list organized by variable, temperature limit, and model, containing
#' SpatRaster objects with layers representing different years. The structure is:
#'   \describe{
#'     \item{\code{variable}}{The degree day variable (e.g., "HDD", "CDD")}
#'     \item{\code{number}}{The temperature limit identifier}
#'     \item{\code{model}}{The climate model name}
#'   }
#'
#' @author Hagen Tockhorn
#'
#' @importFrom terra rast nlyr
#' @importFrom stringr str_extract str_split
#' @importFrom magrittr %>%

.importNoCCData <- function(fileList) {
  # Initialize variable > tlim > model structure
  resultList <- list()

  # Process files
  for (file in fileList) {
    model <- basename(file) %>% str_extract("^[^_]+(?=_)")
    year <- as.numeric(str_extract(file, "\\d{4}(?=\\.tif$)"))

    # Read raster
    rastData <- rast(file)
    layerNames <- names(rastData)

    # Get unique variables and numbers
    varNumPairs <- lapply(layerNames, function(ln) {
      parts <- str_split(ln, "_", simplify = TRUE)
      list(variable = parts[1], number = parts[2])
    })

    variables <- unique(vapply(varNumPairs, function(x) x$variable, character(1)))

    for (variable in variables) {
      if (is.null(resultList[[variable]])) {
        resultList[[variable]] <- list()
      }

      # Filter pairs with matching variable
      matchingPairs <- varNumPairs[vapply(varNumPairs, function(x) x$variable == variable, logical(1))]

      # Extract numbers
      numbers <- unique(vapply(matchingPairs, function(x) x$number, character(1)))

      for (number in numbers) {
        if (is.null(resultList[[variable]][[number]])) {
          resultList[[variable]][[number]] <- list()
        }

        layerName <- paste0(variable, "_", number)
        if (layerName %in% layerNames) {
          singleLayer <- rastData[[layerName]]
          names(singleLayer) <- as.character(year)

          # Add to the model list for this variable/number combination
          if (is.null(resultList[[variable]][[number]][[model]])) {
            resultList[[variable]][[number]][[model]] <- singleLayer
          } else {
            # Check if existing raster has layers
            if (nlyr(resultList[[variable]][[number]][[model]]) > 0) {
              resultList[[variable]][[number]][[model]] <- c(resultList[[variable]][[number]][[model]], singleLayer)
            } else {
              resultList[[variable]][[number]][[model]] <- singleLayer
            }
          }
        }
      }
    }
  }

  return(resultList)
}
