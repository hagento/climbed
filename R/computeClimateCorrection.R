#' Calculate Climate Correction Factors for REMIND Buildings
#'
#' This function implements a climate correction mechanism for building energy demand
#' in the REMIND integrated assessment model. It processes CO2 concentration trajectories
#' from climate assessment data, compares them with reference scenarios, and calculates
#' correction factors that adjust building energy demand based on climate-driven changes.
#'
#' The function expects to be run from a REMIND output directory containing the following files:
#' \itemize{
#'   \item{\code{cfg.txt}}: REMIND configuration file with SSP and RCP scenario information
#'   \item{\code{fulldata_postsolve.gdx}}: GDX file containing baseline building energy demands
#'   \item{\code{p15_climate.gdx}}: GDX file with CO2 concentration projections (p15_co2_conc parameter)
#' }
#'
#' The correction factors are applied to \code{vm_cesIO} variables in REMIND's buildings
#' sector to account for climate-driven changes in heating and cooling demands. The function
#' outputs a GDX file (\code{pm_ClimateCorrection.gdx}) with correction factors by region,
#' time period, and building energy carrier.
#'
#' All data is given for REMIND regions. "EDGE-B resolution" refers to EDGE-B
#' end-use/carrier variables (e.g. \code{space_heating_elecHP_fe}), "REMIND resolution"
#' to REMIND building items (\code{ppfen_buildings_dyn36}, e.g. \code{feelhpb}).
#'
#' @return Invisibly returns \code{NULL}. The function is called for its side effect of
#' writing the \code{pm_ClimateCorrection.gdx} file to the current working directory.
#'
#' @author Hagen Tockhorn
#'
#' @importFrom utils read.csv
#' @importFrom dplyr %>% left_join group_by summarise filter mutate .data
#' @importFrom madrat toolGetMapping
#' @import mrremind
#'
#' @export

computeClimateCorrection <- function() {

  endOfHistory <- 2025

  thermalVariables <- c("space_cooling_elec",
                        "space_heating_biomod",
                        "space_heating_biotrad",
                        "space_heating_coal",
                        "space_heating_elecHP",
                        "space_heating_elecRH",
                        "space_heating_natgas",
                        "space_heating_petrol",
                        "space_heating_heat")

  # LOAD MAPPINGS --------------------------------------------------------------

  # Load RCP scenario mapping
  rcpMap <- read.csv(system.file("extdata", "mappings", "rcpMapping.csv", package = "climbed"))

  # Load REMIND names mapping
  remindNamesMapping <- read.csv(system.file("extdata", "mappings", "remindNamesMapping.csv", package = "climbed"))

  # Load structure mapping from mrcommons
  structureMap <- toolGetMapping("mappingEDGEBuildingsToREMIND.csv",
                                 type = "sectoral",
                                 where = "mrremind")


  # Create column mapping once for reuse
  columnMap <- setNames(remindNamesMapping$remindName, remindNamesMapping$targetName)



  # FILE NAMES -----------------------------------------------------------------

  # Check file availability and get absolute paths
  filePaths <- .checkFileAvailability()

  # Output file name
  outputFile <- "pm_ClimateCorrection.gdx"



  # READ-IN DATA ---------------------------------------------------------------

  # Extract current CO2 concentrations from REMIND/MAGICC output
  currentConc <- .extractCO2Concentrations(filePaths$climateDataFile)

  # Load reference concentration trajectories
  referenceConc <- .loadReferenceConcentrations(rcpMap)

  # Load baseline demand trajectories from GDX
  baselineDemand <- .loadBaselineDemands(filePaths$gdxInputFile, columnMap)

  # Load config data
  configData <- .extractConfigData(filePaths$configFile)



  # PROCESS DATA ---------------------------------------------------------------

  # Filter thermal baseline demands (EDGE-B resolution)
  thermalBaselines <- baselineDemand %>%
    filter(.data$variable %in% paste0(thermalVariables, "_fe"),
           .data$ssp == configData$ssp) %>%
    mutate(period = as.numeric(as.character(.data$period)))

  # Interpolate thermal demands based on CO2 concentrations (EDGE-B resolution)
  interpolatedThermalDemand <- .interpolateDemandFromConcentrations(currentConc,
                                                                    referenceConc,
                                                                    thermalBaselines,
                                                                    endOfHistory)

  # Calculate thermal correction factors as ratio of interpolated to baseline demand (EDGE-B resolution)
  thermalCorrectionFactors <- .calculateClimateCorrectionFactor(interpolatedThermalDemand, thermalBaselines)

  # Calculate weighted correction factors per REMIND item (REMIND resolution)
  weightedCorrectionFactors <- .calculateWeightedFactors(structureMap,
                                                         baselineDemand,
                                                         interpolatedThermalDemand,
                                                         thermalCorrectionFactors,
                                                         configData)

  # Complete missing entries (e.g., feh2b not in structure mapping) with 1.0 (REMIND resolution)
  completeCorrectionFactors <- .completeClimateCorrection(weightedCorrectionFactors, filePaths$gdxInputFile)



  # OUTPUT ---------------------------------------------------------------------

  # Write to REMIND-compatible GDX format
  .writeToREMIND(completeCorrectionFactors, outputFile, columnMap)

}



#' Check Availability of Required Input Files
#'
#' This function verifies that all required input files exist in the current working
#' directory. It checks for the REMIND configuration file, GDX input file, and
#' climate assessment Excel file. The function stops execution if required files
#' are missing.
#'
#' @returns A named list with absolute paths to required files (\code{configFile},
#' \code{gdxInputFile}, \code{climateDataFile}).
#'
.checkFileAvailability <- function() {
  # Get current working directory
  cwd <- getwd()

  files <- list(
    configFile = file.path(cwd, "cfg.txt"),
    gdxInputFile = file.path(cwd, "fulldata_postsolve.gdx"),
    climateDataFile = file.path(cwd, "p15_climate.gdx")
  )

  for (fileName in names(files)) {
    if (!file.exists(files[[fileName]])) {
      stop(fileName, " not found in current directory: ", files[[fileName]])
    }
  }

  return(files)
}


#' Extract CO2 Concentrations from Climate Data
#'
#' This function extracts CO2 concentration trajectories from a GDX climate data file, which contains
#' global CO2 concentrations over time.
#'
#' @param climateDataFile A string specifying the path to the climate data GDX file.
#'
#' @returns A data frame with CO2 concentrations by year.
#'
#' @importFrom gamstransfer readGDX
#' @importFrom dplyr %>% mutate select .data

.extractCO2Concentrations <- function(climateDataFile) {
  climateDataFile %>%
    readGDX("p15_co2_conc")$p15_co2_conc$records %>%
    mutate(period = as.numeric(as.character(.data$period))) %>%
    select("period", "value")
}


#' Fetch Reference CO2 Concentration Files from ISIMIP3b
#'
#' Downloads annual CO2 concentration files for all supported SSP and historical
#' scenarios from the ISIMIP3b data servers into a local cache directory. Files
#' that already exist in \code{dirPath} are skipped. Historical and SSP1-2.6,
#' SSP3-7.0, SSP5-8.5 are sourced from \code{InputData}; SSP1-1.9, SSP2-4.5,
#' and SSP4-6.0 from \code{SecondaryInputData}.
#'
#' @param dirPath Path to the local directory where files are cached.
#'
#' @returns A character vector of absolute paths to all \code{.txt} files in
#'   \code{dirPath} after downloading.
#'
#' @references Büchner, M., Reyer, C.P.O. (2022): ISIMIP3b atmospheric
#'   composition input data (v1.1). ISIMIP Repository.
#'   \doi{10.48364/ISIMIP.482153.1}
#'
#' @importFrom utils download.file

.fetchReferenceConcentrations <- function(dirPath) {
  inputDataBase <- "https://files.isimip.org/ISIMIP3b/InputData/climate/atmosphere_composition/co2"
  secondaryBase <- "https://files.isimip.org/ISIMIP3b/SecondaryInputData/climate/atmosphere_composition/co2"

  scenarios <- list(
    list(scenario = "historical", baseUrl = inputDataBase, start = "1850", end = "2014"),
    list(scenario = "ssp126",     baseUrl = inputDataBase, start = "2015", end = "2100"),
    list(scenario = "ssp370",     baseUrl = inputDataBase, start = "2015", end = "2100"),
    list(scenario = "ssp585",     baseUrl = inputDataBase, start = "2015", end = "2100"),
    list(scenario = "ssp119",     baseUrl = secondaryBase, start = "2015", end = "2100"),
    list(scenario = "ssp245",     baseUrl = secondaryBase, start = "2015", end = "2100"),
    list(scenario = "ssp460",     baseUrl = secondaryBase, start = "2015", end = "2100")
  )

  for (s in scenarios) {
    filename <- paste0("co2_", s$scenario, "_annual_", s$start, "_", s$end, ".txt")
    destPath <- file.path(dirPath, filename)
    if (!file.exists(destPath)) {
      url <- paste(s$baseUrl, s$scenario, filename, sep = "/")
      message("Downloading ", filename, " from ISIMIP...")
      download.file(url, destPath, quiet = TRUE)
    }
  }

  list.files(dirPath, "\\.txt$", full.names = TRUE)
}


#' Load Reference CO2 Concentration Scenarios
#'
#' Fetches any missing scenario files via \code{.fetchReferenceConcentrations},
#' then reads all \code{.txt} files from the cache directory and combines them
#' into a single data frame. The \code{"none"} scenario (no climate change) is
#' appended as a constant equal to the mean 2025 concentration across all SSP scenarios.
#'
#' @param rcpMap A data frame mapping SSP scenario names to RCP codes.
#'
#' @returns A data frame with columns \code{period}, \code{value}, and \code{rcp}.
#'
#' @importFrom stringr str_split
#' @importFrom utils read.table
#' @importFrom dplyr %>% mutate .data
#' @importFrom piamutils getSystemFile
#' @importFrom stats setNames

.loadReferenceConcentrations <- function(rcpMap) {
  scenMap <- setNames(rcpMap$rcp_code, rcpMap$ssp_scenario)

  dirPath <- getSystemFile("extdata", "reference", "co2_ppm", package = "climbed")
  files   <- .fetchReferenceConcentrations(dirPath)

  conc <- do.call("rbind", lapply(files, function(f) {
    read.table(f, col.names = c("period", "value")) %>%
      mutate(rcp = scenMap[[str_split(basename(f), "_")[[1]][[2]]]])
  }))

  noneValue <- mean(conc$value[conc$period == 2025])
  rbind(conc, data.frame(period = 2025:2100, value = noneValue, rcp = "none"))
}


#' Load Baseline Demand Data from GDX
#'
#' This function loads baseline building energy demand data from a GDX file and
#' renames columns according to a provided mapping. The data is extracted from
#' the \code{f_fedemandBuild} parameter.
#'
#' @param gdxInputFile A string specifying the path to the GDX input file.
#'
#' @param columnMap A named vector mapping REMIND column names to target names.
#'
#' @returns A data frame with baseline building energy demands at EDGE-B resolution
#' for all SSP and RCP scenarios.
#'
#' @importFrom gamstransfer readGDX
#' @importFrom dplyr rename

.loadBaselineDemands <- function(gdxInputFile, columnMap) {
  data <- readGDX(gdxInputFile, "f_fedemandBuild")

  # Extract the data frame from the nested list structure
  data <- data$f_fedemandBuild$records
  # Filter columnMap to only include columns that exist in the data
  existingColumnMap <- columnMap[unname(unlist(columnMap)) %in% colnames(data)]
  data <- rename(data, !!!existingColumnMap)

  return(data)
}


#' Extract Scenario Configuration Data from REMIND Config File
#'
#' This function parses a REMIND configuration file to extract scenario parameters.
#' It searches for \code{cm_demScen} (SSP scenario) and \code{cm_rcp_scen} (RCP scenario)
#' settings in the config file.
#'
#' @param configFile A string specifying the path to the REMIND configuration file (\code{cfg.txt}).
#'
#' @returns A named list with scenario parameters (\code{ssp}, \code{rcp}).
#'
#' @importFrom yaml read_yaml

.extractConfigData <- function(configFile) {
  config <- read_yaml(configFile)

  setNames(
    list(ssp = config$cm_demScen,
         rcp = config$cm_rcp_scen),
    c("ssp", "rcp")
  )
}


#' Interpolate Demand from CO2 Concentrations
#'
#' This function interpolates thermal building energy demands based on current CO2
#' concentrations and reference RCP scenarios. For each future time period, it identifies
#' the two bounding RCP scenarios and calculates interpolation weights based on the
#' current concentration's position between them. Historical data is passed through unchanged.
#'
#' @param currentConc A data frame with current global CO2 concentrations.
#'
#' @param referenceConc A data frame with global reference RCP concentration scenarios.
#'
#' @param thermalBaselines A data frame with baseline thermal demands for different RCP scenarios
#' at EDGE-B resolution.
#'
#' @param endOfHistory An integer specifying the last year of the historical period.
#' Defaults to 2025. Data up to this year is not interpolated.
#'
#' @returns A data frame with interpolated thermal demands at EDGE-B resolution. Historical periods contain
#' mean values across scenarios; future periods contain weighted interpolations based
#' on CO2 concentrations.
#'
#' @importFrom dplyr %>% filter left_join group_by summarise mutate select right_join case_when across all_of .data
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom stats approx

.interpolateDemandFromConcentrations <- function(currentConc,
                                                 referenceConc,
                                                 thermalBaselines,
                                                 endOfHistory = 2025) {

  # Future periods: for each (period, region, variable) the reference RCPs give
  # a set of (concentration, demand) points. The demand at the current
  # concentration is their piecewise-linear interpolation, clamped at the ends.
  future <- thermalBaselines %>%
    filter(.data$period > endOfHistory) %>%
    left_join(referenceConc %>%
                rename("conc" = "value"),
              by = c("period", "rcp")) %>%
    left_join(currentConc %>%
                rename("concCurrent" = "value"),
              by = "period") %>%
    filter(!is.na(.data$conc)) %>%
    group_by(across(all_of(c("period", "region", "variable")))) %>%
    summarise(
      value = approx(x = .data$conc,
                     y = .data$value,
                     xout = .data$concCurrent[1],
                     rule = 2)$y,
      .groups = "drop"
    )

  # Historical periods: scenarios share the same history, so collapse to mean.
  history <- thermalBaselines %>%
    filter(.data$period <= endOfHistory) %>%
    group_by(across(all_of(c("period", "region", "variable")))) %>%
    summarise(value = mean(.data$value), .groups = "drop")

  rbind(history, future)
}


#' Calculate Climate Correction Factors
#'
#' This function calculates climate correction factors as the ratio of climate-driven
#' interpolated demand to no-climate-change (noCC) baseline demand. The noCC baseline
#' represents the demand used for REMIND calibration.
#'
#' @param interpolatedThermalDemand A data frame with interpolated thermal demands based on
#' CO2 concentrations at EDGE-B resolution.
#'
#' @param thermalBaselines A data frame with baseline thermal demands at EDGE-B resolution.
#' Must include a scenario with \code{rcp == "none"} representing no climate change.
#'
#' @returns A data frame with thermal correction factors (ratio of interpolated to noCC demand)
#' at EDGE-B resolution.
#'
#' @importFrom dplyr %>% filter select left_join mutate .data

.calculateClimateCorrectionFactor <- function(interpolatedThermalDemand, thermalBaselines) {

  # extract noCC baseline demands (used for REMIND calibration)
  noCC <- thermalBaselines %>%
    filter(.data$rcp == "none") %>%
    select("region", "period", "variable", "value")

  # calculate climate correction factor as ratio of (climate-driven) interpolated to noCC demand
  interpolatedThermalDemand %>%
    select("region", "period", "variable", "value") %>%
    left_join(noCC, by = c("region", "period", "variable"), suffix = c("Interpolated", "NoCC")) %>%
    mutate(value = .data$valueInterpolated / .data$valueNoCC,
           variable = sub("_ue", "", .data$variable),
           .keep = "unused")
}


#' Calculate Weighted Correction Factors for all REMIND Variables
#'
#' This function aggregates disaggregated thermal correction factors to REMIND's
#' building energy carrier level using baseline demand as weights. Non-thermal
#' end uses are assigned a correction factor of unity (1.0).
#'
#' @param structureMap A data frame containing the structure mapping from EDGE-B to REMIND resolution.
#'
#' @param baselineDemand A data frame with baseline building energy demands at EDGE-B resolution.
#'
#' @param interpolatedThermalDemand A data frame with interpolated thermal demands based on
#' CO2 concentrations at EDGE-B resolution.
#'
#' @param thermalCorrectionFactors A data frame with disaggregated correction factors for thermal
#' end uses at EDGE-B resolution.
#'
#' @param configData A list with scenario configuration data (e.g., SSP, RCP) extracted from the REMIND config file.
#'
#' @returns A data frame with weighted correction factors at REMIND resolution.
#' Values beyond 2100 are extrapolated using the 2100 value.
#'
#' @importFrom dplyr %>% select filter group_by ungroup mutate left_join reframe rename anti_join across all_of .data
#' @importFrom tidyr replace_na
#' @importFrom quitte interpolate_missing_periods

.calculateWeightedFactors <- function(structureMap,
                                      baselineDemand,
                                      interpolatedThermalDemand,
                                      thermalCorrectionFactors,
                                      configData) {
  map <- structureMap %>%
    # filter buildings variables
    select("EDGE_buildings_items", "REMINDitems_out") %>%
    filter(grepl("b$", .data$REMINDitems_out),
           .data$EDGE_buildings_items != "",
           grepl("_fe", .data$EDGE_buildings_items))

  baselineDemand %>%
    # select noCC baseline demands for non-thermal end uses
    filter(.data$rcp == "none",
           .data$ssp == configData$ssp) %>%
    select("region", "period", "variable", "value") %>%
    mutate(period = as.numeric(as.character(.data$period))) %>%

    # remove thermal baseline demands and add interpolated ones
    anti_join(interpolatedThermalDemand, by = c("period", "region", "variable")) %>%
    rbind(interpolatedThermalDemand) %>%

    # map EDGE to REMIND variables
    right_join(map, by = c("variable" = "EDGE_buildings_items"), relationship = "many-to-many") %>%

    # join disaggregated thermal correction factors and set non-thermal to unity
    left_join(thermalCorrectionFactors, by = c("region", "period", "variable"), suffix = c("Demand", "Factor")) %>%
    replace_na(list("valueFactor" = 1)) %>%

    # weighted aggregation of correction factors per REMIND variable
    group_by(across(all_of(c("region", "period", "REMINDitems_out")))) %>%
    summarise(value = sum(.data$valueFactor * .data$valueDemand) / sum(.data$valueDemand)) %>%
    ungroup() %>%

    # fix >2100 correction factors at 2100 values
    mutate(value = ifelse(.data$period > 2100, NA, .data$value)) %>%
    interpolate_missing_periods(expand.values = TRUE) %>%

    filter(!is.na(.data$value))
}


#' Complete Climate Correction Factors with Missing Entries
#'
#' This function ensures all ppfen_buildings_dyn36 entries are present in the
#' correction factors. Missing entries (e.g., hydrogen which doesn't have thermal
#' correction) are filled with 1.0 to preserve REMIND's calibrated values.
#' Entries outside the REMIND domain (\code{ttot} x \code{all_regi} x
#' \code{ppfen_buildings_dyn36}) are dropped.
#'
#' @param correctionFactors A data frame with weighted correction factors at REMIND
#' resolution (REMIND regions, REMIND building items).
#' Must have columns: period, region, REMINDitems_out, value.
#'
#' @param gdxPath A string specifying the path to the REMIND GDX file
#' (\code{fulldata_postsolve.gdx}) from which the domain sets are read.
#'
#' @returns A complete data frame at REMIND resolution with all
#' region-period-carrier combinations, with missing values set to 1.0.
#'
#' @importFrom dplyr %>% left_join mutate .data
#' @importFrom tidyr replace_na expand_grid
#' @importFrom gamstransfer readGDX

.completeClimateCorrection <- function(correctionFactors, gdxPath) {
  # Read domain sets - handle different set types correctly
  # ttot is a subset with domain "tall", read from first column
  ttot <- readGDX(gdxPath, "ttot")$ttot$records[[1]]
  # all_regi has domain "*" (none), read from "uni" column
  all_regi <- readGDX(gdxPath, "all_regi")$all_regi$records$uni
  # ppfen_buildings_dyn36 is a subset with domain "all_in", read from first column
  ppfen_buildings_dyn36 <- readGDX(gdxPath, "ppfen_buildings_dyn36")$ppfen_buildings_dyn36$records[[1]]

  # Create complete grid of all combinations
  completeGrid <- expand_grid(
    period = as.numeric(as.character(ttot)),
    region = as.character(all_regi),
    REMINDitems_out = as.character(ppfen_buildings_dyn36)
  )

  # Merge with calculated corrections and fill missing with 1.0
  completeFactors <- completeGrid %>%
    left_join(correctionFactors, by = c("period", "region", "REMINDitems_out"))

  nCalculated <- sum(!is.na(completeFactors$value))
  nFilled <- nrow(completeFactors) - nCalculated

  completeFactors <- completeFactors %>%
    mutate(value = replace_na(.data$value, 1.0))

  message(sprintf("Completed climate correction: %d total entries (%d calculated, %d filled with 1.0)",
                  nrow(completeFactors),
                  nCalculated,
                  nFilled))

  return(completeFactors)
}


#' Write Correction Factors to REMIND GDX Format
#'
#' This function transforms correction factors into REMIND-compatible GDX format and
#' writes them to a file. The output contains the parameter \code{pm_climateCorrection}
#' with dimensions \code{ttot}, \code{all_regi}, and \code{ppfen_buildings_dyn36}.
#' If GDX writing fails, a CSV fallback file is created.
#'
#' @param correctionFactors A data frame with correction factors at REMIND resolution.
#'
#' @param outputPath A string specifying the path for the output GDX file.
#' Defaults to \code{"pmClimateCorrection.gdx"}.
#'
#' @param columnMap A named vector mapping REMIND column names to target names.
#'
#' @returns Invisibly returns \code{NULL}. The function is called for its side effect
#' of writing a GDX file (or CSV fallback) to disk.
#'
#' @importFrom dplyr %>% rename mutate select
#' @importFrom gamstransfer Container Parameter
#' @importFrom utils write.csv

.writeToREMIND <- function(correctionFactors, outputPath = "pm_ClimateCorrection.gdx", columnMap) {

  # Transform to GDX-compatible format for REMIND parameter pm_climateCorrection(ttot,all_regi,ppfen_buildings_dyn36)
  existingColumnMap <- columnMap[unname(unlist(columnMap)) %in% names(correctionFactors)]
  gdxCorrectionFactors <- correctionFactors %>%
    rename(!!!existingColumnMap) %>%
    mutate("ttot" = as.character(.data$ttot),
           "all_regi" = as.character(.data$all_regi),
           "ppfen_buildings_dyn36" = as.character(.data$ppfen_buildings_dyn36)) %>%
    select("ttot", "all_regi", "ppfen_buildings_dyn36", "value")

  # Write to GDX file
  tryCatch({
    # Using gamstransfer for GDX output
    container <- Container$new()
    param <- Parameter$new( # nolint
      container,
      name = "pm_climateCorrection",
      domain = c("ttot", "all_regi", "ppfen_buildings_dyn36"),
      records = gdxCorrectionFactors,
      description = "Climate correction factors for building energy demand"
    )
    container$write(outputPath)
    message(sprintf("Climate correction factors written to GDX: %s", outputPath))
  }, error = function(e) {
    warning(sprintf("Error writing GDX file: %s", e$message))
    # Write to CSV as fallback
    csvPath <- gsub("\\.gdx$", ".csv", outputPath)
    write.csv(gdxCorrectionFactors, csvPath, row.names = FALSE)
    message(sprintf("Correction factors written to CSV fallback: %s", csvPath))
  })
}
