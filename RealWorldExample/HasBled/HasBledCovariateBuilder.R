# Implements HAS-BLED score as a custom covariate builder for FeatureExtraction

library(dplyr)

getHasBledCovariateData <- function(connection,
                                    oracleTempSchema = NULL,
                                    cdmDatabaseSchema,
                                    cdmVersion = "5",
                                    cohortTable = "#cohort_person",
                                    rowIdField = "row_id",
                                    aggregated,
                                    cohortIds,
                                    covariateSettings,
                                    minCharacterizationMean) {
  message("Constructing HAS-BLED covariates")
  if (covariateSettings$useHasBled == FALSE) {
    return(NULL)
  }
  if (aggregated) {
    stop("Aggregation not supported")
  }
  
  # Construct covariates:
  hasBledCohorts <- readr::read_csv("RealWorldExample/HasBled/HasBledCohorts.csv", show_col_types = FALSE)
  sql <- "
  SELECT row_id,
    COUNT(*) AS cohort_count
  FROM (
    SELECT DISTINCT cohort.@row_id_field AS row_id,
      has_bled_cohort.cohort_definition_id
    FROM @cohort_table cohort
    INNER JOIN @has_bled_cohort_database_schema.@has_bled_cohort_table has_bled_cohort
      ON cohort.subject_id = has_bled_cohort.subject_id
        AND cohort.cohort_start_date >= has_bled_cohort.cohort_start_date
    WHERE has_bled_cohort.cohort_definition_id IN (@has_bled_cohort_ids)  
  {@cohort_ids != -1} ? {
      AND cohort.cohort_definition_id IN (@cohort_ids)
  }
  ) tmp
  GROUP BY row_id;
  "
  cohortBasedScore <- DatabaseConnector::renderTranslateQuerySql(
    connection = connection,
    sql = sql,
    cohort_table = cohortTable,
    row_id_field = rowIdField,
    cohort_ids = cohortIds,
    has_bled_cohort_database_schema = covariateSettings$hasBledCohortDatabaseSchema,
    has_bled_cohort_table = covariateSettings$hasBledCohortTable,
    has_bled_cohort_ids = hasBledCohorts$cohortId,
    snakeCaseToCamelCase = TRUE
  )
  sql <- "
  SELECT @row_id_field AS row_id
  FROM @cohort_table cohort
  INNER JOIN @cdm_database_schema.person
    ON cohort.subject_id = person.person_id
  WHERE YEAR(cohort_start_date) - year_of_birth > 65
  {@cohort_ids != -1} ? {
      AND cohort.cohort_definition_id IN (@cohort_ids)
  }
  "
  ageBasedScore <- DatabaseConnector::renderTranslateQuerySql(
    connection = connection,
    sql = sql,
    cohort_table = cohortTable,
    row_id_field = rowIdField,
    cohort_ids = cohortIds,
    cdm_database_schema = cdmDatabaseSchema,
    snakeCaseToCamelCase = TRUE
  )
  covariates <- cohortBasedScore %>%
    select("rowId", "cohortCount") %>%
    full_join(ageBasedScore %>%
                mutate(ageScore = 1), 
              by = join_by(rowId)) %>%
    mutate(covariateValue = if_else(is.na(cohortCount), 0, cohortCount) + if_else(is.na(ageScore), 0, ageScore),
           covariateId = 999) %>%
    select("rowId", "covariateId", "covariateValue")
  
  # Construct covariate reference:
  covariateRef <- data.frame(
    covariateId = 999,
    covariateName = "HAS-BLED",
    analysisId = 999,
    conceptId = 0
  )
  
  # Construct analysis reference:
  analysisRef <- data.frame(
    analysisId = 999,
    analysisName = "HAS-BLED",
    domainId = "Condition",
    startDay = -9999,
    endDay = 0,
    isBinary = "N",
    missingMeansZero = "Y"
  )
  
  # Construct analysis reference:
  metaData <- list(call = match.call())
  result <- Andromeda::andromeda(
    covariates = covariates,
    covariateRef = covariateRef,
    analysisRef = analysisRef
  )
  attr(result, "metaData") <- metaData
  class(result) <- "CovariateData"
  return(result)
}

createHasBledCovariateSettings <- function(useHasBled = TRUE,
                                           hasBledCohortDatabaseSchema,
                                           hasBledCohortTable) {
  covariateSettings <- list(useHasBled = useHasBled,
                            hasBledCohortDatabaseSchema = hasBledCohortDatabaseSchema,
                            hasBledCohortTable = hasBledCohortTable)
  attr(covariateSettings, "fun") <- "getHasBledCovariateData"
  class(covariateSettings) <- "covariateSettings"
  return(covariateSettings)
}