# Code for setting the connection details and schemas for the various databases,
# as well as some folders in the local file system.
options(andromedaTempFolder = "d:/andromedaTemp")
databases <- list()
rootFolder <- "e:/mediatorAnalysis"

# IBM_CCAE ------------------------------------------------------------
databases[[length(databases) + 1]] <- list(
  databaseId = "CCAE",
  connectionDetails = DatabaseConnector::createConnectionDetails(
    dbms = "redshift",
    connectionString = keyring::key_get("redShiftConnectionStringOhdaCcae"),
    user = keyring::key_get("redShiftUserName"),
    password = keyring::key_get("redShiftPassword")
  ),
  cohortDatabaseSchema = "scratch_mschuemi",
  cdmDatabaseSchema = "cdm_truven_ccae_v2756"
)

# OPTUM_Extended_DOD --------------------------------------------------
databases[[length(databases) + 1]] <- list(
  databaseId = "OptumDod",
  connectionDetails = DatabaseConnector::createConnectionDetails(
    dbms = "redshift",
    connectionString = keyring::key_get("redShiftConnectionStringOhdaOptumDod"),
    user = keyring::key_get("redShiftUserName"),
    password = keyring::key_get("redShiftPassword")
  ),
  cohortDatabaseSchema = "scratch_mschuemi",
  cdmDatabaseSchema = "cdm_optum_extended_dod_v2753"
)

# Optum_EHR -----------------------------------------------------------
databases[[length(databases) + 1]] <- list(
  databaseId = "OptumEhr",
  connectionDetails = DatabaseConnector::createConnectionDetails(
    dbms = "redshift",
    connectionString = keyring::key_get("redShiftConnectionStringOhdaOptumEhr"),
    user = keyring::key_get("temp_user"),
    password = keyring::key_get("temp_password")
  ),
  cohortDatabaseSchema = "scratch_mschuemi",
  cdmDatabaseSchema = "cdm_optum_ehr_v2779"
)

# PharMetrics ---------------------------------------------------------
databases[[length(databases) + 1]] <- list(
  databaseId = "PharMetrics",
  connectionDetails = DatabaseConnector::createConnectionDetails(
    dbms = "redshift",
    connectionString = keyring::key_get("redShiftConnectionStringOhdaPharmetrics"),
    user = keyring::key_get("redShiftUserName"),
    password = keyring::key_get("redShiftPassword")
  ),
  cohortDatabaseSchema = "scratch_mschuemi",
  cdmDatabaseSchema = "cdm_iqvia_pharmetrics_plus_v2737"
)

# IBM_MDCR ------------------------------------------------------------
databases[[length(databases) + 1]] <- list(
  databaseId = "MDCR",
  connectionDetails = DatabaseConnector::createConnectionDetails(
    dbms = "redshift",
    connectionString = keyring::key_get("redShiftConnectionStringOhdaMdcr"),
    user = keyring::key_get("redShiftUserName"),
    password = keyring::key_get("redShiftPassword")
  ),
  cohortDatabaseSchema = "scratch_mschuemi",
  cdmDatabaseSchema = "cdm_truven_mdcr_v2755"
)

# Set cohort table and folders -------------------------------------------------
for (i in seq_along(databases)) {
  databases[[i]]$cohortTable <- sprintf(
    "mediation_cohort_%s", 
    databases[[i]]$databaseId
  )
  databases[[i]]$outputFolder <- file.path(
    rootFolder,
    databases[[i]]$databaseId
  )
}