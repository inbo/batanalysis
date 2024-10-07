library(batanalysis)
sc <- c(
  driver = "ODBC Driver 18 for SQL Server", server = "inbo-sql08-prd.inbo.be",
  database = "S0008_00_Meetnetten", port = 1433, uid = "bmkreader",
  pwd = keyring::key_get("meetnetten", username = "bmkreader"), Encrypt = "no"
)
sprintf("%s=%s;", names(sc), sc) |>
  paste(collapse = "") -> constring
origin <- odbc::dbConnect(odbc::odbc(), .connection_string = constring)
target <- "../batanalysis_data"
if (file_test("-d", target)) {
  target <- git2rdata::repository(target)
} else {
  dir.create(target)
  target <- git2r::init(target)
}
import_raw_data(origin = origin, target = target)
DBI::dbDisconnect(origin)

analysis_data <- "../batanalysis_analysis"
if (file_test("-d", analysis_data)) {
  analysis_data <- git2rdata::repository(analysis_data)
} else {
  dir.create(target)
  analysis_data <- git2r::init(analysis_data)
}
prepare_analysis_data(
  raw_data = target, analysis_data = analysis_data, start = 2007,
  n_winter = 4, n_present = 3, strict = FALSE
)

