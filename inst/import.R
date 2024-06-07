sc <- c(
  driver = "ODBC Driver 18 for SQL Server",
  server = "inbo-sql08-prd.inbo.be",
  database = "S0008_00_Meetnetten", port = 1433, uid = "bmkreader",
  pwd = keyring::key_get("meetnetten", username = "bmkreader"), Encrypt = "no"
)
sprintf("%s=%s;", names(sc), sc) |>
  paste(collapse = "") -> constring
origin <- odbc::dbConnect(odbc::odbc(), .connection_string = constring)
