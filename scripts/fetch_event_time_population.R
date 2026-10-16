suppressPackageStartupMessages(library(tidycensus))
dir.create("output/event_time", recursive = TRUE, showWarnings = FALSE)
# Cache within this project; never write the user's Census key to disk or logs.
for (yr in 2009:2022) {
  dest <- sprintf("output/event_time/acs_%s.rds", yr)
  if (file.exists(dest)) next
  z <- tryCatch(get_acs(geography = "zcta", variables = "B01003_001",
                        year = yr, survey = "acs5", geometry = FALSE,
                        cache_table = FALSE),
                error = function(e) NULL)
  if (is.null(z)) {
    message("Population unavailable for ", yr)
  } else {
    saveRDS(z, dest)
    message("Saved population for ", yr)
  }
}
stopifnot(file.exists("output/event_time/acs_2022.rds"))
