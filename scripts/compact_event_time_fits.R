# Compact fit files produced by an earlier version of the diagnostic script.
suppressPackageStartupMessages(library(fixest))
for (f in list.files("output/event_time",pattern="^models_[0-9]+[.]rds$",full.names=TRUE)) {
  z <- readRDS(f)
  if (is.null(z$common)) next
  fit_summary <- lapply(z$models,function(m) list(
    logit=list(coefficients=coef(m$a),fixed_effects=fixef(m$a),
               converged=m$a$convStatus,nobs=nobs(m$a)),
    poisson=list(coefficients=coef(m$b),fixed_effects=fixef(m$b),
                 converged=m$b$convStatus,nobs=nobs(m$b))))
  saveRDS(list(models=fit_summary,standardization_rows=sum(z$common)),f)
  message("Compacted ",f)
}
