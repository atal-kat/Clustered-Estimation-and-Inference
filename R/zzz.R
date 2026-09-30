# The nobs generic is needed while the namespace's S3 table is registered,
# including in sessions started with only the base package attached.
#' @importFrom stats nobs
NULL

# Conditional registration of the tidy()/glance() methods with the
# `generics` package. `generics` sits in Suggests, never Imports: when it is
# absent the package loads and works normally.
.onLoad <- function(libname, pkgname) {
  if (requireNamespace("generics", quietly = TRUE)) {
    ns <- asNamespace(pkgname)
    for (cls in c("cjive", "cjar", "cjscore", "iv_infer")) {
      registerS3method("tidy", cls, get(paste0("tidy.", cls), envir = ns),
                       envir = asNamespace("generics"))
      registerS3method("glance", cls, get(paste0("glance.", cls), envir = ns),
                       envir = asNamespace("generics"))
    }
  }
  invisible()
}
