# DSMHelper - Distance Sampling and Density Surface Modelling Utilities
# Combined from functions.R and ds_utils_0.6.R (Atlantic DSM project)
#
# TODO - functions that reference project-specific globals and need updating
# before they will work outside the Atlantic DSM project:
#   create.survey.data()        - references segProj, study.area
#   create.segdata()            - references segProj, study.area, the.data
#   do.generic.render()         - references ResultsDir, RDir
#   do_oneoff_render()          - sources 'analysis settings.r' directly
#   do.full.prelim()            - references spec.grps and other globals
#   create_subproject_folders() - references GenericRDataDir, GenericShapeDir,
#                                 predLayerDir, GISDir, SubProject
#   get.seas.mean()             - references predgrid
#   various env covar fns       - reference segProj, study.area, RasterDir
# ============================================================================

#' Operators imported for use throughout DSMHelper
#'
#' The pipe operators are used pervasively in this package (\code{\%>\%} alone
#' appears in over 200 places).  They cannot be written with a \code{::}
#' prefix in a pipeline, so unlike the rest of the package's dependencies they
#' have to be imported into the namespace rather than qualified at the call
#' site.
#'
#' Without this, every pipeline relies on the \emph{user} having attached
#' magrittr (or the tidyverse) before calling a DSMHelper function, because an
#' unqualified name inside package code resolves via
#' namespace -> imports -> base -> search path, and only the search path would
#' have it.
#'
#' @importFrom magrittr %>% %<>% %T>%
#' @name DSMHelper-imports
#' @keywords internal
NULL


# Print message when user executes "library(DSMHelper).
# Shamelessly borrowed from mgcv.
.onAttach <- function(...) {
  library(help=DSMHelper)$info[[1]] -> version

  if (!is.null(version)) {
    version <- version[pmatch("Version",version)]
    um <- strsplit(version," ")[[1]]
    version <- um[nchar(um)>0][2]
  } else {
    version <- "Unknown version"
  }

  hello <- paste0("This is DSMHelper ", version, ".")
  packageStartupMessage(hello)
}
