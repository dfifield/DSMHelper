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
