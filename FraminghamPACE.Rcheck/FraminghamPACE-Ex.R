pkgname <- "FraminghamPACE"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('FraminghamPACE')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("getFraminghamPACE")
### * getFraminghamPACE

flush(stderr()); flush(stdout())

### Name: getFraminghamPACE
### Title: Compute Framingham Pace of Aging Calculated from the Epigenome
###   (FraminghamPACE) scores
### Aliases: getFraminghamPACE

### ** Examples

## Not run: 
##D # Basic usage
##D scores <- getFraminghamPACE(betas)
##D 
##D # With a different model, coverage outputs, and a relaxed per-sample coverage threshold
##D res <- getFraminghamPACE(
##D   betas,
##D   model_type = "Alternate1",
##D   sample_coverage_min = 0.70,
##D   return_coverage = TRUE,
##D   return_normalized_betas = FALSE
##D )
##D res$pace_score            # data.frame of scores
##D head(res$coverage$sample) # per-sample predictive coverage (%)
##D res$coverage$threshold    # threshold used
## End(Not run)




cleanEx()
nameEx("getFraminghamPACEProbes")
### * getFraminghamPACEProbes

flush(stderr()); flush(stdout())

### Name: getFraminghamPACEProbes
### Title: List probes used by FraminghamPACE models
### Aliases: getFraminghamPACEProbes

### ** Examples

## Not run: 
##D probes <- getFraminghamPACEProbes(backgroundList = FALSE, model_type = "FraminghamPACE")
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
