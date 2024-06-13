################################################################################
### Prepare quartzbio clinical tracker for processing by preserving the Clinical, Markers, and Response tabs and combine with the Genomics, SNIPDX tabs from the Repare tracker
################################################################################

################################################################################
### Set up
################################################################################

library(tidyverse)
library(openxlsx)
library(readxl)
library(lubridate)
library(stringi)
library(dplyr, warn.conflicts=FALSE)
library(reshape2)
library(argparser)

################################################################################
### Parse arguments
################################################################################

# === Initialize
p <- arg_parser("Process tracker data and output a study object")
p <- add_argument(p, arg = "--quartz_tracker",
                  help = "quartz tracker data path")
p <- add_argument(p, arg = "--repare_tracker",
                  help = "repare tracker data path")
p <- add_argument(p, arg = "--outdir",
                  help = "output directory")
p <- add_argument(p, arg = "--functions",
                  help = "path to functions extension folder (e.g clinicalbioinfoextensions)")
p <- add_argument(p, arg = "--trial",
                  help = "trial number (e.g. RP6306-01")
p <- add_argument(p, arg = "--date",
                  help = "date (format: yyyy-mm-dd)")

argv <- parse_args(p)

date <- as.Date(argv$date)
trial <- argv$trial
quartz_tracker <- argv$quartz_tracker
repare_tracker <- argv$repare_tracker
outdir <- argv$outdir
functions <- argv$functions


################################################################################
### Set testing arguments
################################################################################

interactive <- FALSE
if(interactive){
  quartz_tracker <- '/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-01/clinical/vendor-data/quartzbio/RP3500-01 Patient_Tracker_EDC-2022-10-06.xlsx' 
  repare_tracker <- '/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-01/clinical/vendor-data/2022-12-20_RP3500-01 Patient Tracker.xlsx'
  outdir <- '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/clinical/test/' 
  functions <-'/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoextensions/' 
  trial <- 'RP3500-01' 
  date <- as.Date('2022-12-12') 
}

################################################################################
### Read tracker, keep relevant tabs, clean up columns and data types across sheets
################################################################################
### quartz tracker
quartz_list <- readxlAllSheets(xlsxFile = quartz_tracker)

### Keep only relevant tabs
tabs_to_keep <- c("^Clinical$", "^Response$", "Markers")
tabs_to_rename <- c("Clinical", "Response", "Markers")

quartz_list <- quartz_list[grep(paste0(tabs_to_keep, collapse = "|"),
                                    names(quartz_list),
                                    ignore.case = T)]


### repare tracker
repare_list <- readxlAllSheets(xlsxFile = repare_tracker)

### Keep only relevant tabs
tabs_to_keep <- c("^Dose.*Transfusions", "^Genomics$", "^SNIPDx$")
tabs_to_rename <- c("^Dose.*Transfusions", "^Genomics$", "^SNIPDx$")

repare_list <- repare_list[grep(paste0(tabs_to_keep, collapse = "|"),
                                names(repare_list),
                                ignore.case = T)]

#remove duplicate analysis tabs, applicable to RP3500-03
repare_list <- repare_list[grep("delete|old|Do not use", names(repare_list),invert = T)]

# combine relevant tabs into single list
combined_list <- c(quartz_list, repare_list)

dir.create(outdir, showWarnings = FALSE, recursive = T)
write.xlsx(combined_list,file = paste0(outdir,date,"_", trial, "_", "Patient_Tracker_QuartzBio.xlsx"), overwrite = T)

