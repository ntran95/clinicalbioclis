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
p <- add_argument(p, arg = "--tracker",
                  help = "tracker data path")
p <- add_argument(p, arg = "--outdir",
                  help = "output directory")
p <- add_argument(p, arg = "--functions",
                  help = "functions file")
p <- add_argument(p, arg = "--trial",
                  help = "trial number (e.g. RP6306-01")
p <- add_argument(p, arg = "--date",
                  help = "date (format: yyyy-mm-dd)")
p <- add_argument(p, arg = "--markers",
                  help = "logical to indicate if marker data is available, set to T or F")
p <- add_argument(p, arg = "--response",
                  help = "logical to indicate if response data is available, set to T or F")
p <- add_argument(p, arg = "--snipdx",
                  help = "logical to indicate if SNIPDx data is available, set to T or F")
p <- add_argument(p, arg = "--key",
                  help = "specify the path to the master key labelling each program, defaulted to FALSE if variable is not supplied")
p <- add_argument(p, arg = "--keyAvail",
                  help = "logical to indicate whether path to program key is should be read in, defaulted to FALSE if variable is not supplied",default = FALSE)
p <- add_argument(p, arg = "--inctdna",
                  help = "path to ctDNA vaf dataframe, read in best mVAFR")
p <- add_argument(p, arg = "--genomicstats",
                  help = "logical, indicate whether monthly report of Genomics and best ctDNA response be exported to Stats?")


argv <- parse_args(p)

date <- as.Date(argv$date)
trial <- argv$trial
tracker <- argv$tracker
outdir <- argv$outdir
functions <- argv$functions
markers <- as.logical(argv$markers)
response <- as.logical(argv$response)
snipdx <- as.logical(argv$snipdx)
keyAvail <- as.logical(argv$keyAvail)
inProgramKey <- argv$key
inctdna <- argv$inctdna
genomicstats <- argv$genomicstats



#source(functions)
devtools::load_all(functions)


################################################################################
### Set testing arguments
################################################################################

# date <- ymd("2022-02-28")
# trial <- "RP6306-01"
# outdir <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-01/clinical/"
# tracker <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-01/clinical/20220228_RP6306-01 Patient Tracker.xlsx"
# markers <- T
# response <- T
# snipdx <- T

################################################################################
### Set testing arguments
################################################################################

interactive <- FALSE
if(interactive){
  # tracker <- '/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/clinical/vendor-data/2022-09-29_RP-3500-03 ATTACC Patient Status Tracker_Working Copy.xlsx' 
  tracker <- '/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/clinical/vendor-data/2022-10-29_RP-3500-03 ATTACC Patient Status Tracker_Working Copy.xlsx' 
  
  outdir <- '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/' 
  functions <-'/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/src/functions_ctDNA.R' 
  trial <- 'RP3500-03' 
  date <- as.Date('2022-10-31') 
  markers <- T
  response <- T 
  snipdx <- T
  keyAvail <- T
  inProgramKey <- "/ClinBio/SP-ClinicalBioinformatics/shared/Program-Key/master_key_across_programs.csv"
  genomicstats <- T
  inctdna <- "/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-01/ctDNA/processing/Combined/results/2022-08-11_RP3500-01_processed_ctDNA_Combined_results/2022-08-11_RP3500-01_processed_ctDNA_Combined_Summary.xlsx"
  
  #source(functions)
  devtools::load_all(functions)
  
  
}


################################################################################
### Read tracker, keep relevant tabs, clean up columns and data types across sheets
################################################################################

clinical_list <- readxlAllSheets(xlsxFile = tracker)

### Keep only relevant tabs
tabs_to_keep <- c("^Clinical$", "^Response$", "Markers","^Dose.*Transfusions", "^Genomics$", "^SNIPDx$")
tabs_to_rename <- c("Clinical", "Response", "Markers", "Dose","Genomics","SNIPDx")

clinical_list <- clinical_list[grep(paste0(tabs_to_keep, collapse = "|"),
                                    names(clinical_list),
                                    ignore.case = T)]

#remove duplicate analysis tabs, applicable to RP3500-03
clinical_list <- clinical_list[grep("delete|old|Do not use", 
                                    names(clinical_list),invert = T)]

# match sheets names in to_keep with to_rename, regardless of position of sheet names in to_rename
rename<- sapply(tabs_to_rename, function(tabs) names(clinical_list)[grep(tabs, tabs_to_keep,ignore.case = T)])
clinical_list <- clinical_list[match(rename, names(clinical_list))] 
names(clinical_list) <- names(rename)
                
### Clean columns names, drop rows with all NAs, convert column to correct data type
#clinical_list <- lapply(clinical_list, function(x){x <- x %>% filter(!across(everything(), is.na))}) ##across() is deprecated, use if_any() instead
clinical_list <- lapply(clinical_list, function(x){x <- x %>% filter(if_any(everything(), ~ !is.na(.x)))})
clinical_list <- lapply(clinical_list, function(x){ x <- readr::type_convert(x) })

# iterate and modify each excel sheet
clinical_list <- lapply(clinical_list, function(x){
  # remove columns containing "ignore"
  x <- x[,grep("Ignore", colnames(x), invert = T)]
  #remove empty columns created from readxl
  x <-x[,grep('[...][0-9]', colnames(x),invert = T)]
  
  # remove row entires if Patient ID is NA
  x <- x[!is.na(x$`Patient ID`), ]
})

################################################################################
### Read in program key, switch column names from current trial to the master key that
### meets the mimimal requirement to run through pipeline
################################################################################
#load key
if (keyAvail){
  print(class(inProgramKey))
  programKeyTbl  <- read.csv(inProgramKey,
                             header = T,
                             sep = ",",
                             check.names = F,
                             na.strings = "",
                             fileEncoding="UTF-8-BOM")
  
  print(programKeyTbl)
  
  clinical_list <- lapply(clinical_list, function(x){ 
    x <- programToMasterKey(trial = trial,
                            programKeyTbl = programKeyTbl, 
                            data = x)})
  
}else{
  print("No program conversion file inputed")
}
 
Clinical <- clinical_list$Clinical %>% select(any_of(clinical_cols))
Genomics <- clinical_list$Genomics %>% select(any_of(genomics_cols))
if(markers) { Markers <-  clinical_list$Markers %>% select(any_of(markers_cols)) }
if(response) { Response <- clinical_list$Response %>% select(any_of(response_cols)) }
if(snipdx) { SNIPDx <- clinical_list$SNIPDx %>% select(any_of(snipdx_cols)) }


################################################################################
### Clean clinical sheet
################################################################################

### Set dates to Date class
Clinical$`Tx Discont.` <- as.Date(as.character(Clinical$`Tx Discont.`))
Clinical$`First dose` <- as.Date(as.character(Clinical$`First dose`))

### Set time-stamps to date of data copy if necessary

Clinical$`Off-Tx Flag`[ Clinical$`Tx Discont.` > date ] <- NA

# Modify Days, weeks, cycles in accordance with the data extraction date
Clinical$Days <- as.numeric(Clinical$`Tx Discont.`- Clinical$`First dose`)
Clinical$Weeks <- round(Clinical$Days / 7, 1)
Clinical$Cycles <-  round(Clinical$Weeks / 3, 1)

################################################################################
### Clean response sheet, calculate BOR
################################################################################

if(response) {

# ensure only RECIST levels are converted, ignore everything after levels
response_lvls <- c("CR","PR", "SD",  "PD")
regex <- paste0("(", paste(response_lvls,collapse = "|"), ").*")
pattern <- paste(regex, collapse = "|")
Response$`Response Type` <- gsub(regex, "\\1",Response$`Response Type`)
Response$responseTypeLong <- ResponseShort2Long(Response$`Response Type`)

### Calculate time in days from baseline

Response$`% Change from Baseline`[ Response$Assessment == "Baseline"] <- 0
# if cases of duplicate Assessments, choose , the latest date 
Response <- Response %>% group_by(`Patient ID`, Assessment) %>% 
  slice_max(`Response Date`) %>%
  ungroup()
# ensure NAs are removed from assessment levels
assessments <- unique(Response$Assessment[ Response$Assessment != "Baseline" & !is.na(Response$Assessment) ])
timepoints <- Response %>% select(`Patient ID`, `Assessment`, `Response Date`) %>% filter(!is.na(Assessment))
timepoints <- pivot_wider(timepoints, id_cols = `Patient ID`, names_from =Assessment, values_from = `Response Date`)
timepoints <- pivot_longer(timepoints, cols = all_of(assessments), names_to = "Assessment", values_to = "Response Date")
timepoints <- timepoints %>% mutate(`Response Days` = as.numeric(`Response Date` - Baseline))
Response <- merge(Response, select(timepoints, `Patient ID`, Assessment, `Response Days`), all.x = T)

### Calculate BOR

possible_responses <- c("CR", "PR", "SD", "PD")
BestResponse <- Response %>% distinct(`Patient ID`,`Response Type`,.keep_all = T) %>% 
  pivot_wider(., id_cols = `Patient ID`, names_from = `Response Type`, values_from = `% Change from Baseline`, values_fill = NA)
missing <- setdiff(possible_responses, colnames(BestResponse))
BestResponse[missing] <- NA
BestResponse$BOR <- ifelse( !is.na(BestResponse$CR), "CR",
                    ifelse( !is.na(BestResponse$PR), "PR",
                    ifelse( !is.na(BestResponse$SD), "SD",
                    ifelse( !is.na(BestResponse$PD), "PD", NA))))

### Calculate best percent change in tumor volume from baseline

BestPercentChange <- clinical_list$Response %>% select(`Patient ID`, `% Change from Baseline`) %>%  group_by(`Patient ID`) %>% summarise_all(min, na.rm=T)
colnames(BestPercentChange)[2] <- "Best % Change From Baseline"
BestResponse <- merge(BestPercentChange,BestResponse)
BestResponse <- select(BestResponse, `Patient ID`, BOR, `Best % Change From Baseline`)

# if `Best % Change From Baseline is infinite`, convert to NAs
BestResponse[is.infinite(BestResponse$`Best % Change From Baseline`),
             "Best % Change From Baseline"] <- NA

}

################################################################################
### Clean marker sheet, calculate marker response
################################################################################

if(markers) {
  
MarkerResponse <- select(Markers, `Patient ID`, 
                         `Tumor marker`,
                         BestMarkerChange = `Tumor Marker % Change From Baseline`) %>%
  #remove NAs from Tumor Marker levels
  filter(!is.na( `Tumor marker`))

MarkerResponse <- MarkerResponse %>% 
  group_by(`Patient ID`, `Tumor marker`) %>% 
  summarise_all(min, na.rm =T) %>% 
  ungroup()

MarkerResponse <- MarkerResponse %>% mutate( BestMarkerResponse = case_when(BestMarkerChange <= 50 ~ "Response",
                                                                            BestMarkerChange > 50 ~ "Non-Response",
                                                                            is.na(BestMarkerChange) ~ "Unavailable"))
#pivot wider by Tumor marker, ensure unique rowwise patient id before merging with clinical
MarkerResponseWide <- MarkerResponse %>%
  pivot_wider(id_cols = c(`Patient ID`), names_from = `Tumor marker`, values_from = BestMarkerResponse, values_fill = NA)

#derive COR from Marker first, then tumor after
MarkerResponseWide <- MarkerResponseWide %>% 
  rowwise() %>%
  mutate(COR = ifelse(any(grepl('^Response$', c_across())), 'PR', NA)) 
}


################################################################################
### Merge all data into one table
################################################################################

Combined <- Clinical
Combined <- merge(Combined, Genomics, all.x = T)
if(snipdx) {Combined <- merge(Combined, SNIPDx, all.x = T)}
if(response) { Combined <- merge(Combined, BestResponse, all.x = T) }
if(markers) { Combined <-  merge(Combined, MarkerResponseWide, all.x = T) }

################################################################################
### Calculate combined overall response
################################################################################

#### Derive Combined Response

if(markers & response) {
  
  #tumor_markers <-  names(MarkerResponseWide)[names(MarkerResponseWide) != "Patient ID"]
  
  Combined <- Combined %>% mutate(COR = case_when(BOR == "PR"~"PR",
                                                  BOR == "SD"~"SD",
                                                  BOR == "PD"~"PD"))
  
  
} else if (response) {
  
  Combined <- Combined %>% mutate(COR = BOR)
  
}

################################################################################
#### Genomic Stats Output ####
################################################################################
#### Export monthly report derived from Genomics tab and ctDNA best reponses
if(genomicstats){
  ctdna_dfs <- readxlAllSheets(inctdna)
  
  dfmVAFRBest <- ctdna_dfs$dfmVAFRBest
  
  dfmVAFRBest <- dfmVAFRBest %>% 
    mutate(Molecular.Response = case_when(`Molecular Response` == "Response" ~ "MR",`Molecular Response` == "Stable" ~ "Non-MR",`Molecular Response` == "Progression" ~ "Non-MR",TRUE~"")) %>%
    select(Patient_ID,`Molecular.Response`,mVAFR) %>% 
    dplyr::rename(c( `Patient ID` = Patient_ID ,`Molecular Response` = Molecular.Response, `Best mVAFR` = mVAFR))
  
  # preserve all patient ids present in Genomics not present in Clinical
  stats_out <- Combined %>%
    full_join(Genomics) %>%
    select(`Patient ID`,`Enrollment Gene`,TMB,MSI,Reversion,`Enrollment Test Type`,`germline status`,`Type of loss`,`ATM IHC`)%>% 
    left_join(dfmVAFRBest,by="Patient ID")
  
  
}

################################################################################
#### Plot Helper Variables ####
################################################################################
#### Generate Short and Long IDs
Combined$`Patient No` <- gsub("[0-9]+-0","", Combined$`Patient ID`)

Combined$ID.disease <- paste(Combined$`Tissue Type`, Combined$`Patient ID`)

Combined$ID.gene <- paste(Combined$`Enrollment Gene`, Combined$`Patient ID`)

Combined$ID.disease.short <- paste(Combined$`Tissue Type`, Combined$`Patient No`)

Combined$ID.gene.short <- paste(Combined$`Enrollment Gene`, Combined$`Patient No`)

################################################################################
### Save individual clean tables to clinical object
################################################################################

study <- new("Study", clinical = Clinical, genomics = Genomics, combined = Combined,dose = clinical_list$Dose, trial = trial, date = date)
if(response) { study@response = Response }
if(markers) { study@marker = Markers }
if(snipdx) { study@snipdx = SNIPDx}

dir.create(outdir, showWarnings = FALSE, recursive = T)
saveRDS(study, paste0(outdir, "/", date, "_", trial, "_processedStudy.rds"))
write.xlsx(stats_out,paste(outdir,"/",date,"_", trial, "_Genomics_Summary.xlsx",sep=""))

################################################################################

