################################################################################
####  Initilize environment ####
################################################################################

library(tidyverse)
library(tidyr)
library(openxlsx)
library(readxl)
library(stringi)
library(data.table)
library(dplyr, warn.conflicts=FALSE)   # Manipulating the dataframes
library(reshape2) # Reformmating dataframes a
library(argparser)

print(paste0("The current working directory is:"))
getwd()

## === Initialize argparse =====
p <- arg_parser("Process clinical data and output into a list of dataframes")
p <- add_argument(p, arg = "--file", 
                  help = "specify the input Patient tracker data, use relative or absolute paths")
#omit, use regex
#p <- add_argument(p, arg = "--date", 
#                  help = "specify the date that data is process, use the date on Patient Tracker")
p <- add_argument(p, arg = "--trial",
                  help = "specify the trial, i.e) RP3500-01")
p <- add_argument(p, arg = "--output_dir", 
                  help = "specify the path of the output directory, can be relative or absolute paths. Include '/'")
p <- add_argument(p, arg = "--git_extension", 
                  help = "specify the path of the gitextension, can be relative or absolute paths. Include '/'")
p <- add_argument(p, arg = "--markers",
                  help = "logical to indicate if marker data is available, set to T or F")
p <- add_argument(p, arg = "--response",
                  help = "logical to indicate if response data is available, set to T or F")
p <- add_argument(p, arg = "--key",
                  help = "specify the path to the master key labelling each program")
p <- add_argument(p, arg = "--ctdna", 
                  help = "specify the ctDNA file")

argv <- parse_args(p)

#date <- argv$date

#trial <- argv$trial

git_extension <- argv$git_extension
output_dir <- argv$output_dir
inClinical <- argv$file
markers <- as.logical(argv$markers)
response <- as.logical(argv$response)
inProgramKey <- argv$key
trial <- argv$trial
inctDNA <- argv$ctdna


correctClinicalTimestamp2 <- function(data, date){
  tx_dis <- data$`Tx Discont.`
  first_dose <- data$`First dose`
  
  if(any(tx_dis > as.POSIXct(date)  & 
         !is.na(tx_dis))) {
    
    data[which(tx_dis > as.POSIXct(date)),]$`Tx Discont.` <-
      as.POSIXlt(paste(date, "00:00:00"),"UTC")
    
    
    data[which(tx_dis > date),"Off-Tx Flag"] <- NA
    
  }
    #Modify Days, weeks, cycles in accordance with the data extraction date
    data$Days <- as.numeric(tx_dis - first_dose)
    data$Weeks = as.numeric(data$Days / 7 )
    data$Cycles = as.numeric(data$Weeks / 3)
  
  
  return(data)
}

################################################################################
#### Interactive Session ####
################################################################################
interactive <- F
if(interactive){
  git_extension <- here::here("/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoextensions/")
  
  inClinical <- here::here("/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/clinical/vendor-data/2022-11-11_RP-3500-03 ATTACC Patient Status Tracker_Working Copy.xlsx")
  
  inProgramKey <- here::here("/ClinBio/SP-ClinicalBioinformatics/shared/Program-Key/master_key_across_programs.csv")
  
  markers <- F
  
  response <- T
  
  output_dir <- "/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/clinical/processing/data/"
  
  inctDNA <- "/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/ctDNA/processing/Combined/results/2022-10-26_RP3500-03_processed_ctDNA_Combined_Batch7_Revised_results/2022-10-26_RP3500-03_processed_ctDNA_Combined_Batch7_Revised_Summary.xlsx"
  
  trial <- "RP3500-03"
}

################################################################################
### Load dataset
################################################################################
date <- str_extract(string = inClinical, "\\d{4}-\\d{2}-\\d{2}")

# trial <- str_extract(string = inClinical, "RP\\d{4}-\\d{2}")

# at some point we have to install the extension as a R system library using devtools::install()
devtools::load_all(git_extension)

#load key
programKeyTbl  <- read.csv(inProgramKey,
                          header = T,
                          sep = ",",
                          check.names = F,
                          na.strings = "",
                          fileEncoding="UTF-8-BOM")

#print msg loading in 
message(paste0("Loading in Patient Tracker file from: ", inClinical))

message(paste("The extraction date of the Patient Tracker file is:", date, 
              "from the following trial:", trial))

# from the patient tracker spread sheet, load in every sheet as a df in a named list
# convert "NA" and "" into NA
  suppressMessages(
  clinical_list <- readxlAllSheets(xlsxFile = inClinical,
                                               type.convert = FALSE)
)

################################################################################
#### Clean each tabs ####
################################################################################  
  
# the names of each df in the clinical_list will inherit the sheet names from the excel sheet
print("Extracting the following excel sheets:")
names(clinical_list)

# rename list elements to match isilverman's variables 
to_rename <- c("Clinical","Response", "Markers", 
               "DoseTransfusions",  "Genomic", "SNiPDx")

#ensure tab selection is compatible with rp6306
to_keep <- c("Clinical", "Response", "markers", 
             "Dose", "Genomic", "snipdx")

names(to_keep) <- to_rename

clinical_list <- clinical_list[grep(paste0(to_keep, collapse = "|"),
                          names(clinical_list),
                          ignore.case = T)]

if(trial == "RP3500-03"){
  #remove duplicate analysis tabs
  clinical_list <- clinical_list[grep("delete|old|Do not use|Sample Level", 
                                      names(clinical_list),invert = T)]
}

# match sheets names in to_keep with to_rename, regardless of position of sheet names in to_rename
names(clinical_list) <- unlist(lapply(to_keep, 
                                      grep,
                                      to_rename, 
                                      ignore.case = TRUE,
                                      value = T))


# iterate and modify each excel sheet
clinical_list <- lapply(clinical_list, function(x){
  # remove columns containing "ignore"
  x <- x[,grep("Ignore", colnames(x), invert = T)]
  #remove empty columns created from readxl
  x <-x[,grep('[...][0-9]', colnames(x),invert = T)]
  
  # remove row entires if Patient ID is NA
  x <- x[!is.na(x$`Patient ID`), ]
  
  # convert trial columns to master key
  x <- programToMasterKey(trial = trial,
                          programKeyTbl = programKeyTbl, 
                          data = x)
})


### Add Response in Long Form
clinical_list$Response$responseTypeLong <- 
  ResponseShort2Long(clinical_list$Response$`Response Type`)

### Convert column class
# validate that all numeric cols have been converted
num_cols <- c("Dose",
          "Daily Dose",
          "% Change from Baseline",
          "Tumor Marker % Change From Baseline",
          "Tumor Marker Weeks",
          "Tumor Marker Value")

clinical_list <- lapply(clinical_list, function(df){
  mutate(df, across(any_of(num_cols), as.numeric))
})

################################################################################
#### Correct Timestamp ####
################################################################################

### Change Timestamps,
# if Patient has an Off-Tx flag but empty Tx.Discont, replace with EOT Visit
Off.Tx.Flag <- clinical_list$Clinical$`Off-Tx Flag`
Tx.Discont <- clinical_list$Clinical$`Tx Discont.`

clinical_list$Clinical[which(Off.Tx.Flag %in% "Off-Tx" &
                               is.na(Tx.Discont)), "Tx Discont."] <- 
  clinical_list$Clinical[which(Off.Tx.Flag %in% "Off-Tx" &
                                 is.na(Tx.Discont)), "EOT visit"] 


# double check pat: 1008-0006, 1008-0015, 1008-0018, 1015-0030 are really off-tx
# labelled as off-tx but no tx discount. date
clinical_list$Clinical <- correctClinicalTimestamp2(data = clinical_list$Clinical,
                          date = date)

################################################################################
#### Merge Genomics, SNiPDx to Clinical ####
################################################################################

# As of 10-22-2021, merge SNiPDx with Genomics prior to merging Genomics with Clinical
if(any(grepl("SNIPDx", names(clinical_list),ignore.case = T))){
  print("Merging Genomics data with SNIPDX data")
  
  SNiPDx <- clinical_list[[grep("SNIPDx",
                               names(clinical_list),
                               ignore.case = T)]]
  
  clinical_list$Genomic <- merge(clinical_list$Genomic,
                                 SNiPDx, all.x = T)
  
  colsGenomic <- colnames(clinical_list$Genomic)
  #sometimes Column orders shift in the Patient Tracker so we have to specify by name
  to_keep <- colsGenomic[colsGenomic %in% c("Patient ID", 
                                            "Enrollment Gene",
                                            "cDNA Change",
                                            "AA Change", 
                                            "Other Alteration",
                                            "Alteration (S)" ,
                                            "Alteration Classification", 
                                            "Enrollment VAF", 
                                            "TP53 Status", 
                                            "Other Drivers", 
                                            "TMB", 
                                            "MSI","Reversion",
                                            "Enrollment Test Provider",
                                            "Enrollment Test", 
                                            "Enrollment Test Type", 
                                            "germline status",
                                            "Germline call reason", 
                                            "Liquid biopsy Guardant Baseline",
                                            "Liquid Biopsy Tempus Baseline", 
                                            "ATM IHC",
                                            "RNASEH2 IHC",
                                            "baseline gH2AX",
                                            "Paired Biopsy (Y/N)", 
                                            "SNiPDx PBMC",
                                            "SNiPDx Tumor Mutation",
                                            "SNiPDx Gene",
                                            "SNIPDx cDNA Change",
                                            "SNiPDx AA Change", 
                                            "SNiPDx Other",
                                            "SNiPDx Alteration",
                                            "SNiPDx Alteration Type", 
                                            "SNiPDx Tumor VAF", "SNiPDx Tumor CN",
                                            "Type of loss")]

  print("Merging Clinical data with Genomics data")
  
  clinical_list$Clinical <- clinical_list$Clinical %>%
    select(-any_of("Notes")) %>%
    left_join(select(clinical_list$Genomic, any_of(to_keep)))
  
}else{
  print("Merging Clinical data with Genomics data")
  
  ## Merge Clincial Genomic and Duration of Treatment Data Tables
  clinical_list$Clinical <- merge(clinical_list$Clinical[,c(1:ncol(clinical_list$Clinical) - 1)],
                                  clinical_list$Genomic,all.x=T)
}


### Set GeneGroups and Merge "other STEP2" alterations
clinical_list$Clinical <- setGeneGroup(df = clinical_list$Clinical)

################################################################################
#### Response ####
################################################################################
if(response){
suppressWarnings(
  ### Derive BOR and Best %Change
  BestResponse <- dcast(clinical_list$Response,
                        `Patient ID`~`Response Type`,
                        value.var="% Change from Baseline")
  
) #end warning msg

# PR takes precedences over SD, SD takes precedence over PD
BestResponse[BestResponse$PD > 0,"BOR"] <- "PD"
BestResponse[BestResponse$SD >0,"BOR"] <- "SD"
BestResponse[BestResponse$PR >0,"BOR"] <- "PR"
BestResponse[BestResponse$CR >0,"BOR"] <- "CR"

#convert chr "NA" to NA
BestResponse$BOR[BestResponse$BOR == "NA"] <- NA

suppressWarnings(
  BestPercentChange <- dcast(clinical_list$Response,
                             `Patient ID`~`Response Type`,
                             value.var="% Change from Baseline",
                             function(x) min(x,na.rm = T))
) #end warning msg

#BestResponse <- merge(BestPercentChange,BestResponse[,c(1,6)])

BestResponse <- BestResponse %>%
  select(`Patient ID`, BOR) %>%
inner_join(BestPercentChange)

BestResponse$`Best % Change From Baseline` <- NA
BestResponse[which(BestResponse$BOR == "PD"),"Best % Change From Baseline"] <- BestResponse[which(BestResponse$BOR == "PD"),"PD"]
BestResponse[which(BestResponse$BOR == "SD"),"Best % Change From Baseline"] <- BestResponse[which(BestResponse$BOR == "SD"),"SD"]
BestResponse[which(BestResponse$BOR == "PR"),"Best % Change From Baseline"] <- BestResponse[which(BestResponse$BOR == "PR"),"PR"]
BestResponse[which(BestResponse$BOR == "CR"),"Best % Change From Baseline"] <- BestResponse[which(BestResponse$BOR == "CR"),"CR"]

# if `Best % Change From Baseline is infinite`, convert to NAs
BestResponse[is.infinite(BestResponse$`Best % Change From Baseline`),
             "Best % Change From Baseline"] <- NA

# merge `BOR` and `Best % Change From Baseline to Clinical`
# clinical_list$Clinical <- merge(clinical_list$Clinical,
#                                 BestResponse[,c(1,6:7)],
#                                 by="Patient ID",all.x=T)

clinical_list$Clinical <- clinical_list$Clinical %>%
  left_join(select(BestResponse, `Patient ID`,BOR, `Best % Change From Baseline`))
}

################################################################################
#### Markers ####
################################################################################
### Derive PSA Response ###
# coerce string "NA" to NA in  `Marker Response` levels
#clinical_list$Markers[which(clinical_list$Markers$`Marker Response` %in% "NA"),]$`Marker Response` <- NA
if(markers){
  suppressWarnings(
    clinical_list$BestMarkerResponseTime <- dcast(clinical_list$Markers,
                                                  `Patient ID`+`Tumor marker`~`Marker Response`,
                                                  value.var="Tumor Marker Weeks",
                                                  function(x) min(x,na.rm = T))
  ) #end warning msg
  
  
  clinical_list$BestMarkerResponseTime <- 
    clinical_list$BestMarkerResponseTime %>%
    select(any_of(c("Patient ID", "Tumor marker", "Response")))
  
  clinical_list$BestMarkerResponseTime <- 
    clinical_list$BestMarkerResponseTime[clinical_list$BestMarkerResponseTime$Response != Inf,]
  
  suppressWarnings(
    BestMarkerResponse <- dcast(clinical_list$Markers,
                                `Patient ID`~`Tumor marker`,
                                value.var="Tumor Marker % Change From Baseline",
                                function(x) min(x,na.rm = T))
  ) #end warning msg
  
  suppressWarnings(
    BestMarkerChange <- dcast(clinical_list$Markers,
                              `Patient ID`~`Tumor marker`,
                              value.var="Tumor Marker % Change From Baseline",
                              function(x) min(x,na.rm = T))
  ) #end warning msg
  
  BestMarkerResponse$`NA` <- NULL
  BestMarkerChange$`NA` <- NULL
  
  BestMarkerResponse[BestMarkerChange <= 50] <- "Response"
  BestMarkerResponse[BestMarkerChange > 50] <- "Non-Response"
  BestMarkerResponse[BestMarkerChange == Inf] <- NA
  BestMarkerResponse$`Patient ID` <- BestMarkerChange$`Patient ID`
  
  print("Merging Clinical data with BestMarkerResponse data")
  
  clinical_list$Clinical <- merge(clinical_list$Clinical,
                                  BestMarkerResponse,
                                  by="Patient ID",all.x=T)
  
  #### Derive Combined Response, subjct to change to "response, stable, progression"
  # ensure response columns exists before mutating COR
  
  # CA15-3 present in RP3500-03
  #return all marker levels
  tumor.markers <- unique(na.omit(clinical_list$Markers$`Tumor marker`))
  
  tumor.markers <- tumor.markers[tumor.markers !="-"]
    
  clinical_list$Clinical <- clinical_list$Clinical %>%
  mutate(COR = case_when(
    BOR == "CR"~"CR",
    BOR == "PR"~"PR",
    # mutate if the tumor is present (diff markers each trial)
     across(any_of(tumor.markers)) %in% "Response" ~"PR",
    #across(any_of("CA-125")) == "Response" ~"PR",
    BOR == "SD"~"SD",
    BOR == "PD"~"PD"))
}
################################################################################
#### Dose ####
################################################################################

### Derive Actualy Dose
Dose <- clinical_list$DoseTransfusions[clinical_list$DoseTransfusions$`Event Type` == "Dose",]
#was this suppose to be clinical_list$Clinical <- merge(clinical_list$Clinical, Dose) instead?
#Dose df is not exported to .Rda
Dose <- merge(clinical_list$Clinical,Dose)

################################################################################
#### Plot Helper Variables ####
################################################################################
#### Generate Short and Long IDs
clinical_list$Clinical$`Patient No` <- gsub("[0-9]+-0","",
                                            clinical_list$Clinical$`Patient ID`)

clinical_list$Clinical$ID.disease <- paste(clinical_list$Clinical$`Tissue Type`,
                                           clinical_list$Clinical$`Patient ID`)

clinical_list$Clinical$ID.gene <- paste(clinical_list$Clinical$`Enrollment Gene`,
                                        clinical_list$Clinical$`Patient ID`)

clinical_list$Clinical$ID.disease.short <- paste(clinical_list$Clinical$`Tissue Type`,
                                                 clinical_list$Clinical$`Patient No`)

clinical_list$Clinical$ID.gene.short <- paste(clinical_list$Clinical$`Enrollment Gene`,
                                              clinical_list$Clinical$`Patient No`)

################################################################################
#### Spider Plots ####
################################################################################
### Response Spider
if(response){
  DataSpidResponse <- merge(clinical_list$Clinical,
                            clinical_list$Response)
  
  # revalue all Baseline % change to 0 unless it contains an NA (omit downstream)
  DataSpidResponse[which(DataSpidResponse$Assessment == "Baseline" & !is.na(DataSpidResponse$`% Change from Baseline`)),"% Change from Baseline"] <- 0
  
  DataSpidResponse$AssessmentTimePoint <- DataSpidResponse$`Response Weeks`
  DataSpidResponse <- DataSpidResponse[order(DataSpidResponse$AssessmentTimePoint,decreasing = T),]
  DataSpidResponse <- DataSpidResponse[!is.na(DataSpidResponse$`% Change from Baseline`),]
  
  DataSpidResponse[DataSpidResponse$Assessment == "Baseline","AssessmentTimePoint"] <- 0
  DataSpidResponse <- DataSpidResponse[which(DataSpidResponse$`Size (mm)` != "Missing"),]
  
  clinical_list$DataSpidResponse <- DataSpidResponse
}  

if(markers){
  ### Marker Spider

  DataSpidMarkers <- merge(clinical_list$Clinical,
                           clinical_list$Markers)
  DataSpidMarkers <- DataSpidMarkers[which(!is.na(DataSpidMarkers$`Tumor Marker Value`)),]
  DataSpidMarkers <- DataSpidMarkers[which(!is.na(DataSpidMarkers$`Tumor Marker Weeks`)),]
  DataSpidMarkers[which(is.na(DataSpidMarkers$`Tumor Marker % Change From Baseline`)),"Tumor Marker % Change From Baseline"] <- 100
  
  clinical_list$DataSpidMarkers <- DataSpidMarkers
}

### Make output tables for biostats
dfmVAFRBest <- read_xlsx(inctDNA,sheet=9)

dfmVAFRBest <- dfmVAFRBest %>% 
  mutate(Molecular.Response = case_when(`Molecular Response` == "Response" ~ "MR",`Molecular Response` == "Stable" ~ "Non-MR",`Molecular Response` == "Progression" ~ "Non-MR",TRUE~"")) %>%
  select(`Patient ID` = Patient_ID,
         `Molecular Response` = `Molecular.Response`,
         `Best mVAFR` = mVAFR) 

# no ctDNA response data outputted to stats before approval from translational
stats_out <- clinical_list$Genomic %>% 
  select(`Patient ID`,`Enrollment Gene`,TMB,MSI,Reversion,`Enrollment Test Type`,`germline status`,`Type of loss`,`ATM IHC`) %>% 
  #full_join(clinical_list$Clinical) %>%
  left_join(dfmVAFRBest,by="Patient ID")


##### Write Analysis Folder

outFolder <- 
  paste0(output_dir, date)

# dir.create(output_dir, showWarnings = FALSE,recursive = T)

dir.create(outFolder, showWarnings = FALSE, recursive = T)
outFolderHeader <- paste0(outFolder, "/",
                          date,"_",trial, "_processed_clinical_data")

# === add readme

clinical_list$readme <- paste("Clinical data was processed using input data: ",
                              inClinical,"on",date, "and the output directory is: ",
                              outFolder)

# === Save
save(clinical_list,
     file = paste(outFolderHeader, ".Rda", sep = ""))

print(paste("Clinical data is saved as an .Rda object in the following directory:",
            paste0(outFolderHeader, ".Rda")))

# === Write
write.xlsx(stats_out,paste(outFolder,"/",date,"_", trial, "_Genomics_Summary.xlsx",sep=""))

