################################################################################
### Initilize environment
################################################################################

library(tidyverse)
library(openxlsx)
library(readxl)
library(stringi)
library(dplyr, warn.conflicts=FALSE)   # Manipulating the dataframes
library(reshape2) # Reformmating dataframes a
library(argparser)

print(paste0("The current working directory is:"))
getwd()

## === Initialize argparse =====
p <- arg_parser("Process clinical data and output into a list of dataframes")
p <- add_argument(p, arg = "--file", 
                  help = "specify the input Patient tracker data, use relative or absolute paths")
p <- add_argument(p, arg = "--ctdna", 
                  help = "specify the ctDNA file")
#omit, use regex
#p <- add_argument(p, arg = "--date", 
#                  help = "specify the date that data is process, use the date on Patient Tracker")
# p <- add_argument(p, arg = "--trial", 
#                   help = "specify the trial, i.e) RP3500-01")
p <- add_argument(p, arg = "--output_dir", 
                  help = "specify the path of the output directory, can be relative or absolute paths. Include '/'")
p <- add_argument(p, arg = "--git_extension", 
                  help = "specify the path of the gitextension, can be relative or absolute paths. Include '/'")

argv <- parse_args(p)

#date <- argv$date

#trial <- argv$trial

git_extension <- argv$git_extension

output_dir <- argv$output_dir

# at some point we have to install the extension as a R system library using devtools::install()
devtools::load_all(git_extension)


inClinical <- argv$file
inctDNA <- argv$ctdna
date <- str_extract(string = inClinical, "\\d{4}-\\d{2}-\\d{2}")

trial <- str_extract(string = inClinical, "RP\\d{4}-\\d{2}")


################################################################################
### Interactive Session
################################################################################
if(FALSE){
  inClinical <- "/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-01/clinical/vendor-data/2022-10-26_RP3500-01 Patient Tracker.xlsx"
  inctDNA <- "/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-01/ctDNA/processing/Combined/results/2022-08-11_RP3500-01_processed_ctDNA_Combined_results/2022-08-11_RP3500-01_processed_ctDNA_Combined_Summary.xlsx"

  git_extension <- "/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoextensions/"
  devtools::load_all(git_extension)
  output_dir <- "/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-01/clinical/processing/data/"
  
  date <- str_extract(string = inClinical, "\\d{4}-\\d{2}-\\d{2}")
  
  trial <- str_extract(string = inClinical, "RP\\d{4}-\\d{2}")
}

################################################################################
### Load dataset
################################################################################
#print msg loading in 
message(paste0("Loading in Patient Tracker file from: ", inClinical))

message(paste("The extraction date of the Patient Tracker file is:", date, 
              "from the following trial:", trial))

# from the patient tracker spread sheet, load in every sheet as a df in a named list
  suppressMessages(
  clinical_list <- readxlAllSheets(xlsxFile = inClinical,
                                               type.convert = FALSE)
)

# the names of each df in the clinical_list will inherit the sheet names from the excel sheet
print("Extracting the following excel sheets:")
names(clinical_list)

# rename list elements to match isilverman's variables 
to_rename <- c("Clinical", "Response", "Markers",
               "DoseTransfusions",  "Genomic")

to_keep <- c("Clinical", "Response", "Tumor markers", 
             "Dose Mods_Transfusions", "Genomics")

if("SNIPDx" %in% names(clinical_list)){
  print("adding SNIPDx tab to clinical data")
  
  to_rename <- c(to_rename, "SNIPDx")
  
  to_keep <- c(to_keep, "SNIPDx")
  
  print(paste("Keeping the following clinical data from the excel sheets:",
              to_keep))

    clinical_list <- clinical_list[to_keep]
  names(clinical_list) <- to_rename
  
}else{
  print(paste("Keeping the following clinical data from the excel sheets:", 
              to_keep))
  
  clinical_list <- clinical_list[to_keep]
  names(clinical_list) <- to_rename
  
}

# remove columns containing "ignore"
clinical_list <- lapply(clinical_list, function(x){
  x <- x[,grep("Ignore", colnames(x), invert = T)]
  #remove empty columns created from readxl
  x <-x[,grep('[...][0-9]', colnames(x),invert = T)]
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

### Change Timestamps,
# ensure dates are set to the date Patient Tracker was copied to vendor-data
# internal function in excel sheet that updates Tx Discont, Days, Weeks, Cycles. to the current date upon opening 
if(any(clinical_list$Clinical$`Tx Discont.` > as.POSIXct(date) & 
       !is.na(clinical_list$Clinical$`Tx Discont.`))) {
  print(paste("Tx. Discont. column contains timestamps dated after",
              "the date copied to vendor-data direcory on", date, ".",
              "Correcting dates back to,", date, "..."))
  
  clinical_list$Clinical <- 
    clinicalbioinfoextensions::correctClinicalTimestamp(data = clinical_list$Clinical,
                                                        date = date)
}

# As of 10-22-2021, merge SNiPDx with Genomics prior to merging Genomics with Clinical
if("SNIPDx" %in% names(clinical_list)){
  print("Merging Genomics data with SNIPDX data")
  
  clinical_list$Genomic <- merge(clinical_list$Genomic,
                                 clinical_list$SNIPDx)
  
  colsGenomic <- colnames(clinical_list$Genomic)
  #sometimes Column orders shift in the Patient Tracker so we have to specify by name
  to_keep <- colsGenomic[colsGenomic %in% c("Patient ID", "Enrollment Gene","cDNA Change",
                                            "AA Change",  "Other Alteration", "Alteration (S)" ,
                                            "Alteration Classification",  "Enrollment VAF", 
                                            "TP53 Status", "Other Drivers", "TMB", 
                                            "MSI","Reversion", "Enrollment Test Provider", "Enrollment Test", "Enrollment Test Type", 
                                            "germline status", "Germline call reason", "Liquid biopsy Guardant Baseline",
                                            "Liquid Biopsy Tempus Baseline", "ATM IHC","RNASEH2 IHC",
                                            "baseline gH2AX", "Paired Biopsy (Y/N)", "SNiPDx PBMC",
                                            "SNiPDx Tumor Mutation", "SNiPDx Gene","SNIPDx cDNA Change",
                                            "SNiPDx AA Change", "SNiPDx Other", "SNiPDx Alteration",
                                            "SNiPDx Alteration Type",  "SNiPDx Tumor VAF", "SNiPDx Tumor CN",
                                            "Type of loss")]
  # to_keep <-colnames(clinical_list$Genomic)[c(1:34)]
  
  print("Merging Clinical data with Genomics data")
  
  # clinical_list$Clinical <- merge(clinical_list$Clinical[,c(1:ncol(clinical_list$Clinical) - 1)],
  #                                 clinical_list$Genomic[,to_keep],all.x=T)
  
  clinical_list$Clinical <- clinical_list$Clinical %>%
    select(-any_of("Notes")) %>%
    left_join(select(clinical_list$Genomic, any_of(to_keep)))
}else{
  print("Merging Clinical data with Genomics data")
  
  ## Merge Clincial Genomic and Duration of Treatment Data Tables
  clinical_list$Clinical <- merge(clinical_list$Clinical[,c(1:ncol(clinical_list$Clinical) - 1)],
                                  clinical_list$Genomic[,c(1:16)],all.x=T)
}


### Set GeneGroups and Merge "other STEP2" alterations
clinical_list$Clinical <- setGeneGroup(df = clinical_list$Clinical)

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


### Derive PSA Response 
# coerce string "NA" to NA in  `Marker Response` levels
clinical_list$Markers[clinical_list$Markers$`Marker Response` %in% "NA",]$`Marker Response` <- NA

suppressWarnings(
  clinical_list$BestMarkerResponseTime <- dcast(clinical_list$Markers,
                                                `Patient ID`+`Tumor marker`~`Marker Response`,
                                                value.var="Tumor Marker Weeks",
                                                function(x) min(x,na.rm = T))
) #end warning msg


clinical_list$BestMarkerResponseTime <- 
  clinical_list$BestMarkerResponseTime %>%
  select(`Patient ID`, `Tumor marker`, Response)

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
clinical_list$Clinical <- clinical_list$Clinical %>% mutate(COR = case_when(
  BOR == "CR"~"CR",
  BOR == "PR"~"PR",
  PSA == "Response" ~"PR",
  `CA-125` == "Response" ~"PR",
  BOR == "SD"~"SD",
  BOR == "PD"~"PD"))

### Derive Actualy Dose
Dose <- clinical_list$DoseTransfusions[clinical_list$DoseTransfusions$`Event Type` == "Dose",]
#was this suppose to be clinical_list$Clinical <- merge(clinical_list$Clinical, Dose) instead?
#Dose df is not exported to .Rda
Dose <- merge(clinical_list$Clinical,Dose)

#### Generate Short and Long IDs, group Modules
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
# consolidate panels
clinical_list$Clinical <- clinical_list$Clinical %>% 
  mutate(Consolidated.Modules  = case_when(grepl("^1", Module) ~ "M1",
                                           grepl("^2_Arm", Module) ~ "M2",
                                           grepl("^3", Module) ~ "M3",
                                           grepl("^4", Module) ~ "M4"))

###  Spider Plots
DataSpidResponse <- merge(clinical_list$Clinical,clinical_list$Response)

# revalue all Baseline % change to 0 unless it contains an NA (omit downstream)
DataSpidResponse[which(DataSpidResponse$Assessment == "Baseline" & !is.na(DataSpidResponse$`% Change from Baseline`)),"% Change from Baseline"] <- 0

DataSpidResponse$AssessmentTimePoint <- DataSpidResponse$`Response Weeks`
DataSpidResponse <- DataSpidResponse[order(DataSpidResponse$AssessmentTimePoint,decreasing = T),]
DataSpidResponse <- DataSpidResponse[!is.na(DataSpidResponse$`% Change from Baseline`),]

DataSpidResponse[DataSpidResponse$Assessment == "Baseline","AssessmentTimePoint"] <- 0
DataSpidResponse <- DataSpidResponse[which(DataSpidResponse$`Size (mm)` != "Missing"),]

clinical_list$DataSpidResponse <- DataSpidResponse

### Marker Spider
print("Merging Clinical data with BestSpidMarker data")

DataSpidMarkers <- merge(clinical_list$Clinical,
                         clinical_list$Markers)
DataSpidMarkers <- DataSpidMarkers[which(!is.na(DataSpidMarkers$`Tumor Marker Value`)),]
DataSpidMarkers <- DataSpidMarkers[which(!is.na(DataSpidMarkers$`Tumor Marker Weeks`)),]
DataSpidMarkers[which(is.na(DataSpidMarkers$`Tumor Marker % Change From Baseline`)),"Tumor Marker % Change From Baseline"] <- 100

clinical_list$DataSpidMarkers <- DataSpidMarkers

### Make output tables for biostats
dfmVAFRBest <- read_xlsx(inctDNA,sheet=9)

dfmVAFRBest <- dfmVAFRBest %>% mutate(Molecular.Response = case_when(`Molecular Response` == "Response" ~ "MR",`Molecular Response` == "Stable" ~ "Non-MR",`Molecular Response` == "Progression" ~ "Non-MR",TRUE~"")) %>%
  select(Patient_ID,`Molecular.Response`,mVAFR) %>% rename(c( `Patient ID` = Patient_ID ,`Molecular Response` = Molecular.Response, `Best mVAFR` = mVAFR))

# no ctDNA reponse outputs to stats until approved by Translational
stats_out <- clinical_list$Clinical %>% select(`Patient ID`,`Enrollment Gene`,TMB,MSI,Reversion,`Enrollment Test Type`,`germline status`,`Type of loss`,`ATM IHC`) %>% left_join(dfmVAFRBest,by="Patient ID")

##### Write Analysis Folder

outFolder <- 
  paste0(output_dir, date)

dir.create(outFolder, showWarnings = FALSE)

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
            paste(outFolderHeader, ".Rda")))

# === Write
write.xlsx(stats_out,paste(outFolder,"/",date,"_RP3500-01_Genomics_Summary.xlsx",sep=""))


# }) # end warning msg

