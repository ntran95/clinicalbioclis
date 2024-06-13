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
library(cowplot)
library(data.table)

################################################################################
### Parse arguments
################################################################################

# === Initialize
p <- arg_parser("Process ctDNA and output a ctDNA object")
p <- add_argument(p, arg = "--manifest",
                  help = "manifest directory path")
p <- add_argument(p, arg = "--qc",
                  help = "qc file path")
p <- add_argument(p, arg = "--mmf",
                  help = "mmf file path")
p <- add_argument(p, arg = "--cnv",
                  help = "cnv file path")
p <- add_argument(p, arg = "--monitoring",
                  help = "monitoring file path")
p <- add_argument(p, arg = "--monitoring_avail",
                  help = "logical to indicate if monitoring file is available, set to T or F")
p <- add_argument(p, arg = "--outdir",
                  help = "output directory")
p <- add_argument(p, arg = "--functions",
                  help = "functions file")
p <- add_argument(p, arg = "--trial",
                  help = "trial number (e.g. RP6306-01")
p <- add_argument(p, arg = "--date",
                  help = "date (format: yyyy-mm-dd)")
p <- add_argument(p, arg = "--ctfe",
                  help = "path to ctfe file")
p <- add_argument(p, arg = "--ctfe_avail",
                  help = "logical, indicating whether ctfe file is available")

argv <- parse_args(p)


date <- as.Date(argv$date)
trial <- argv$trial
manifest_path <- argv$manifest
qc_path <- argv$qc
mmf_path <- argv$mmf
cnv_path <- argv$cnv
monitoring_avail <- argv$monitoring_avail
if(monitoring_avail) {monitoring_path <- argv$monitoring}
outdir <- argv$outdir
functions <- argv$functions
inctfe <- argv$ctfe
ctfe_avail <- argv$ctfe_avail

# source(functions)
devtools::load_all(functions)

################################################################################
### Set testing arguments/ Interactive session
################################################################################
interactive <- FALSE
if(interactive){
  manifest <-  '/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/ctDNA/vendor-data/Tempus/Manifests/'
  qc <-'/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/ctDNA/vendor-data/Tempus/Batch7_2022.09.23/Repare_3500-03_xFSeq_QC_2022_09_22.xlsx'
  mmf <- '/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/ctDNA/vendor-data/Tempus/Batch7_2022.09.23/g_molecular_master_file_filtered.csv'
  cnv <- '/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/ctDNA/vendor-data/Tempus/Batch7_2022.09.23/g_cnv_segments_cnvkit.csv'
  monitoring <-  '/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/ctDNA/processing/Tempus/data/monitoring/2022-10-26_RP3500-03_processed_ctDNA_Tempus_Batch7_Revised_MonitorIn_2022-10-26.xlsx'
  monitoring_avail <- T
  outdir <- '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/data/'
  functions <-  '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/src/functions_ctDNA.R'
  trial <- 'RP3500-03'
  date <- '2022-11-07'
  inctfe <- "/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/ctDNA/vendor-data/Tempus/Batch7_2022.09.23/Repare_3500-03_xFSeq_ctFE.xlsx"
  ctfe_avail <- T

  
  date <- as.Date(date)
  manifest_path <- manifest
  qc_path <- qc
  mmf_path <- mmf
  cnv_path <- cnv
  if(monitoring_avail) {monitoring_path <- monitoring}

  #source(functions)
  devtools::load_all(functions)
  
}

################################################################################
### Set testing arguments
################################################################################
#
# date <- ymd("2022-03-05")
# trial <- "RP6306-01"
# outdir <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-01/ctDNA/processing/Tempus/data/"
#
# manifest_path <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-01/ctDNA/vendor-data/Tempus/Manifests/Tempus Manifest/"
# mmf_path <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-01/ctDNA/vendor-data/Tempus/Batch5_20220303/g_molecular_master_file_filtered.csv"
# qc_path <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-01/ctDNA/vendor-data/Tempus/Quality Reports/RP6306-01_Tempus_ctDNA_Batch5_Quality-Report_03Mar22.xlsx"
# cnv_path <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-01/ctDNA/vendor-data/Tempus/Batch5_20220303/g_cnv_segments_cnvkit.csv"
# monitoring_path <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-01/ctDNA/processing/Tempus/monitor/2022-02-28_RP6306-01_monitoring_completed.xlsx"
# monitoring_avail <- T
#
################################################################################
### Read and combine manifest data from all available manifests
################################################################################

manifest <- list.files(manifest_path, full.names = T,pattern = "Manifest",recursive = T)
manifest <- lapply(manifest, FUN= function(x) {read_excel(x, skip =14,na = c(""))}) %>% rbindlist()
manifest <- readr::type_convert(manifest)
manifest$`Date of Collection` <- as.Date(manifest$`Date of Collection`, tryFormats = date_formats) # different date formats across trials
manifest$`Patient ID` <- as.character(manifest$`Patient ID`)
manifest$`Sample ID` <- as.character(manifest$`Sample ID`)
#manifest$`Sample ID` <- substr(manifest$`Sample ID`, 1, 11) # not applicable to RP3500-03
manifest <- unique(manifest)

################################################################################
### Read and clean QC file
################################################################################

qclist <- readxlAllSheets(xlsxFile = qc_path,
                          check.names = F,
                          skip = 4)
qclist <- qclist[ grepl("Delivery", names(qclist))]
qclist <- lapply(qclist, function(x){ x <- readr::type_convert(x) })
qclist <- rev(qclist)
# label batches 
for(currBatch in names(qclist)){
  qclist[[currBatch]]$Batch <- currBatch
}

qc <- qclist %>% purrr::reduce(merge, all = T)

qc <- filter(qc, `Path Review Decision` == "Approved")
qc$cfDNA_ng.ml <- as.numeric(qc$`DNA QC1 Yield (ng)`) / as.numeric(qc$`Plasma Volume (mL)`)

qc$`Study Sample ID` <- as.character(qc$`Study Sample ID`)
qc$`Study Patient ID` <- as.character(qc$`Study Patient ID`)
qc$`Sample Collection Date` <- as.Date(qc$`Sample Collection Date`, tryFormats = date_formats) # different date formats across trials

### correct Prescreening sample IDs, match mmf
screeningSamples <- grep("pre|screening|\\d{4}-\\d{4}", 
                         qc$`Study Sample ID`,
                         ignore.case = T )

if(length(screeningSamples) != 0){
  qc[screeningSamples, "Study Sample ID"] <-  paste(qc[screeningSamples, ]$`Study Patient ID`, "Pre-Screening",sep = "_")
  
}


#Clean up of manual errors - this cannot be built into code, need a better solution here
# qc <- separate(qc, `Study Sample ID`, into = c("Study Sample ID", "unit"), sep = "_", extra = "merge")
# qc <- separate_rows(qc, unit, sep = "_")
# qc$`Study Sample ID` <- paste0(qc$`Study Sample ID`, qc$unit)
# qc <- select(qc, -unit)

################################################################################
### Read and clean molecular master file
################################################################################
mmf <- read.csv(mmf_path)
#mmf$partner_sample_id <- substr(mmf$partner_sample_id, 1, 11) ### not applicable to rp3500-03
mmf <- select(mmf, all_of(variant_cols))
mmf <- filter(mmf, variant_type == "Short Variant" & ! functional_impact %in% c("B", "LB"))
mmf$VAF <- mmf$variant_allele_freq * 100
mmf$amino_acid_change <- AAA2A(mmf$amino_acid_change)
mmf <- mmf %>% mutate(amino_acid_change = case_when(is.na(amino_acid_change) | amino_acid_change == "" ~ "",
                                                    !grepl("^p.", amino_acid_change) ~ paste0("p.", amino_acid_change),
                                                    TRUE ~ amino_acid_change))

mmf$genomic_change <- paste0("chr", mmf$chromosome, ":", mmf$position_1,":", mmf$reference, ">", mmf$alternative)


mmf <- mmf %>% mutate( variant_effect = case_when(grepl("fs$", amino_acid_change) ~ "Frameshift",
                                                  grepl("del|ins", amino_acid_change) ~ "In-frame indel",
                                                  grepl("dup", amino_acid_change) ~"In-frame indel",
                                                  grepl("\\*", amino_acid_change) ~ "Nonsense",
                                                  grepl("\\?", amino_acid_change) ~ "Start Loss",
                                                  gsub("p\\.([A-Z]+)[0-9].*", "\\1", amino_acid_change) == gsub(".*[0-9]([A-Z]+)", "\\1", amino_acid_change) ~ "Synonymous",
                                                  gsub("p\\.([A-Z]+)[0-9].*", "\\1", amino_acid_change) != gsub(".*[0-9]([A-Z]+)", "\\1", amino_acid_change) ~ "Missense",
                                                  grepl("c\\.[-|\\+]", nucleotide_change) ~ "Promoter",
                                                  grepl("c\\.[0-9]+[-|\\+][0-9]", nucleotide_change) ~ "Splice"
  ))


### correct Patient ID
mmf$partner_patient_id <- gsub("^(.{4})(.*)$","\\1-\\2",gsub("-","",mmf$partner_patient_id))

### correct Prescreening sample IDs, match Manifest/QC
# match screening samples with ManifestQC
# screening samples in mmfk contain either "pre", "screening" or just the Patient ID (XXXX-XXXX), standardize screening sample ids to merge downstream with TempusQC and Manifest
screeningSamples <- grep("pre|screening|\\d{4}-\\d{4}", 
                         mmf$partner_sample_id,
                         ignore.case = T )

if(length(screeningSamples) != 0){
  mmf[screeningSamples, "partner_sample_id"] <-  paste(mmf[screeningSamples, ]$partner_patient_id, "Pre-Screening",sep = "_")
  
}

# ensure empty cells are converted to NAs in amino_acid_change column
# prevents errors in 1:1 merging by variants
mmf <-  mmf %>% mutate_all(na_if,"")


#Clean up of manual errors - this cannot be built into code, need a better solution here
# mmf$partner_sample_id <- gsub("-", "", mmf$partner_sample_id)
# mmf$partner_sample_id <- gsub("03_", "_03_", mmf$partner_sample_id)
# mmf$partner_sample_id <- gsub("__", "_", mmf$partner_sample_id)
# mmf <- separate(mmf, partner_sample_id, into = c("partner_sample_id", "unit"), sep = "_", extra = "merge")
# mmf <- separate_rows(mmf, unit, sep = "_")
# mmf$partner_sample_id <- paste0(mmf$partner_sample_id, mmf$unit)
# mmf <- select(mmf, -unit)

################################################################################
### Read in CNV segments file
################################################################################

cnv <- read.csv(cnv_path)
#cnv$partner_sample_id <- substr(cnv$partner_sample_id, 1, 11) ## not applicable to RP3500-03
cnv <- select(cnv, all_of(cnv_cols))
cnv <- separate_rows(cnv, gene, sep = ",")
### correct Patient ID
cnv$partner_patient_id <- gsub("^(.{4})(.*)$","\\1-\\2",gsub("-","",cnv$partner_patient_id))

#Clean up of manual errors - this cannot be built into code, need a better solution here
# cnv$partner_sample_id <- gsub("-", "", cnv$partner_sample_id)
# cnv$partner_sample_id <- gsub("03_", "_03_", cnv$partner_sample_id)
# cnv$partner_sample_id <- gsub("__", "_", cnv$partner_sample_id)
# cnv <- separate(cnv, partner_sample_id, into = c("partner_sample_id", "unit"), sep = "_", extra = "merge")
# cnv <- separate_rows(cnv, unit, sep = "_")
# cnv$partner_sample_id <- paste0(cnv$partner_sample_id, cnv$unit)
# cnv <- select(cnv, -unit)
# cnv <- unique(cnv)

################################################################################
### Read monitoring file if available
################################################################################

if(monitoring_avail) {

  monitoring <- read_excel(monitoring_path)

}

################################################################################
### Merge manifest and qc information
################################################################################

assay <- merge(manifest, qc, by.x=c("Patient ID","Sample ID"), by.y=c("Study Patient ID","Study Sample ID"), all.y=T)
assay <- select(assay, any_of(assay_cols))
### correct Patient ID
assay$`Patient ID` <- gsub("^(.{4})(.*)$","\\1-\\2",gsub("-","",assay$`Patient ID`))
### if Date of Collection is NA, inherit Sample Collection Date (applicable to RP3500-03)
assay$`Date of Collection` <- as.Date(ifelse(is.na(assay$`Date of Collection`),yes = as.character(assay$`Sample Collection Date`),
                                     no = as.character(assay$`Date of Collection`)))


################################################################################
### Read ctfe file
################################################################################
if(ctfe_avail){
  ctfe <- read_excel(inctfe,skip = 5)
  ctfe$`Study Patient ID` <- gsub("^(.{4})(.*)$","\\1-\\2",gsub("-","",ctfe$`Study Patient ID`))
  ctfe$`Study Sample ID` <- as.character(ctfe$`Study Sample ID`)
  
  ### correct Prescreening sample IDs, match Manifest/QC
  # screening samples in mmfk contain either "pre", "screening" or just the Patient ID (XXXX-XXXX), standardize screening sample ids to merge downstream with TempusQC and Manifest
  screeningSamples <- grep("pre|screening|\\d{4}-\\d{4}", 
                           ctfe$`Study Sample ID`,
                           ignore.case = T )
  
  if(length(screeningSamples) != 0){
    ctfe[screeningSamples, "Study Sample ID"] <-  paste(ctfe[screeningSamples, ]$`Study Patient ID`, "Pre-Screening",sep = "_")
    
  }
}

################################################################################
### Create ctDNA object and write to file
################################################################################

ctdna <- new("ctDNA",
              assay = assay,
              variants = mmf,
              cnv = cnv)
if(monitoring_avail) { ctdna@monitoring = monitoring }
if(ctfe_avail) { ctdna@ctfe = ctfe }

dir.create(outdir, showWarnings = F,recursive = T)

saveRDS(ctdna, paste0(outdir, "/", date, "_", trial, "_processedCtDNA.rds"))

################################################################################



