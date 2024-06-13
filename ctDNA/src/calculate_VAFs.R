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
p <- arg_parser("Calculate mAVAFs and generate a monitoring")
p <- add_argument(p, arg = "--study",
                  help = "study object path")
p <- add_argument(p, arg = "--ctdna",
                  help = "ctdna object path")
p <- add_argument(p, arg = "--mondir",
                  help = "monitoring directory path")
p <- add_argument(p, arg = "--datadir",
                  help = "output directory")
p <- add_argument(p, arg = "--monitoring_avail",
                  help = "logical to indicate if monitoring file is available, set to T or F")
p <- add_argument(p, arg = "--functions",
                  help = "path to functions extension folder (e.g clinicalbioinfoextensions)")
p <- add_argument(p, arg = "--trial",
                  help = "trial number (e.g. RP6306-01")
p <- add_argument(p, arg = "--date",
                  help = "date (format: yyyy-mm-dd)")
p <- add_argument(p, arg = "--rmarkdown_report",
                  help = "logical indicating whether to generate unfiltered VAF rmarkdown report",default = FALSE)
p <- add_argument(p, arg = "--resultsdir",
                  help = "results directory, applicable when rmarkdown_report is TRUE")
p <- add_argument(p, arg = "--rmd_path",
                  help = "path to unfiltered rmd script, applicable when rmarkdown_report is TRUE")

argv <- parse_args(p)


date <- as.Date(argv$date)
trial <- argv$trial
study_path <- argv$study
ctdna_path <- argv$ctdna
mondir <- argv$mondir
datadir <- argv$datadir
monitoring_avail <- argv$monitoring_avail
functions <- argv$functions
rmarkdown_report <- argv$rmarkdown_report
if(rmarkdown_report){
  resultsdir <- argv$resultsdir
  rmd_path <- argv$rmd_path
}

devtools::load_all(functions)

################################################################################
### Set testing arguments
################################################################################

# date <- ymd("2022-03-09")
# trial <- "RP6306-02"
# study_path <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-02/clinical/2022-03-09_RP6306-02_processedStudy.rds"
# ctdna_path <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-02/ctDNA/processing/Tempus/data/2022-03-09_RP6306-02_processedCtDNA.rds"
# datadir <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-02/ctDNA/processing/Tempus/data/"
# mondir <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-02/ctDNA/processing/Tempus/monitor/"
# monitoring_avail <- F

################################################################################
### Set testing arguments
################################################################################

# date <- ymd("2022-03-05")
# trial <- "RP6306-01"
# study_path <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-01/clinical/2022-02-28_RP6306-01_processedStudy.rds"
# ctdna_path <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-01/ctDNA/processing/Tempus/data/2022-03-05_RP6306-01_processedCtDNA.rds"
# datadir <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-01/ctDNA/processing/Tempus/data/"
# mondir <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-01/ctDNA/processing/Tempus/monitor/"
# monitoring_avail <- T

################################################################################
### Set testing arguments
################################################################################
interactive <- FALSE
if(interactive){
  date <- ymd("2022-12-06")
  trial <- "RP3500-03"
  # study_path <- "/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/2022-10-31_RP3500-03_processedStudy.rds"
  study_path <- "/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/clinical/processing/data/2022-12-06/2022-12-06_RP3500-03_processedStudy.rds"
  
  # ctdna_path <- "/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/data/2022-11-07_RP3500-03_processedCtDNA.rds"
  ctdna_path <- "/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/ctDNA/processing/Tempus/data/2022-11-22_RP3500-03_processedCtDNA.rds"
  
  datadir <- "/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/data/"
  mondir <- "/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/data/monitoring/"
  monitoring_avail <- T
  functions <-  '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoextensions/'
  rmarkdown_report <- T
  resultsdir <- "/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/results/"
  rmd_path <- "/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/src/generate_rmardown_report_unfiltered_ctDNA.Rmd"
  
  devtools::load_all(functions)
  
  
  
}


################################################################################
### Read in study and ctDNA objects and add study metadata to ctDNA
################################################################################

study <- readRDS(study_path)
ctdna <- readRDS(ctdna_path)

ctdna@metadata <- study@combined

################################################################################
### Combine mutation data with metadata and filter screen fail samples
################################################################################

### Copy date of collection from assay data
variants <- merge(ctdna@variants, select(ctdna@assay, `Patient ID`, `Sample ID`, `Date of Collection`, `cfDNA_ng.ml`, Batch, `Study Timepoint`, `Sample Collection Date`),
                  by.x = c("partner_patient_id", "partner_sample_id"),
                  by.y = c("Patient ID", "Sample ID"),
                  all.x = T)

### Merge with clinical meta data
variants <- merge(ctdna@metadata, variants, by.x = c("Patient ID"), by.y = c("partner_patient_id"), all.y = T)


  

### Save screen fail and filter from variants table
screenfail <- variants %>% filter(is.na(`First dose`))
screenfail <- screenfail[ , colSums(is.na(screenfail)) < nrow(screenfail)] 


################################################################################
### Combine CNV data with metadata 
################################################################################

### Copy date of collection from assay data
cnv <- merge(ctdna@cnv, select(ctdna@assay, `Patient ID`, `Sample ID`, `Date of Collection`, Batch, `Study Timepoint`, `Sample Collection Date`),
             by.x = c("partner_patient_id", "partner_sample_id"),
             by.y = c("Patient ID", "Sample ID"),
             all.x = T)

### Merge with clinical meta data
cnv <- merge(ctdna@metadata, cnv, by.x = c("Patient ID"), by.y = c("partner_patient_id"), all.y = T)

### Filter screenfail samples

cnv <- cnv %>% filter(!is.na(`First dose`))

################################################################################
### Combine ctFe estimates 
################################################################################
# merge ctfe estimates if df is not empty
if(nrow(ctdna@ctfe) != 0){
  ctfe <- ctdna@ctfe
  
  # reomve dup sample IDs
  ctfe <-ctfe %>%
    distinct(`Study Sample ID`, .keep_all = T) %>%
    mutate(ctFE = ctFE * 100)
  
  ctfe$`Study Sample ID` <- as.character(ctfe$`Study Sample ID`)
  
  variants <- variants %>%
    left_join(ctfe, by = c("Patient ID" = "Study Patient ID",
                            "partner_sample_id" = "Study Sample ID"))
}

################################################################################
### Calculate visit ID
################################################################################
### ensure all samples missing first dose have a pseudo first dose, inherit Date of Collection from Tempus
variants$`First dose` <- as.Date(ifelse(is.na(variants$`First dose`),
                                        as.character.Date(variants$`Date of Collection`),
                                        as.character.Date(variants$`First dose`)))
# 
# ### ensure all samples with missing Date of Collection (start date) to calculate TimeFromC1D1, inherit Sample Collection Date from TempusQC if missing
# variants$`Date of Collection` <- as.Date(ifelse(is.na(variants$`Date of Collection`),
#                                                 as.character.Date(variants$`Sample Collection Date`),
#                                                 as.character.Date(variants$`Date of Collection`)))

variants <- variants %>% mutate( TimeFromC1D1 = as.numeric(`Date of Collection` - `First dose`))

# baselines labelled between -21 and 10 days
variants <- variants %>% mutate( VisitID = calcAdjVisitID(TimeFromC1D1))

#### Baseline Corrections ####
# revalue baseline samples to pre-tx if there are dups in each patient
samplesAnalyzed <- variants %>%
  select(`Patient ID`, TimeFromC1D1, VisitID, partner_sample_id) %>% # select SamplesAnalyzed column
  distinct() %>% # remove dups form variants, sample level only
  group_by(`Patient ID`, VisitID) %>% 
  mutate(n.baseline = sum(VisitID %in% "Baseline"))  %>% # calculate number of baselines per patient
  mutate(n.negative = sum(VisitID %in% "Baseline" & TimeFromC1D1 < 0)) %>% # count number of positive baselines
  mutate(n.positive = sum(VisitID %in% "Baseline" & TimeFromC1D1 > 0)) %>% # count number of negative baselines
  ungroup()

############################### test baseline logic ############################
if(FALSE){
samplesAnalyzed<-  data.frame("Patient ID" = c("001","001", "001", "002", "002"),
                              "TimeFromC1D1" = c(-114, -1, 0, -1,-5),
                              "VisitID" = c("Pre-Tx", "Baseline", "Baseline", "Baseline", "Cycle 1"),
                              "partner_sample_id" = seq(1,5),check.names = F) %>%
  distinct() %>% # remove dups form variants, sample level only
  group_by(`Patient ID`, VisitID) %>% 
  mutate(n.baseline = sum(VisitID %in% "Baseline")) %>%
  mutate(n.negative = sum(VisitID %in% "Baseline" & TimeFromC1D1 < 0)) %>%
  mutate(n.positive = sum(VisitID %in% "Baseline" & TimeFromC1D1 > 0)) %>%
  ungroup()
}

# correct VisitIDs with more than 2 baselines
baselineCorrections <- samplesAnalyzed %>%
  filter(n.baseline >=2) %>% # filter samples with >=2 baselines per patient 
  group_by(`Patient ID`) %>%
  mutate(VisitID = case_when(n.positive >= 2 & abs(TimeFromC1D1 - 0) == min(abs(TimeFromC1D1 - 0)) ~ "Baseline", # if both BLs are positive days from C1D1, choose the value closest to 0
                             n.negative >= 2 & abs(TimeFromC1D1 - 0) == min(abs(TimeFromC1D1 - 0)) ~ "Baseline", # if both BLs are negative days from C1D1, choose the value closest to 0
                             n.positive >= 1 & n.negative >=1 & TimeFromC1D1 == min(TimeFromC1D1)  ~ "Baseline",# if sample has 2 or more at baseline, revalue the min TimeFromC1D1 to Baseline
                             n.positive >= 0 & n.negative >=1 & TimeFromC1D1 == min(TimeFromC1D1)  ~ "Baseline",# if sample has 2 or more at baseline, revalue the min TimeFromC1D1 to Baseline
                             TRUE ~ as.character("Pre-Tx"))) %>% # label others as Pre-tx
  ungroup()

# merge back with samples analyzed
samplesAnalyzed <- filter(samplesAnalyzed, partner_sample_id %notin% baselineCorrections$partner_sample_id) %>% # remove dup samples present in baselineCorrections
  full_join(baselineCorrections) %>% #merge 
  arrange(`Patient ID`, TimeFromC1D1) %>% # order by patient ID and ascending timepoints 
  select(-n.baseline, -n.negative, -n.positive) %>% # remove baseline freq (no longer needed)
  ungroup()

# merge back corrected VisitIDs
variants <- samplesAnalyzed %>%
  inner_join(select(variants, -TimeFromC1D1, -VisitID))
  

cnv <- cnv %>% mutate( TimeFromC1D1 = as.numeric(`Date of Collection` - `First dose`))
cnv <- cnv %>% mutate( VisitID = calcAdjVisitID(TimeFromC1D1))

# merge back corrected VisitIDs
cnv <- samplesAnalyzed %>%
  inner_join(select(cnv, -TimeFromC1D1, -VisitID))

### relabelling visit ids to EOT not implemented in RP6306, may be a different label, check back
### relabelling EOT timepoints, some samples are labelled as Early termination or EOT
variants <- variants %>% mutate(VisitID = case_when(`Study Timepoint` %in% c("Early Termination", "EOT") ~ "EOT",
                                        TRUE ~ as.character(VisitID))) 

cnv <- cnv %>% mutate(VisitID = case_when(`Study Timepoint` %in% c("Early Termination", "EOT") ~ "EOT",
                                                    TRUE ~ as.character(VisitID))) 




################################################################################
### Count visit per patient and  alterations per visit
################################################################################

variants <- variants %>% mutate( alteration = case_when(!is.na(amino_acid_change) & amino_acid_change !="" ~  amino_acid_change,
                                                        TRUE ~ nucleotide_change))

visits_per_patient <- variants %>% select(`Patient ID`, TimeFromC1D1) %>% unique() %>% group_by(`Patient ID`) %>% tally() %>% ungroup()
# this calculates the occurance of unique alts per patient and sample, causes diff NumVisits per alt, patient, and timepoint
# i.e) n = (alt1 + alt2 + alt3) = 3 for patient 001 at baseline
alts_per_visit <- variants %>% select(`Patient ID`, TimeFromC1D1, gene, alteration) %>% unique() %>% group_by(`Patient ID`, TimeFromC1D1) %>% tally() %>% ungroup()

#calculate the num visits per patient and alt
# i.e) n = baseline, C1, C2, C4 = 4 for patient 001 with alt p.L19F
visits_per_patient_alt <- variants %>%
  group_by(`Patient ID`,gene,alteration) %>%
  distinct( TimeFromC1D1) %>% #use TimeFromC1D1 in edge cases where Visit.ID is the same despite haveing 2 diff TimeFromC1D1
  dplyr::summarize(NumVisitsPerAlt = n()) %>%
  ungroup()

variants <- merge(variants, visits_per_patient) %>% dplyr::rename(NumVisits = n)
variants <- merge(variants, alts_per_visit) %>% dplyr::rename(NumAltsPerVisit = n)
variants <- merge(variants, visits_per_patient_alt) 

################################################################################
### Handle duplicate variants within the same patient/sample
################################################################################
# if alterations within same patient, alt, and sample contain duplicate alt with diff nucleotide change, concatenate p. and c., else keep alt as same
variants <- variants %>%
  group_by(`Patient ID`, alteration, partner_sample_id) %>% 
  #select(`Patient ID`, gene, alteration, nucleotide_change, partner_sample_id ) %>%
  # within grouped pat, sample and alt, if there are more than one occurance of the same alt, concatenate p. and c. to make unique alts
  mutate(alteration = case_when(n() > 1 ~ paste(alteration, nucleotide_change, sep = "_"), 
                                TRUE ~ as.character(alteration))) %>%
  ungroup()

### test case
if(FALSE){
  data.frame(`Patient ID` = c(001,001,002,003),
             gene = c("RET", "RET", "BAP1", "TP53"),
             alteration = c("p.L730V", "p.L730V", "p.Q28*", "c.559+1G>A"),
             nucleotide_change = c("c.2188C>G", "c.2188C>G", "c.82C>T", "c.559+1G>A"),
             partner_sample_id =  c(001,001,002,003),check.names = F)  %>%
    group_by(`Patient ID`, alteration, partner_sample_id) %>%
    #mutate(n.alt.pat.dups = n()) %>%
    select(`Patient ID`, gene, alteration, nucleotide_change, partner_sample_id ) %>%
    mutate(alteration = case_when(n() > 1 ~ paste(alteration, nucleotide_change, sep = "_"),
                                  TRUE ~ as.character(alteration))) %>%
    View()
}


################################################################################
### Calculate mean alteration VAF for each patient
################################################################################

mavaf <- variants %>% select(`Patient ID`, gene, alteration, VAF) %>% 
                      unique() %>%
                      group_by(`Patient ID`, gene, alteration) %>% 
                      summarize(mAVAF = round(mean(VAF, na.rm = T),2),  
                                sdAVAF = round(sd(VAF, na.rm = T),2))

variants <- merge(variants, mavaf)

################################################################################
### Calculate mean coverage for each sample, and mean coverage for each alteration across visits
################################################################################

# Calculate alteration min coverage across visits
altCov <- variants %>% select(`Patient ID`,`TimeFromC1D1`, `gene`, `alteration`, coverage, mAVAF, sdAVAF) %>%
  group_by(`Patient ID`, `gene`, `alteration`) %>%
  summarize(minAltCov = min(coverage, na.rm = T)) %>% 
  ungroup()

variants <- merge(variants, altCov, all.x = T)

# Calculate mean sample coverage
samcov <- variants %>% select(`Patient ID`,`TimeFromC1D1`, `gene`, `alteration`, coverage, partner_sample_id) %>% 
  unique() %>% 
  group_by(partner_sample_id) %>% 
  summarize( meanSampleCov = mean(coverage, na.rm = T)) %>%
  ungroup()

variants <- merge(variants, samcov, all.x = T)

################################################################################
### Clean
################################################################################
### organize patient, alterations, and timepoints in ascending order
variants <- variants %>%
  arrange(`Patient ID`, gene, TimeFromC1D1)

### relocate significant columns to before first column
variants <- variants %>%
  relocate(any_of(vafsummary_cols))

################################################################################
### Append or generate new monitoring file
################################################################################

if(monitoring_avail) {
  
  monitoring <- ctdna@monitoring
  # genomic_change previously labelled as GenomicPositions in RP3500-03
  monitoring <- monitoring %>%select(., any_of(monitoringIn_cols)) %>%
    # rename old monitoring columns if neccessary
    dplyr::rename(any_of(c(`Patient ID` = "Patient_ID",
                  `gene` = "Gene",
                  `alteration` = "Alteration")))

  monitoringOut <- variants %>% select(any_of(monitoring_cols)) %>%
    distinct(`Patient ID`, gene,  nucleotide_change,  amino_acid_change, .keep_all = T)
  
   monitoring <- left_join(monitoringOut, monitoring ) %>%
     unique() %>%
     relocate(any_of(c("Status", "Enrollment","Notes", "Annotation_Date","Annotater")),
              .after = last_col()) %>%
     arrange(`Patient ID`, gene, alteration)
  
} else {
  
  keep_cols <- intersect(colnames(variants), monitoring_cols)
  add_cols <- setdiff(monitoring_cols, colnames(variants)) 
  monitoring <- select(variants, all_of(keep_cols)) %>% unique()
  monitoring[add_cols] <- NA
  
}

### relabel Status if alteration level mean cov cutoff is "low coverage"
monitoring <- monitoring %>%
  mutate(Status = case_when(minAltCov <= 500 ~ "low variant mean coverage <= 500",
                            #else leave Status as is
                            TRUE ~ as.character(Status)))

################################################################################
### Create VAF plots
################################################################################

VAFPlots <- list()
VAFPlotsWatermark <- list()

for(patient in unique(variants$`Patient ID`)){
  
  print(patient)
  
  ### Filter unique alterations for given patient
  
  plotdf <- filter(variants, `Patient ID` == patient & VAF > 0)
  
  plotdf <- plotdf %>% 
    mutate( ordalteration = paste0(gene, " ", alteration)) %>% 
    arrange(-mAVAF) %>%
    mutate(ordalteration = factor(ordalteration, levels = unique(ordalteration)))
    
  title <- unique(paste(plotdf$`Enrollment Gene`, patient , plotdf$`Tissue Type`))
  
  ### if pre-screen sample or if patient is not entered in clinical 
  clinical_avail <- unique(plotdf$ID.gene)
  if(is.na(clinical_avail)){
    title <- paste(patient, unique(plotdf$partner_sample_id))
  }
  
  ### ggplot              
  
  p <- ggplot(plotdf, aes(x=TimeFromC1D1, y=VAF/100, color=ordalteration, group=ordalteration))
  p <- p + geom_point()
  p <- p + geom_line()
  p <- p + scale_color_discrete(name = "Gene Alteration")
  p <- p + geom_vline(xintercept = 0,color="black",linetype = 2,alpha=0.2)
  p <- p + scale_y_continuous(labels = scales::percent, limits = c(0,1))
  p <- p + ggtitle(title) + ylab("Allele Frequency") + xlab ("Time on Treatment (Days)") 
  p <- p + theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 12))
  
  w <- drawWatermark(g = p, date = date)
  
  ### Append to lists
  
  VAFPlots[[patient]] <- p
  VAFPlotsWatermark[[patient]] <- w
  
}

################################################################################
### Add clean tables to ctDNA object
################################################################################

ctdna@vafsummary <- variants
ctdna@cnvsummary <- cnv
ctdna@monitoring <- monitoring
ctdna@screenfail <- screenfail
ctdna@VAFplots <- VAFPlots
ctdna@VAFplots_wm <- VAFPlotsWatermark

dir.create(datadir, showWarnings = FALSE, recursive = T)
dir.create(mondir, showWarnings = FALSE, recursive = T)
saveRDS(ctdna, paste0(datadir, "/", date, "_", trial, "_unfilteredCtDNA.rds"))
write.xlsx(monitoring, paste0(mondir, "/", date, "_", trial, "_monitoring_required.xlsx"), rowNames = F,overwrite = T)

################################################################################
### Generate vaf markdown
################################################################################

if(rmarkdown_report){
  knit_root_dir <- getwd()
  script_name <- paste(date, trial, "unfiltered", "ctDNA", "Tempus", sep = "_")
  ### knit rmd
  rmarkdown::render(rmd_path, 
                    output_dir = resultsdir,
                    output_file = paste0(resultsdir,
                                         script_name,
                                         '.html'),
                    knit_root_dir = knit_root_dir,
                    params=list(
                      ctdna = paste0(datadir, "/", date, "_", trial, "_unfilteredCtDNA.rds"),
                      date = date,
                      trial = trial,
                      functions = functions,
                      resultsdir =  resultsdir)
  )
  
  ### convert html knit out to aspx for compatibility with Sharepoint
  system(paste("cp", paste0(resultsdir,script_name,'.html'), paste0(resultsdir, script_name,'.aspx')))
}



