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
library(ggpubr)
library(survminer)
library(survival)
library(data.table)

################################################################################
### Parse arguments
################################################################################

# === Initialize
p <- arg_parser("Calculate mVAFs and generate plots")
p <- add_argument(p, arg = "--ctdna",
                  help = "ctdna object path")
p <- add_argument(p, arg = "--monitoring",
                  help = "annotated monitoring file path")
p <- add_argument(p, arg = "--outdir",
                  help = "output directory")
p <- add_argument(p, arg = "--functions",
                  help = "path to functions extension folder (e.g clinicalbioinfoextensions)")
p <- add_argument(p, arg = "--cnvgene",
                  help = "gene of interest for CNV analysis, separate multiple genes by comma (i.e BRCA1,BRCA2)")
# p <- add_argument(p, arg = "--response",
#                   help = "logical to indicate if response data is available, set to T or F")
p <- add_argument(p, arg = "--trial",
                  help = "trial number (e.g. RP6306-01")
p <- add_argument(p, arg = "--date",
                  help = "date (format: yyyy-mm-dd)")


argv <- parse_args(p)


date <- as.Date(argv$date)
trial <- argv$trial
ctdna_path <- argv$ctdna
monitoring_path <- argv$monitoring
outdir <- argv$outdir
cgoi <- argv$cnvgene
functions <- argv$functions

cgoi <- unlist(strsplit(cgoi,split = ","))

devtools::load_all(functions)

################################################################################
### Set testing arguments
################################################################################
# 
# date <- ymd("2022-03-05")
# trial <- "RP6306-01"
# ctdna_path <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-01/ctDNA/processing/Tempus/data/2022-03-05_RP6306-01_unfilteredCtDNA.rds"
# monitoring_path <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-01/ctDNA/processing/Tempus/monitor/2022-03-05_RP6306-01_monitoring_completed.xlsx"
# outdir <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-01/ctDNA/processing/Tempus/data/"
# cgoi <- "CCNE1"
# response <- F

################################################################################
### Set testing arguments RP3500-03
################################################################################
interactive <- FALSE
data_location <- "local"
if(interactive){
  
  if(data_location == "rserver"){
    date <- ymd("2022-12-06")
    trial <- "RP3500-03"
    ctdna_path <- "/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/ctDNA/processing/Tempus/data/2022-12-06_RP3500-03_unfilteredCtDNA.rds"
    monitoring_path <- "/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/ctDNA/processing/Tempus/data/monitoring/2022-11-22_RP3500-03_monitoring_completed.xlsx"
    outdir <- "/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/data/"
    cgoi <- c("BRCA1","BRCA2")
    functions <-  '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoextensions/'
    
    #source(functions)
    devtools::load_all(functions)
  }
  if(data_location == "local"){
    date <- ymd("2022-11-07")
    trial <- "RP3500-03"
    ctdna_path <- "./../test/2022-11-07_RP3500-03_unfilteredCtDNA.rds"
    monitoring_path <- "../test/data/monitoring/2022-11-07_RP3500-03_monitoring_completed.xlsx"
    outdir <- "../test/data/"
    cgoi <- c("BRCA1","BRCA2")
    functions <-  '../../../clinicalbioinfoextensions/'
    
    #source(functions)
    devtools::load_all(functions)
  }
  
}


################################################################################
### Read ctDNA object and monitoring files
################################################################################

ctdna <- readRDS(ctdna_path)
vafsummary <- ctdna@vafsummary
cnvsummary <- ctdna@cnvsummary
clinical <- ctdna@metadata

monitoring <- readxl::read_excel(monitoring_path)
ctdna@monitoring <- monitoring

################################################################################
### Calculate mean coverage for each sample, and mean coverage for each alteration across visits - moved to calculate_VAFs
################################################################################

# # Calculate alteration min coverage across visits

altCov <- vafsummary %>% select(`Patient ID`,`TimeFromC1D1`, `gene`, `alteration`, coverage) %>%
  group_by(`Patient ID`, `gene`, `alteration`) %>%
  summarize(minAltCov = min(coverage, na.rm = T))

vafsummary <- merge(vafsummary, altCov, all.x = T)

# Calculate mean sample coverage

samcov <- vafsummary %>% select(`Patient ID`,`TimeFromC1D1`, `gene`, `alteration`, coverage) %>% unique() %>% group_by(`Patient ID`,`TimeFromC1D1`) %>% summarize( meanSampleCov = mean(coverage, na.rm = T))
vafsummary <- merge(vafsummary, samcov, all.x = T)


################################################################################
### Calculate max VAF for each alteration and max VAF for each patient
################################################################################

vafsummary <- merge(vafsummary, select(monitoring, any_of(monitoringIn_cols)), all = TRUE)

# Calculate alteration max VAF across visits
altvaf <- vafsummary %>% select(`Patient ID`,`TimeFromC1D1`, `gene`, `alteration`, VAF) %>%
                         group_by(`Patient ID`, `gene`, `alteration`) %>%
                         summarize(maxAltVAF = max(VAF, na.rm = T))

vafsummary <- merge(vafsummary, altvaf, all.x = T)

# Calculate patient max VAF across variants

ptvaf <- vafsummary %>% 
  filter(Status == "somatic") %>% 
  select(`Patient ID`, VAF) %>% unique() %>% 
  group_by(`Patient ID`) %>%
  summarize( maxPatientVAF = max(VAF, na.rm = T))

vafsummary <- merge(vafsummary, ptvaf, all.x = T)

################################################################################
### Flag variants to omit from further analysis and filter
################################################################################
orig_vafsummary <- vafsummary

# Flag variants 

vafsummary$Inclusion <- ifelse(vafsummary$somatic_germline %in% "S" &
                               vafsummary$Status %in% "somatic" &
                               vafsummary$maxAltVAF > 0.5 & # not implemented in Rp3500, Use this instead
                               # vafsummary$VAF > 0.5 & # currently implemented in RP3500
                               vafsummary$maxPatientVAF > 1 & # max Patient VAF is not filtered out until responses, confirmed 
                               vafsummary$minAltCov >= 500 & # this is 
                               vafsummary$meanSampleCov > 750, 
                              "include", "exclude")

table(vafsummary$Inclusion)
# Save omitted variants
omitted <- filter(vafsummary, Inclusion == "exclude")
omitted <- omitted %>% mutate( Reason = case_when(Status != "somatic" | somatic_germline != "S" ~ "not somatic",
                                                  VAF < 0.5 ~ "somatic variant <0.5%",
                                                  maxPatientVAF <1 ~ "patient vaf < 1",
                                                  maxAltVAF < 0.5 ~ "variant vaf < 0.5",
                                                  meanSampleCov < 750 ~ "sample coverage < 750",
                                                  minAltCov <= 500 ~ "variant coverage <= 500")) 
###TODO: relocate omitted columns, add columns to functions

# preserve unfiltered timepoints and corresponding sample level estimates
unique_visits <- vafsummary %>% 
   distinct(`Patient ID`, TimeFromC1D1, VisitID, partner_sample_id,`Tempus Sample ID`,ctFE, cfDNA_ng.ml) 


# Filter somatic only
vafsummary <- filter(vafsummary, Inclusion %in% "include")

################################################################################
### Calculate patient-level mean VAF (mVAF) - to edit 
################################################################################

# mean VAF with filling 0.1 for missing timepoints

# distinct filtered alts
unique_alts <- vafsummary %>% 
  select(any_of(c(alteration_summary_cols))) %>% 
  unique()

# report all somatic patients and total somatic num alts per patient
unique_patients <- vafsummary %>%
    group_by(`Patient ID`) %>%
    distinct(alteration) %>%
    mutate(SomaticNumAltsPerVisit = n()) %>%
    dplyr::select(-alteration) %>%
    unique()

# report all distinct somatic alts
alts_vafs <- vafsummary %>%
  select(`Patient ID`, TimeFromC1D1, VisitID, partner_sample_id,`Tempus Sample ID`,ctFE, cfDNA_ng.ml, gene, alteration, VAF, Inclusion) %>% 
  # # select all vaf summary columns for cleaner merge downstream
  # select(any_of(c(vafsummary_cols))) %>%
  unique()

# pivot wider and merge with all samples analyzed
alts_vafs <- dcast(setDT(alts_vafs), ... ~ gene + alteration + Inclusion, value.var = "VAF",sep = "|")

# merge alts with samples analyzed, preserve all timepoints
visits_alts <- alts_vafs %>% right_join(unique_visits) 

# merge alts with patient level df, perserve distinct patients with total somatic num alts per patient
visits_alts <- visits_alts %>% left_join(unique_patients)

# calculate VAFsum
# visits_alts$VAFsum <- rowSums(visits_alts[,grep("include",colnames(visits_alts),value = T)[1]],na.rm = T)

# pivot back to longer
visits_alts <- melt(setDT(visits_alts),
     id.vars=c("Patient ID",
               "TimeFromC1D1",
               "VisitID",
               "partner_sample_id",
               "Tempus Sample ID",
               "cfDNA_ng.ml",
               "ctFE",
               "SomaticNumAltsPerVisit"),
      variable.name = "Alteration",
      value.name = "VAF")

# separate alterations back into gene and mutation
visits_alts <- separate(visits_alts,Alteration,
                       into = c ("gene","alteration","Inclusion"),
                       sep = "\\|",remove = T)

visits_alts <- visits_alts %>% inner_join(unique_alts)


# add pseudoVAFs of 0.1 to alterations with missing timepoints
visits_alts$VAF <- ifelse( is.na(visits_alts$VAF), 0.1, visits_alts$VAF)

# calculate mVAF
mvaf <- visits_alts %>% group_by(`Patient ID`, TimeFromC1D1) %>% mutate(mVAF = mean(VAF))

# merge new vafsummary with clinical cols
vafsummary <-inner_join(mvaf, clinical) %>%
  ungroup()


################################################################################
### Calculate mVAF ratio and mVAF difference from baseline
################################################################################

baseline_mvaf <- vafsummary %>% filter(VisitID == "Baseline") %>% select(`Patient ID`, mVAF) %>% unique() %>% rename(mVAF_Baseline = mVAF)
other_mvaf <- vafsummary %>% select(`Patient ID`, VisitID, TimeFromC1D1, mVAF) %>%  unique()

mvafrd <- merge(baseline_mvaf, other_mvaf, by = "Patient ID")
mvafrd <- mvafrd %>% mutate(mVAFR = round(mVAF/mVAF_Baseline, 2) - 1)
## Set maximum to 100
mvafrd[which(mvafrd$mVAFR > 1),"mVAFR"] <- 1
mvafrd <- mvafrd %>% mutate(mVAFD = round(mVAF - mVAF_Baseline, 2))
mvafrd <- mvafrd %>% select(`Patient ID`, VisitID, TimeFromC1D1, mVAF, mVAFR, mVAFD)
# only evaluate on-tx timepoints for best mVAFR, lowest mVAFR,if patient has 2 best mVAFR @ same cycle, choose lowest by TimeFromC1D1
# exclude timepoints past 80 days from C1D1, exclude baseline 
# mvafrd <- mvafrd %>%
#   group_by(`Patient ID`) %>%
#   mutate(Best_Visit = case_when(mVAFR == min(mVAFR[TimeFromC1D1 <80 & VisitID %notin% c("Baseline", "Pre-Tx", "EOT")]) ~ "Yes",
#                                   TRUE ~ "No"))

### label best mvafr
# only evaluate on-tx timepoints for best mVAFR, lowest mVAFR,if patient has 2 best mVAFR @ same cycle, choose lowest by TimeFromC1D1
# exclude timepoints past 80 days from C1D1, exclude baseline
best_mvaf <- mvafrd %>%
  # arrange in order of low->high patient, mvafr, then timepoint
  arrange(`Patient ID`, mVAFR,TimeFromC1D1) %>%
  group_by(`Patient ID`) %>%
  filter(TimeFromC1D1 <80,
         VisitID %notin% c("Baseline", "Pre-Tx", "EOT")) %>% 
  # choose the best (lowest) mvafr at earliest timepoint
  slice(1) %>%
  mutate(Best_Visit = "Yes") %>% 
  ungroup()

# merge back with mvafrd
mvafrd <- best_mvaf %>% right_join(mvafrd) %>%
  mutate(Best_Visit = case_when(Best_Visit %notin% "Yes" ~ "No",
                                TRUE ~ as.character(Best_Visit))) %>%
  arrange(`Patient ID`,TimeFromC1D1, mVAFR)


vafsummary <- merge(vafsummary, mvafrd, all.x = T)

################################################################################
### Build patient, sample, and alteration level summaries
################################################################################

patient_summary <- vafsummary %>% select(`Patient ID`, VisitID, mVAF, Inclusion) %>% 
                          filter(VisitID == "Baseline") %>% 
                          unique() %>%
                          rename( `Baseline mVAF` = mVAF ) %>%
                          select( - VisitID)

# sample_summary <- vafsummary %>% select(`Patient ID`, VisitID, TimeFromC1D1, mVAF, Inclusion) %>% 
#                          unique()

sample_summary <- vafsummary %>%
  select(any_of(c(sample_summary_cols))) %>%
  distinct()

# alteration_summary <- vafsummary %>% select(`Patient ID`, gene, alteration, Status, Inclusion) %>% unique()

alteration_summary <-  vafsummary %>%
  select(any_of(c(alteration_summary_cols))) %>%
  distinct()

best_mvaf <- vafsummary %>% 
  filter(Best_Visit == "Yes") %>%
  dplyr::select(any_of(c(mvaf_summary_cols, clinical_cols))) %>%
  unique()

################################################################################
### Define molecular response and increase/decrease groups
################################################################################

vafsummary <- vafsummary %>% mutate( `Molecular Response` = case_when( mVAFR > 0.5 ~ "progression",
                                                                       mVAFR <= -0.5 ~ "response",
                                                                       !is.na(mVAFR) ~ "stable",
                                                                       TRUE ~ "NA"))
vafsummary <- vafsummary %>% mutate( `ctDNA trend` = case_when(mVAFD >= 0 ~ "increasing",
                                                               mVAFD < 0 ~ "decreasing"))                                     

################################################################################
### Plot patient-wise VAF and mVAF line plots 
################################################################################

# VAF plots

VAFPlots <- list()
VAFPlotsWatermark <- list()

for(patient in unique(vafsummary$`Patient ID`)) {
  
  # plotdf <- filter(vafsummary, `Patient ID` == patient)
  # plotdf$ordalteration <- paste0(plotdf$gene, " ", plotdf$alteration)
  # plotdf <- plotdf %>% arrange(-mAVAF)
  
  print(patient)
  
  ### Filter unique alterations for given patient
  
  plotdf <- filter(vafsummary, `Patient ID` == patient)
  
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
  
  #title <- unique(paste0(patient, " ", plotdf$`Enrollment Gene`, " ", plotdf$`Tissue Type`))
  
  p <- ggplot(plotdf, aes(x=TimeFromC1D1, y=VAF/100, color=ordalteration))
  p <- p + geom_point()
  p <- p + geom_line()
  p <- p + ggtitle(title)
  p <- p + ylab("Allele Frequency")
  p <- p + xlab ("Time on Treatment (Days)")
  p <- p + geom_vline(xintercept = 0,color="black",linetype = 2,alpha=0.2)
  p <- p + scale_y_continuous(labels = scales::percent, limits = c(0,1.25 * max(plotdf$VAF/100)))
  p <- p + theme_minimal()
  p <- p + labs(color="Gene Alteration")
  p <- p + theme(legend.position = "right", plot.title = element_text(hjust = 0.5,size = 18))
  p <- p + theme(text = element_text(size = 14))
  
  VAFPlots[[patient]] <- p
  VAFPlotsWatermark[[patient]] <- drawWatermark(g = p, date = date)
  
}

# mVAF plots

mVAFPlots <- list()
mVAFPlotsWatermark <- list()

for(patient in unique(vafsummary$`Patient ID`)) {
  
  plotdf <- filter(vafsummary, `Patient ID` == patient)

  title <- unique(paste(plotdf$`Enrollment Gene`, patient , plotdf$`Tissue Type`))
  
  ### if pre-screen sample or if patient is not entered in clinical 
  clinical_avail <- unique(plotdf$ID.gene)
  if(is.na(clinical_avail)){
    title <- paste(patient, unique(plotdf$partner_sample_id))
  }
  
  p <- ggplot(plotdf, aes(x=TimeFromC1D1, y=mVAF/100))
  p <- p + geom_point()
  p <- p + geom_line()
  p <- p + ggtitle(title)
  p <- p + ylab("Mean Allele Frequency")
  p <- p + xlab ("Time on Treatment (Days)")
  p <- p + geom_vline(xintercept = 0,color="black",linetype = 2,alpha=0.2)
  p <- p + scale_y_continuous(labels = scales::percent, limits = c(0,1.25 * max(plotdf$VAF/100)))
  p <- p + theme_minimal()
  p <- p + theme(legend.position = "right", plot.title = element_text(hjust = 0.5,size = 18))
  p <- p + theme(text = element_text(size = 14))
  
  mVAFPlots[[patient]] <- p
  mVAFPlotsWatermark[[patient]] <- drawWatermark(g = p, date = date)
  
}

################################################################################
### Filter for CNVs that meet the minimum log2 threshold 
################################################################################

cnvlog2 <- cnvsummary %>% select(`Patient ID`, TimeFromC1D1, gene, `log2`) %>%
                          group_by(`Patient ID`, gene) %>%
                          summarize(maxCNVlog2 = max(log2, na.rm=T))

cnvsummary <- merge(cnvsummary, cnvlog2, all.x = T)

################################################################################
### Flag CNVs to omit from further analysis
################################################################################

# Flag variants 

cnvsummary$Inclusion <- ifelse(cnvsummary$maxCNVlog2 > 0.5, "include", "exclude")

################################################################################
### Plot CNV line plots
################################################################################

cgoi <- cnvsummary %>% filter(gene %in% cgoi & Inclusion == "include")

CNVPlots <- list()
CNVPlotsWatermark <- list()

for(patient in unique(cgoi$`Patient ID`)) {
  
  plotdf <- filter(cgoi, `Patient ID` == patient)

  title <- unique(paste0(patient, " ", plotdf$`Enrollment Gene`, " ", plotdf$`Tissue Type Consolidated`))
  p <- ggplot(plotdf, aes(x=TimeFromC1D1, y= log2, color = gene))
  p <- p + geom_point()
  p <- p + geom_line()
  p <- p + ggtitle(title)
  p <- p + ylab("Log Ratio")
  p <- p + xlab ("Time on Treatment (Days)")
  p <- p + geom_vline(xintercept = 0,color="black",linetype = 2,alpha=0.2)
  #p <- p + scale_y_continuous(labels = scales::percent, limits = c(0,1.25 * max(plotdf$VAF/100)))
  p <- p + theme_minimal()
  p <- p + labs(color="Gene")
  p <- p + theme(legend.position = "right", plot.title = element_text(hjust = 0.5,size = 24))
  p <- p + theme(text = element_text(size = 20))
  
  CNVPlots[[patient]] <- p
  CNVPlotsWatermark[[patient]] <- drawWatermark(g = p, date = date)
  
}

################################################################################
### Save updated tables and plots to objects
################################################################################

ctdna@patient_summary <- patient_summary
ctdna@alteration_summary <- alteration_summary
ctdna@sample_summary <- sample_summary
ctdna@addtl_summaries[["mVAFR"]] <- mvafrd
ctdna@best_mvaf <- best_mvaf
ctdna@omitted <- omitted
ctdna@vafsummary <- vafsummary
ctdna@cnvsummary <- cnvsummary
ctdna@VAFplots <- VAFPlots
ctdna@VAFplots_wm <- VAFPlotsWatermark
ctdna@CNVplots <- CNVPlots
ctdna@CNVplots_wm <- CNVPlotsWatermark
ctdna@mVAFplots <- mVAFPlots
ctdna@mVAFplots_wm <- mVAFPlotsWatermark


saveRDS(ctdna, paste0(outdir, "/", date, "_", trial, "_filteredCtDNA.rds"))





