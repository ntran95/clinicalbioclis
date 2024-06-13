################################################################################
### Set up
################################################################################

library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(RColorBrewer)
library(openxlsx)
library(readxl)
library(ggrepel)
library(argparser)

################################################################################
### Parse arguments
################################################################################

# === Initialize
p <- arg_parser("process pd for baseline and paired biopsy")
p <- add_argument(p, arg = "--functions",
                  help = "path to functions extension folder (e.g clinicalbioinfoextensions)")
p <- add_argument(p, arg = "--trial",
                  help = "trial number (e.g. RP6306-01")
p <- add_argument(p, arg = "--date",
                  help = "date (format: yyyy-mm-dd)")
p <- add_argument(p, arg = "--baseline_biopsy",
                  help = "path to baseline biopsy data ")
p <- add_argument(p, arg = "--paired_biopsy",
                  help = "path to paired biopsy data")
p <- add_argument(p, arg = "--study",
                  help = "path to processed clinical study object, RDS file")
p <- add_argument(p, arg = "--rmd_path",
                  help = "path to pd rmd script to render report")
p <- add_argument(p, arg = "--resultsdir",
                  help = "results directory to store report")
p <- add_argument(p, arg = "--outdir",
                  help = "output directory to store pd rds")


argv <- parse_args(p)

date <- as.Date(argv$date)
trial <- argv$trial
study_path <- argv$study
outdir <- argv$outdir
functions <- argv$functions
baseline_biopsy <- argv$baseline_biopsy
paired_biopsy <- argv$paired_biopsy
rmd_path <- argv$rmd_path
resultsdir <- argv$resultsdir

devtools::load_all(functions)

################################################################################
### Set testing arguments
################################################################################

interactive <- FALSE
if(interactive){
  
  outdir <- '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/pd/test/data/' 
  functions <-'/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoextensions/' 
  trial <- 'RP3500-01' 
  date <- as.Date('2022-12-30') 
  study_path <- "/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-01/clinical/processing/data/2022-12-20/2022-12-20_RP3500-01_processedStudy.rds"
  baseline_biopsy <- "/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-01/pd/vendor-data/Histowiz15468_ScoreSheet.xlsx"
  paired_biopsy <- "/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-01/pd/vendor-data/Cumulative_RP-3500-01_paired biopsy scoresheet_10NOV2021.xlsx"
  rmd_path <- "/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/pd/src/generate_rmarkdown_report_pd.Rmd"
  resultsdir <- "/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/pd/test/results/"
  
  #source(functions)
  devtools::load_all(functions)
  
}

################################################################################
### Read in clinical and pd data
################################################################################

clinical <- readRDS(study_path)

Clinical <- clinical@combined

Genomic <- clinical@genomics

SNIPDX <- clinical@snipdx

# preserve type of loss from snipdx
Genomic <- merge(Genomic, SNIPDX)

PairedBiopsy <- read_excel(paired_biopsy,sheet = 1)

BaselineBiopsy <- read_excel(baseline_biopsy,sheet = 2) 

################################################################################
### Clean Paired Biopsy
################################################################################

#### Remove Duplicates
PairedBiopsy <- PairedBiopsy[!duplicated(PairedBiopsy[,c(1:14)]),]

#### fix gamma notation
PairedBiopsy[PairedBiopsy$Stain == "p-histone (gammaH2AX)","Stain"] <- "pH2AX"

#### Remove QC Fail
PairedBiopsy <- PairedBiopsy[is.na(PairedBiopsy$`Exclude from dataset for QC reasons`),]

### Change to Baseline and C1D10
PairedBiopsy$VISIT <- gsub("Screening","Pre-Tx",PairedBiopsy$VISIT)
PairedBiopsy$VISIT <- gsub("Cycle 2 Day 10|Unplanned","On-Tx",PairedBiopsy$VISIT)

### Paired Analysis
### If Score is 0 set to 1
PairedBiopsy$`H-Score` <- ifelse(PairedBiopsy$`H-Score`==0, 
                                 1, 
                                 PairedBiopsy$`H-Score`)

PairedBiopsy$`%3+` <- ifelse(PairedBiopsy$`%3+`==0,  
                             1, 
                             PairedBiopsy$`%3+`)


###

PairedBiopsy <- PairedBiopsy %>% group_by(Stain,SUBJID) %>% #group by patient and strain 
  mutate(`H-Score % Change` = (`H-Score`[VISIT == "On-Tx"] /`H-Score`[VISIT == "Pre-Tx"] ) -1) %>% #modify the H-Score % Change entries by dividing the H-score of On-tx by pre-tx of each grouping
  ungroup() 

#why do we change pre-tx H-score % of change to zero?
PairedBiopsy[PairedBiopsy$VISIT == "Pre-Tx","H-Score % Change"] <- 0

# levels for plottings
PairedBiopsy$VISIT <- factor(PairedBiopsy$VISIT,levels = c('Pre-Tx', 'On-Tx'))

################################################################################
### Merge Paired Biopsy
################################################################################

PairedBiopsy<- merge(Clinical,PairedBiopsy,by.x="Patient ID",by.y="SUBJID")

PairedBiopsy <- setGeneGroup(PairedBiopsy, colname = "Enrollment Gene")


################################################################################
### pKAP1 and H2AX strain subsets in paired biopsies
################################################################################

PairedBiopsyStain <- split(PairedBiopsy, PairedBiopsy$Stain)
names(PairedBiopsyStain) <- c("H2AX", "pKAP1")

# PairedBiopsyStainTest <- lapply(PairedBiopsyStain, function(x){
#   x <- dcast(x,`Patient ID`~VISIT+Stain,value.var="H-Score")
#   #
# })

################################################################################
### Paired Biopsy plots
################################################################################
#### paired line plots
paired_h_score_line_p <- lapply(PairedBiopsyStain, paired_h_score_line)

paired_h_score_change_line_p <- lapply(PairedBiopsyStain, paired_h_score_change_line)

#### paired bar plots
for (i in seq_along(PairedBiopsyStain)){
  PairedBiopsyStain[[i]]$VISIT <- factor(PairedBiopsyStain[[i]]$VISIT,
                                         levels = c('On-Tx','Pre-Tx'))
}
levels(PairedBiopsyStain$H2AX$VISIT) #check

paired_stain_barplt_p <- lapply(PairedBiopsyStain,paired_stain_barplt)

################################################################################
### Clean Baseline Biopsy
################################################################################

###Remove tumors with no scores
BaselineBiopsy <- BaselineBiopsy[!is.na(BaselineBiopsy$`Score 1`),] #46 obs
Remove <- c("1004-0041_s")
BaselineBiopsy <- BaselineBiopsy[!(BaselineBiopsy$`Histowiz Sample ID` %in% Remove),] #45 obs
BaselineBiopsy$Source <- "Baseline"

### Pull Baseline From Paired, "screening" = "Pre-tx" = baseline
# this df is used to rbind baseline df with baseline paired entries
PairedBiopsyPreTx <- PairedBiopsy[PairedBiopsy$VISIT == "Pre-Tx" & PairedBiopsy$Stain != "pKAP1",] #22 obs

colnames(PairedBiopsyPreTx) <- gsub("H-Score","H Score",colnames(PairedBiopsyPreTx))
colnames(PairedBiopsyPreTx) <- gsub("Patient ID","Subject ID",colnames(PairedBiopsyPreTx))
PairedBiopsyPreTx$Source <- "Paired"

################################################################################
### Merge Baseline Biopsy
################################################################################
BaselineBiopsyCombined <- rbind(BaselineBiopsy[,c("Subject ID","H Score","Source")],PairedBiopsyPreTx[,c("Subject ID","H Score","Source")]) # 69 obs

# baseline biopsy has patients from Genomics not present in Clinical
BaselineBiopsyCombined<- merge(Genomic,BaselineBiopsyCombined,by.x="Patient ID",by.y="Subject ID")

BaselineBiopsyCombined[is.na(BaselineBiopsyCombined$`Enrollment Gene`),"Enrollment Gene"] <- "Missing"
BaselineBiopsyCombined[is.na(BaselineBiopsyCombined$`Type of loss`),"Type of loss"] <- "missing"

# consolidate gene groups
BaselineBiopsyCombined <- setGeneGroup(BaselineBiopsyCombined,colname = "Enrollment Gene")

BaselineBiopsyCombined$`Type of loss` <- factor(BaselineBiopsyCombined$`Type of loss`)


#### add to pd slots

################################################################################
### Baseline Biopsy plots
################################################################################

baseline_h_score_boxplot_p <- baseline_h_score_boxplot(BaselineBiopsyCombined)

baseline_h_score_facet_sources_p <- baseline_h_score_facet_sources(BaselineBiopsyCombined)

################################################################################
### Visit subsets in baseline biopsies
################################################################################

#split by stain and visit, include only baseline
PairedBiopsyStainVisit <- split(PairedBiopsy, list(PairedBiopsy$Stain, PairedBiopsy$VISIT))

PairedBiopsyStainVisit <- PairedBiopsyStainVisit[c(grep("Pre-Tx", names(PairedBiopsyStainVisit)))]
names(PairedBiopsyStainVisit) <- gsub(".Pre-Tx","",  names(PairedBiopsyStainVisit))

################################################################################
### Baseline H-score facet enrollment gene by stain
################################################################################

stain_bl_barplt_p <- lapply(PairedBiopsyStainVisit, stain_bl_barplt)


################################################################################
### Save to slots
################################################################################

BaselineBiopsyPlots <- list("Baseline Boxplot" = baseline_h_score_boxplot_p,
                            "Baseline Boxplot by Sources" = baseline_h_score_facet_sources_p,
                            "Baseline Stain Barplot" = stain_bl_barplt_p)

PairedBiopsyPlots <- list("Paired Stain Line Plot" = paired_h_score_line_p,
                          "Paired Stain % Change Line Plot" = paired_h_score_change_line_p,
                          "Paired Stain Barplot" = paired_stain_barplt_p)


pd <- new("pd",
          PairedBiopsy = PairedBiopsy,
          BaselineBiopsy = BaselineBiopsy,
          BaselineBiopsyCombined = BaselineBiopsyCombined,
          PairedBiopsyStain = PairedBiopsyStain,
          BaselineBiopsyPlots = BaselineBiopsyPlots,
          PairedBiopsyPlots = PairedBiopsyPlots,
          trial = trial,
          date = date)

dir.create(outdir, showWarnings = FALSE, recursive = T)
saveRDS(pd, paste0(outdir, "/", date, "_", trial, "_processedPD.rds"))


################################################################################
### Knit out report
################################################################################

knit_root_dir <- getwd()
script_name <- paste(date, trial, "pd", "analysis", sep = "_")

### knit rmd
rmarkdown::render(rmd_path, 
                  output_dir = resultsdir,
                  output_file = paste0(resultsdir,
                                       script_name,
                                       '.html'),
                  knit_root_dir = knit_root_dir,
                  params=list(
                    pd = paste0(outdir, "/", date, "_", trial, "_processedPD.rds"),
                    date = date,
                    trial = trial,
                    functions = functions,
                    resultsdir =  resultsdir)
)

### convert html knit out to aspx for compatibility with Sharepoint
system(paste("cp", paste0(resultsdir,script_name,'.html'), paste0(resultsdir, script_name,'.aspx')))

