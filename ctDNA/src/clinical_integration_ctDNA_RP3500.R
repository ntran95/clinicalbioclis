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
library(patchwork)
library(ggrepel)

################################################################################
### Parse arguments
################################################################################

# === Initialize
p <- arg_parser("generate trial specific response plots")
p <- add_argument(p, arg = "--ctdna",
                  help = "filtered ctdna object path, ")
p <- add_argument(p, arg = "--datadir",
                  help = "output directory")
p <- add_argument(p, arg = "--functions",
                  help = "path to functions extension folder (e.g clinicalbioinfoextensions) ")
p <- add_argument(p, arg = "--trial",
                  help = "trial number (e.g. RP6306-01")
p <- add_argument(p, arg = "--date",
                  help = "date (format: yyyy-mm-dd)")
p <- add_argument(p, arg = "--rmarkdown_report",
                  help = "logical indicating whether to generate filtered markdown reports with MR reponse plots, set to 'T' or 'F'")
p <- add_argument(p, arg = "--resultsdir",
                  help = "results directory, applicable when rmarkdown_report is TRUE")
p <- add_argument(p, arg = "--rmd_path",
                  help = "path to filtered rmd script, applicable when rmarkdown_report is TRUE")


argv <- parse_args(p)


date <- as.Date(argv$date)
trial <- argv$trial
ctdna_path <- argv$ctdna
datadir <- argv$datadir
functions <- argv$functions
rmarkdown_report <- as.logical(argv$rmarkdown_report)
if(rmarkdown_report){
  resultsdir <- argv$resultsdir
  rmd_path <- argv$rmd_path
}

devtools::load_all(functions)


################################################################################
### Set testing arguments RP3500-03
################################################################################
interactive <- FALSE
if(interactive){
  date <- ymd("2022-12-02")
  trial <- "RP3500-03"
  ctdna_path <- "/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/ctDNA/processing/Combined/data/2022-12-02_RP3500-03_filteredCtDNA.rds"
  datadir <- "/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/data/"
  functions <-  '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoextensions/'
  rmarkdown_report <- T
  resultsdir <- "/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/ctDNA/processing/Combined/results/"
  rmd_path <- "/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/src/generate_rmarkdown_report_filtered_ctDNA.Rmd"
  
  devtools::load_all(functions)
  
  
}


################################################################################
### Read ctDNA object and monitoring files
################################################################################

ctdna <- readRDS(ctdna_path)
vafsummary <- ctdna@vafsummary

################################################################################
### plotting utils
################################################################################
# consolidate enrollment gene groups
vafsummary <- setGeneGroup(df = vafsummary, 
                           colname = "Enrollment Gene", 
                           nonStep2Genes = "BRCA1|BRCA2|ATM|CDK12")

clinical_cols <- c(colnames(ctdna@metadata), "GeneGroup", "GeneGroup2")

# specify clinical levels for shape legend
scaleShapeBreaks <- c("biallelic","suspect biallelic",
                      "monoallelic","CHIP",
                      "germline","somatic",
                      "TP53 Mut","TP53 WT","Missing")

################################################################################
### plotting utils
################################################################################
best_mvaf <- vafsummary %>% 
  filter(Best_Visit == "Yes") %>%
  dplyr::select(any_of(c(mvaf_summary_cols, clinical_cols))) %>%
  unique()
  

################################################################################
#### compile all response plot to list ####
################################################################################

generate_response_plots <- function(best_mvaf, scaleShapeBreaks, date){
  ### plot_best_mVAFR ###
  plot_best_mVAFR_p <- plot_best_mVAFR(best_mvaf, gene_col = "GeneGroup") 
  
  ### add RP3500 specific shape legend 
  plot_best_mVAFR_p <- plot_best_mVAFR_p + 
    geom_point(aes(y=-1.3,shape=`TP53 Status`),size=3) +
    geom_point(aes(y=-1.1,shape=`Type of loss`),size=3) +
    geom_point(aes(y=-1.2,shape=`germline status`),size=3) +
    geom_point(aes(y=-1.3,shape=`TP53 Status`),size=3) +
    scale_shape_manual(breaks=scaleShapeBreaks,
                       values=scaleShapeBreaksID[scaleShapeBreaks %in%
                                                   names(scaleShapeBreaksID)])
  
  ### facet by gene group
  plot_best_mVAFR_genegroup <-plot_best_mVAFR_p + 
    facet_grid(.~GeneGroup2,scales="free",space="free",switch="y")
  
  ### simple, no COR labels or shape legend
  plot_best_mVAFR_simple_p <- plot_best_mVAFR_simple(best_mvaf,gene_col = "GeneGroup")
  
  ### simple, no COR labels or shape legend, faceted by gene groups
  plot_best_mVAFR_simple_genegroup <- plot_best_mVAFR_simple_p +
    facet_grid(.~GeneGroup2,scales="free",space="free",switch="y")
  
  ### ENA simple, no COR labels or shape legend, replace Patient No with Tissue Type on y axis ###
  plot_best_mVAFR_simple_ENA <- plot_best_mVAFR_simple_p +
    scale_x_discrete(breaks =plot_best_mVAFR_simple_p$data$`Patient No`,
                     labels = plot_best_mVAFR_simple_p$data$`Tissue Type` )
  
  ### ENA simple, no COR labels or shape legend, replace Patient No with Tissue Type on y axis, faceted by gene group
  plot_best_mVAFR_simple_ENA_genegroup <- plot_best_mVAFR_simple_ENA +
    facet_grid(.~GeneGroup2,scales="free",space="free",switch="y")
  

  ### plot_mVAFR_vs_TumorVol
  plot_mVAFR_vs_TumorVol_p <- plot_mVAFR_vs_TumorVol(best_mvaf, gene_col = "GeneGroup")
  
  # regression lines grouped by consolidated genes
  plot_mVAFR_vs_TumorVol_byEG_p <- plot_mVAFR_vs_TumorVol_byEG(best_mvaf, gene_col = "GeneGroup2")
  

  ### plot_mVAFR_vs_BOR: Box plot comparing mVAFR across different BOR groups
  plot_mVAFR_vs_BOR_p <- plot_mVAFR_vs_BOR(best_mvaf, gene_col = "GeneGroup")
  

  ### Kaplan-Meier curve grouped by molecular response and ctdna trend ###
  plot_MR_onTx_KMP_p <- plot_MR_onTx_KM(best_mvaf)
  
  plot_trend_onTx_KM_p <- plot_trend_onTx_KM(best_mvaf)
  
  responsePlots <- list("Response Waterfall" = plot_best_mVAFR_p,
                        "Response Waterfall by Gene" = plot_best_mVAFR_genegroup,
                        "Response Waterfall Simple" = plot_best_mVAFR_simple_p,
                        "Response Waterfall Simple by Gene" = plot_best_mVAFR_simple_genegroup,
                        "Response Waterfall ENA" = plot_best_mVAFR_simple_ENA,
                        "Response Waterfall ENA by Gene" = plot_best_mVAFR_simple_ENA_genegroup,
                        "Response vs TumorVol" = plot_mVAFR_vs_TumorVol_p,
                        "Response vs TumorVol by Gene" = plot_mVAFR_vs_TumorVol_byEG_p,
                        "Response vs BOR" = plot_mVAFR_vs_BOR_p,
                        "KM MR on-Tx" = plot_MR_onTx_KMP_p,
                        "KM ctDNA Trend on-Tx" = plot_trend_onTx_KM_p)
  
  ### Watermark plots ###
  
  w1 <- drawWatermark(plot_best_mVAFR_p, date)
  w2 <- drawWatermark(plot_best_mVAFR_genegroup, date)
  w3 <- drawWatermark(plot_best_mVAFR_simple_p, date)
  w4 <- drawWatermark(plot_best_mVAFR_simple_genegroup, date)
  w4 <- drawWatermark(plot_best_mVAFR_simple_ENA, date)
  w5 <- drawWatermark(plot_best_mVAFR_simple_ENA_genegroup, date)
  w6 <- drawWatermark(plot_mVAFR_vs_TumorVol_p, date)
  w7 <- drawWatermark(plot_mVAFR_vs_TumorVol_byEG_p, date)
  w8 <- drawWatermark(plot_mVAFR_vs_BOR_p, date)
  
  # plot grid specific to survival plot objs
  w9 <-  cowplot::plot_grid(plot_MR_onTx_KMP_p$plot,
                            plot_MR_onTx_KMP_p$table,
                            nrow = 2, rel_heights = c(2,1))
  w9 <- drawWatermark(w9, date = date)
  
  w10 <-  cowplot::plot_grid(plot_trend_onTx_KM_p$plot,
                             plot_trend_onTx_KM_p$table,
                            nrow = 2, rel_heights = c(2,1))
  
  w10 <- drawWatermark(w10, date = date)
  
  responsePlotsWatermark <- list("Response Waterfall" = w1,
                                 "Response Waterfall by Gene" = w2,
                                 "Response Waterfall Simple" = w3,
                                 "Response Waterfall Simple by Gene" = w4,
                                 "Response Waterfall ENA" = plot_best_mVAFR_simple_ENA,
                                 "Response Waterfall ENA by Gene" = w5,
                                 "Response vs TumorVol" = w6,
                                 "Response vs TumorVol by Gene" = w7,
                                 "Response vs BOR" = w8,
                                 "KM MR on-Tx" = w9,
                                 "KM ctDNA Trend on-Tx" = w10)
  
  return(list(responsePlots, responsePlotsWatermark))
}

################################################################################
### Call plotting functions for different time windows (needs testing) - move to trial specific responses trial
################################################################################

responsePlots <- list()
responsePlotsWatermark <- list()

### All

responsePlots <- generate_response_plots(best_mvaf = best_mvaf,
                                         scaleShapeBreaks =scaleShapeBreaks,
                                         date = date )[[1]]

responsePlotsWatermark <- generate_response_plots(best_mvaf = best_mvaf,
                                                  scaleShapeBreaks = scaleShapeBreaks,
                                                  date = date)[[2]]

### TODO:
# ### Cycles 1-2
# 
# responsePlots <- generate_response_plots(vafsummary, "cycles", min = 1, max = 3, "Cycles1to3", responsePlots, responsePlotsWatermark)[[1]]
# responsePlotsWatermark <- generate_response_plots(vafsummary, "cycles", min = 1, max = 3, "Cycles1to3", responsePlots, responsePlotsWatermark)[[2]]
# 
# 
# ### Weeks 1-5
# 
# responsePlots <- generate_response_plots(vafsummary, "weeks", min = 1, max = 5, "Weeks1to5", responsePlots, responsePlotsWatermark)[[1]]
# responsePlotsWatermark <- generate_response_plots(vafsummary, "weeks", min = 1, max = 5, "Weeks1to5", responsePlots, responsePlotsWatermark)[[2]]


################################################################################
### Save updated tables and plots to objects
################################################################################


ctdna@best_mvaf <- best_mvaf
ctdna@responsePlots <- responsePlots
ctdna@responsePlots_wm <- responsePlotsWatermark

saveRDS(ctdna, paste0(datadir, "/", date, "_", trial, "_filtered_mVAFR_responses_CtDNA.rds"))


if(rmarkdown_report){
  knit_root_dir <- getwd()
  script_name <- paste(date, trial, "filtered", "ctDNA", sep = "_")
  
  ### knit rmd
  rmarkdown::render(rmd_path, 
                    output_dir = resultsdir,
                    output_file = paste0(resultsdir,
                                         script_name,
                                         '.html'),
                    knit_root_dir = knit_root_dir,
                    params=list(
                      ctdna = paste0(datadir, "/", date, "_", trial, "_filtered_mVAFR_responses_CtDNA.rds"),
                      date = date,
                      trial = trial,
                      functions = functions,
                      resultsdir =  resultsdir,
                      response = "T")
  )
  
  ### convert html knit out to aspx for compatibility with Sharepoint
  system(paste("cp", paste0(resultsdir,script_name,'.html'), paste0(resultsdir, script_name,'.aspx')))
}


