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

#### ==== Initialize argparse ==== ####
p <- arg_parser("Process local data from SolvBio")
p <- add_argument(p, arg = "--file", 
                  help = "specify the path to raw local data")
p <- add_argument(p, arg = "--output_dir", 
                  help = "specify the path of the output directory, can be relative or absolute paths. Include '/'")
p <- add_argument(p, arg = "--git_extension", 
                  help = "path to functions extension folder (e.g clinicalbioinfoextensions)")
p <- add_argument(p, arg = "--local_reference_path", 
                  help = "specify the path of the consolidation file for local alterations, can be relative or absolute paths.")
p <- add_argument(p, arg = "--snipdx_reference_path", 
                  help = "specify the path of the reference file for snipdx genes, can be relative or absolute paths.")
p <- add_argument(p, arg = "--trial", 
                  help = "specify trial")

argv <- parse_args(p)

#### ==== Interactive Session ==== ####
interactive <- FALSE
if(interactive){
  file <- "/ClinBio/SP-ClinicalBioinformatics/shared/RP6306/RP6306-01/SolveBio/vendor-data/RP-6306-01\ Retrospective\ Local\ Reports\ 15APR2022.xlsx"
  
  output_dir <- "/ClinBio/SP-ClinicalBioinformatics/shared/RP6306/RP6306-01/SolveBio/processing/data/"
  
  git_extension <- "/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoextensions/"
  
  local_reference_path <- "/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/snipdx_local_concordance/ref/SolveBio_Local_Test_Consolidation.xlsx"
  
  snipdx_reference_path <- "/ClinBio/SP-ClinicalBioinformatics/shared/SNiPDx/Reference/panel_genes_annotated_V2.tsv"
  
  local_alts <- read.xlsx(file,
                          sheet = "Variants")
  
  trial = "RP6306-01"
  
  #date = "2022-05-04"
}

git_extension <- argv$git_extension

output_dir <- argv$output_dir

file <- argv$file

local_reference_path <- argv$local_reference_path

snipdx_reference_path <- argv$snipdx_reference_path

trial <- argv$trial

date <- Sys.Date()

################################################################################
### Read in File
################################################################################

devtools::load_all(git_extension)

local_alts <- read_xlsx(file, sheet = "Variants")

local_ref <- read_xlsx(local_reference_path)

snipdx_reference <- read_tsv(snipdx_reference_path)

################################################################################
### Clean and Partition by Types
################################################################################

# label genes 
# local_alts <- local_alts %>%
#   mutate(SNiPDx_Gene = case_when(Gene %in% snipdx_reference$hgnc_symbol ~ "TRUE",
#                                   TRUE ~ "FALSE")) 

# #### ---- convert numeric dates to dates ---- ####
### note: date format is different across trials, omit
# multiple_dates_values <- local_alts$`Collection Date`[grep(";", local_alts$`Collection Date`)]
# 
# multiple_dates_idx <- grep(";", local_alts$`Collection Date`)
# 
# local_alts$`Collection Date` <- as.character(as.Date(as.numeric(local_alts$`Collection Date`),
#                                         origin = "1899-12-30"))
# insert multiple dates back to column
#local_alts[multiple_dates_idx, "Collection Date"] <- multiple_dates_values

# local_alts <- inner_join(local_alts,local_ref)
### adjust splice site to c.
local_alts_cleaned <- local_alts %>%
  ### TODO: may need to do a left join (local, ref) instead to preserve all entries?
  inner_join(local_ref) %>%
  mutate(AAChange = case_when(grepl("splice",`SV Protein Change`)~`SV CDS Change`,
                              is.na(`SV Protein Change`)~`SV CDS Change`, 
                              TRUE~`SV Protein Change`)) %>%
  mutate(short_class = case_when(`SV Change Type` %in% "snv" & 
                                   `Variant Type` %in% "short" ~"SNV",
                                   `Variant Type`%in% "short" & 
                                   `SV Change Type` %notin% "snv"~"InDel")) %>%
  # change `NGS; MLPA` to MLPA only
  mutate(`Assay Type` = case_when(`Assay Type` %in% "NGS; MLPA" ~ "MLPA",
                                  TRUE ~ as.character(`Assay Type`))) %>%
  mutate(`Test Type` = case_when(`Analyte Type` %in% "cfDNA"~"ctDNA",
                                 `Analyte Type` %in% "Protein" ~ "IHC",
                                 `Tissue Type`%in% "Normal"~"germline",
                                 TRUE~"tissue NGS")) %>%
  select(-any_of(c("Row ID",
                   "Report PDF",
                   "Accession Date",
                   "Report Date",
                   "Diagnosis Verbatim",
                   "Specimen Site Verbatim",
                   "Sample Type",
                   "Tissue Type",
                   #"Analyte Type", #preserve to partition by "Protein" types
                   "Variant Assessment Verbatim",
                   "Transcript",
                   # "Coverage",
                   # "Depth",
                   "SB Variant ID",
                   "CNV Verbatim",
                   "RE Verbatim",
                   "IHC Verbatim",
                   "Notes"))) %>%
  mutate(SNiPDx_Gene = case_when(Gene %in% snipdx_reference$hgnc_symbol ~ "TRUE",
                                 TRUE ~ "FALSE")) 

  

#### all separate tables are  filtered on short variants, 
### non-intron variants or pathogenic assessments,
### to do so, filter by "short_class"

#### ==== PARTITION DF BY BY ASSAY TYPES ==== #### 
assay.types <- split(local_alts_cleaned,local_alts_cleaned$`Assay Type`)

# ngs.snv.indel <- assay.types$NGS %>%
#   filter(`SV Change Type` == "indel" | `SV Change Type` == "snv")

ngs.snv.indel <- assay.types$NGS %>%
  filter(short_class == "InDel" | short_class == "SNV")

ngs.snv.indel.test.type <- split(ngs.snv.indel, ngs.snv.indel$`Test Type`)

# names(ngs.snv.indel.test.type) <- paste0("snv.indel.", 
#                                    gsub(" ", ".", names(ngs.snv.indel.test.type)))

ngs.cnv <- assay.types$NGS %>%
  filter(`Variant Type` %in% "CNV")

ngs.cnv.test.type <- split(ngs.cnv, ngs.cnv$`Test Type`)

# names(ngs.cnv.test.type) <- paste0("cnv.", 
#                                    gsub(" ", ".", names(ngs.cnv.test.type)))

ngs.re <- assay.types$NGS %>%
  filter(`RE Type` %in% "rearrangement")

ngs.re.test.type <-  split(ngs.re, ngs.re$`Test Type`)

# names(ngs.re.test.type) <-  paste0("re.", 
#        gsub(" ", ".", names(ngs.re.test.type)))

if(!is.null(assay.types[["Array"]])){
  array.cnv <- assay.types$Array %>%
    filter(`Variant Type` %in% "CNV")   
}else{
  array.cnv <- NULL
}


if(!is.null(array.cnv)){
  array.cnv.test.type <- split(array.cnv, array.cnv$`Test Type`)
}else{
  array.cnv.test.type <- NULL
}

# names(array.cnv.test.type) <- paste0("cnv.", 
#                                    gsub(" ", ".", names(array.cnv.test.type)))

if(!is.null(assay.types[["MLPA"]])){
   mlpa.re <- assay.types$MLPA %>%
    filter(`RE Type` %in% "rearrangement")
   
   mlpa.re.test.type <- split(mlpa.re, mlpa.re$`Test Type`)
   
}else{
  mlpa.re <- NULL
  
  mlpa.re.test.type <- NULL
}



# sanger.sequencing.snv.indel <- assay.types$`Sanger sequencing` %>%
#   filter(`SV Change Type` == "indel" | `SV Change Type` == "snv") 
if(!is.null(assay.types[["Sanger sequencing"]])){
  sanger.sequencing.snv.indel <- assay.types$`Sanger sequencing` %>%
    filter(short_class == "InDel" | short_class == "SNV")
  
  sanger.sequencing.snv.indel.test.type <- 
    split(sanger.sequencing.snv.indel,
          sanger.sequencing.snv.indel$`Test Type`)
}else{
  sanger.sequencing.snv.indel <- NULL
  
  sanger.sequencing.snv.indel.test.type <- NULL
}



# names(sanger.sequencing.snv.indel.test.type) <- paste0("snv.indel.", 
#                                      gsub(" ", ".", names(sanger.sequencing.snv.indel.test.type)))

#### ==== STORE ==== ####
ngs.list <- list(ngs = assay.types$NGS,
                 snv.indel = c(list(snv.indel = ngs.snv.indel),
                                      ngs.snv.indel.test.type),
                 cnv = c(list(cnv = ngs.cnv),
                             ngs.cnv.test.type),
                 re = c(list(re = ngs.re),
                            ngs.re.test.type)
                 )

array.list <- list(array = assay.types$Array,
                   cnv = c(list(cnv = array.cnv),
                                array.cnv.test.type)
                   )

mlpa.list <- list(mlpa = assay.types$MLPA,
                  re = c(list(re = mlpa.re),
                         mlpa.re.test.type)
                  )

sanger.list <- list(sanger = assay.types$`Sanger sequencing`,
                    snv.indel = c(list(snv.indel = sanger.sequencing.snv.indel),
                                  sanger.sequencing.snv.indel.test.type)
                    )

################################################################################
### Initialize S4 object
################################################################################

## move to clinicalbioinfoextensions repo
# setClass("local", representation(orig_local_alts = "data.frame",
#                                  local_alts_cleaned = "data.frame",
#                                  ngs = "list",
#                                  ihc = "list",
#                                  fish = "list",
#                                  array = "list",
#                                  mlpa = "list",
#                                  sanger = "list",
#                                  trial = "character",
#                                  date = "Date"))
 
local <- new("local",
             orig_local_alts = local_alts,
             local_alts_cleaned = local_alts_cleaned,
             ngs = ngs.list,
             ihc = list(ihc = assay.types$IHC),
             fish = list(fish = assay.types$FISH),
             array = array.list,
             mlpa = mlpa.list,
             sanger = sanger.list,
             trial = trial,
             date = as.Date(date))
           


################################################################################
### Export
################################################################################
outFolder <- paste0(output_dir,"/","local-processed-data", "/", local@date, "/")
dir.create(outFolder, recursive = T,showWarnings = F)

saveRDS(local, 
     file = paste0(outFolder, 
                   local@date, "_", 
                   local@trial, "_",
                   "processed_local_data.rds"))
