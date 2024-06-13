# SNiPDx and Local Concordance

This pipeline processes the local data and merges with the process snipdx object in order to quantify the concordance of short variants and CNVs between the two platforms.

______________________

# Table of Contents

- [Workflow overview](#workflow-overview)
- [Prerequisites](#prerequisites)
  * [Box credentials](#box-credentials)
  * [SolveBio Local test consolidation](#solvebio-local-test-consolidation)
  * [SNiPDx gene panel](#snipdx-gene-panel)
  * [Processed SNiPDx object](#processed-snipdx-object)
- [Usage](#usage)
  * [1) Process Local data](#1--process-local-data)
    + [Command line parameters](#command-line-parameters)
    + [Running the pipeline](#running-the-pipeline)
  * [2) Quantify SNiPDx and Local Concordance](#2--quantify-snipdx-and-local-concordance)
    + [Command line parameters](#command-line-parameters-1)
    + [Running the pipeline](#running-the-pipeline-1)

## Workflow overview

![](../workflows/snipdx_local_concordance.PNG?raw=true)

## Prerequisites

### Box credentials 

Raw local data is pulled from the [Box](https://solvebio.account.box.com/login) software provided by the vendor ***SolveBio***. To extract the vendor data, users need to request access to Box.

### SolveBio Local test consolidation

This file can be found here: [./ref/SolveBio_Local_Test_Consolidation.xlsx](./ref/SolveBio_Local_Test_Consolidation.xlsx)

### SNiPDx gene panel 

This file can be found here: [../snipdx/ref/panel_genes_annotated_V2.tsv](../snipdx/ref/panel_genes_annotated_V2.tsv)

### Processed SNiPDx object

The concordance analysis executed following the snipdx project level pipeline on (nf-snipdx)[https://bitbucket.org/reparecompbio/nf-snipdx/src/master/README.md]. The integrated snipdx inventory from the processed snipdx s4 object will be the input parameter used to determine concordance by variant types.



## Usage

The stepwise snipdx and local concordance pipeline should be executed in the following order:

### 1) Process Local data

This processing command line script performs the following:

* reads in the raw local data and merges with local test consolidatio file
* cleans and partitions by assay and variant types to their respective local s4 object slots
* reads in the SNiPDx gene panel file and labels whether genes from local data is a gene on the SNiPDx panel
* exports the cleaned local data as `[date]_[tria]_processed_local_data.rds` 


#### Command line parameters

| Parameters            | Description                                                            |
|-----------------------|------------------------------------------------------------------------|
| trial                 | specify trial                                                          |
| file                  | specify the path to raw local data                                     |
| output_dir            | specify the path of the output directory                               |
| git_extension         | path to functions extension folder (e.g clinicalbioinfoextensions)     |
| local_reference_path  | specify the path of the consolidation file for local alterations (e.g  |
| snipdx_reference_path | specify the path of the reference file for snipdx genes                |

#### Running the pipeline

```
Rscript ../git-repos/reparecompbio/clinicalbioinfoCLIs/snipdx_local_concordance/src/local_data_processing.R \
--file ../vendor-data/RP-3500-01\ Retrospective\ Local\ Reports\ 18DEC2022.xlsx \
--output_dir ../data/ \
--git_extension ../git-repos/reparecompbio/clinicalbioinfoextensions/ \
--local_reference_path ../git-repos/reparecompbio/clinicalbioinfoCLIs/snipdx_local_concordance/ref/SolveBio_Local_Test_Consolidation.xlsx \
--snipdx_reference_path ../SNiPDx/Reference/panel_genes_annotated_V2.tsv \
--trial "RP3500-01" 
```

### 2) Quantify SNiPDx and Local Concordance

This processing command line script performs the following:

* read in processed local and snipdx object, pulls ngs assay, detects shared or unique snv/indel or cnv calls between both platforms 
* builds an html report and exports the concordance s4 object labelled as `[date]_[trial]_local_and_snipdx_concordance.rds`

#### Command line parameters

| Parameters          | Description                                                                                                                          |
|---------------------|--------------------------------------------------------------------------------------------------------------------------------------|
| currTrial           | specify trial                                                                                                                        |
| inExtensions        | path to functions extension folder (e.g clinicalbioinfoextensions)                                                                   |
| inSNiPDx            | specify the path to the processed SNiPDx object                                                                                      |
| inLocal             | specify the path to the processed local object                                                                                       |
| knit_root_dir       | specify the path to the root directory to output the concordance report `./results/[date]_[trial]_local_and_snipdx_concordance.html` |
| blacklistedPatients | specify a vector of blacklisted patients to ignore during the concordance analysis                                                   |

#### Running the pipeline

```
currTrial="RP3500-01"

inExtensions="/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoextensions/"

inSNiPDx="/ClinBio/SP-ClinicalBioinformatics/shared/SNiPDx/RP3500-01/processing/data/2022-07-07_RP3500-01_SNiPDx_Complete_Inventory_Processed_Data.rds"

inLocal="/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-01/SolveBio/processing/data/local-processed-data/2022-07-04/2022-07-04_RP3500-01_processed_local_data.rds"

knit_root_dir="/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-01/SolveBio/processing/"

blacklistedPatients="3003-0211"

 Rscript -e \
  "rmarkdown::render('/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/snipdx_local_concordance/src/generate_local_SNiPDx_concordance_report.Rmd', 
  knit_root_dir = '$knit_root_dir',
  output_dir = paste0('$knit_root_dir/results/'),
  output_file = paste0(Sys.Date(),'_', '$currTrial','_local_snipdx_concordance',
                       '.html'),
  params=list(
  inExtensions = '$inExtensions',
  inSNiPDx = '$inSNiPDx',
  inLocal = '$inLocal',
  currTrial = '$currTrial',
  blacklistedPatients = '$blacklistedPatients'))"
```


