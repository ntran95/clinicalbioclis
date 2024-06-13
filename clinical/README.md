# Clinical analysis

This 2 step pipeline outputs the following:

1) Processes the raw clinical tracker and stores all summary tables into their respective slots in the s4 object
2) Reads in the processed clinical s4 object and generates contingency tables (i.e patient enrollment gene by status, type of loss by enrollment gene) and plots and knits outputs into an html report.

______________________

# Table of Contents

- [Workflow overview](#workflow-overview)
- [Prerequisites](#prerequisites)
- [Future implementation](#future-implementation)
- [Usage](#usage)
  * [0.5) Prepare Quarzbio tracker](#05--prepare-quarzbio-tracker)
    + [Command line parameters](#command-line-parameters)
  * [1) Process clinical](#1--process-clinical)
    + [Command line parameters](#command-line-parameters-1)
    + [Running the pipeline](#running-the-pipeline)
  * [2) Clinical reports (deprecated)](#2--clinical-reports--deprecated-)


## Workflow overview

![](../workflows/clinical.PNG?raw=true)

## Prerequisites

Column names within the raw trackers differ between trials. For example, the column labeling `Module` in RP3500-01 is labelled as `Arm` in RP3500-03--making these columns inaccessible in the processing scripts. Harmonizing the differences in clinical column names between trials is a prerequisite to processing the clinical pipeline. To do so, variations in column names across trials are tracked in the [master_key_across_programs.csv](ref/master_key_across_programs.csv) file. 

The `master_key_across_programs.csv` file contains a list of column names from the clinical tracker that meets the minimum requirement to execute the clinical pipeline. Based on the user specified `trial` parameter, the file will take the current trial and map column names back to the uniform names from the `Master_Key` column. This action is executed in the [`clinicalbioinfoextensions::programToMasterKey()`](https://bitbucket.org/reparecompbio/clinicalbioinfoextensions/src/master/R/utils.R) function. 

## Future implementation

An implementation yet to be in fully in place is a tracker sourced from Quartzbio. The purpose of Quartzbio tracker is provide a consistent version of the clinical tracker across trials. The Quartzbio tracker will contain the Clinical, Response, and Tumor marker tabs so we will still need to extract Genomic and SNiPDx data from the Repare tracker. An auxiliary script was developed to reformat the Quartzbio tracker prior the clinical processing step, see [./src/process_quartz_clinical_tracker.R](./src/process_quartz_clinical_tracker.R)

## Usage

To show the list of parameter descriptions for each processing command line script use the following command:
```
Rscript ./src/process_ClinicalAndGenomics.R--help
```

The stepwise clinical pipeline should be executed in the following order:

### 0.5) Prepare Quarzbio tracker

If the tracker is sourced from Quartzbio, use this step to reformat the file so that it is compatible for processing in [step #1](#1--process-clinical)

#### Command line parameters
| Parameters     | Description                                                        |
|----------------|--------------------------------------------------------------------|
| quartz_tracker | quartz tracker data path                                           |
| repare_tracker | repare tracker data path                                           |
| outdir         | output directory                                                   |
| functions      | path to functions extension folder (e.g clinicalbioinfoextensions) |
| trial          | trial number (e.g. RP6306-01"                                      |
| date           | date (format: yyyy-mm-dd)                                          |

### 1) Process clinical 

This processing command line script performs the following:

* reads in the raw clinical tracker
* cleans and preserves the Clinical, Response, Markers, Dose Transfusion, Genomics and SNiPDX tabs
* maps trial specific differences in column names within each clinical tab and relabels to uniform names found the `Master_Key` column 
* derives BOR from response tab
* derives COR from markers tab
* merges all clinical tabs into a single dataframe in the clinical@Combined slot
* if `genomicstats` is `T`, merge ctDNA molecular response with genomic info for biostats
* exports the cleaned clinical s4 object labelled as `[date]_[trial]_processedStudy.rds

#### Command line parameters

| Parameter    | Description                                                                                          |
|--------------|------------------------------------------------------------------------------------------------------|
| date         | date (format: yyyy-mm-dd), labels all output files by specific date                                  |
| trial        | trial number (e.g. RP6306-01), labels all output files by specific trials                            |
| tracker      | tracker data path                                                                                    |
| outdir       | output directory                                                                                     |
| functions    | path to functions extension folder (e.g clinicalbioinfoextensions)                                   |
| markers      | logical to indicate if marker data is available, set to T or F                                       |
| response     | logical to indicate if response data is available, set to T or F                                     |
| snipdx       | logical to indicate if SNIPDx data is available, set to T or F                                       |
| keyAvail     | logical to indicate whether path to program key is should be read in                                 |
| key          | specify the path to the master key labelling each program (e.g ./ref/master_key_across_programs.csv) |
| genomicstats | logical, indicate whether monthly report of Genomics and best ctDNA response be exported to Stats    |
| inctdna      | path to ctDNA vaf dataframe, read in best mVAFR                                                      |

#### Running the pipeline

```
tracker="../vendor-data/2022-12-12_RP3500-01 Patient Tracker.xlsx"
date="2022-12-20"
Rscript '../git-repos/reparecompbio/clinicalbioinfoCLIs/clinical/src/process_ClinicalAndGenomics_v1.0.1.R' \
		--tracker "$tracker" \
		--outdir "../data/"$date"/" \
		--functions '../git-repos/reparecompbio/clinicalbioinfoextensions/' \
		--trial 'RP3500-01' \
		--date "$date" \
		--markers 'T' \
		--response 'T' \
		--snipdx 'T'\
		--keyAvail 'T' \
		--key '../Program-Key/master_key_across_programs.csv' \
		--genomicstats 'T' \
		--inctdna '../ctDNA/processing/Combined/results/2022-12-09_RP3500-01_filtered_ctDNA_results/2022-12-09_RP3500-01_filtered_ctDNA_Summary.xlsx'

```
### 2) Clinical reports (deprecated)

This pipeline is ***deprecated*** and superseded by the clinical shiny app maintained by BioStats [https://rshiny.reparerx.net/RP3500-1/](https://rshiny.reparerx.net/RP3500-1/). The output in this secondary step in the pipeline is a rmarkdown rendered html report containing contingency tables and patient tumor/response plots. The historical  code is stored in the clinical analysis [link to clinical] subfolder for reference however requires previous version the clinical objet outputted from the old processing script found [here](./src/processing_clinical_v1.0.4.R) 


