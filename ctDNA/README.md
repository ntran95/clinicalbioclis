# ctDNA analysis
______________________

# Table of Contents

- [Workflow overview](#workflow-overview)
- [Prerequisites](#prerequisites)
- [Usage](#usage)
  * [1) Process ctDNA (unfiltered)](#1--process-ctdna--unfiltered-)
    + [Command line parameters](#command-line-parameters)
    + [Running the pipeline](#running-the-pipeline)
  * [2) Calculate VAFs (unfiltered)](#2--calculate-vafs--unfiltered-)
    + [Command line parameter](#command-line-parameter)
    + [Running the pipeline](#running-the-pipeline-1)
  * [3) Calculate mVAFs (filtered)](#3--calculate-mvafs--filtered-)
    + [Command line parameters](#command-line-parameters-1)
    + [Running the pipeline](#running-the-pipeline-2)
  * [4) Trial specific clinical integration (filtered)](#4--trial-specific-clinical-integration--filtered-)
    + [Command line parameters](#command-line-parameters-2)
    + [Running the pipeline](#running-the-pipeline-3)


## Workflow overview

![](../workflows/ctdna.png?raw=true)

## Prerequisites

The processed clinical s4 object is a prerequisite to initializing the calculating VAFs (step 2) of the ctDNA pipeline. The processed clinical table (`clinical@combined`) will be integrated with the VAF summary table at this step

## Usage

A stepwise demo to execute the complete ctDNA pipeline, see the [test.sh](test.sh) example.

To show the list of parameter descriptions for each processing command line script use the following command:
```
Rscript ./src/[script_name.R] --help
```

The stepwise ctDNA pipeline should be executed in the following order:

### 1) Process ctDNA (unfiltered)

This processing command line script performs the following:

* reads in the external tempus data, clean manifest, qc, and variant. 
* integrates qc and manifest. 
* stores all data into corresponding slots in ctDNA S4 obj and exports RDS as `[date]_[trial]_processedCtDNA.rds`


#### Command line parameters

| Parameter        | Description                                                                                                                 |
|------------------|-----------------------------------------------------------------------------------------------------------------------------|
| date             | date (format: yyyy-mm-dd), labels all output files by specific date                                                         |
| trial            | trial number (e.g. RP6306-01), labels all output files by specific trials                                                   |
| manifest         | manifest directory path (e.g /ClinBio/SP-ClinicalBioinformatics/shared/RP6306/RP6306-01/ctDNA/vendor-data/Tempus/Manifests/ |
| qc               | qc file path from latest batch                                                                                              |
| mmf              | mmf file path from latest batch                                                                                             |
| cnv              | cnv file path from latest batch                                                                                             |
| monitoring_avail | logical to indicate if monitoring file is available, set to T or F                                                          |
| monitoring       | file path to previous, annotated monitoring file if monitoring_avail is T                                                   |
| ctfe_avail       | logical, indicating whether ctfe file is available, set to T or F                                                           |
| ctfe             | path to ctfe file                                                                                                           |
| outdir           | output directory                                                                                                            |
| functions        | path to functions extension folder (e.g clinicalbioinfoextensions)                                                          |

#### Running the pipeline

```
Rscript ../clinicalbioinfoCLIs/ctDNA/src/process_ctDNA.R \
		--manifest '../vendor-data/Tempus/Manifests/' \
		--qc '../vendor-data/Tempus/Batch7_2022.09.23/Repare_3500-03_xFSeq_QC_2022_09_22.xlsx' \
		--mmf '../vendor-data/Tempus/Batch7_2022.09.23/g_molecular_master_file_filtered.csv' \
		--cnv '../vendor-data/Tempus/Batch7_2022.09.23/g_cnv_segments_cnvkit.csv' \
		--monitoring '../data/monitoring/2022-10-26_RP3500-03_processed_ctDNA_Tempus_Batch7_Revised_MonitorIn_2022-11-21.xlsx' \
		--monitoring_avail 'T' \
		--outdir '../data/' \
		--functions '..//git-repos/reparecompbio/clinicalbioinfoextensions/' \
		--trial 'RP3500-03' \
		--date '2022-11-22' \
		--ctfe '../vendor-data/Tempus/Batch7_2022.09.23/Repare_3500-03_xFSeq_ctFE.xlsx' \
		--ctfe_avail 'T'

```

### 2) Calculate VAFs (unfiltered)

This processing command line script performs the following: 

* reads in ctnda obj
* filters by SV, integrate manifestQC with variants, annotate Visit IDs, integrate monitoring In with new variants
* read in clinical s4 obj and integrates with the integrated vaf summary table
* build VAF plots with unfiltered variants for each patient and stores into the `ctdn@VAFplots` slot
* exports new monitoring Out file with new variants from latest batches as `[date]_[trial]_monitoring_required.xlsx`
* exports an updated ctdna s4 object as `[date]_[trial]_unfilteredCtDNA.rds`
* generates a unfiltered ctdna html report rendered by rmarkdown if the `--rmarkdown_report` is set to `T`, labelled as `[date]_[trial]_unfiltered_ctDNA-Tempus.html` and as `[date]_[trial]_unfiltered_ctDNA-Tempus.aspx` (specifc to browse report via a Sharepoint link)

#### Command line parameter

| Parameter        | Description                                                                                                                                                                                                                                                                                                                  |
|------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| date             | date (format: yyyy-mm-dd), labels all output files by specific date                                                                                                                                                                                                                                                          |
| trial            | trial number (e.g. RP6306-01), labels all output files by specific trials                                                                                                                                                                                                                                                    |
| study            | study object path matching trial (e.g [date]_[trial]_processedStudy.rds)                                                                                                                                                                                                                                                     |
| ctdna            | ctdna object path matching trial (e.g [date]_trial]_processedCtDNA.rds)                                                                                                                                                                                                                                                      |
| datadir          | [output directory                                                                                                                                                                                                                                                                                                            |
| monitoring_avail | logical to indicate if monitoring file is available, set to T or F if monitoring_avail is set to F, a new monitoring file is generated if monitoring_avail is set to T, pass the previous monitoring file from `ctdna@monitoring` and integrate with new variants from latest batch before exporting new monitoring out file |
| mondir           | monitoring directory path                                                                                                                                                                                                                                                                                                    |
| functions        | path to functions extension folder (e.g clinicalbioinfoextensions)                                                                                                                                                                                                                                                           |
| rmarkdown_report | logical indicating whether to generate unfiltered VAF rmarkdown report                                                                                                                                                                                                                                                       |
| resultsdir       | results directory, applicable when rmarkdown_report is TRUE                                                                                                                                                                                                                                                                  |
| rmd_path         | path to unfiltered rmd script, applicable when rmarkdown_report is TRUE (e.g ./src/generate_rmarkdown_report_unfiltered)ctDNA.Rmd)                                                                                                                                                                                           |

#### Running the pipeline

```
### change variables below
date="2022-12-09"
trial="RP3500-03"
study="../clinical/processing/data/2022-12-06/2022-12-06_RP3500-03_processedStudy.rds"
ctdna="../data/2022-11-22_RP3500-03_processedCtDNA.rds"
rmarkdown_report="T"
Rscript /ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/src/calculate_VAFs.R \
		--date "$date" \
		--trial "$trial" \
		--study "$study" \
		--ctdna "$ctdna" \
		--datadir "../data/" \
		--mondir "../data/monitoring/" \
		--monitoring_avail 'T'\
		--functions "../git-repos/reparecompbio/clinicalbioinfoextensions/" \
		--rmarkdown_report "$rmarkdown_report" \
		--resultsdir "../results/" \
		--rmd_path "..//git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/src/generate_rmarkdown_report_unfiltered_ctDNA.Rmd"
```

### 3) Calculate mVAFs (filtered)

This processing command line script performs the following: 

* reads in the unfiltered ctDNA object and annotated monitoring file, integrate
* filters for somatic variants, calculate mVAF, mVAFR, mVAFD, and label MR
* create variant, patient, and sample level table summaries and store into their respective slots in the ctDNA object
* build VAF plots with filtered somatic variants for each patient and update the `ctdn@VAFplots` slot
* build mVAF plots with filtered somatic variants for each patient and update the `ctdn@mVAFplots` slot
* builds CNV plots based on genes listed by the `cnvgene` parameters
* exports an updated ctdna s4 object labelled as `[date]_[trial]_filteredCtDNA.rds`

#### Command line parameters

| Parameter  | Description                                                               |
|------------|---------------------------------------------------------------------------|
| date       | date (format: yyyy-mm-dd), labels all output files by specific date       |
| trial      | trial number (e.g. RP6306-01), labels all output files by specific trials |
| ctdna      | unfiltered ctdna object path matching trial (e.g [date]_trial]_unfilteredCtDNA.rds)  |
| monitoring | annotated monitoring file path (e.g [date]_[trial]_completed.xlsx)        |
| outdir     | output directory for rds file                                             |
| functions  | path to functions extension folder (e.g clinicalbioinfoextensions)        |

#### Running the pipeline

```
### replace variable here ###
date="2022-12-11"
trial="RP3500-03"
ctdna="../Tempus/data/2022-12-09_RP3500-03_unfilteredCtDNA.rds"
monitoring="../Tempus/data/monitoring/2022-11-22_RP3500-03_monitoring_completed.xlsx"
Rscript ../git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/src/calculate_mVAFR.R \
                --date "$date" \
                --trial "$trial" \
                --ctdna "$ctdna" \
                --monitoring "$monitoring" \
                --outdir "../data/" \
                --cnvgene "BRCA1","BRCA2" \
                --functions "../../../../../../git-repos/reparecompbio/clinicalbioinfoextensions/"
                
                
#### optional: stand alone command to knit unfiltered report after generating mVAF response rds ####
rmarkdown_report=true
if($rmarkdown_report); then
### change variables below
ctdna="../data/2022-12-11_RP3500-03_filteredCtDNA.rds"
date="2022-12-11"
trial="RP3500-03"
knit_root_dir="$PWD"
  Rscript -e \
    "rmarkdown::render('../git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/src/generate_rmarkdown_report_filtered_ctDNA.Rmd', 
    output_dir = '../results/',
    output_file = paste0('../results/','$date','_', '$trial', '_filtered_ctDNA', '.html'),
    knit_root_dir = '$knit_root_dir',
    params=list(
    ctdna = '$ctdna',
    date = '$date',
    trial = '$trial',
    functions = '../git-repos/reparecompbio/clinicalbioinfoextensions/',
    resultsdir = '../results/',
    response = 'T'
    ))"
fi

```

### 4) Trial specific clinical integration (filtered)

This processing command line script performs the following: 

* read in filtered ctdna object
* build best ctDNA response plots, and clinical associations, store in ctDNA obj
* generates a filtered ctdna html report rendered by rmarkdown if the `--rmarkdown_report` is set to `T`, labelled as `[date]_[trial]_filtered_ctDNA-Tempus.html` and as `[date]_[trial]_filtered_ctDNA-Tempus.aspx` (specifc to browse report via a Sharepoint link)

The command script can be trial specific with varying formats for response plots.

#### Command line parameters

| Parameter        | Description                                                                                                             |
|------------------|-------------------------------------------------------------------------------------------------------------------------|
| date             | date (format: yyyy-mm-dd), labels all output files by specific date                                                     |
| trial            | trial number (e.g. RP6306-01), labels all output files by specific trials                                               |
| ctdna            | filtered ctdna object path matching trial (e.g [date]_trial]_filteredCtDNA.rds)                                         |
| datadir          | output directory to store reponse rds file                                                                              |
| functions        | path to functions extension folder (e.g clinicalbioinfoextensions)                                                      |
| rmarkdown_report | logical indicating whether to generate filtered markdown reports with MR reponse plots, set to 'T' or 'F'               |
| resultsdir       | results directory to store html report raw pngs, applicable when rmarkdown_report is TRUE                               |
| rmd_path         | path to filtered rmd script, applicable when rmarkdown_report is TRUE (e.g generate_rmarkdown_report_filtered_ctDNA.Rmd |

#### Running the pipeline

```
### replace variable here ###
date="2022-12-11"
trial="RP3500-03"
ctdna="../data/2022-12-11_RP3500-03_filteredCtDNA.rds"
Rscript ../git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/src/clinical_integration_ctDNA_RP3500.R \
                --date "$date" \
                --trial "$trial" \
                --ctdna "$ctdna" \
                --datadir "../data/" \
		            --functions "../git-repos/reparecompbio/clinicalbioinfoextensions/" \
		            --rmarkdown_report "T" \
		            --resultsdir "../results/" \
		            --rmd_path "../git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/src/generate_rmarkdown_report_filtered_ctDNA.Rmd"
		            
```



