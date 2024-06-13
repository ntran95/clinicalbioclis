# PD analysis

This pipeline outputs the following:

1) Cleans and integrates the baseline and paired biopsy with clinical data
2) Builds H-score plots for baseline and paired biopsy
3) Renders a html report showing tables, significance tests, and plots for baseline and paired biopsy


## Workflow overview

![](../workflows/pd.PNG?raw=true)

## Usage

### Command line parameters

| Parameter       | Description                                                        |
|-----------------|--------------------------------------------------------------------|
| functions       | path to functions extension folder (e.g clinicalbioinfoextensions) |
| trial           | trial number (e.g. RP6306-01)                                      |
| date            | date (format: yyyy-mm-dd)                                          |
| baseline_biopsy | path to baseline biopsy data                                       |
| paired_biopsy   | path to paired biopsy data                                         |
| study           | path to processed clinical study object, RDS file                  |
| rmd_path        | path to pd rmd script to render report                             |
| resultsdir      | results directory to store report                                  |

#### Running the pipeline

```
conda activate R-4.0.5

outdir="/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/pd/test/data/"
functions="/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoextensions/"
trial="RP3500-01"
date="2022-12-30"
study="/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-01/clinical/processing/data/2022-12-20/2022-12-20_RP3500-01_processedStudy.rds"
baseline_biopsy="/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-01/pd/vendor-data/Histowiz15468_ScoreSheet.xlsx"
paired_biopsy="/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-01/pd/vendor-data/Cumulative_RP-3500-01_paired biopsy scoresheet_10NOV2021.xlsx"
rmd_path="./src/generate_rmarkdown_report_pd.Rmd"
resultsdir="/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/pd/test/results/"
Rscript ./src/process_pd.R \
		--outdir "$outdir" \
		--functions "$functions" \
		--trial "$trial" \
		--date "$date" \
		--study "$study" \
		--baseline_biopsy "$baseline_biopsy" \
		--paired_biopsy "$paired_biopsy" \
		--rmd_path "$rmd_path" \
		--resultsdir "$resultsdir"
```		
