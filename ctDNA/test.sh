
Rscript ./src/process_ClinicalAndGenomics.R \
		--tracker '/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/clinical/vendor-data/2022-10-29_RP-3500-03 ATTACC Patient Status Tracker_Working Copy.xlsx' \
		--outdir '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/data/' \
		--functions '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoextensions/' \
		--trial 'RP3500-03' \
		--date '2022-10-31' \
		--markers 'T' \
		--response 'T' \
		--snipdx 'T'\
		--keyAvail 'T' \
		--key '/ClinBio/SP-ClinicalBioinformatics/shared/Program-Key/master_key_across_programs.csv'\
		--genomicstats 'T' \
		--inctdna '/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/ctDNA/processing/Combined/results/2022-08-30_RP3500-03_processed_ctDNA_Combined_results/2022-08-30_RP3500-03_processed_ctDNA_Combined_Summary.xlsx'



Rscript ./src/process_ctDNA.R \
		--manifest '/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/ctDNA/vendor-data/Tempus/Manifests/' \
		--qc '/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/ctDNA/vendor-data/Tempus/Batch7_2022.09.23/Repare_3500-03_xFSeq_QC_2022_09_22.xlsx' \
		--mmf '/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/ctDNA/vendor-data/Tempus/Batch7_2022.09.23/g_molecular_master_file_filtered.csv' \
		--cnv '/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/ctDNA/vendor-data/Tempus/Batch7_2022.09.23/g_cnv_segments_cnvkit.csv' \
		--monitoring '/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/ctDNA/processing/Tempus/data/monitoring/2022-10-26_RP3500-03_processed_ctDNA_Tempus_Batch7_Revised_MonitorIn_2022-10-26.xlsx' \
		--monitoring_avail 'T' \
		--outdir '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/data/' \
		--functions '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoextensions/' \
		--trial 'RP3500-03' \
		--date '2022-11-07' \
		--ctfe '/ClinBio/SP-ClinicalBioinformatics/shared/RP3500/RP3500-03/ctDNA/vendor-data/Tempus/Batch7_2022.09.23/Repare_3500-03_xFSeq_ctFE.xlsx' \
		--ctfe_avail 'T'


outdir="/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/data/"
Rscript ./src/calculate_VAFs.R \
		--date '2022-11-07' \
		--trial 'RP3500-03' \
		--study '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/data/2022-10-31_RP3500-03_processedStudy.rds' \
		--ctdna '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/data/2022-11-07_RP3500-03_processedCtDNA.rds' \
		--datadir "$outdir" \
		--mondir '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/data/monitoring/' \
		--monitoring_avail 'T'\
		--functions '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoextensions/' \
		--rmarkdown_report 'F' \
		--resultsdir '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/results/' \
		--rmd_path '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/src/generate_rmarkdown_report_unfiltered_ctDNA.Rmd'


### execute if calculate_vaf rds is already generated and rendering report is to follow		
rmarkdown_report=false
if($rmarkdown_report); then
  ctdna="/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/data/2022-11-10_RP3500-03_filtered_responses_CtDNA.rds"
  rmd_path="/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/src/generate_rmardown_report_unfiltered_ctDNA.Rmd"
  knit_root_dir="$PWD"
Rscript -e \
  "rmarkdown::render('$rmd_path', 
  output_dir = paste0('/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/results/'),
  output_file = paste0('2022-11-03','_', 'RP3500-03','_unfiltered_ctDNA_response',
                       '.html'),
  knit_root_dir = '$knit_root_dir',
  params=list(
  ctdna = '$ctdna',
  date= '2022-11-10',
  trial= 'RP3500-03',
  functions= '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoextensions/',
  resultsdir= '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/results/',
  ))"
fi
	

Rscript ./src/calculate_mVAFR.R \
                --date '2022-11-07' \
                --trial 'RP3500-03' \
                --ctdna '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/data/2022-11-07_RP3500-03_unfilteredCtDNA.rds' \
                --monitoring '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/data/monitoring/2022-11-07_RP3500-03_monitoring_completed.xlsx' \
                --outdir '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/data/' \
                --cnvgene 'BRCA1','BRCA2' \
		--functions '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoextensions/'


### trial specific repsonse plot processing script (optional)
Rscript ./src/clinical_integration_ctDNA_RP3500.R \
                --date '2022-11-10' \
                --trial 'RP3500-03' \
                --ctdna '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/data/2022-11-07_RP3500-03_filteredCtDNA.rds' \
                --datadir '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/data/' \
		--functions '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoextensions/' \
		--rmd_path '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/src/generate_rmarkdown_report_filtered_ctDNA.Rmd' \
		--rmarkdown_report 'F' \
		--resultsdir '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/results/' 



		
# generate filtered Rmd report
# if response plots have been generated and stored in R obj, the reponse section will populate in report
ctdna="/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/data/2022-11-10_RP3500-03_filtered_responses_CtDNA.rds"
rmd_path="/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/src//generate_rmarkdown_report_filtered_ctDNA.Rmd"
knit_root_dir="$PWD"
Rscript -e \
  "rmarkdown::render('$rmd_path', 
  output_dir = paste0('/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/results/'),
  output_file = paste0('2022-11-07','_', 'RP3500-03','_filtered_ctDNA_response',
                       '.html'),
  knit_root_dir = '$knit_root_dir',
  params=list(
  ctdna = '$ctdna',
  date= '2022-11-10',
  trial= 'RP3500-03',
  functions= '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoextensions/',
  resultsdir= '/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoCLIs/ctDNA/test/results/',
  response ='T'))"
  
  


