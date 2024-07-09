# Clinical Bioinformatics CLIs

A repository of command line processing scripts to generate clinical, ctDNA, snipdx, local/snipdx concordance analysis

______________________

# Table of Contents

- [Prerequisites](#prerequisites)
  * [Importing Conda enviroment](#importing-conda-enviroment)
    + [Conda installation](#conda-installation)
    + [Updating existing environment](#updating-existing-environment)
  * [Importing clinicalbioinfoextensions package](#importing-clinicalbioinfoextensions-package)
- [ctDNA analysis](#ctdna-analysis)
  * [Workflow overview](#workflow-overview)
- [SNiPDx analysis](#snipdx-analysis)
  * [Workflow overview](#workflow-overview-1)
- [SNiPDx and Local Concordance analysis](#snipdx-and-local-concordance-analysis)
  * [Workflow overview](#workflow-overview-2)
- [Clinical analysis](#clinical-analysis)
  * [Workflow overview](#workflow-overview-3)



## Prerequisites 

### Importing Conda enviroment

For reproducibility, users can copy the exact conda environment exported to the [environment.yml](./environment.yml) file. 

#### Conda installation

```
conda env create --name R-4.0.5 -f environment.yml
```

#### Updating existing environment 

```
conda activate R-4.0.5

conda install -c conda-forge r-[package_name] 

conda env export > environment.yml

```


### Importing clinicalbioinfoextensions package

Each analysis usually contains a parameter that reads in the path to the `clinicalbioinfoextensions` package. This package  contains supplemental functions aiding in the specific processes as well as the initializing the s4 class for each analysis. 

In an interactive R session, the `clinicalbioinfoextensions` package can be loaded using the following command:

```
repo_path <- "/ClinBio/SP-ClinicalBioinformatics/shared/git-repos/reparecompbio/clinicalbioinfoextensions/"
devtools::load_all(repo_path)
```

______________________

## ctDNA analysis

The stepwise ctDNA pipeline should be executed in the following order:

1) Process ctDNA (unfiltered)
2) Calculate VAFs (unfiltered)
3) Calculate mVAFs (filtered)
4) Trial specific clinical integration (filtered)

For full details on the parameters and output structure from this pipeline is found [here](./ctDNA/README.md)

### Workflow overview

![](workflows/ctdna.png?raw=true)


______________________

## SNiPDx analysis

This analysis is deprecated and moved to the [nf-snipdx](https://bitbucket.org/reparecompbio/nf-snipdx/src/master/) pipeline

For full details on the parameters and output structure from this pipeline is found [here](./snipdx/README.md)

### Workflow overview

![](workflows/snipdx.PNG?raw=true)

______________________

## SNiPDx and Local Concordance analysis

The stepwise snipdx and local concordance pipeline should be executed in the following order:

1) Process Local data
2) Quantify SNiPDx and Local Concordance

For full details on the parameters and output structure from this pipeline is found [here](./snipdx_local_concordance/README.html)

### Workflow overview

![](workflows/snipdx_local_concordance.PNG?raw=true)

______________________

## Clinical analysis 

This 2 step pipeline outputs the following:

1) Processes the raw clinical tracker and stores all summary tables into their respective slots in the s4 object
2) Reads in the processed clinical s4 object and generates contingency tables (i.e patient enrollment gene by status, type of loss by enrollment gene) and plots and knits outputs into an html report.

For full details on the parameters and output structure from this pipeline is found [here](./clinical/README.md)

### Workflow overview

![](workflows/clinical.PNG?raw=true)



