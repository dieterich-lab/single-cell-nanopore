# Singlesquare

Singlesquare is designed for cell barcode assignment from Nanopore sequencing data.
It requires bam files from both Nanopore and Illumina reads, then builds a model based on the parameters estimated from the two libraries.

# Installation:

1. Most software dependecies are managed using **`conda`**. Please install as described at  <br> [https://docs.conda.io/projects/conda/en/latest/user-guide/install/](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).
```
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    bash
```
2. Download the files into folder named `single-cell-nanopore`. The program requires the **`NanoSim`** which should be installed first through **`conda`**
```
    conda install -c bioconda nanosim
    git clone https://github.com/nanoporetech/single-cell-nanopore.git single-cell-nanopore
```
3. Change working directory into the new `single-cell-nanopore` folder 
```
    cd single-cell-nanopore
```

4. Use these commands for building **`SingleCellPipe`**: 
```
    cmake .
    make
```
5. Install conda software dependencies with
```
    conda env create --name single-cell-nanopore --file environment.yaml
```
6. Initialise conda environment with 
```
    conda activate single-cell-nanopore
```

# Quick run

1. *optional* edit the provided **`config.yaml`** file to match your own sequence files, reference genome and annotation. Update the adapter and polyT length that fit your libraries.

2. Run the **`snakemake`** command under the conda environment to perform the bioinformatics analysis on the specified sequence files. Several analysis steps can benefit from multiple computer cores; use the `-j` parameter to parallise the analysis (this document is assuming that 8 cores are available).
```
    snakemake -j 8 all
```
You can also submit the job via job schedulers. We have provided an example on the SLURM. Change the account settings in cluster.json before using it.
```
    snakemake -j 8 --cluster-config cluster.json --cluster "sbatch -A {cluster.account} --mem={cluster.mem} -t {cluster.time} -c {cluster.threads}"
```

Input files:

* [barcodes.tsv.gz](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices)  and [possorted_genome_bam.bam](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/overview) from Cell Ranger.

* Nanopore reads ([sample](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/data/)) in FATSQ format

* Reference genome and transcriptome in FATSA format

* Annotation file in ensemble GTF/GFF formats

Output files:

* Barcode assignments above the threshold (0-99, 30 in the example we provided).

* Demultiplexed FATSQ file

# Workflow

[![](https://mermaid.ink/img/eyJjb2RlIjoiZ3JhcGggVERcbiAgICBBW0NlbGxSYW5nZXI8YnI-Y291bnRdIC0tPnxidWlsZF9nZW5vbWUuc2h8QihBcnRpZmljYWwgR2Vub21lKVxuICAgIEEgLS0-fGdldF9jYmZyZXEuc2h8QjIoQ2VsbCBCYXJjb2RlPGJyPkZyZXF1ZW5jaWVzKVxuICAgIEEgLS0-fGdldF9jYmMuc2h8QzIoQ2VsbCBCYXJjb2RlPGJyPldoaXRlbGlzdClcbiAgICBDW05hbm9wb3JlPGJyPlJlYWRzXSAtLT58YnVpbGRfbmFub3NpbS5zaHxEKE5hbm9TaW08YnI-TW9kZWxzKVxuICAgIEMyLS0-IHxmaW5kX2Rpc3QucnxDMyhFcnJvci1jb250YWluaW5nPGJyPldoaXRlbGlzdClcbiAgICBDMy0tPiBIXG4gICAgQiAtLT4gfHNpbV9yZWFkcy5zaHxFKFNpbXVsYXRlZDxicj5OYW5vcG9yZTxicj5SZWFkcylcbiAgICBEIC0tPiBFXG4gICAgRSAtLT4gfGJ1aWxkX2FsaWduLnNofEYoU2ltdWxhdGVkPGJyPk5hbm9wb3JlPGJyPkFsaWdubWVudHMpXG4gICAgQiAtLT4gfGdldF9iYXJjb2Rlcy5zaHxHKEdyb3VuZCBUcnV0aDxicj5CYXJjb2RlcylcbiAgICBFIC0tPiBHXG4gICAgRiAtLT4gfHJ1bl9waXBlLnNofEgoTGFiZWxlZDxicj5GZWF0dXJlcylcbiAgICBNIC0tPiB8cHJlZC5yfExcbiAgICBIIC0tPiBMKFByb2JhYmlsaXRpZXM8YnI-b2YgU2ltdWxhdGVkPGJyPkJhcmNvZGUgTWF0Y2hlcylcbiAgICBMIC0tPiB8c2ltX2JlbmNoLnIsIHN0YXQuc2h8TihDbGFzc2lmaWNhdGlvbjxicj5QZXJmb3JtYW5jZSlcbiAgICBIIC0tPiB8ZmVhdF9zdGF0LnJ8SShGZWF0dXJlczxicj5JbXBvcnRhbmNlKVxuICAgIEggLS0-IHxmZWF0X3Zhci5yfE8oRmVhdHVyZXM8YnI-VmFyaWFuY2UpXG4gICAgRyAtLT4gSFxuICAgIEggLS0-IHxidWlsZF9tb2RlbC5yfE0oTmFpdmUgQmF5ZXM8YnI-TW9kZWwpXG4gICAgc3R5bGUgSSBmaWxsOiNmZmNcbiAgICBzdHlsZSBOIGZpbGw6I2ZmY1xuICAgIHN0eWxlIE8gZmlsbDojZmZjXG4gICAgY2xpY2sgQiBcImh0dHBzOi8vZ2l0aHViLmNvbS9kaWV0ZXJpY2gtbGFiL3NpbmdsZS1jZWxsLW5hbm9wb3JlL2Jsb2IvbWFzdGVyL3BpcGVsaW5lcy9SRUFETUUubWQjYnVpbGRfZ2Vub21lc2hcIiIsIm1lcm1haWQiOnsidGhlbWUiOiJkZWZhdWx0In0sInVwZGF0ZUVkaXRvciI6ZmFsc2V9)](https://mermaid-js.github.io/mermaid-live-editor/#/edit/eyJjb2RlIjoiZ3JhcGggVERcbiAgICBBW0NlbGxSYW5nZXI8YnI-Y291bnRdIC0tPnxidWlsZF9nZW5vbWUuc2h8QihBcnRpZmljYWwgR2Vub21lKVxuICAgIEEgLS0-fGdldF9jYmZyZXEuc2h8QjIoQ2VsbCBCYXJjb2RlPGJyPkZyZXF1ZW5jaWVzKVxuICAgIEEgLS0-fGdldF9jYmMuc2h8QzIoQ2VsbCBCYXJjb2RlPGJyPldoaXRlbGlzdClcbiAgICBDW05hbm9wb3JlPGJyPlJlYWRzXSAtLT58YnVpbGRfbmFub3NpbS5zaHxEKE5hbm9TaW08YnI-TW9kZWxzKVxuICAgIEMyLS0-IHxmaW5kX2Rpc3QucnxDMyhFcnJvci1jb250YWluaW5nPGJyPldoaXRlbGlzdClcbiAgICBDMy0tPiBIXG4gICAgQiAtLT4gfHNpbV9yZWFkcy5zaHxFKFNpbXVsYXRlZDxicj5OYW5vcG9yZTxicj5SZWFkcylcbiAgICBEIC0tPiBFXG4gICAgRSAtLT4gfGJ1aWxkX2FsaWduLnNofEYoU2ltdWxhdGVkPGJyPk5hbm9wb3JlPGJyPkFsaWdubWVudHMpXG4gICAgQiAtLT4gfGdldF9iYXJjb2Rlcy5zaHxHKEdyb3VuZCBUcnV0aDxicj5CYXJjb2RlcylcbiAgICBFIC0tPiBHXG4gICAgRiAtLT4gfHJ1bl9waXBlLnNofEgoTGFiZWxlZDxicj5GZWF0dXJlcylcbiAgICBNIC0tPiB8cHJlZC5yfExcbiAgICBIIC0tPiBMKFByb2JhYmlsaXRpZXM8YnI-b2YgU2ltdWxhdGVkPGJyPkJhcmNvZGUgTWF0Y2hlcylcbiAgICBMIC0tPiB8c2ltX2JlbmNoLnIsIHN0YXQuc2h8TihDbGFzc2lmaWNhdGlvbjxicj5QZXJmb3JtYW5jZSlcbiAgICBIIC0tPiB8ZmVhdF9zdGF0LnJ8SShGZWF0dXJlczxicj5JbXBvcnRhbmNlKVxuICAgIEggLS0-IHxmZWF0X3Zhci5yfE8oRmVhdHVyZXM8YnI-VmFyaWFuY2UpXG4gICAgRyAtLT4gSFxuICAgIEggLS0-IHxidWlsZF9tb2RlbC5yfE0oTmFpdmUgQmF5ZXM8YnI-TW9kZWwpXG4gICAgc3R5bGUgSSBmaWxsOiNmZmNcbiAgICBzdHlsZSBOIGZpbGw6I2ZmY1xuICAgIHN0eWxlIE8gZmlsbDojZmZjXG4gICAgY2xpY2sgQiBcImh0dHBzOi8vZ2l0aHViLmNvbS9kaWV0ZXJpY2gtbGFiL3NpbmdsZS1jZWxsLW5hbm9wb3JlL2Jsb2IvbWFzdGVyL3BpcGVsaW5lcy9SRUFETUUubWQjYnVpbGRfZ2Vub21lc2hcIiIsIm1lcm1haWQiOnsidGhlbWUiOiJkZWZhdWx0In0sInVwZGF0ZUVkaXRvciI6ZmFsc2V9)

# Step-by-step instructions

## Model estimation

[1) Parameter estimation from Illumina data](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/pipelines/README.md#illumina-library)

We use the filtered cellular barcodes ([barcodes.tsv.gz](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices)) detected from Cell Ranger pipeline, and produce a ([fasta file](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/pipelines/README.md#get_cbcsh)) as barcode whitelist for the following barcode alignment. 

In addtion to the filtered cellular barcodes, we also need to add some more unfiltered cellular barcodes into our whitelist, to make sure we do not align the read to sub-optimal barcodes due to the absence of the real barcode sequences. This step is done by a [R script](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/pipelines/README.md#find_distr) which retrieves all the other barcodes within two edit-distances from the filtered cellular barcodes.

Next we need the read counts for each barcode as prior knowledge. The method is documented in the [10xgenomics website](https://kb.10xgenomics.com/hc/en-us/articles/360007068611-How-do-I-get-the-read-counts-for-each-barcode-) and our [pipeline](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/pipelines/README.md#get_cbfreqsh).

[2) Parameter estimation from Nanopore data](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/pipelines/README.md#nanopore-library)

The Nanopore reads should be aligned to the reference genome before the analyses. We recommend using Minimap2 with [long-reads settings](https://github.com/lh3/minimap2#map-long-mrnacdna-reads), also check out an example from the [pipelines](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/pipelines/README.md#align_longreadssh).

[3) Build predictive model](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/pipelines/README.md#nanopore-library)

We generated an artifical "genome" which contains only the [cDNA primer from 10x Chromium Single Cell V3](https://kb.10xgenomics.com/hc/en-us/articles/217268786-How-do-I-design-a-custom-targeted-assay-for-Single-Cell-3-), cellular barcode and UMI sequences as the same counts as the Illumina library, followed by 20bp oligo-dT and 32bp cDNA sequences in our [pipeline](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/pipelines/README.md#build_genomesh), in order to estimate the likelihood of barcode mismatches and indels in our model.

The Nanopore error profile is produced using the ["read_analysis.py"](https://github.com/bcgsc/NanoSim#1-characterization-stage) from NanoSim. The following [step](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/pipelines/README.md#build_nanosimsh) in our pipeline creates a directory containing the error profile. Next we generate [one million Nanopore reads](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/pipelines/README.md#sim_readssh) based on the artifical we built previously using NanoSim. The generated Nanopore reads were [trimmed to 100bp and wrote to a bam file](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/pipelines/README.md#build_alignsh), and the ground truth barcode sequences can be known by looking into the corresponding genomic locations from the artifical genome in the [pipeline](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/pipelines/README.md#get_barcodessh). 

In the next step, we run our feature extraction [pipeline](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/pipelines/README.md#run_pipesh) on the bam file of simulated reads. By comparing to our ground truth barcode sequences, we assign either 0 or 1 as labels to each barcode alignment, indicating whether the corresponding alignment is correct or not. Then we build a [naive bayesian model](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/pipelines/README.md#build_modelr) based on the label and previously extracted barcode alignment features. 

Using the naive bayesian model built in the last step, we can predict the likelihood given the alignment features from the other Nanopore reads sequenced from the same cDNA library. Then we use the bayesian theorem to calculated the posterior probabilities that the barcode alignment is correct among all potential barcodes. The [predicted probabilities](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/pipelines/README.md#predr) allow us benchmarking our predictions with the other simulated reads, or do barcode assignment with the real Nanopore reads. User may [set a cutoff and select the highest from these probabilities](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/pipelines/README.md#filter_predsh)

## Benchmarking
[1) Simulated Nanopore reads](https://github.com/dieterich-lab/single-cell-nanopore)

We have provided an [R script](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/pipelines/sim_bench.r) that allows user to inspect the various metrics (accuracy, precision, sensitivity, specificity, and F1) along different cut-off thresholds on the predicted probabilities. User may also want to plot the [feature variances](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/pipelines/README.md#feat_varr) and [feature importance](https://github.com/dieterich-lab/single-cell-nanopore/tree/master/pipelines#feat_statr) for a better understanding of the correct and false barcode assignments.

[2) Real Nanopore reads](https://github.com/dieterich-lab/single-cell-nanopore)


## Authors

* Qi Wang <[qwang.big@gmail.com](mailto:qwang.big@gmail.com)>

* Sven Bönigk <[sven.boenigk@fu-berlin.de](mailto:sven.boenigk@fu-berlin.de)>

* Christoph Dieterich
