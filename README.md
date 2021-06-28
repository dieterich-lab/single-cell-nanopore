# ScNapBar

ScNapBar is designed for cell barcode assignment from Nanopore sequencing data.
It requires bam files from both Nanopore and Illumina reads, then builds a model based on the parameters estimated from the two libraries.

## Quick run

1. *optional* edit the provided **`config.yaml`** file to match your own sequence files, reference genome and annotation. Update the adapter and polyT length that fit your libraries.

2. Run the **`snakemake`** command under the conda environment to perform the bioinformatics analysis on the specified sequence files. Several analysis steps can benefit from multiple computer cores; use the `-j` parameter to parallelize the analysis (this document is assuming that 8 cores are available).
```
snakemake -j 8 all
```
You can also submit the job via job schedulers. We have provided an example on the SLURM. Change the account settings in cluster.json before using it.
```
snakemake -j 8 --cluster-config cluster.json --cluster "sbatch -A {cluster.account} --mem={cluster.mem} -t {cluster.time} -c {cluster.threads}"
```

## Two modes for the cell barcode assignment

1. ScNapBar `(option 1)` uses a probabilistic model for barcode assignment by default (see command below), which performs very well in cases of low sequencing saturation. 

```
snakemake -j 8
```

2. ScNapBar `(option 2)` assigns barcode based on matched Illumina UMIs without additional probabilistic modeling. Use the following command to start this mode: 

```
snakemake -j 8 --until run_umi_seq
```

## Installation:

1. Most software dependencies are managed using **`conda`**. Please install as described at  <br> [https://docs.conda.io/projects/conda/en/latest/user-guide/install/](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
bash
```
2. Download the files into a folder named `single-cell-nanopore`. The program requires the **`NanoSim`** and **`tbb`** which should be installed first through **`conda`**
```
conda install -c bioconda nanosim
conda install -c conda-forge tbb tbb-devel 
git clone https://github.com/nanoporetech/single-cell-nanopore.git single-cell-nanopore
```
3. Change working directory into the new `single-cell-nanopore` folder 
```
cd single-cell-nanopore
```
4. Put **`seqan`** library files into the empty `seqan` folder
```
wget https://github.com/seqan/seqan/releases/download/seqan-v2.4.0/seqan-library-2.4.0.tar.xz
tar xJf seqan-library-2.4.0.tar.xz
cp -r seqan-library-2.4.0/include single-cell-nanopore/seqan/
```

5. Use these commands for building **`ScNapBar`**: 
```
cmake .
make
```
6. Install conda software dependencies with
```
conda env create --name single-cell-nanopore --file environment.yaml
```
7. Initialise conda environment with 
```
conda activate single-cell-nanopore
```

## Input files:

For testing purposes, we suggest downloading the [GSE130708](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130708) dataset or [PRJNA722142](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA722142)(which we have already extracted _chr17_ as the demo dataset in the `data` folder), and rename the files according to the following naming convention.

* [`possorted_genome_bam.bam`](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/overview) from Cell Ranger, **`or`** [`barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz`](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices) in the `filtered_feature_bc_matrix` folder of Cell Ranger.

* Background cell barcodes named as [`raw.tsv.gz`](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/data/), could be renamed from the `barcodes.tsv.gz` in the `raw_feature_bc_matrix` folder of Cell Ranger. 

* Nanopore reads ([sample](https://github.com/dieterich-lab/single-cell-nanopore/blob/master/data/)) in FASTQ format

* Reference genome and transcriptome in FASTA format

* Annotation file in [`refFlat`](http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/refFlat.txt.gz) formats and `.refFlat` filename extension.

## Output files:

The output files are put in the `analysis` folder of the pipeline. The main target output files are `sim.label` and `real.label` (three columns: read_id, barcode, score). In addition, you will find the following files:

* `sim.label`: the barcode assignment from the simulated reads with scores. The scores range from 0-99, and larger scores indicate higher confidence for the assignment. Reads assigned to multiple barcodes only have the assignment with the highest score retained. 

* `real.label`: the barcode assignment from the real Nanopore reads with scores. The scores range from the cutoff set in `config.yaml` to 99. Reads assigned to multiple barcodes are removed if both are above the score cutoff. 

* `real.umi`: The barcodes assignment of the Nanopore reads with matched Illumina UMIs from the same cell and the same gene.

* `sim.prob` and `real.prob`: The complete feature tables used to generate the probability scores of the simulated reads and the real Nanopore reads, respectively. 

Other intermediate files which contain the useful information in the model estimation and benchmarking includes:

* `sim_barcodes.txt`: The ground-truth of cell barcodes in the simulated reads. 

* `sim.model.rda`: The naive bayes model trained from the simulated reads. 

* `genome.fa`: The artificial genome inputted into NanoSim for generating the simulated reads. 

* `Nanopore.bam`: The mapping of the real Nanopore reads generated by minimap2.

* `sim_umi.fasta` and `real_umi.fasta`: The rest DNA sequences after removing the adapter and barcode sequences. 

* `sim_umi.tab`: The barcodes assignment of the simulated reads with matched Illumina UMIs from the same cell and the same gene.

## Parameters:

### Parameters of **`ScNapBar`**

The parameters are set in the `config.yaml`. If the entry is a file, then it must be placed under the `data` folder. 

* `reference_genome`: The reference genome sequences in FASTA format. 

* `nanopore_fq`: Nanopore reads you want to process in FASTQ format. 

* `adapter`: 10x genomics P1 adapter sequences. 
 
* `polyTlength`: Number of poly-Ts you want to simulate. 

* `cdnalength`: Number of nucleotides used to append to each entry in the artifical genome. 

* `umilength`: Number of nucleotides of the UMI sequences. 

* `barcodelength`: Number of nucleotides of the cell barcode sequences. 

* `numSimReads`: Number of Nanopore reads to simulate. 

* `numSimReads`: Number of Illumina reads to sample from the Illumina sequencing. 

* `cutoff`: Score cutoff of the barcode assignment of the real Nanopore reads. 

* `percent_raw`: A fraction number representing the percentage of additional simulated reads you want to use as true negatives. These reads contain the cell barcodes from the background rather than the whitelist. From our experience, there are about 20% reads do not contain the cell barcodes from the whitelist in the 10x genomics library. 

* `threads`: Number of CPUs for the multiple-threaded jobs. 

* `cdnaseq`: DNA sequences used to append to each entry in the artificial genome. 

### Parameters of **`singleCellPipe`**

singleCellPipe is modified from [flexbar](https://github.com/seqan/flexbar) which performs the adapter and barcode alignments. Most of the flexbar parameters are available in this program, yet there are a few distinctive parameters as follows:

* `-ul`: Number of nucleotides of barcode and UMI sequences. E.g., the parameter should be 26 for a 16bp barcode and 10bp UMI library. 

* `-kb`: Number of additional nucleotides to search after the adapter alignment. 

* `-fl`: Number of nucleotides from both ends for searching the barcode sequences. 

## Step-by-step instructions

In this paragraph, we explain the use of each major job in the pipeline. 

### Parameter estimation from Illumina data

* `build_illumina`: If Illumina bam file is not provided, we use the filtered cellular barcodes ([barcodes.tsv.gz](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices)) detected from Cell Ranger pipeline, and produce an Illumina bam file based on the frequencies in the matrices.  

* `find_dist`: add some more background cellular barcodes into the cell barcode whitelist, to make sure we do not align the read to sub-optimal barcodes due to the absence of the real barcode sequences. It retrieves all the other barcodes within two edit-distances from the filtered cellular barcodes.

* `get_cbfreq`: use the read counts for each barcode as prior knowledge in the Bayesian model. 

* `align_longreads`: We align the Nanopore reads to the reference genome using Minimap2 with [long-reads settings](https://github.com/lh3/minimap2#map-long-mrnacdna-reads).

### Build predictive model

* `build_genome`: We generated an artificial "genome" which contains only the [cDNA primer from 10x Chromium Single Cell V3](https://kb.10xgenomics.com/hc/en-us/articles/217268786-How-do-I-design-a-custom-targeted-assay-for-Single-Cell-3-), cellular barcode and UMI sequences as the same counts as the Illumina library, followed by 20bp oligo-dT and 32bp cDNA sequences in our pipeline, in order to estimate the likelihood of barcode mismatches and indels in our model.

* `build_nanosim`: The Nanopore error profile is produced using the ["read_analysis.py"](https://github.com/bcgsc/NanoSim#1-characterization-stage) from NanoSim, and creates a directory under the `analysis` folder that contains the error profile of the Nanopore reads. 

* `sim_reads`: We generate a number of Nanopore reads based on the artificial we built previously using NanoSim. 

* `build_test`: The generated Nanopore reads were trimmed to 100bp and wrote to a bam file, and the ground truth barcode sequences can be known by looking into the corresponding genomic locations from the artificial genome in the pipeline. 

* `run_pipe_sim` and `run_pipe_real`: run the feature extraction pipeline on the bam file of the simulated and the real Nanopore reads, respectively. 

* `add_label`: By comparing to our ground truth barcode sequences, we assign either 0 or 1 as labels to each barcode alignment, indicating whether the corresponding alignment is correct or not. 

* `build_model`: build a naive bayesian model based on the label and previously extracted barcode alignment features. 

* `pred`: Use the naive bayesian model previously built and predict the likelihood given the alignment features from the other Nanopore reads sequenced from the same cDNA library. Then we use the Bayesian theorem to calculate the posterior probabilities that the barcode alignment is correct among all potential barcodes. The predicted probabilities allow benchmarking our predictions with the other simulated reads, or do barcode assignment with the real Nanopore reads. 

* `get_gene_umi`: create UMI whitelist for each gene based on the Illumina data. 

* `run_umi_sim`: perform UMI alignment against the UMI whitelist of the same gene on the sequences that have the barcode trimmed.

* `filter_sim` and `filter_pred`: output the reads which passed the cutoff of the predicted probabilities. 

## Installation troubleshooting

Q: `fatal error: tbb/pipeline.h: No such file or directory` when compiling **`singleCellPipe`**.

A: Please run `conda install tbb tbb-devel` to install the required **`TBB`** library. 

Q: `fatal error: seqan/basic.h: No such file or directory` when compiling **`singleCellPipe`**.

A: Please download [SeqAn](https://github.com/seqan/seqan/releases/download/seqan-v2.4.0/seqan-library-2.4.0.tar.xz) first and move the **`SeqAn`** include folder to **seqan**:

## Authors

* Qi Wang <[qwang.big@gmail.com](mailto:qwang.big@gmail.com)>

* Sven Bönigk <[sven.boenigk@fu-berlin.de](mailto:sven.boenigk@fu-berlin.de)>

* Christoph Dieterich
