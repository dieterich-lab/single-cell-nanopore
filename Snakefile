import os
import re
import sys
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

configfile: "config.yaml"

usableThreads = 8
dir_in = "data/"
dir_out = "analysis/"
ref_genome = config["reference_genome"]
ref_cdna = config["reference_cdna"]
ref_gff = config["reference_geneset"]
barcode = config["barcode"]
bam_illumina = config["illumina_bam"]
fq_nanopore = config["nanopore_fq"]
bam_nanopore = os.path.splitext(os.path.basename(fq_nanopore))[0] + '.bam'
fa_barcode = os.path.splitext(os.path.basename(barcode))[0] + '.fa'

rule all:
  input:
    dir_out + ".label"

rule get_cbc:
  input:
    barcode = dir_in + barcode
  output:
    fa_barcode = dir_out + fa_barcode
  shell:
    """
    zcat {input} | perl -ne 'print ">$1\n$1\n" if /^(\w+)/' > {output}
    """

rule get_cbfreq:
  input:
    bam = dir_in + bam_illumina
  output:
    reads_per_barcode = dir_out + "reads_per_barcode"
  shell:
    """
    samtools view {input} | perl -ne 'print "$1\n" if /GN:Z:.*CB:Z:([ACGT]+)/' | sort | uniq -c > {output}
    """

rule find_dist:
  input:
    reads_per_barcode = dir_out + "reads_per_barcode",
    fa_barcode = dir_out + fa_barcode
  output:
    barcode = dir_out + 'whitelist.fa'
  shell:
    """
    Rscript pipelines/find_dist.r {input.reads_per_barcode} {input.fa_barcode} {output}.tmp
    sort {output}.tmp|uniq|perl -ne 'print ">$_$_"' > {output}
    """

rule align_longreads:
  input:
    fq = dir_in + fq_nanopore,
    ref_genome = dir_in + ref_genome
  output:
    bam = dir_out + bam_nanopore
  shell:
    """
    minimap2 -o {output} -t 5 -ax splice --MD -ub {input.ref_genome} {input.fq}
    """

rule build_nanosim:
  input:
    fq = dir_in + fq_nanopore,
    ref_genome = dir_in + ref_genome,
    ref_cdna = dir_in + ref_cdna,
    ref_gff = dir_in + ref_cdna
  output:
    model = dir_out + "nanosim_model"
  shell:
    """
    read_analysis.py transcriptome -i {input.fq} -rg {input.ref_genome} -rt {input.ref_cdna} -annot {input.ref_gff} -t 10 -o {output}/sim
    """

rule build_genome:
  input:
    bam = dir_in + bam_illumina
  output:
    fa_sim = dir_out + "genome.fa"
  params:
    adapter = config["adapter"]
    polyTlength = config["polyTlength"]
  shell:
    """
    samtools view {input} |perl -ne '@t=split(/\t/);print ">",++$j,"\n" if $i++%25e5==0;print "{params.adapter}$3$4","T"x{params.polyTlength},substr($t[9],0,32),"\n" if /(TX|AN):Z:(\w+).*CB:Z:([ACGT]+).*UB:Z:([ACGT]+)/' > {output}
    samtools faidx {output}
    """

rule sim_reads:
  input:
    model = dir_out + "nanosim_model",
    fa_sim = dir_out + "genome.fa"
  output:
    sim = dir_out + "sim_reads.fasta"
  shell:
    """
    simulator.py genome -rg {input.fa_sim} -c {input.model}/sim -o {output} -n 1000000
    """

rule build_align:
  input:
    sim = dir_out + "sim_reads.fasta"
  output:
    fa = dir_out + "sim_test.fa",
    bam = dir_out + "sim_test.bam"
  shell:
    """
    perl pipelines/get_firstread.pl {input} > {output.fa}
    perl pipelines/fa2sam.pl {output.fa}|samtools view -bS > {output.bam}
    """

rule get_barcodes:
  input:
    sim = dir_out + "sim_reads.fasta",
    fa_sim = dir_out + "genome.fa"
  output:
    barcode = dir_out + "_barcodes.txt"
  shell:
    """
    perl -ne '$L=100;next unless /^>/;$_=substr($_,1);@t=split(/_/);$d=$t[1]+$L-$t[1]%$L;print "$t[0]\t$d\t",$d+$L,"\t$_"' {input.sim}|sort -k1,1 -k2,2n > {input.sim}.bed
    bedtools getfasta -fi {input.fa_sim} -bed {input.sim}.bed -name > {input.sim}.fa
    perl -ne 'print "$1\t" if /^>(\w+):/;print substr($_,22,16),"\n" unless /^>/' {input.sim}.fa > {output}
    """

rule run_pipe:
  input:
    barcode = dir_out + 'whitelist.fa',
    bam = dir_out + "sim_test.bam"
  output:
    tab = dir_out + "sim.tab"
  params:
    adapter = config["adapter"],
    prefix = "sim"
  threads: config["threads"]
  shell:
    """
    bin/singleCellPipe -n 20 -r {input.bam} -t {params.prefix} -w {input.barcode} -as {params.adapter} -ao 10 -ae 0.3 -ag -2 -hr T -hi 10 -he 0.3 -bo 5 -be 0.2 -bg -2 -ul 26 -kb 3 -fl 100
    awk '$2!="NA" || NR==1' {params.prefix}.tab > {output}
    rm {params.prefix}.tab {params.prefix}.fasta {params.prefix}parameterLog.log
    """

rule add_label:
  input:
    tab = dir_out + "sim.tab",
    barcode = dir_out + "_barcodes.txt"
  output:
    tab = dir_out + "sim.tab1"
  shell:
    """
    Rscript pipelines/add_label.r {input.tab} {input.barcode} {output}
    """

rule build_model:
  input:
    tab = dir_out + "sim.tab1"
  output:
    model = dir_out + ".model.rda"
  shell:
    """
    Rscript pipelines/build_model.r {input} {output}
    """

rule pred:
  input:
    model = dir_out + ".model.rda",
    tab = dir_out + "sim.tab",
    reads_per_barcode = dir_out + "reads_per_barcode"
  output:
    prob = dir_out + ".prob"
  shell:
    """
    Rscript pipelines/build_model.r {input.model} {input.tab} {input.reads_per_barcode} {output}
    """

rule filter_pred:
  input:
    prob = dir_out + ".prob"
  output:
    label = dir_out + ".label"
  shell:
    """
    awk '$15>30' {input} | sed 's/_end[1|2]//' | awk '{a[$1]++;b[$1]=$0}END{for(i in a){if(a[i]==1)print b[i]}}' | cut -f1-2 > {output}
    """
