import os
import re
import sys
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

configfile: "config.yaml"

dir_in = "data/"
dir_out = "analysis/"
ref_genome = config["reference_genome"]
barcode = config["barcode"]
bam_illumina = config["illumina_bam"]
fq_nanopore = config["nanopore_fq"]
simreadlength = len(config["adapter"]) + config["barcodelength"] + config["umilength"] + config["polyTlength"] + config["cdnalength"]
_nanopore = os.path.splitext(os.path.basename(fq_nanopore))[0]
fa_barcode = os.path.splitext(os.path.basename(barcode))[0] + '.fa'

localrules: all, unzip_fq, get_cbc, build_genome, build_align

rule all:
  input:
    dir_out + "real.label"

rule unzip_fq:
  input:
    file = dir_in + fq_nanopore + '.gz'
  output:
    file = dir_in + fq_nanopore
  run:
    shell("gunzip -c {input} > {output}")

rule get_cbc:
  input:
    barcode = dir_in + barcode
  output:
    fa_barcode = dir_out + fa_barcode
  shell:
    """
    zcat {input} | perl -ne 'print ">$1\\n$1\\n" if /^(\\w+)/' > {output}
    """

rule get_cbfreq:
  input:
    bam = dir_in + bam_illumina
  output:
    reads_per_barcode = dir_out + "reads_per_barcode"
  shell:
    """
    samtools view {input} | perl -ne 'print "$1\\n" if /GN:Z:.*CB:Z:([ACGT]+)/' | sort | uniq -c > {output}
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
    sam = dir_out + _nanopore + '.sam',
    bam = dir_out + _nanopore + '.bam'
  shell:
    """
    minimap2 -v1 -t {threads} -ax splice --MD -ub {input.ref_genome} {input.fq} > {output.sam}.tmp
    grep '^@' {output.sam}.tmp |sort|uniq > {output.sam}.head
    grep -v '^@' {output.sam}.tmp |sort -snk3 -k4|uniq > {output.sam}.body
    cat {output.sam}.head {output.sam}.body > {output.sam}
    samtools view -bS -o {output.bam} {output.sam}
    rm {output.sam}.head {output.sam}.body {output.sam}.tmp
    """

rule build_nanosim:
  input:
    fq = dir_in + fq_nanopore,
    genome_alignment = dir_out + _nanopore + '.sam'
  output:
    model = dir_out + "nanosim_model/sim_model_profile"
  shell:
    """
    read_analysis.py genome -i {input.fq} -ga {input.genome_alignment} -t {threads} -o {dir_out}nanosim_model/sim
    """

rule build_illumina:
  input:
    barcode = dir_in + config["barcode"],
    feature = dir_in + config["feature"],
    matrix  = dir_in + config["matrix"],
  output:
    bam = dir_in + bam_illumina
  params:
    cdnaseq = config["cdnaseq"],
    umilength = config["umilength"],
    nreads = 1000000
  shell:
    """
    python pipelines/simreads.py {input.barcode} {input.feature} {input.matrix} {params.nreads} {params.umilength} {params.cdnaseq} | samtools view -bS > {output}
    """

rule build_genome:
  input:
    bam = dir_in + bam_illumina
  output:
    fa_sim = dir_out + "genome.fa"
  params:
    adapter = config["adapter"],
    cdnalength = config["cdnalength"],
    cdnaseq = config["cdnaseq"][:config["cdnalength"]],
    polyTlength = config["polyTlength"]
  shell:
    """
    samtools view {input} | perl -ne '@t=split(/\\t/);print ">",++$j,"\\n" if $i%25e5==25e5-1 or $j==0;if (/GN:Z:(\\w+).*CB:Z:([ACGT]+).*UB:Z:([ACGT]+)/){{print "{params.adapter}$2$3","T"x{params.polyTlength},{params.cdnaseq},"\\n";$i++}}' > {output}
    samtools faidx {output}
    """

rule sim_reads:
  input:
    model = dir_out + "nanosim_model/sim_model_profile",
    fa_sim = dir_out + "genome.fa"
  output:
    sim = dir_out + "sim_reads.fasta"
  params:
    num = config["numSimReads"],
  shell:
    """
    simulator.py genome -rg {input.fa_sim} -c {dir_out}nanosim_model/sim -o {dir_out}sim -n {params.num}
    cat {dir_out}sim_aligned_reads.fasta {dir_out}sim_unaligned_reads.fasta > {output}
    """

rule build_align:
  input:
    sim = dir_out + "sim_reads.fasta"
  output:
    fa = dir_out + "sim_test.fa",
    fq = dir_out + "sim_test.fq",
    bam = dir_out + "sim_test.bam"
  params:
    readlen = simreadlength,
    cdnalength = config["cdnalength"]
  shell:
    """
    perl -ne '$L={params.readlen};chomp;if (/^>/){{$id=$_}}else{{@t=split(/_/,$id);$s=$_;$d=$t[5];if ($t[4] eq "R"){{$s=reverse $s;$s=~tr/ATGCatgc/TACGtacg/;$d=$t[7]}}$s=substr($s,$d,$t[6]);print $id,"\\n",$t[6]>$L?substr($s,$L-$t[1]%$L,$L):$s,"\\n"}}' {input} > {output.fa}
    perl -ne 'chomp;print "\\@SQ\\tSN:15\\tLN:101991189\\n" unless defined $b;$b=/^>/;$l=length($_);$q="I"x$l;$s=substr($_,1) if $b;$l2=$l-{params.cdnalength};print "$s\\t0\\t15\\t69452800\\t31\\t${{l2}}S{params.cdnalength}M\\t*\\t0\\t0\\t$_\\t$q\\n" if $l=={params.readlen} and !$b' {output.fa} | samtools view -bS > {output.bam}
    perl -ne 'if(/^>/){{print "@",substr($_,1)}}else{{chomp;print "$_\\n+\\n","I"x length($_),"\\n"}}' {output.fa} > {output.fq}
    """

rule get_barcodes:
  input:
    sim = dir_out + "sim_reads.fasta",
    fa_sim = dir_out + "genome.fa"
  output:
    barcode = dir_out + "sim_barcodes.txt",
    umi = dir_out + "sim_umi.txt"
  params:
    readlen = simreadlength,
    adapterlength = len(config["adapter"]),
    umilength = config["umilength"],
    barcodelength = config["barcodelength"]
  shell:
    """
    perl -ne '$L={params.readlen};next unless /^>/;$_=substr($_,1);@t=split(/_/);$d=$t[1]+$L-$t[1]%$L;print "$t[0]\\t$d\\t",$d+$L,"\\t$_"' {input.sim}|sort -k1,1 -k2,2n > {input.sim}.bed
    bedtools getfasta -fi {input.fa_sim} -bed {input.sim}.bed -name > {input.sim}.fa
    perl -ne 'print "$1\\t" if /^>(\\w+):/;print substr($_,{params.adapterlength},{params.barcodelength}),"\\n" unless /^>/' {input.sim}.fa > {output.barcode}
    perl -ne 'print "$1\\t" if /^>(\\w+):/;print substr($_,{params.adapterlength}+{params.barcodelength},{params.umilength}),"\\n" unless /^>/' {input.sim}.fa > {output.umi}
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
  shell:
    """
    bin/singleCellPipe -n {threads} -r {input.bam} -t {params.prefix} -w {input.barcode} -as {params.adapter} -ao 10 -ae 0.3 -ag -2 -hr T -hi 10 -he 0.3 -bo 5 -be 0.2 -bg -2 -ul 26 -kb 3 -fl 100
    awk '$2!="NA" || NR==1' {params.prefix}.tab > {output}
    rm {params.prefix}.tab {params.prefix}.fasta {params.prefix}parameterLog.log
    """

rule run_pipe2:
  input:
    barcode = dir_out + 'whitelist.fa',
    bam = dir_out + _nanopore + '.bam'
  output:
    tab = dir_out + "real.tab"
  params:
    adapter = config["adapter"],
    prefix = "real"
  shell:
    """
    bin/singleCellPipe -n {threads} -r {input.bam} -t {params.prefix} -w {input.barcode} -as {params.adapter} -ao 10 -ae 0.3 -ag -2 -hr T -hi 10 -he 0.3 -bo 5 -be 0.2 -bg -2 -ul 26 -kb 3 -fl 100
    awk '$2!="NA" || NR==1' {params.prefix}.tab > {output}
    rm {params.prefix}.tab {params.prefix}.fasta {params.prefix}parameterLog.log
    """

rule add_label:
  input:
    tab = dir_out + "sim.tab",
    barcode = dir_out + "sim_barcodes.txt"
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
    model = dir_out + "sim.model.rda"
  shell:
    """
    Rscript pipelines/build_model.r {input} {output}
    """

rule run_pred:
  input:
    model = dir_out + "sim.model.rda",
    tab = dir_out + "real.tab",
    reads_per_barcode = dir_out + "reads_per_barcode"
  output:
    prob = dir_out + "real.prob"
  shell:
    """
    Rscript pipelines/pred.r {input.model} {input.tab} {input.reads_per_barcode} {output}
    """

rule filter_pred:
  input:
    prob = dir_out + "real.prob"
  output:
    label = dir_out + "real.label"
  params:
    cutoff = config["cutoff"],
  shell:
    """
    awk 'NR>1 && $15>{params.cutoff}' {input} | sed 's/_end[1|2]//' | awk '{{a[$1]++;b[$1]=$0}}END{{for(i in a){{if(a[i]==1)print b[i]}}}}' | cut -f1-2 > {output}
    """
