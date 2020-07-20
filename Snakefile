import os
import re
import sys
import gzip
import scipy.io

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

def get_most_abundant_gene(mtx_file, feature_file, output_file):
    mat = (scipy.io.mmread(gzip.open(mtx_file, 'r')))
    features = gzip.open(feature_file, 'rt').read().splitlines()
    with open(output_file, 'w') as f:
        f.write(features[mat.sum(axis=1).argmax()].split('\t')[1])
    return ""

localrules: all, unzip_fq, get_cbc, build_genome, build_align

rule all:
  input:
    sim  = dir_out + "sim.label",
    real = dir_out + "real.label"

rule unzip_fq:
  input:
    file = dir_in + fq_nanopore + '.gz'
  output:
    file = dir_in + fq_nanopore
  run:
    shell("gunzip -c {input} > {output}")

rule get_a_gene:
  input:
    feature = dir_in + config["feature"],
    matrix  = dir_in + config["matrix"]
  output:
    file = dir_out + 'gene'
  run:
    get_most_abundant_gene(input.matrix, input.feature, output.file)

rule get_cbc:
  input:
    barcode = dir_in + barcode
  output:
    fa_barcode = dir_out + fa_barcode
  shell:
    """
    gzip -dc {input} | perl -ne 'print ">$1\\n$1\\n" if /^(\\w+)/' > {output}
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

rule get_gene_umi:
  input:
    bam = dir_in + bam_illumina
  output:
    umi = dir_out + "umis/umi"
  shell:
    """
    if [ ! -d {dir_out}umis ]; then mkdir {dir_out}umis; fi
    samtools view {input} | perl -ne 'print "$1\\t$2\\t$3\\n" if /GN:Z:(\S+).*CB:Z:([ACGT]+).*UB:Z:([ACGT]+)/' | sort -k1,1 | perl -F"\\t" -ane 'if ($F[0] ne $x){{close O;open O, ">{dir_out}umis/$F[0].txt"}}print O $F[1],"\\t",$F[2];$x=$F[0];END{{close O}}' > {output}
    for f in {dir_out}umis/*.txt;do sort $f|uniq|perl -F"\\t" -ane 'print ">$F[0]\\n$F[1]"' > ${{f%.*}}.fa;done
    rm {dir_out}umis/*.txt
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
    num = config["numIlmReads"]
  shell:
    """
    python pipelines/simreads.py {input.barcode} {input.feature} {input.matrix} {params.num} {params.umilength} {params.cdnaseq} | samtools view -bS > {output}
    """

rule build_genome:
  input:
    bam = dir_in + bam_illumina,
    bk_barcode = dir_out + "bk_barcodes.tsv.gz"
  output:
    fa_sim = dir_out + "genome.fa"
  params:
    adapter = config["adapter"],
    cdnalength = config["cdnalength"],
    umilength = config["umilength"],
    cdnaseq = config["cdnaseq"][:config["cdnalength"]],
    polyTlength = config["polyTlength"],
    num = int(config["numIlmReads"]*config["percent_raw"])
  shell:
    """
    samtools view {input.bam} | perl -ne 'print ">",++$j,"\\n" if $i%25e5==25e5-1 or $j==0;if (/CB:Z:([ACGT]+).*UB:Z:([ACGT]+)/){{print "{params.adapter}$1$2","T"x{params.polyTlength},{params.cdnaseq},"\\n";$i++}}' > {output}
    gzip -dc {input.bk_barcode} | shuf -r -n {params.num} | perl -ne 'print ">chr",++$j,"\\n" if $i%25e5==25e5-1 or $j==0;if (/([ACGT]+)/){{$u="";$u.=[A,T,G,C]->[rand 4]for 1..{params.umilength};print "{params.adapter}$1$u","T"x{params.polyTlength},{params.cdnaseq},"\\n";$i++}}' >> {output}
    samtools faidx {output}
    """

rule bk_barcodes:
  input:
    barcode = dir_in + config["barcode"],
    raw_barcode = dir_in + config["barcode_raw"]
  output:
    bk_barcode = dir_out + "bk_barcodes.tsv.gz"
  shell:
    """
    awk 'NR==FNR{{a[$0];next}}!($0 in a)' <(gzip -dc {input.barcode}) <(gzip -dc {input.raw_barcode}) | gzip > {output.bk_barcode}
    """

rule sim_reads:
  input:
    model = dir_out + "nanosim_model/sim_model_profile",
    fa_sim = dir_out + "genome.fa"
  output:
    sim = dir_out + "sim_reads.fasta"
  params:
    num = config["numSimReads"],
    num2 = int(config["numSimReads"]*1.35)
  shell:
    """
    simulator.py genome -rg {input.fa_sim} -c {dir_out}nanosim_model/sim -o {dir_out}sim -n {params.num2}
    perl -ne '$i++ if /^>/;print if $i<={params.num}' {dir_out}sim_aligned_reads.fasta > {output}
    """

rule build_test:
  input:
    sim = dir_out + "sim_reads.fasta"
  output:
    fa = dir_out + "sim_test.fa"
  params:
    readlen = simreadlength,
    cdnalength = config["cdnalength"]
  shell:
    """
    perl -ne '$L={params.readlen};chomp;if (/^>/){{$id=$_}}else{{@t=split(/_/,$id);$s=$_;$d=$t[5];if ($t[4] eq "R"){{$s=reverse $s;$s=~tr/ATGCatgc/TACGtacg/;$d=$t[7]}}$s=substr($s,$d,$t[6]);print $id,"\\n",$t[6]>$L?substr($s,$L-$t[1]%$L,$L):$s,"\\n"}}' {input} > {output.fa}
    """

rule build_align:
  input:
    sim = dir_out + "sim_test.fa",
    barcode = dir_out + "sim_barcodes.txt"
  output:
    bam = dir_out + "sim_test.bam"
  params:
    readlen = simreadlength,
    cdnalength = config["cdnalength"]
  shell:
    """
    perl pipelines/fa2sam.pl {input.sim} {input.barcode} {params.readlen} {params.cdnalength} | samtools view -bS > {output}
    """

rule get_sim_umi_reads:
  input:
    umi = dir_out + "sim_umi.fasta",
    genes = dir_out + "sim_genes.txt"
  output:
    bam = dir_out + "umis/bam"
  params:
    adapter = 'AAA'
  shell:
    """
    perl pipelines/fa2sams.pl {input.umi} {input.genes} {params.adapter} {dir_out}umis > {output}
    for f in {dir_out}umis/*.txt;do samtools view -bS $f > ${{f%.*}}.bam;done
    rm {dir_out}umis/*.txt
    """

rule get_barcodes:
  input:
    bam = dir_in + bam_illumina,
    sim = dir_out + "sim_reads.fasta",
    fa_sim = dir_out + "genome.fa",
    gene = dir_out + 'gene'
  output:
    genes = dir_out + "sim_genes.txt",
    barcode = dir_out + "sim_barcodes.txt"
  params:
    readlen = simreadlength,
    adapterlength = len(config["adapter"]),
    barcodelength = config["barcodelength"],
    umilength = config["umilength"]
  shell:
    """
    perl -ne '$L={params.readlen};next unless /^>/;$_=substr($_,1);@t=split(/_/);$d=$t[1]+$L-$t[1]%$L;print "$t[0]\\t$d\\t",$d+$L,"\\t$_"' {input.sim}|sort -k1,1 -k2,2n > {input.sim}.bed
    bedtools getfasta -fi {input.fa_sim} -bed {input.sim}.bed -name > {input.sim}.fa
    samtools view {input.bam} | perl pipelines/gtruth.pl {input.sim}.fa {params.adapterlength} {params.barcodelength} {params.umilength} {input.gene} > {output.barcode}
    cut -f1-2 {output.barcode} > {output.genes}
    """

rule run_pipe_sim:
  input:
    barcode = dir_out + 'whitelist.fa',
    bam = dir_out + "sim_test.bam"
  output:
    tab = dir_out + "sim.tab",
    umi = dir_out + "sim_umi.fasta"
  params:
    adapter = config["adapter"],
    prefix = "sim"
  shell:
    """
    if [ -f umi.fasta ]; then rm umi.fasta; fi
    bin/singleCellPipe -n {threads} -r {input.bam} -t {params.prefix} -w {input.barcode} -as {params.adapter} -ao 10 -ae 0.3 -ag -2 -hr T -hi 10 -he 0.3 -bo 5 -be 0.2 -bg -2 -ul 26 -kb 3 -fl 100
    awk '$2!="NA" || NR==1' {params.prefix}.tab > {output.tab}
    perl -ne 'if(/^>/){{print}}else{{chomp;print substr($_,3,18),"\\n"}}' umi.fasta > {output.umi}
    rm {params.prefix}.tab {params.prefix}.fasta {params.prefix}parameterLog.log umi.fasta
    """

rule run_pipe_real:
  input:
    barcode = dir_out + 'whitelist.fa',
    bam = dir_out + _nanopore + '.bam'
  output:
    tab = dir_out + "real.tab",
    log = dir_out + "real.log",
    umi = dir_out + "real_umi.fasta"
  params:
    adapter = config["adapter"],
    barumilength = config["barcodelength"] + config["umilength"],
    prefix = "real"
  shell:
    """
    if [ -f umi.fasta ]; then rm umi.fasta; fi
    bin/singleCellPipe -n {threads} -r {input.bam} -t {params.prefix} -w {input.barcode} -as {params.adapter} -ao 10 -ae 0.3 -ag -2 -hr T -hi 10 -he 0.3 -bo 5 -be 0.2 -bg -2 -ul {params.barumilength} -kb 3 -fl 100
    awk '$2!="NA" || NR==1' {params.prefix}.tab > {output.tab}
    printf "Aligned to genome\t" > {output.log}
    cut -f1 {params.prefix}.tab|perl -npe 's/_end[1|2]//'|sort|uniq|wc -l >> {output.log}
    printf "Aligned to adapter\t" >> {output.log}
    awk '$3=="no"' {params.prefix}.tab|cut -f1|perl -npe 's/_end[1|2]//'|sort|uniq|wc -l >> {output.log}
    printf "Aligned to barcode\t" >> {output.log}
    awk '$2!="NA"' {params.prefix}.tab|cut -f1|perl -npe 's/_end[1|2]//'|sort|uniq|wc -l >> {output.log}
    perl -ne 'if(/^>/){{print}}else{{chomp;print substr($_,3,18),"\\n"}}' umi.fasta > {output.umi}
    rm {params.prefix}.tab {params.prefix}.fasta {params.prefix}parameterLog.log umi.fasta
    """

rule run_umi_sim:
  input:
    bam = dir_out + "umis/bam",
    umi = dir_out + "umis/umi"
  output:
    tab = dir_out + "sim_umi.tab"
  params:
    adapter = 'AAA'
  shell:
    """
    for f in $(ls -1 analysis/umis/*.bam);do bin/singleCellPipe -n {threads} -r $f -t $(basename ${{f%.*}}) -w ${{f%.*}}.fa -as {params.adapter} -ao 3 -ae 0 -hr T -hi 10 -he 0.3 -bo 5 -be 0.2 -bg -2 -ul 26 -kb 3 -fl 100;done
    cat *.tab|grep -v '^read_id'|perl -F"\\t|\\.\\." -ane 'print if $F[1] eq $F[2]' > {output}
    rm *.tab *.fasta *parameterLog.log
    """

rule add_umiflag:
  input:
    tab = dir_out + "sim.tab",
    umi = dir_out + "sim_umi.tab"
  output:
    tab = dir_out + "sim.tab0"
  shell:
    """
    perl pipelines/add_umi.pl {input.tab} {input.umi}|cut -f1-13 > {output}
    """

rule add_label:
  input:
    tab = dir_out + "sim.tab0",
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

rule sim_pred:
  input:
    model = dir_out + "sim.model.rda",
    tab = dir_out + "sim.tab1",
    reads_per_barcode = dir_out + "reads_per_barcode"
  output:
    prob = dir_out + "sim.prob"
  shell:
    """
    Rscript pipelines/pred.r {input.model} {input.tab} {input.reads_per_barcode} {output}
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

rule filter_sim:
  input:
    prob = dir_out + "sim.prob"
  output:
    label = dir_out + "sim.label"
  params:
    cutoff = config["cutoff"]
  shell:
    """
    awk 'NR>1' {input} | sed 's/_end[1|2]//' | cut -f1-2,15 | sort -k3,3rn | awk '!_[$1]++' > {output}
    """

rule filter_pred:
  input:
    prob = dir_out + "real.prob"
  output:
    label = dir_out + "real.label"
  params:
    cutoff = config["cutoff"]
  shell:
    """
    awk 'NR>1 && $15>{params.cutoff}' {input} | sed 's/_end[1|2]//' | awk '{{a[$1]++;b[$1]=$0}}END{{for(i in a){{if(a[i]==1)print b[i]}}}}' | cut -f1-2,15 > {output}
    """

rule report:
  input:
    slabel= dir_out + "sim.label",
    rlabel= dir_out + "real.label",
    sprob = dir_out + "sim.prob",
    rprob = dir_out + "real.prob",
    rlog  = dir_out + "real.log",
    barcode = dir_out + "sim_barcodes.txt"
  output:
    file = dir_out + "report.pdf"
  params:
    cutoff = config["cutoff"],
  shell:
    """
    Rscript -e 'rmarkdown::render("report.rmd",output_file="{output}")' {input.barcode} {input.slabel} {input.sprob} {input.rlabel} {input.rprob} {input.rlog} {params.cutoff} null
    """
