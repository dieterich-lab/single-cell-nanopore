# Overview
## var.sh Define bash variables
```
export samtools=`which samtools`
export minimap2="$HOME/bin/minimap2-master/minimap2"
export single_cell_pipe="$HOME/bin/singleCellPipe"
export nanosim_path="$HOME/bin/NanoSim-2.4-beta/src"
export ref_genome="GRCh38_90.fa"
export ref_cdna="Homo_sapiens.GRCh38.cdna.all.fa"
export ref_gff="Homo_sapiens.GRCh38.90.gff3"
export input_nanosim="Nanopore.fq"
export output_sim="sim"
export nanosim_model="NanosimModel"
export sim_genome="genome.fa"
```
## main.sh
```
sh src/get_cbc.sh
sh src/get_cbfreq.sh
Rscript src/find_dist.r
sort MoreBarcodes.txt|uniq|perl -ne 'print ">$_$_"' > MoreBarcodes.fasta
sh src/align_longreads.sh
mkdir $nanosim_model 
sh src/build_nanosim.sh
sh src/build_genome.sh
sh src/sim_reads.sh
sh src/build_align.sh
sh src/get_barcodes.sh
sh src/run_pipe.sh
Rscript src/add_label.r
Rscript src/build_model.r ${output_sim}_model.rda ${output_sim}.tab1
Rscript src/pred.r ${output_sim}_model.rda ${output_sim}.tab1
sh src/filter_pred.sh ${output_sim}.tab1.prob
```
## install_packages.r
```
install.packages('stringdist')
install.packages('caret')
install.packages('e1071')
```
# Illumina library
## get_cbc.sh
```
zcat barcodes.tsv.gz | perl -ne 'print ">$1\n$1\n" if /^(\w+)/' > CellBarcodes.fasta
```
## get_gfpcbc.sh
```
$samtools view possorted_genome_bam.bam pcDNA5:1-6000 | perl -ne 'print "$1\n" if /CB:Z:([ACGT]+)/' > gfp.txt
```
## find_dist.r
```
library(stringdist)
options(stringsAsFactors = FALSE)
x=read.table('reads_per_barcode')
y=read.table('CellBarcodes.fasta')
r=unlist(lapply(y[seq(2,nrow(y),2),1],function(d) x[stringdist(d,x[,2],method='dl')<2,2]))
write.table(r,file='MoreBarcodes.txt',sep="\t",quote=F,row.names=F,col.names=F)
```
## get_cbfreq.sh
```
$samtools view possorted_genome_bam.bam | grep GN:Z:.*CB:Z: | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > reads_per_barcode
```
# Nanopore library
## align_longreads.sh
```
$minimap2 -o Nanopore.bam -t 5 -ax splice --MD -ub $ref_genome $input_nanosim
```
## assign_genenames.sh
```
java -jar -Xmx64g Jar/Sicelore-1.0.jar AddGeneNameTag I=FC1.bam O=FC1.GE.bam REFFLAT=refFlat.txt GENETAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT
$samtools view FC1.GE.bam|perl -ne 'print "$1\t$2\n" if /^(\S+).*GE:Z:(\S+)/' > FC1.ge
```
# Build artifical genome
## build_nanosim.sh
```
$nanosim_path/read_analysis.py transcriptome -i $input_nanosim -rg $ref_genome -rt $ref_cdna -annot $ref_gff -t 10 -o $nanosim_model/sim
```
## build_genome.sh
```
$samtools view possorted_genome_bam.bam |perl -ne '@t=split(/\t/);print ">",++$j,"\n" if $i++%25e5==0;print "CTACACGACGCTCTTCCGATCT$3$4TTTTTTTTTTTTTTTTTTTT",substr($t[9],0,32),"\n" if /(TX|AN):Z:(\w+).*CB:Z:([ACGT]+).*UB:Z:([ACGT]+)/' > $sim_genome
$samtools faidx $sim_genome 
```
## sim_reads.sh
```
$nanosim_path/simulator.py genome -rg $sim_genome -c $nanosim_model/sim -o $output_sim -n 1000000
```
## fa2sam.pl
```
print '@SQ	SN:1	LN:100',"\n";
while(<>){chomp;
$b=/^>/;
$l=length($_);
$q='I'x$l;
print substr($_,1),"\t" if $b;
print "0\t1\t21\t31\t${l}S\t*\t0\t0\t$_\t$q\n" unless $b;
}
```
## build_align.sh
```
perl -ne '$L=100;chomp;if (/^>/){$id=$_}else{@t=split(/_/,$id);$s=$_;$d=$t[5];if ($t[4] eq 'R'){$s=reverse $s;$s=~tr/ATGCatgc/TACGtacg/;$d=$t[7]}$s=substr($s,$d,$t[6]);print $id,"\n",$t[6]>$L?substr($s,$L-$t[1]%$L,68):$s,"\n"}' ${output_sim}_reads.fasta > ${output_sim}_test.fa
perl src/fa2sam.pl ${output_sim}_test.fa|$samtools view -bS > ${output_sim}.bam
```
## get_barcodes.sh
```
perl -ne '$L=100;next unless /^>/;$_=substr($_,1);@t=split(/_/);$d=$t[1]+$L-$t[1]%$L;print "$t[0]\t$d\t",$d+$L,"\t$_"' ${output_sim}_reads.fasta|sort -k1,1 -k2,2n > ${output_sim}_real.bed
bedtools getfasta -fi $sim_genome -bed ${output_sim}_real.bed -name > ${output_sim}_real.fa
perl -ne 'print "$1\t" if /^>(\w+):/;print substr($_,22,16),"\n" unless /^>/' ${output_sim}_real.fa > ${output_sim}_barcodes.txt
```
# Features extraction and build model
## run_pipe.sh
```
$single_cell_pipe -n 20 -r ${output_sim}.bam -t $output_sim -w MoreBarcodes.fasta -as CTACACGACGCTCTTCCGATCT -ao 10 -ae 0.3 -ag -2 -hr T -hi 10 -he 0.3 -bo 5 -be 0.2 -bg -2 -ul 26 -kb 3 -fl 100
awk '$2!="NA" || NR==1' ${output_sim}.tab > ${output_sim}.tab1
```
## add_label.r
```
options(stringsAsFactors = FALSE)
dat=Sys.getenv('output_sim')
x=read.table(paste0(dat,'.tab1'),sep="\t",header=TRUE)
y=read.table(paste0(dat,'_barcodes.txt'),sep="\t",header=FALSE,row.names=1)
x[,ncol(x)] = 1-as.integer(y[x[,1],1]==x[,2])
colnames(x)[ncol(x)]='label'
write.table(x,file=paste0(dat,'.tab1'),sep="\t",quote=F,row.names=FALSE)
```
## build_model.r
```
library(caret)
library(e1071)
Sys.setlocale("LC_NUMERIC","C")
args = commandArgs(trailingOnly=TRUE)
mdat =args[1]
dat = args[2]
options(stringsAsFactors = FALSE)
j=c(5,7,8,9,10,11) # We use feature 3-9
x=read.table(dat,sep="\t",header=TRUE)
d=preProcess(x[,j], method=c("YeoJohnson","range"))
x[,j]=predict(d,x[,j])
model = naiveBayes(label~., data = x[,c(j,ncol(x))])
save(model,d,file=mdat)
```
## pred.r
```
library(caret)
library(e1071)
Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
j=c(5,7,8,9,10,11) # We use feature 3-9
args = commandArgs(trailingOnly=TRUE)
mdat =args[1]
dat = args[2]
load(mdat)
y=read.table(dat,sep="\t",header=TRUE)
y2=y
y[,j]=predict(d,y[,j])
pred = predict(model, y[,j], "raw")[,1]
# Prior knowledge from Illumina sequencing
cnt = read.table('reads_per_barcode')
n = cnt[match(y[,2],cnt[,2]),1]
n[is.na(n)]=0
# Pseudocount for the barcodes which are not in the whitelist
p=1
y1 = split(cbind(y2,n,pred),y2[,1])
# Do the bayesian calculation
m=do.call(rbind,lapply(y1,function(d){
d[,ncol(d)]=unlist(lapply(seq_len(nrow(d)),function(i){
prob=d[i,'n']*d[i,'pred']/(sum(d[,'n']*d[,'pred'])+p)
as.integer(round(min(d[i,'pred'],prob)*100))
}))
d
}))
write.table(m,file=paste0(dat,'.prob'),sep="\t",quote=F,row.names=F,col.names=F)
```
## filter_pred.sh
```
awk '$15>30' $1 | sed 's/_end[1|2]//' | awk '{a[$1]++;b[$1]=$0}END{for(i in a){if(a[i]==1)print b[i]}}' | cut -f1-2 > $1.label
```
## feat_stat.r
```
mat = cor(x[,j])
findCorrelation(mat, cutoff=0.5)
corrplot(mat, method="circle")
control = trainControl(method="repeatedcv", number=10, repeats=3)
model = train(label~., data=x[,c(j,14)], method="lvq", preProcess=c("YeoJohnson","range"), trControl=control)
importance = varImp(model, scale=FALSE)
```
## feat_var.r
```
x=read.table('FC1.prob',header=TRUE)
df=data.frame(assignment=c('unassigned','assigned')[as.integer(x[,ncol(x)]>30)+1],melt(x[,5:11],id.vars=NULL))
p=ggplot(df,aes(value,fill=assignment))+geom_histogram()+facet_wrap(~variable,scale="free")
ggsave("fc1-fp.pdf",width=6,height=4.5,units="in")
```
