# Illumina library
## get_cbc.sh
```
zcat outs/filtered_feature_bc_matrix/barcodes.tsv.gz | perl -ne 'print ">$1\n$1\n" if /^(\w+)/' > CellBarcodes.fasta
```
## get_gfpcbc.sh
```
samtools view possorted_genome_bam.bam pcDNA5:1-6000 | perl -ne 'print "$1\n" if /CB:Z:([ACGT]+)/' > gfp.txt
```
## find_dist.r
```
library(stringdist)
options(stringsAsFactors = FALSE)
x=read.table('../reads_per_barcode')
y=read.table('CellBarcodes.fasta')
r=unlist(lapply(y[seq(2,nrow(y),2),1],function(d) x[stringdist(d,x[,2],method='dl')<2,2]))
write.table(r,file='10k.barcodes.fa',sep="\t",quote=F,row.names=F,col.names=F)
```
## get_cbfreq.sh
```
samtools view possorted_genome_bam.bam | grep GN:Z:.*CB:Z: | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > reads_per_barcode
```
# Nanopore library
## align_longreads.sh
```
~/minimap2-2.17/minimap2 -o FC1 -t 5 -ax splice --MD -ub GRCh38_90.mmi pass_reads_guppy333.fastq.gz
```
## assign_genenames.sh
```
java -jar -Xmx64g Jar/Sicelore-1.0.jar AddGeneNameTag I=FC1.bam O=FC1.GE.bam REFFLAT=refFlat.txt GENETAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT
samtools view FC1.GE.bam|perl -ne 'print "$1\t$2\n" if /^(\S+).*GE:Z:(\S+)/' > FC1.ge
```
# Build artifical genome
## build_nanosim.sh
```
~/bin/NanoSim-2.4-beta/src/read_analysis.py transcriptome -i FC1.fq -rg GRCh38_90.fa -rt Homo_sapiens.GRCh38.cdna.all.fa -annot Homo_sapiens.GRCh38.90.gff3 -t 10 -o ../NanosimModel/FC1
```
## build_genome.sh
```
samtools view possorted_genome_bam.bam |perl -ne '@t=split(/\t/);print ">",++$j,"\n" if $i++%25e5==0;print "CTACACGACGCTCTTCCGATCT$3$4TTTTTTTTTTTTTTTTTTTT",substr($t[9],0,32),"\n" if /(TX|AN):Z:(\w+).*CB:Z:([ACGT]+).*UB:Z:([ACGT]+)/' > genome.fa
samtools faidx genome.fa
```
## sim_reads.sh
```
~/bin/NanoSim-2.4-beta/src/simulator.py genome -rg genome.fa -c ../NanosimModel/FC1 -o fc1 -n 1000000
```
## build_align.sh
```
perl -ne '$L=100;chomp;if (/^>/){$id=$_}else{@t=split(/_/,$id);$s=$_;$d=$t[5];if ($t[4] eq 'R'){$s=reverse $s;$s=~tr/ATGCatgc/TACGtacg/;$d=$t[7]}$s=substr($s,$d,$t[6]);print $id,"\n",$t[6]>$L?substr($s,$L-$t[1]%$L,68):$s,"\n"}' fc1_reads.fasta > fc1.test.fa
perl fa2sam.pl fc1.test.fa|samtools view -bS > fc1.bam
```
## get_barcodes.sh
```
perl -ne '$L=100;next unless /^>/;$_=substr($_,1);@t=split(/_/);$d=$t[1]+$L-$t[1]%$L;print "$t[0]\t$d\t",$d+$L,"\t$_"' fc1_reads.fasta|sort -k1,1 -k2,2n > fc1.real.bed
bedtools getfasta -fi genome.fa -bed fc1.real.bed -name > fc1.real.fa
perl -ne 'print "$1\t" if /^>(\w+):/;print substr($_,22,16),"\n" unless /^>/' fc1.real.fa > fc1.barcodes.txt
```
# Features extraction and build model
## run_pipe.sh
```
singleCellPipe -n 20 -r fc1.bam -t FC1 -w 10k.barcodes.fa -as CTACACGACGCTCTTCCGATCT -ao 10 -ae 0.3 -ag -2 -hr T -hi 10 -he 0.3 -bo 5 -be 0.2 -bg -2 -ul 26 -kb 3 -fl 100
awk '$2!="NA" || NR==1' fc1.tab > fc1.tab1
```
## add_label.r
```
options(stringsAsFactors = FALSE)
x=read.table('fc1.tab1',sep="\t",header=TRUE)
y=read.table('fc1.barcodes.txt',sep="\t",header=FALSE,row.names=1)
x[,ncol(x)] = 1-as.integer(y[x[,1],1]==x[,2])
colnames(x)[ncol(x)]='label'
write.table(x,file='fc1.tab1',sep="\t",quote=F,row.names=FALSE)
```
## build_model.r
```
library(caret)
library(e1071)
Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
args = commandArgs(trailingOnly=TRUE)
j=c(5,7,8,9,10,11) # We use feature 3-9
x=read.table('fc1.tab1',sep="\t",header=TRUE)
d=preProcess(x[,j], method=c("YeoJohnson","range"))
x[,j]=predict(d,x[,j])
model = naiveBayes(label~., data = x[,c(j,ncol(x))])
save(model,d,file='fc1.model.rda')
```
## pred.r
```
y=read.table('FC1.tab1',sep="\t",header=TRUE)
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
write.table(m,file='FC1.prob',sep="\t",quote=F,row.names=F,col.names=F)
```
## filter_pred.sh
```
awk '$15>30' FC1.prob | sed 's/_end[1|2]//' | sort -k15,15rn |sort -u -k1,1 | cut -f1-2 > FC1.label
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
