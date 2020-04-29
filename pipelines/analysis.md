## gene_names.sh
```
samtools view possorted_genome_bam.bam|perl -ne 'print "$2\t$1\n" if /GN:Z:(\S+).*CB:Z:([ACGT]+)/' > bc.gene
awk -v OFS='\t' '{print $2,$1}' bc.gene > gene.bc
```
## count_umi.sh
```
perl -F"\t" -ane '$h{$F[0]}{$F[1]}++; END { print "$_\t".(keys %{$h{$_}})."\n" for keys %h }' fc1.barumi > fc1.barumi.cnt
```
## count_isoform.sh
```
perl -ne 'print unless /cov \"[0|1]\./' FC1new.gtf|perl -ne 'print "$1\n" if /gene_name \"(\S+)\"/'|sort|uniq -c > fc1.uipg
```
## ts_cov.sh
```
grep transcript Illu.gtf |perl -F"\t" -ane 'print $F[4]-$F[3],"\t$1\n" if /cov "(\S+)\"/' > illu.cov
```
## cor_gene.r
```
Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
y=read.table('FC1.label',col.names=c('id','barcode'))
nm=read.table('FC1.ge')
i=match(y$id,nm[,1])
df=data.frame(y$barcode,nm[i,2])
colnames(df)=c('barcode','gene')
df=df[!is.na(df[,2]),]
x=reshape2::dcast(df, gene~barcode, length)
rownames(x)=x[,1]
x=x[,-1]
lw=function(d) length(which(d))
i=apply(x,1,function(d) lw(d>0))
j=apply(x,2,function(d) lw(d>0))
x=x[i>=3,j>=200]
cb=read.table(gzfile("barcodes.tsv.gz"),sep='-')[,1]
x=x[,colnames(x) %in% cb]
y=read.table('gene.bc',col.names=c('gene','barcode'))
y=y[y$barcode %in% colnames(x),]
y=y[y$gene %in% rownames(x),]
y=reshape2::dcast(y, gene~barcode, length)
rownames(y)=y[,1]
y=y[,-1]
save(x,y,file="ge.rda")
ba=intersect(colnames(x),colnames(y))
ge=intersect(rownames(x),rownames(y))
sb=sample(ba)
same=unlist(lapply(seq_len(length(ba)),function(i)cor(x[ge,ba[i]],y[ge,ba[i]])))
diff=unlist(lapply(seq_len(length(ba)),function(i)cor(x[ge,ba[i]],y[ge,sb[i]])))
library(ggplot2)
library(ggbeeswarm)
df=data.frame(V1=same,V2='same cell')
df=rbind(df,data.frame(V1=diff,V2='different cell'))
p=ggplot(df, aes(x=V2, y=V1)) + geom_quasirandom(dodge.width=.8,cex=.2) + labs(x='',y='Correlation Illumina / Nanopore')
ggsave(p,file='fc1-same.pdf',height=6,width=6)
```
## count_umi.r
```
Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
x=read.table('fc1.barumi.cnt')
y=read.table('illu.barumi.cnt')
cb=read.table(gzfile("barcodes.tsv.gz"),sep='-')[,1]
x=x[x[,1] %in% cb,]
y=y[y[,1] %in% cb,]
xy=merge(x,y,by='V1')
png('FC1-UMIcount.png')
plot(xy[,2:3],main="Number of UMI sequences",xlab="Nanopore",ylab="Illumina")
text(100,1e5,paste0("r=",round(cor(xy[,2],xy[,3]), digits=2)))
dev.off()
```
## count_upig.r
```
Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
y=read.table('FC2.label',col.names=c('id','barcode'))
nm=read.table('FC2.ge')
i=match(y$id,nm[,1])
df=data.frame(y$barcode,nm[i,2])
colnames(df)=c('barcode','gene')
df=df[!is.na(df[,2]),]
x=reshape2::dcast(df, gene~barcode, length)
rownames(x)=x[,1]
x=x[,-1]
lw=function(d) length(which(d))
i=apply(x,1,function(d) lw(d>0))
j=apply(x,2,function(d) lw(d>0))
x=x[i>=3,j>=200]
cb=read.table(gzfile("barcodes.tsv.gz"),sep='-')[,1]
x=x[,colnames(x) %in% cb]
y=read.table('fc2.uipg',row.names=2)
df=data.frame(V1=i,V2=y[names(i),])
df=df[!is.na(df[,2]),]
png('FC2-UIPG.png')
plot(df,xlab="Cells expressing gene",ylab="Unique isoforms per gene",ylim=c(0,200),col='darkgrey')
points(df[grepl('^RPS',rownames(df)),],col='red')
points(df[grepl('^RPL',rownames(df)),],col='red')
dev.off()
```
