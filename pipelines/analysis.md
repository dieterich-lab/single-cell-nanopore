## gene_names.sh
```
samtools view possorted_genome_bam.bam|perl -ne 'print "$2\t$1\t$3\n" if /GN:Z:(\S+).*CB:Z:([ACGT]+).*UB:Z:([ACGT]+)/' > illu.gene
cut -f1,2 illu.gene|sort|uniq|cut -f1|sort|uniq -c|perl -npe 's/^ +//;s/ +/\t/' > illu.gene.cnt
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
## test_barcode.sh
```
perl  gt1.pl fc1sv.prob > fc1sv.prob.lbl
perl sis1.pl fc1sv.prob.lbl > fc1sv.prob.lbl.sis
sort -k15,15rn fc1sv.prob.lbl.sis | perl -ne 'chomp;@t=split(/\t/);print "$_\t";print defined $h{$t[0]}?1:0;print "\n";$h{$t[0]}=1' > fc1sv.prob.lbl.sis2
awk '$17==0 && $18==1' fc1sv.prob.lbl.sis2|cut -f19|sort|uniq -c
```
## ts_cov.sh
```
grep transcript Illu.gtf |perl -F"\t" -ane 'print $F[4]-$F[3],"\t$1\n" if /cov "(\S+)\"/' > illu.cov
```
## gb_cov.sh
```
~/bin/ExpressionAnalysis-ea-utils-bd148d4/clipper/gtf2bed GRCh38.90.gtf > 1.bed
perl -F"\t" -ane '$s=$F[2]-$F[1];print "$s\t$_" if $s>100' 1.bed|sort -k1,1n|cut -f2-|perl -ne 'print if $i++%40==0' > 2.bed
geneBody_coverage.py -i FC1.bam,FC2.bam -r 2.bed -o coverage
```
## stat.sh
```
x="fc1.tab"
cut -f1 $x|perl -npe 's/_end[1|2]//'|sort|uniq|wc -l
awk '$3=="no"' $x|cut -f1|perl -ne 's/_end[1|2]//;print'|sort|uniq|wc -l
awk '$2!="NA"' $x|cut -f1|perl -ne 's/_end[1|2]//;print'|sort|uniq|wc -l
```
## sim_ed.sh
```
samtools view fc1.GEUS10xAttributes.bam|perl -F"\t" -ane '$F[0]=~s/_REV.*//;print "$F[0]\t$2\t$1\n" if /B1:i:(\d+).*BC:Z:(\w+)/' > fc1.sis
perl gt1feat.pl fc1.prob > fc1.prob.ed
perl gt1feat.sis.pl fc1.sis > fc1.sis.ed
cat fc1.sis.ed fc1.prob.ed|perl -F"\t" -ane '$h{"$F[0]\t$F[1]"}="$F[0]\t$F[1]\t$F[2]\t3\n"; END {print $h{$_} foreach keys %h}' > fc1.all.ed
cat fc1.sis.ed fc1.prob.ed fc1.all.ed > fc1.tab.ed
```
## gfp_isoform.sh
```
samtools view -h FC1new.bam 17:76734115-76737374|awk '/^@/||$2==0||$2==16' > SRSF2.FC1.sam
grep -v '^@' SRSF2.FC1.sam|awk '{print "@"$1"\n"$10"\n+\n"$11}' > SRSF2.FC1.fq
cat ../FC1.label ../FC2.label > label.pl
perl label.pl SRSF2.FC1.sam
for f in $(ls -1 *.sam);do samtools view -bS $f > ${f%.*}.bam;samtools index ${f%.*}.bam;done
for f in *.bam;do regtools junctions extract $f -s 0|intersectBed -a stdin -b as.bed -F 0.5;done
```
## fisher_test.r
```
Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
x=read.table('fc1sv.prob.lbl.sis2')
(m=table(x[x$V17==0,c(19,18)]))
fisher.test(matrix(matrix(m)[3:6,1],ncol=2))
x$V20=abs(x$V8)>=3
(m=table(x[x$V17==0,c(20,18)]))
#ours
(m=table(x[x$V16>30,c(20,17)]))
fisher.test(matrix(matrix(m)[,1],ncol=2))
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
## diff_isoform.r
```
options(stringsAsFactors = FALSE)
Sys.setlocale("LC_NUMERIC","C")
library(ballgown)
bg=ballgown(dataDir="/scratch/qw/nanopore/", samplePattern = "GFP", meas='all')
pData(bg)=data.frame(id=sampleNames(bg), group=substr(sampleNames(bg),1,6), batch=substr(sampleNames(bg),8,10))
df = stattest(bg, feature="transcript", covariate="batch", getFC=TRUE, meas="FPKM")
df$significant=df$pval < 0.01
library(ggplot2)
p=ggplot(df) + geom_point(aes(x=log2(fc), y=-log10(pval), colour=significant))+scale_x_continuous(limits = c(-3,3))
ggsave(p,file='fc1-dffiso.png',height=6,width=6)
```
## gb_cov.r
```
df=melt(data_matrix)
colnames(df)[1]='Library'
p=ggplot(data=df, aes(x=Var2, y=value, group=Library)) +
  geom_line(aes(color=Library))+labs(title='Gene body coverage',x="Gene body percentile (5'->3')", y="Coverage")
ggsave(p,file='gbcov.pdf',height=6,width=6)
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
colnames(xy)=c('barcode','Nano','Illu')
library(ggplot2)
p=ggplot(xy[,2:3], aes(x=Nano, y=Illu)) + geom_point() + labs(title="Number of UMI sequences",x="Nanopore", y="Illumina") + geom_text(x=100,y=1e5,label=paste0("r=",round(cor(xy[,2],xy[,3]), digits=2)))
ggsave(p,file='fc1-corumi.pdf',height=6,width=6)
```
## count_gene.r
```
Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
y=read.table('FC1.label',col.names=c('id','barcode'))
nm=read.table('FC1.ge')
i=match(y$id,nm[,1])
df=data.frame(y$barcode,nm[i,2])
colnames(df)=c('barcode','gene')
df=df[!is.na(df[,2]),]
write.table(df,sep='\t',quote=F,row.names=F,col.names=F,file='FC1.gene')
#sort FC1.gene |uniq|cut -f1|sort|uniq -c|perl -npe 's/^ +//;s/ +/\t/' > FC1.gene.cnt
Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
x=read.table('FC1.gene.cnt')
y=read.table('illu.gene.cnt')
cb=read.table(gzfile("barcodes.tsv.gz"),sep='-')[,1]
x=x[x[,2] %in% cb,]
y=y[y[,2] %in% cb,]
xy=merge(x,y,by='V2')
colnames(xy)=c('gene','Nano','Illu')
library(ggplot2)
p=ggplot(xy[,2:3], aes(x=Nano, y=Illu)) + geom_point() + labs(title="Number of Genes",x="Nanopore", y="Illumina") + geom_text(x=100,y=8e3,label=paste0("r=",round(cor(xy[,2],xy[,3]), digits=2)))
ggsave(p,file='fc1-corgene.pdf',height=6,width=6)
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
df$RP=grepl('^RP[S|L]',rownames(df))
p=ggplot(df, aes(x=V1, y=V2, color=RP)) + geom_point() + labs(title="Number of unique isoforms",x="Cells expressing gene", y="Unique isoforms per gene")+scale_y_continuous(limits = c(0,200))
ggsave(p,file='fc1-numiso.pdf',height=6,width=6)
```
## exp_ratio.r
```
Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
x=read.table('illu.cov')
y=read.table('fc1.cov')
breaks=seq(0,5700,300)
tags=paste(breaks,breaks+300,sep='-')
breaks=c(breaks,6000)
x$V1 = cut(x$V1,breaks=breaks,include.lowest=TRUE,right=FALSE,labels=tags)
y$V1 = cut(y$V1,breaks=breaks,include.lowest=TRUE,right=FALSE,labels=tags)
x=do.call(rbind,lapply(split(x,x$V1),function(d)d[sample(seq_len(nrow(d)),1000),]))
y=do.call(rbind,lapply(split(y,y$V1),function(d)d[sample(seq_len(nrow(d)),1000),]))
df=x
df$V2 = y$V2/x$V2
df=df[is.finite(df$V2),]
df=df[df$V2<20,]
library(ggplot2)
p=ggplot(df, aes(x = V1, y = V2)) + geom_boxplot(outlier.shape = NA) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + labs(x='Transcript length (bp)',y="Expression level ratio \n Nanopore / Illumina") + coord_cartesian(ylim = c(0,10)) 
ggsave(p,file='fc1-ratio.pdf')
```
## runtime.r
```
#https://gist.github.com/qwang-big/658d3ca9ad775503bb20123151e37b55
perm_without_replacement = function(n, r) factorial(n)/factorial(n-r)
d  <- do.call(rbind,lapply(0:4,function(i){
l=2*i+1
j=2000
data.frame(pos=i,Sicelore=l*perm_without_replacement(i+16,3)*j,SingleCellPipe=(l+16)*16*j)
}))
x  <- 'pos'
y1 <- 'SingleCellPipe'
y2 <- 'Sicelore'
a            <- range(d[[y1]])
b            <- range(d[[y2]])
scale_factor <- diff(a)/diff(b)
d[[y2]]      <- ((d[[y2]] - b[1]) * scale_factor) + a[1]
trans <- ~ ((. - a[1]) / scale_factor) + b[1]
ggplot(d, aes(x = pos)) +
  geom_line(aes(y = Sicelore, colour = "Sicelore")) + 
  geom_line(aes(y = SingleCellPipe, colour = "SingleCellPipe"))  + scale_colour_manual(values = c("blue", "red"))+ 
  scale_y_continuous(sec.axis = sec_axis(trans=trans, name=y2)) +theme(legend.position="top") + labs(title="Runtime comparison",x="Barcode positions", y=y1)
```
## tsne.r
```
library(Rtsne)
lbl=read.table("gfp.txt")[,1]
lbl=names(which(table(lbl)>20))
lw=function(d) length(which(d))
i=apply(x,1,function(d) lw(d>0))
j=apply(x,2,function(d) lw(d>0))
x=t(x[i>=3,j>=200])
cb=read.table(gzfile("barcodes.tsv.gz"),sep='-')[,1]
x=x[rownames(x) %in% cb,]
tsne <- Rtsne(x, dims=2, perplexity=30, initial_dims=10, num_threads=0)
df=data.frame(tsne$Y)
df$GFP = rownames(x) %in% lbl
p=ggplot(df, aes(x=X1, y=X2, color=GFP)) + geom_point() + labs(title="Gene expression Nanopore",x="tSNE_1", y="tSNE_2")
ggsave(p,file='fc1-tsne.pdf',height=6,width=6)
df=read.csv('projection.csv')
df$Barcode=substr(df$Barcode,1,16)
df$GFP = df$Barcode %in% lbl
p=ggplot(df, aes(x=TSNE.1, y=TSNE.2, color=GFP)) + geom_point() + labs(title="Gene expression Illumina",x="tSNE_1", y="tSNE_2")
ggsave(p,file='illu-tsne.pdf',height=6,width=6)
```
## tsne_isoform.r
```
df=data.frame(do.call(cbind,lapply(ff,function(f){x=read.table(f)
as.numeric(x[seq(2,nrow(x),2),1])})))
rownames(df)=x[seq(1,nrow(x),2),1]
write.table(df,sep='\t',quote=F,file='SRSF2.isoform.txt')

cnt=read.table('SRSF2.isoform.txt')
cname=gsub('\\.','-',colnames(cnt))
exp=cnt[match(rownames(x),rownames(cnt)),]
exp[is.na(exp)]=0
exp2=exp[,c(1,3)]+exp[,c(2,4)]
expression=exp2[,1]
p=ggplot(data.frame(tsne$Y), aes(x=X1, y=X2, color=expression)) + geom_point() +theme(legend.position="top")+ labs(title='SRSF2-WT',x="tSNE_1", y="tSNE_2")+scale_color_gradient2(low="black", mid="yellow", high="yellow", midpoint = 4.5)
ggsave(p,file='SRSF2-WT.pdf',height=6,width=5)
```
## seurat.r
```
#https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html
pbmc = CreateSeuratObject(counts = x, project = "FC1", min.cells = 3, min.features = 200)
p = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
p=VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
p=DimPlot(pbmc, reduction = "pca")
p=DimHeatmap(pbmc, dims = 1, cells = 200, balanced = TRUE)
pbmc <- RunUMAP(pbmc, dims = 1:10)
p=DimPlot(pbmc, reduction = "umap")
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10)
p=DimPlot(pbmc, reduction = "tsne")
ggsave(p, file='FC1-tsne.pdf',height=6,width=6)
```
## stat.r
```
x=read.csv(text="Steps,FC1,FC2
Aligned to adapter,815243,807878
Aligned to barcode,751903,741537
Assigned to barcode,484458,464685
Correctly assigned,458145,440433")
y=read.csv(text="Steps,FC1,FC2
Total reads,13126013,11923896
Aligned to genome,11158994,10164820
Aligned to adapter,8223382,7571238
Aligned to barcode,7457842,6882174
Assigned to barcode,4692389,4359238")
x=read.csv(text="event,WT,KD
NMD_ex,556,950
NMD_in,1113,1890
A5SS,441,788")
x=melt(x)
colnames(x)[2]="RunId"
x$Steps=factor(x$Steps,levels=rev(levels(as.factor(x$Steps))))
p=ggplot(x, aes(x=Steps, y=value, fill=RunId)) + geom_bar(stat="identity", position=position_dodge())+theme(legend.position="top")+coord_flip()+ labs(title='Number of simulated reads per step',x="Processing steps", y="Number of simulated reads")
```
## feat_stat.r
```
mat = cor(x[,j])
#findCorrelation(mat, cutoff=0.5)
#corrplot(mat, method="circle")
p=ggcorrplot(mat, method = "circle")
ggsave(p,file='fc1-cor.pdf',height=5,width=5)

control = trainControl(method="repeatedcv", number=10, repeats=3)
model = train(label~., data=x[,c(j,14)], method="lvq", preProcess=c("YeoJohnson","range"), trControl=control)
importance = varImp(model, scale=FALSE)
df = data.frame(importance$importance)
df[,1] = rownames(df)
p=ggplot(df, aes(x=reorder(X0, -X1), y=X1, fill='blue')) + geom_bar(stat="identity", position=position_dodge())+theme(legend.position="none")+coord_flip()+ labs(title='Feature importances',y="Importance", x="")
ggsave(p,file='fc1-varimp.pdf',height=4,width=4)
```
## feat_var.r
```
x=read.table('FC1.prob',header=TRUE)
df=data.frame(assignment=c('unassigned','assigned')[as.integer(x[,ncol(x)]>30)+1],melt(x[,5:11],id.vars=NULL))
p=ggplot(df,aes(value,fill=assignment))+geom_histogram()+facet_wrap(~variable,scale="free")
ggsave("fc1-fp.pdf",width=6,height=4.5,units="in")
```
## den_pred.r
```
Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
x=read.table('fc1.prob.den',header=T,col.names=c('Score','Assignment'))
x[,2]=c('correct','false')[x[,2]+1]
p=ggplot(x, aes(x = Score, group = Assignment)) + geom_density(aes(color = Assignment))+ labs(title="Density of barcode scores",x="Scores", y='Density')
ggsave(p,file='fc1-pden.pdf',height=6,width=6)
```
## edit_dist.r
```
Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
x=read.table('fc1.tab.ed')
table1=function(x){x=data.frame(table(x));x[,1]=as.integer(as.character(x[,1]));x[x[,1]<10,]}
df=data.frame(table1(x[x$V4==3,3]),method="ground_truth")
df=rbind(df,data.frame(table1(x[x$V4==2,3]),method="Sicelore"))
df=rbind(df,data.frame(table1(x[x$V4==0,3]),method="Single^2"))
df[nrow(df)+1,]=c(4,0,'Sicelore')
df[nrow(df)+1,]=c(5,0,'Sicelore')
df[nrow(df)+1,]=c(5,0,'Single^2')
df$Freq=as.integer(df$Freq)
p=ggplot(data=df, aes(x=x, y=Freq, fill=method)) + geom_bar(stat="identity", position=position_dodge())+theme(legend.position="top")+ labs(title='Edit distance of simulated reads',x="Edit distance", y="Percentage of simulated reads")+coord_cartesian(xlim = c(0,6)) 
ggsave(p,file='fc1-ed.pdf',height=6,width=6)
```
## sis1.pl
```
while(<DATA>){chomp;
@t=split(/\t/);
$h{$t[0]}=$t[1];
}
while(<>){chomp;
@t=split(/\t/);
$t[0]=~s/_end\d//;
print $_,"\t";
print defined $h{$t[0]} ? ($h{$t[0]} eq $t[1] ? 0 : 1) : -1;
print "\n"
}
__DATA__
```
## label.pl
```
while(<DATA>){chomp;
@t=split(/\t/);
$h{$t[0]}=$t[1]
}
while(<>){chomp;
next if /^@/;
@t=split(/\t/);
$h1{$h{$t[0]}}.=$_
}
foreach(keys %h1){
mkdir($_);
open(F,">$_/seq.sam");
print F "\@SQ\tSN:17\tLN:83257441\n";
print F $h1{$_};
close F;
}
__DATA__
```
