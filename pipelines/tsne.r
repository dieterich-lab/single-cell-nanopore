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
png("FC1-gfp.png")
plot(tsne$Y,col=c('black','red')[as.integer(rownames(x) %in% lbl)+1])
dev.off()
