if(!require(caret)){
    install.packages("caret",repos="https://cloud.r-project.org")
    library(caret)
}
if(!require(e1071)){
    install.packages("e1071",repos="https://cloud.r-project.org")
    library(e1071)
}
Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
j=c(5,7,8,9,10,11) # We use feature 3-9
args = commandArgs(trailingOnly=TRUE)
load(args[1])
y=read.table(args[2],sep="\t",header=TRUE)
y2=y
y[,j]=predict(d,y[,j])
pred = predict(model, y[,j], "raw")[,1]
# Prior knowledge from Illumina sequencing
cnt = read.table(args[3])
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
write.table(m,file=args[4],sep="\t",quote=F,row.names=FALSE,col.names=TRUE)
