options(stringsAsFactors = FALSE)
Sys.setlocale("LC_NUMERIC","C")
args = commandArgs(trailingOnly=TRUE)
x=read.table(args[1],sep="\t",header=TRUE)
y=read.table(args[2],sep="\t",header=FALSE,row.names=1)
x[,ncol(x)] = 1-as.integer(y[x[,1],2]==x[,2])
colnames(x)[ncol(x)]='label'
write.table(x,file=args[3],sep="\t",quote=F,row.names=FALSE)
