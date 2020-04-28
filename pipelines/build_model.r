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
args = commandArgs(trailingOnly=TRUE)
j=c(5,7,8,9,10,11) # We use feature 3-9
x=read.table(args[1],sep="\t",header=TRUE)
d=preProcess(x[,j], method=c("center","scale","BoxCox"))
x[,j]=predict(d,x[,j])
model = naiveBayes(label~., data = x[,c(j,ncol(x))])
save(model,d,file=args[2])
