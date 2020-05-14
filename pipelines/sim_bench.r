Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
library(reshape2)
library(ggplot2)
t1=read.table('gtruth.txt')
t2=read.table('sicelore.txt')
i=match(t1[,1], t2[,1])
pred=as.integer(is.na(t2[i,2]))
fal=t2[i,2]!=t1[,2]
fal[is.na(fal)]=FALSE
unal=grepl('_unal',t1[,1])
gtruth=as.integer(unal | fal)
res=caret::confusionMatrix(factor(pred), factor(gtruth))
ss=c(res$overall,res$byClass)[c(1,8,9,12,13,14)]
v = seq(0,100)
m= read.table('fc1.prob')
df=do.call(rbind,lapply(v,function(d) {
t2=m[m[,16]>d,1:2]
t2[,1]=substr(t2[,1],1,nchar(t2[,1])-5)
i=match(t1[,1], t2[,1])
pred=as.integer(is.na(t2[i,2]))
fal=t2[i,2]!=t1[,2]
fal[is.na(fal)]=FALSE
unal=grepl('_unal',t1[,1])
gtruth=as.integer(unal | fal)
res=caret::confusionMatrix(factor(pred), factor(gtruth))
c(res$overall,res$byClass)[c(1,8,9,12,13,14)]
}))
rownames(df)=v
df=df[!is.na(df[,4]),]
# ROC curve
dat=data.frame(Specificity=c(1-ss[3],1-df[,3]),Sensitivity=c(ss[2],df[,2]),method=c('sicelore',rep('bayesian',nrow(df))))
p=ggplot(dat,aes(Specificity,Sensitivity))+geom_point(aes(color=method))+labs(x='1-Specificity',y='Sensitivity')
ggsave('fc1-roc.pdf')
# all benchmarks
df=melt(df)
df1=melt(data.frame(ss,var=names(ss)))
p=ggplot(df,aes(Var1,value,col=Measures))+geom_line()+geom_vline(xintercept=35, linetype=4)+ annotate("text", x=38,y=0.2,label="35")+geom_hline(data=df1,aes(yintercept=value,colour=var),linetype="dashed")+
labs(title='Performance of the predicted probabilities',x='Predicted probability cutoff',y='Score')
ggsave('fc1-t.pdf')
