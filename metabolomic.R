####heterogeneity between replicates####

#RNA-Seq
log.fc.base=data.frame(base=apply(cpm(y)[,c(11,12,23,33,34)],1,mean))
cv.base.gene=data.frame(base=apply(cpm(y)[,c(11,12,23,33,34)],1,function(x) var(x)/(mean(x))^2))
cv.base.gene.sqrt=sqrt(cv.base.gene)
dat.pre=data.frame(log2(cpm(y)[,-c(11,12,23,33,34)]/log.fc.base$base))
group_rep=paste(substr(colnames(dat.pre),1,1),substr(colnames(dat.pre),8,8),sep = "_")

logfc.emp=data.frame(rep1=apply(dat.pre[,c(1,11,21)],1,mean),
                     rep2=apply(dat.pre[,c(2,12,22)],1,mean),
                     rep3=apply(dat.pre[,c(3,13,23)],1,mean),
                     rep4=apply(dat.pre[,c(4,14,24)],1,mean),
                     rep5=apply(dat.pre[,c(5,15,25)],1,mean),
                     rep6=apply(dat.pre[,c(6,16,26)],1,mean),
                     rep7=apply(dat.pre[,c(7,17,27)],1,mean),
                     rep8=apply(dat.pre[,c(8,18,30)],1,mean),
                     rep9=apply(dat.pre[,c(9,19,28)],1,mean),
                     rep10=apply(dat.pre[,c(10,20,29)],1,mean))

cor.logfc.emp=cor(logfc.emp)[upper.tri(cor(logfc.emp))]

cor_display=function(x,y,...) {
  text(0,0,labels=round(cor(x,y,method = "pearson"),2),cex=1,...)
}
scatter.plot=function(x,y,...){
  points(x,y,asp=1,pch=19)
  abline(v=0,col="grey")
  abline(h=0,col="grey")
}
png("../../gene_logfc/pairwise_scatter_plot_logFC_gene_pearson.png",height = 12.5,width = 12.5,units = "cm",res=600,pointsize = 9)
pairs(logfc.emp,upper.panel = scatter.plot,
      lower.panel = cor_display,xlim=c(-4,4),ylim=c(-4,4),asp=1)
dev.off()

dat.use=as.matrix(dat.pre)
geno=group_rep
res=substr(colnames(dat.use),9,9)
meta.table=data.frame(geno=geno,tech=res,row.names = colnames(dat.use))
annot=data.frame(labelDescription=c("Factor levels","Factor levels"))
annot_factors=AnnotatedDataFrame(meta.table,annot)
expr.set=ExpressionSet(dat.use,annot_factors)
pvca_res=pvcAnaly(expr.set, 0.75,c("geno","tech"))
barplot(c(0.5246,0.4729),col="blue")

bp=barplot(c(0.52,0.47),  xlab = "Effects",
           ylab = "Weighted average proportion variance", ylim= c(0,1.1),
           col = c("darkblue"), las=2, main="PVCA estimation bar chart")
axis(1, at = bp, labels = c("replicate","residual"), xlab = "Effects", cex.axis = 1, las=1)
fig_label("B",cex = 2)
text(bp,c(0.52,0.47),c(0.52,0.47),pos = 3,cex = 0.75)


cor.logfc.permut=list()
for (i in 1:100) {
  ind=sample(1:30,30,replace = F)
  logfc.permut=data.frame(rep1=apply(dat.pre[,ind[1:3]],1,mean),
                       rep2=apply(dat.pre[,ind[4:6]],1,mean),
                       rep3=apply(dat.pre[,ind[7:9]],1,mean),
                       rep4=apply(dat.pre[,ind[10:12]],1,mean),
                       rep5=apply(dat.pre[,ind[13:15]],1,mean),
                       rep6=apply(dat.pre[,ind[16:18]],1,mean),
                       rep7=apply(dat.pre[,ind[19:21]],1,mean),
                       rep8=apply(dat.pre[,ind[22:24]],1,mean),
                       rep9=apply(dat.pre[,ind[25:27]],1,mean),
                       rep10=apply(dat.pre[,ind[28:30]],1,mean))
  cor.logfc.permut[[i]]=cor(logfc.permut)[upper.tri(cor(logfc.permut))]
  print(paste0(i," time"))
}

plot(density(cor.logfc.permut[[1]]),xlim=c(0.4,1),
     main="Distribution of permuted correlation across replicates",lty=2,col="grey",xlab="correlation coefficient",ylim=c(0,15))
for (i in 2:100) {
  lines(density(cor.logfc.permut[[i]]),add=T,lty=2,col="grey")
}
lines(density(cor.logfc.emp),add=T,lty=2,col="red",lwd=1.5)
legend("topright",bty = "n",legend = c("oberved","permuted"),col = c("red","grey"),lty=2,lwd = 2)



####replicate-specific effect####
test.rep=matrix(NA,nrow=10780,ncol = 3)
row.names(test.rep)=row.names(dat.pre)
for (i in 1:10780) {
  temp=data.frame(obs=as.numeric(dat.pre[i,]),trt=group_rep)
  fit=anova(lm(temp$obs~temp$trt))
  test.rep[i,1]=fit$`F value`[1]
  test.rep[i,2]=fit$`Pr(>F)`[1]
}
test.rep[,3]=p.adjust(test.rep[,2],method = "BH")


#permuted#
hetero.permut=c()
for (j in 1:100) {
  test.permut=c()
  for (i in 1:10780) {
    temp=data.frame(obs=as.numeric(dat.pre[i,]),trt=sample(group_rep,30,replace = F))
    fit=anova(lm(temp$obs~temp$trt))
    test.permut[i]=fit$`F value`[1]
  }
  hetero.permut=rbind(hetero.permut,test.permut)
  print(j)
}

plot(density(log(hetero.permut[1,])),xlim=c(-4,6),
     main="Distribution of F-value across replicates",lty=2,col="grey",xlab="F-value")
for (i in 2:100) {
  lines(density(log(hetero.permut[i,])),add=T,lty=2,col="grey")
}
lines(density(log(test.rep[,1])),add=T,lty=2,col="red",lwd=1.5)
legend("topright",bty = "n",legend = c("oberved","permuted"),col = c("red","grey"),lty=2,lwd = 2)




#metabolites
setwd("/Volumes/cluster/Wei-Yun/metabolomic_data/")
dat.met=read.csv(file = "metabolome.csv")
dat.met.bhqc=dat.met[,c(1:8,24:41,50:57)]
log.fc.base=data.frame(base=apply(dat.met.bhqc[,c(4:8)],1,mean))
cv.base.met=data.frame(base=apply(dat.met.bhqc[,c(4:8)],1,function(x) var(x)/(mean(x))^2))
cv.base.met.sqrt=sqrt(cv.base.met)
dat.pre=data.frame(log2(dat.met.bhqc[,c(9:26)]/log.fc.base$base))
group_rep=substr(colnames(dat.pre),46,49)
dat.logfc.met.obs=data.frame(rep4=apply(dat.pre[,1:3],1,mean),
                         rep5=apply(dat.pre[,4:6],1,mean),
                         rep6=apply(dat.pre[,7:9],1,mean),
                         rep7=apply(dat.pre[,10:12],1,mean),
                         rep8=apply(dat.pre[,13:15],1,mean),
                         rep9=apply(dat.pre[,16:18],1,mean))
cor.logfc.met.obs=cor(dat.logfc.met.obs)[upper.tri(cor(dat.logfc.met.obs))]

cor_display=function(x,y,...) {
  text(0,0,labels=round(cor(x,y,method = "pearson"),2),cex=1.5,...)
}
scatter.plot=function(x,y,...){
  points(x,y,asp=1,pch=19)
  abline(v=0,col="grey")
  abline(h=0,col="grey")
}
png("../../meta_logfc/pairwise_scatter_plot_logFC_met_pearson.png",height = 12.5,width = 12.5,units = "cm",res=600)
pairs(dat.logfc.met.obs,upper.panel = scatter.plot,
      lower.panel = cor_display,xlim=c(-3,3),ylim=c(-3,3),asp=1)
dev.off()





cor.logfc.met.permut=list()
for (i in 1:100) {
  ind=sample(1:18,18,replace = F)
  dat.permut=data.frame(rep4=apply(dat.pre[,ind[1:3]],1,mean),
                        rep5=apply(dat.pre[,ind[4:6]],1,mean),
                        rep6=apply(dat.pre[,ind[7:9]],1,mean),
                        rep7=apply(dat.pre[,ind[10:12]],1,mean),
                        rep8=apply(dat.pre[,ind[13:15]],1,mean),
                        rep9=apply(dat.pre[,ind[16:18]],1,mean))
  cor.logfc.met.permut[[i]]=cor(dat.permut)[upper.tri(cor(dat.permut))]
  print(paste0(i," time"))
}

plot(density(cor.logfc.met.permut[[1]]),xlim=c(0.4,1),
     main="Distribution of permuted correlation across replicates",lty=2,col="grey",xlab="correlation coefficient",ylim=c(0,15))
for (i in 2:100) {
  lines(density(cor.logfc.met.permut[[i]]),add=T,lty=2,col="grey")
}
lines(density(cor.logfc.met.obs),add=T,lty=2,col="red",lwd=1.5)
legend("topright",bty = "n",legend = c("oberved","permuted"),col = c("red","grey"),lty=2,lwd = 2)



####replicate specific effects####
test.met.rep=matrix(NA,nrow=940,ncol = 3)
row.names(test.met.rep)=dat.met.bhqc[,1]
for (i in 1:940) {
  temp=data.frame(obs=as.numeric(dat.pre[i,]),trt=group_rep)
  fit=anova(lm(temp$obs~temp$trt))
  test.met.rep[i,1]=fit$`F value`[1]
  test.met.rep[i,2]=fit$`Pr(>F)`[1]
}
test.met.rep[,3]=p.adjust(test.met.rep[,2],method = "BH")



hetero.permut=c()
for (j in 1:100) {
  test.permut=c()
  for (i in 1:940) {
    temp=data.frame(obs=as.numeric(dat.pre[i,]),trt=sample(group_rep,18,replace = F))
    fit=anova(lm(temp$obs~temp$trt))
    test.permut[i]=fit$`F value`[1]
  }
  hetero.permut=rbind(hetero.permut,test.permut)
  print(j)
}

plot(density(log(hetero.permut[1,])),xlim=c(-4,6),ylim=c(0,0.55),
     main="Distribution of F-value across replicates",lty=2,col="grey",xlab="F-value")
for (i in 2:100) {
  lines(density(log(hetero.permut[i,])),lty=2,col="grey")
}
lines(density(log(test.met.rep[,1])),lty=2,col="red",lwd=1.5)
legend("topright",bty = "n",legend = c("oberved","permuted"),col = c("red","grey"),lty=2,lwd = 2)

lines(density(log(as.numeric(paste(test.rep[,1])))),lty=2,col="red",lwd=1.5)

t.test(log(as.numeric(paste(test.rep[,1]))),log(test.met.rep[,1]))
