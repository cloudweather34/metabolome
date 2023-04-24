#permutation-DE genes

####JI_empirical observation####
dat_B_H=read.csv(file = "~/Dropbox (PopGen)/backup/Wei-Yun/monster_fullset_readcounts.csv",sep = ";")
colnames(dat_B_H)[1]="name"
dat_B_H_filtered=dat_B_H[apply(cpm(dat_B_H[,-1]),1,function(x){!any(x<=1)}),]
group=c("H1.1","H2.1","H3.1","H4.1","H5.1","H6.1","H7.1","H8.1","H9.1","H10.1","B27","B30",
        "H1.2","H2.2","H3.2","H4.2","H5.2","H6.2","H7.2","H8.2","H9.2","H10.2","B28","H1.3",
        "H2.3","H3.3","H4.3","H5.3","H6.3","H7.3","H9.3","H10.3","B26","B29","H8.3")

row.names(dat_B_H_filtered)=dat_B_H_filtered[,1]
dat_B_H_filtered=dat_B_H_filtered[,-1]
dat_B_H_filtered_log=log(dat_B_H_filtered)
group_rep=paste(substr(colnames(dat_B_H_filtered),1,1),substr(colnames(dat_B_H_filtered),8,8),sep = "_")
group_rep[c(11,12,23,33,34)]="B"
y=DGEList(counts = dat_B_H_filtered,group = group_rep)
ModelDesign_rep=model.matrix(~0+group_rep)
DGE_rep=estimateDisp(y,design = ModelDesign_rep,robust = T)
GLM_rep=glmFit(DGE_rep,design = ModelDesign_rep)
res_table=list()
res_logfc=data.frame(matrix(NA,10780,10))
sig_ID_all=list()
#sig_ID_up=list()
#sig_ID_dn=list()

for (i in 1:10){
  my.contrasts_m=rep(0,length(unique(group_rep)))
  my.contrasts_m[c(1,i+1)]=c(-1,1)
  LRT_m=glmLRT(GLM_rep,contrast = my.contrasts_m)
  res_m=LRT_m$table
  row.names(res_m)=row.names(dat_B_H_filtered)
  res_m$padj=p.adjust(res_m$PValue,method = "BH")
  res_table[[i]]=res_m
#  sig_ID_up[[i]]=row.names(res_m[res_m$padj<0.05&res_m$logFC>0,])
#  sig_ID_dn[[i]]=row.names(res_m[res_m$padj<0.05&res_m$logFC<0,])
  sig_ID_all[[i]]=row.names(res_m[res_m$padj<0.05,])
  res_logfc[,i]=res_m$logFC
  #  write.table(res_m[res_m$padj<0.05,],file = paste0("/Users/weiyun/Dropbox (PopGen)/Wei-Yun (1)/manuscript_targetofselection/table/DE results/DE_",i,".txt"),quote = F,sep = "\t")
}

ind=combinations(10,2,set=TRUE, repeats.allowed=FALSE)
ja_de=c()
for (i in 1:45) {
  inte=length(intersect(unlist(sig_ID_all[ind[i,1]]),unlist(sig_ID_all[ind[i,2]])))
  uni=length(union(unlist(sig_ID_all[ind[i,1]]),unlist(sig_ID_all[ind[i,2]])))
  ja_de[i]=inte/uni
}

cor_de=cor(res_logfc)[upper.tri(cor(res_logfc))]

####JI_permutation####

dat_B_H=read.csv(file = "~/Dropbox (PopGen)/backup/Wei-Yun/monster_fullset_readcounts.csv",sep = ";")
colnames(dat_B_H)[1]="name"
dat_B_H_filtered=dat_B_H[apply(cpm(dat_B_H[,-1]),1,function(x){!any(x<=1)}),]
group=c("H1.1","H2.1","H3.1","H4.1","H5.1","H6.1","H7.1","H8.1","H9.1","H10.1","B27","B30",
        "H1.2","H2.2","H3.2","H4.2","H5.2","H6.2","H7.2","H8.2","H9.2","H10.2","B28","H1.3",
        "H2.3","H3.3","H4.3","H5.3","H6.3","H7.3","H9.3","H10.3","B26","B29","H8.3")
row.names(dat_B_H_filtered)=dat_B_H_filtered[,1]
dat_B_H_filtered=dat_B_H_filtered[,-1]
dat_B_H_filtered_log=log(dat_B_H_filtered)
group_rep=paste(substr(colnames(dat_B_H_filtered),1,1),substr(colnames(dat_B_H_filtered),8,8),sep = "_")
group_rep[c(11,12,23,33,34)]="B"

ja_gene_null=list()
cor_gene_null=list()

for (j in 1:100) {
print(j)

ind1=group_rep[sample(c(1:10,13:22,24:32,35))]
ind2=c(c(1:10,13:22,24:32,35))
for (i in 1:30) {
  group_rep[ind2[i]]=ind1[i]
}

print(group_rep)
print(table(group_rep))

y=DGEList(counts = dat_B_H_filtered,group = group_rep)
ModelDesign_rep=model.matrix(~0+group_rep)
DGE_rep=estimateDisp(y,design = ModelDesign_rep,robust = T)
GLM_rep=glmFit(DGE_rep,design = ModelDesign_rep)
res_logfc=data.frame(matrix(NA,10780,10))
sig_ID_all=list()
#sig_ID_up=list()
#sig_ID_dn=list()

for (i in 1:10){
  my.contrasts_m=rep(0,length(unique(group_rep)))
  my.contrasts_m[c(1,i+1)]=c(-1,1)
  LRT_m=glmLRT(GLM_rep,contrast = my.contrasts_m)
  res_m=LRT_m$table
  row.names(res_m)=row.names(dat_B_H_filtered)
  res_m$padj=p.adjust(res_m$PValue,method = "BH")
  #  sig_ID_up[[i]]=row.names(res_m[res_m$padj<0.05&res_m$logFC>0,])
  #  sig_ID_dn[[i]]=row.names(res_m[res_m$padj<0.05&res_m$logFC<0,])
  sig_ID_all[[i]]=row.names(res_m[res_m$padj<0.05,])
  res_logfc[,i]=res_m$logFC
  #  write.table(res_m[res_m$padj<0.05,],file = paste0("/Users/weiyun/Dropbox (PopGen)/Wei-Yun (1)/manuscript_targetofselection/table/DE results/DE_",i,".txt"),quote = F,sep = "\t")
}

ind=combinations(10,2,set=TRUE, repeats.allowed=FALSE)
ja_de_per=c()
for (i in 1:45) {
  inte=length(intersect(unlist(sig_ID_all[ind[i,1]]),unlist(sig_ID_all[ind[i,2]])))
  uni=length(union(unlist(sig_ID_all[ind[i,1]]),unlist(sig_ID_all[ind[i,2]])))
  ja_de_per[i]=inte/uni
}

ja_gene_null[[j]]=ja_de_per
cor_gene_null[[j]]=cor(res_logfc)[upper.tri(cor(res_logfc))]
}

png("~/Dropbox (PopGen)/Wei-Yun (1)/manuscript_metabolome/figure/fig_s1.png",width = 14,height = 7,units = "cm",pointsize = 6,res = 600)
par(mfrow=c(1,2))
plot(density(sapply(ja_gene_null,mean)),xlim=c(0.2,0.4),col="grey",main="",xlab="Jaccard Index",lty=2)
abline(v=mean(ja_de),col="red")
p.val.ja=sum(mean(ja_de)>sapply(ja_gene_null,mean))/100
legend("topright",legend = c("observed","permuted"),col=c("red","grey"),lty = c(1,2),bty="n")

plot(density(sapply(cor_gene_null,mean)),xlim=c(0.6,0.8),col="grey",main="",xlab="Pearson's correlation coefficient",lty=2)
abline(v=mean(cor_de),col="red")
p.val.cor=sum(mean(cor_de)>sapply(cor_gene_null,mean))/100
legend("topright",legend = c("observed","permuted"),col=c("red","grey"),lty = c(1,2),bty="n")

dev.off()
