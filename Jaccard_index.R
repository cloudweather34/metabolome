setwd("/Users/weiyun/Dropbox (PopGen)/backup/Wei-Yun/metabolomic_data/joint_pathway_anal/")
datbas=readRDS(file = "geneMetIDinCommonSigPathway.rds")
library(pheatmap)


#logFC-DE genes
gene_logfc=data.frame(gene_name=row.names(res_table[[1]])) #res_table from DE_analysis.R
for (i in 1:10) {
  gene_logfc[,(i+1)]=res_table[[i]]$logFC
  colnames(gene_logfc)[i+1]=paste0("evo",i)
}

write.table(gene_logfc,file = "gene_logfc.txt",sep = "\t",quote = F)

#logFC-metabolites
metabolite_logfc=data.frame(name=met_ID_anno$KEGG_ID,
                            evo4=dat.met.bhqc.use.fc[,1],evo5=dat.met.bhqc.use.fc[,2],
                            evo6=dat.met.bhqc.use.fc[,3],evo7=dat.met.bhqc.use.fc[,4],
                            evo8=dat.met.bhqc.use.fc[,5],evo9=dat.met.bhqc.use.fc[,6])

write.table(metabolite_logfc,file = "metabolite_logfc.txt",sep = "\t",quote = F)

#example_pentose phosphate pathway
ppp_gene=datbas$`Pentose phosphate pathway`$genes
ppp_metabolite=datbas$`Pentose phosphate pathway`$metabolites

cor_gene=cor(gene_logfc[(gene_logfc$gene_name%in%ppp_gene),2:11])[upper.tri(cor(gene_logfc[(gene_logfc$gene_name%in%ppp_gene),2:11]))]
cor_metabolite=cor(metabolite_logfc[metabolite_logfc$name%in%ppp_metabolite,2:7])[upper.tri(cor(metabolite_logfc[metabolite_logfc$name%in%ppp_metabolite,2:7]))]

boxplot(cor_gene,cor_metabolite,main="pentose phosphate pathway",names=c("gene","metabolite"),
        ylab="correlation coefficient")

pheatmap(metabolite_logfc[metabolite_logfc$name%in%ppp_metabolite,2:7])
pheatmap(gene_logfc[(gene_logfc$gene_name%in%ppp_gene),2:11])

#example_purine metabolism
pm_gene=sample(datbas$`Purine metabolism`$genes,19)
pm_metabolite=datbas$`Purine metabolism`$metabolites

cor_gene=cor(gene_logfc[(gene_logfc$gene_name%in%pm_gene),2:11])[upper.tri(cor(gene_logfc[(gene_logfc$gene_name%in%pm_gene),2:11]))]
cor_metabolite=cor(metabolite_logfc[metabolite_logfc$name%in%pm_metabolite,2:7])[upper.tri(cor(metabolite_logfc[metabolite_logfc$name%in%pm_metabolite,2:7]))]

boxplot(cor_gene,cor_metabolite,main="purine metabolism",names=c("gene","metabolite"),
        ylab="correlation coefficient")

pheatmap(metabolite_logfc[metabolite_logfc$name%in%pm_metabolite,2:7])
pheatmap(gene_logfc[(gene_logfc$gene_name%in%pm_gene),2:11])

#example_purine metabolism
pm_gene=datbas$`Purine metabolism`$genes
pm_metabolite=datbas$`Purine metabolism`$metabolites

cor_gene=cor(gene_logfc[(gene_logfc$gene_name%in%pm_gene),2:11])[upper.tri(cor(gene_logfc[(gene_logfc$gene_name%in%pm_gene),2:11]))]
cor_metabolite=cor(metabolite_logfc[metabolite_logfc$name%in%pm_metabolite,2:7])[upper.tri(cor(metabolite_logfc[metabolite_logfc$name%in%pm_metabolite,2:7]))]

boxplot(cor_gene,cor_metabolite,main="purine metabolism",names=c("gene","metabolite"),
        ylab="correlation coefficient")

pheatmap(metabolite_logfc[metabolite_logfc$name%in%pm_metabolite,2:7])
pheatmap(gene_logfc[(gene_logfc$gene_name%in%pm_gene),2:11])

#tyrosine metabolism
tm_gene=datbas$`Tyrosine metabolism`$genes
tm_metabolite=datbas$`Tyrosine metabolism`$metabolites

cor_gene=cor(gene_logfc[(gene_logfc$gene_name%in%tm_gene),2:11])[upper.tri(cor(gene_logfc[(gene_logfc$gene_name%in%tm_gene),2:11]))]
cor_metabolite=cor(metabolite_logfc[metabolite_logfc$name%in%tm_metabolite,2:7])[upper.tri(cor(metabolite_logfc[metabolite_logfc$name%in%tm_metabolite,2:7]))]

boxplot(cor_gene,cor_metabolite,main="tyrosine metabolism",names=c("gene","metabolite"),
        ylab="correlation coefficient")

pheatmap(metabolite_logfc[metabolite_logfc$name%in%tm_metabolite,2:7])
pheatmap(gene_logfc[(gene_logfc$gene_name%in%tm_gene),2:11])




####Jaccard Index####
setwd("/Users/weiyun/Dropbox (PopGen)/backup/Wei-Yun/project4/")

#gene
DE_gene=list()
list.files(path = "gene/")

for (i in 3:8) {
  tmp=read.table(file = paste0("gene/","DE_",i,".txt"),header = T,stringsAsFactors = F)
  DE_gene[[i-2]]=row.names(tmp)[tmp$padj<0.05]
}

ind=combinations(6,2,set=TRUE, repeats.allowed=FALSE) # library(gtools)
ja_gene=c()
for (i in 1:15) {
  inte=length(intersect(unlist(DE_gene[ind[i,1]]),unlist(DE_gene[ind[i,2]])))
  uni=length(union(unlist(DE_gene[ind[i,1]]),unlist(DE_gene[ind[i,2]])))
  ja_gene[i]=inte/uni
}

#metabolite
DE_metabolite=list()
list.files(path = "metabolite/")

for (i in 1:6) {
  tmp=read.delim(file = paste0("metabolite/",i,".txt"),stringsAsFactors = F,sep = "\t")
  DE_metabolite[[i]]=which(tmp$p.adj<0.05)
}

ind=combinations(6,2,set=TRUE, repeats.allowed=FALSE)
ja_metabolite=c()
for (i in 1:15) {
  inte=length(intersect(unlist(DE_metabolite[ind[i,2]]),unlist(DE_metabolite[ind[i,2]])))
  uni=length(union(unlist(DE_metabolite[ind[i,1]]),unlist(DE_metabolite[ind[i,2]])))
  ja_metabolite[i]=inte/uni
}

#pathway
DE_pathway=list()
list.files(path = "pathway/")

#ind=names(which(table(unlist(DE_pathway))==6))
for (i in 1:6) {
  tmp=read.delim(file = paste0("pathway/joint_pathway_analysis_R",i+3,".txt"),stringsAsFactors = F,sep = "\t")
  DE_pathway[[i]]=tmp$pathway[tmp$FDR_comb<0.05]
  #DE_pathway[[i]]=tmp[tmp$pathway%in%ind,5]
}

mean_FDR=c()
for (i in 1:29) {
  mean_FDR[i]=mean(DE_pathway[[1]][i],DE_pathway[[2]][i],DE_pathway[[3]][i],DE_pathway[[4]][i],DE_pathway[[5]][i],DE_pathway[[6]][i])
}


ind=combinations(6,2,set=TRUE, repeats.allowed=FALSE)
ja_pathway=c()
for (i in 1:15) {
  inte=length(intersect(unlist(DE_pathway[ind[i,1]]),unlist(DE_pathway[ind[i,2]])))
  uni=length(union(unlist(DE_pathway[ind[i,1]]),unlist(DE_pathway[ind[i,2]])))
  ja_pathway[i]=inte/uni
}


#Fig. 2
dat.met=read.csv(file = "/Volumes/cluster/Wei-Yun/metabolomic_data/metabolome.csv")
dat.met.bh=dat.met[,c(1:8,24:41)]

colnames(dat.met.bh)[-c(1:3)]=c(paste0("B",1:5),paste0("H",rep(4:9,each=3)))
pca=prcomp(t(log(dat.met.bh[,-c(1:3)])))
pca$sdev^2/sum(pca$sdev^2)


png(filename = "Fig2.png",width = 20,height = 10,units = "cm",res = 600,pointsize = 9)
par(mfrow=c(1,2))
plot(pca$x,col=rep(c("forestgreen","salmon"),c(5,18)),pch=19,xlab="PC1 (23%)",ylab="PC2 (11%)")
text(pca$x[,c(1,2)],labels = colnames(dat.met.bh)[-c(1:3)],pos = 1,col = rep(c("forestgreen","salmon"),c(5,18)))
boxplot(ja_gene,ja_metabolite,names=c("gene","metabolite"),ylab="Jaccard Index")
dev.off()


#fig.3
table(table(unlist(DE_pathway)))

png(filename = "FigS5.png",width = 10,height = 10,units = "cm",res = 600,pointsize = 9)
par(mfrow=c(1,1))
barplot(table(table(unlist(DE_pathway))),main = "Replicate frequency specturm",
        ylab = "number of enriched pathway",xlab = "number of population")
dev.off()

#significant pathway across six replicates
pathway_6=names(table(unlist(DE_pathway)))[table(unlist(DE_pathway))==6]
pathway_6.use=pathway_6[pathway_6%in%names(datbas)]

ind=c()
for (i in 1:27) {
  ind[i]=length(datbas[[pathway_6.use[[i]]]]$metabolites)
}



pathway_6_met.use=pathway_6.use[which(ind>5)]


metabolite_logfc=data.frame(name=met_ID_anno$KEGG_ID,dat.met.bhqc.use.fc)

path_cor_gene=list()
path_cor_metabolite=list()
for (i in 1:length(pathway_6_met.use)) {
  tmp_gene=gene_logfc[gene_logfc$gene_name%in%datbas[[pathway_6_met.use[[i]]]]$genes,5:10]
  path_cor_gene[[i]]=cor(tmp_gene)[upper.tri(cor(tmp_gene))]
  
  tmp_metabolite=metabolite_logfc[metabolite_logfc$name%in%datbas[[pathway_6_met.use[[i]]]]$metabolites,2:7]
  path_cor_metabolite[[i]]=cor(tmp_metabolite)[upper.tri(cor(tmp_metabolite))]
}

#lines(density(cor(gene_logfc[,2:11])[upper.tri(cor(gene_logfc[,2:11]))]),lty=2,col="red")
for (i in 1:27) {
  print(length(datbas[[pathway_6.use[[i]]]]$metabolites))
}


#ja_genes

ja_snp=c(0.44,0.37,0.48,0.36,0.36,0.35,0.32,0.37,0.24,0.33,0.32,0.39,0.43,0.38,0.29)



