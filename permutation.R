save.image("permutation.RData")
load("permutation.RData")

head(dat.met.bhqc.use)

ind=seq(6,21,3)
ind2=combinations(6,2,repeats.allowed = F)

PC2.met=c()
for (i in 1:15) {
  x1=ind[ind2[i,1]]
  x2=ind[ind2[i,2]]
  pca=prcomp(t(dat.met.bhqc.use[,c(1:5,x1:(x1+2),x2:(x2+2))]))
  temp=pca$sdev^2/sum(pca$sdev^2)
  PC2.met[i]=sum(temp[-c(1,2)])
}

expression.use=log(cpm(y))
ind3=unique(group_rep)[1:10]
ind4=combinations(10,2,repeats.allowed = F)

PC2.gene=c()
base=which(group_rep=="B")
for (i in 1:45) {
  x1=which(group_rep==ind3[ind4[i,1]])
  x2=which(group_rep==ind3[ind4[i,2]])
  pca=prcomp(t(expression.use[,c(base,x1,x2)]))
  temp=pca$sdev^2/sum(pca$sdev^2)
  PC2.gene[i]=sum(temp[-c(1,2)])
}

mean(PC2.gene)
mean(PC2.met)

#fstatistics-gene

dat_use=log(cpm(y))#y from DE analysis
fstat_gene = data.frame(matrix(NA,ncol = 2,nrow = 10780))
b_idx=group_rep%in%"B"
dat_fc=log2(dat_use[,!b_idx]/apply(dat_use[,b_idx],1,mean))
group_rep_evo=paste(substr(colnames(dat_fc),1,1),substr(colnames(dat_fc),8,8),sep = "_")
test.rep=matrix(NA,nrow=10780,ncol = 2)
row.names(test.rep)=row.names(dat_fc)
for (i in 1:10780) {
  temp=data.frame(obs=as.numeric(dat_fc[i,]),trt=group_rep_evo)
  fit=anova(lm(temp$obs~temp$trt))
  test.rep[i,1]=fit$`Sum Sq`[1]/sum(fit$`Sum Sq`)
  test.rep[i,2]=fit$`Mean Sq`[2]
}

#fstatistics-metabolites

dat_use_met=dat.met.bhqc.use[,1:23]
b_idx=1:5
dat_fc_met=dat_use_met[,6:23]-apply(dat_use_met[,b_idx],1,mean)
group_rep_evo=substr(colnames(dat_fc_met),1,2)
test.rep.met=matrix(NA,nrow=940,ncol = 2)

row.names(test.rep.met)=row.names(dat_fc_met)
for (i in 1:940) {
  temp=data.frame(obs=as.numeric(dat_fc_met[i,]),trt=group_rep_evo)
  fit=anova(lm(temp$obs~temp$trt))
  test.rep.met[i,1]=fit$`Sum Sq`[1]/sum(fit$`Sum Sq`)
  test.rep.met[i,2]=fit$`Mean Sq`[2]
}

#permutation

#gene

dat.use=cpm(y)
B=which(group_rep=="B")
H=which(group_rep!="B")
permuted.gene.cor=c()
for (j in 1:1000) {
  gene_logfc_permuted=data.frame(matrix(data = NA,nrow = 10780,ncol = 10))
  for (i in 1:10) {
    gene_logfc_permuted[,i]=log2(apply(dat.use[,sample(H,3)], 1, mean)/apply(dat.use[,B], 1, mean))
  }
  permuted.gene.cor[j]=mean(cor(gene_logfc_permuted)[upper.tri(cor(gene_logfc_permuted))])
  print(j)
  }


#metabolite
#ja

dat.met.bhqc.use.shuffle.pre=dat.met.bhqc.use[,1:23]
ind_evo=sample(6:23)


ja.per.met=c()
for (c in 1:100) {
  print(c)
  ind_per=sample(ind_evo,3)
  met.test.permute.emp.p=list()
  for (z in 1:6) {
    met.test.permute.evo=list()
    for (i in 1:1000) {
      print(i)
      temp_permute=matrix(NA,nrow = 940,ncol = 3)
      ind_all=c(1:5,ind_per)
      ind_base=sample(ind_all,5,replace = F)
      for (j in 1:940) {
        a=t.test(dat.met.bhqc.use.shuffle.pre[j,ind_base],dat.met.bhqc.use.shuffle.pre[j,ind_all[!ind_all%in%ind_base]])
        temp_permute[j,1]=a$estimate[2]
        temp_permute[j,2]=a$statistic
        temp_permute[j,3]=a$p.value
      }
      met.test.permute.evo[[i]]=temp_permute
    }
    t.permute.evo=sapply(met.test.permute.evo,function(x) x[,2])
    emp.p=c()
    for (j in 1:940){
      emp.p[j]=sum(t.permute.evo[j,]>abs(met.test[[z]][j,2]))/100
    }
    met.test.permute.emp.p[[z]]=emp.p
  }
  
  
  permute.p.adj=lapply(met.test.permute.emp.p, function(x) p.adjust(x,"BH"))
  ja.use=lapply(permute.p.adj, function(x)which(x<0.05))
  
  
  ind=combinations(6,2,set=TRUE, repeats.allowed=FALSE)
  ja_metabolite_pa=c()
  for (b in 1:15) {
    inte=length(intersect(unlist(ja.use[ind[b,2]]),unlist(ja.use[ind[b,2]])))
    uni=length(union(unlist(ja.use[ind[b,1]]),unlist(ja.use[ind[b,2]])))
    ja_metabolite_pa[b]=inte/uni
  }
  
  ja.per.met[c]=mean(ja_metabolite_pa)
}

ja.per.met=ja.per.met-0.2799251



#permutation_cor
dat.use=dat.met.bhqc.use
B=1:5
H=6:23
permuted.cor=c()
for (j in 1:100) {
  meta_logfc_permuted=data.frame(matrix(data = NA,nrow = 940,ncol = 6))
  for (i in 1:6) {
    meta_logfc_permuted[,i]=log2(apply(dat.use[,sample(H,3)], 1, mean)/apply(dat.use[,B], 1, mean))
  }
  permuted.cor[j]=mean(cor(meta_logfc_permuted)[upper.tri(cor(meta_logfc_permuted))])
  print(j)
}


png("~/Dropbox (PopGen)/Wei-Yun (1)/manuscript_metabolome/figure/fig_s2.png",width = 14,height = 7,units = "cm",pointsize = 6,res = 600)
par(mfrow=c(1,2))
plot(density(ja.per.met),xlim=c(0.25,0.45),col="grey",main="",xlab="Jaccard Index",lty=2)
abline(v=mean(ja_metabolite),col="red")
p.val.ja=sum(ja.per.met<mean(ja_metabolite))/100
legend("topright",legend = c("observed","permuted"),col=c("red","grey"),lty = c(1,2),bty="n")

plot(density(permuted.cor),xlim=c(0.7,1),col="grey",main="",xlab="Pearson's correlation coefficient",lty=2)
abline(v=mean(cor.logfc.met.obs),col="red")
p.val.cor=sum(mean(cor.logfc.met.obs)>permuted.cor)/100
legend("topright",legend = c("observed","permuted"),col=c("red","grey"),lty = c(1,2),bty="n")

dev.off()


