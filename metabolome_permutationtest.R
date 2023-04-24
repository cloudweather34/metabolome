setwd("/Users/weiyun/Dropbox (PopGen)/backup/Wei-Yun/project4/")

#metabolome raw data
#rm(list=ls())
dat.met=read.csv(file = "/Users/weiyun/Dropbox (PopGen)/backup/Wei-Yun/metabolomic_data/metabolome.csv")
dat.met.bhqc=dat.met[,c(1:8,24:41,50:57)]
colnames(dat.met.bhqc)[-c(1:3)]=c(paste0("B",1:5),paste0("H",rep(4:9,each=3),rep(1:3,6)),paste0("QC",1:8))

pca=prcomp(t(log(dat.met.bhqc[,-c(1:3)])))
pca$sdev^2/sum(pca$sdev^2)
par(mfrow=c(1,1))
plot(pca$x,col=rep(c("forestgreen","salmon","grey70"),c(5,18,8)),pch=19,xlab="PC1 (23%)",ylab="PC2 (11%)")
text(pca$x[,c(1,2)],labels = colnames(dat.met.bhqc)[-c(1:3)],pos = 1,col = rep(c("forestgreen","salmon","grey70"),c(5,18,8)))

dat.met.bhqc.use=log(dat.met.bhqc[,-c(1:3)])


#metabolome permutation test
dat.met.bhqc.use.shuffle.pre=dat.met.bhqc.use[,1:23]
ind_evo=seq(6,23,3)

set.seed(100)

met=list()
met_diff=list()
for (z in 1:6) {
  print(Sys.time())
  ind_all=c(1:5,ind_evo[z]:(ind_evo[z]+2))
  emp.p=c()
  diff=c()
  
  for (j in 1:940) {
    
    #observed value
    emp_diff=mean(as.numeric(dat.met.bhqc.use.shuffle.pre[j,ind_all[1:5]]))-mean(as.numeric(dat.met.bhqc.use.shuffle.pre[j,ind_all[6:8]]))
    diff[j]=as.numeric(emp_diff)
    
    
    #null  
    null_diff=c()
    for (i in 1:1000) {
      ind_base=sample(ind_all,5,replace = F)
      a=mean(as.numeric(dat.met.bhqc.use.shuffle.pre[j,ind_base]))-mean(as.numeric(dat.met.bhqc.use.shuffle.pre[j,ind_all[!ind_all%in%ind_base]]))
      null_diff[i]=as.numeric(a)
    }
    
    #empirical p-value
    
    emp.p[j]=sum(abs(emp_diff)<abs(null_diff))/1000
  }
  
  emp.p.adj=p.adjust(emp.p,method = "BH")
  
  met_diff[[z]]=diff
  met[[z]]= emp.p.adj
  
  print(z)
  print(Sys.time())
}



#result output
#met_ID_anno=read.csv(file = "~/Dropbox (PopGen)/backup/Wei-Yun/metabolomic_data/metabolome_ID_matching_complete.csv",stringsAsFactors = F,header = T)

#for (i in 1:6) {
#  tmp=data.frame(matrix(ncol=0,nrow=940))
#  tmp[,1:4]=met_ID_anno[1:4]
#  tmp$Avg=apply(dat.met.bhqc.use[,c(1:5,ind_evo[i]:(ind_evo[i]+2))],1,mean)
#  tmp$Diff=apply(dat.met.bhqc.use[,ind_evo[i]:(ind_evo[i]+2)],1,mean)-apply(dat.met.bhqc.use[,1:5],1,mean)
#  tmp$p.adj=met[[i]]
#  write.csv(tmp,paste0("/Users/weiyun/Dropbox (PopGen)/backup/Wei-Yun/project4/",i,".csv"),
#              row.names = F)
#  print(c(sum(tmp$Diff>0&tmp$p.adj<0.05),sum(tmp$Diff<0&tmp$p.adj<0.05)))
#}


