rm(list=ls())
setwd("/Users/weiyun/Dropbox (PopGen)/backup/Wei-Yun/metabolomic_data/")
dat.met=read.csv(file = "metabolome.csv")
dat.met.bhqc=dat.met[,c(1:8,24:41,50:57)]

colnames(dat.met.bhqc)[-c(1:3)]=c(paste0("B",1:5),paste0("H",rep(4:9,each=3),rep(1:3,6)),paste0("QC",1:8))
pca=prcomp(t(log(dat.met.bhqc[,-c(1:3)])))
pca$sdev^2/sum(pca$sdev^2)
par(mfrow=c(1,1))
plot(pca$x,col=rep(c("forestgreen","salmon","grey70"),c(5,18,8)),pch=19,xlab="PC1 (23%)",ylab="PC2 (11%)")
text(pca$x[,c(1,2)],labels = colnames(dat.met.bhqc)[-c(1:3)],pos = 1,col = rep(c("forestgreen","salmon","grey70"),c(5,18,8)))

met.test=list()
ind=seq(6,23,3)
dat.met.bhqc.use=log(dat.met.bhqc[,-c(1:3)])


for (i in 1:6) {
  temp=matrix(NA,nrow = 940,ncol = 4)
    for (j in 1:940) {
    a=t.test(dat.met.bhqc.use[j,1:5],dat.met.bhqc.use[j,c(ind[i]:(ind[i]+2))])
    temp[j,1]=a$estimate[2]
    temp[j,2]=a$statistic
    temp[j,3]=a$p.value
  }
  temp[,4]=p.adjust(temp[,3],method = "BH")
  met.test[[i]]=temp
}

dat.met.bhqc.use.fc=matrix(NA,nrow = 940,ncol = 6)

for (i in 1:6) {
  dat.met.bhqc.use.fc[,i]=apply(dat.met.bhqc.use[,c(ind[i]:(ind[i]+2))], 1, mean)-apply(dat.met.bhqc.use[,1:5],1,mean)
}
cor(dat.met.bhqc.use.fc)





ja.use=lapply(met.test, function(x) dat.met.bhqc[x[,4]<0.05,1])
sapply(ja.use,length)


library(gtools)
ind=combinations(6,2,set=TRUE, repeats.allowed=FALSE)


ja_met=c()
for (j in 1:15) {
    inte=length(intersect(unlist(ja.use[ind[j,1]]),unlist(ja.use[ind[j,2]])))
    uni=length(union(unlist(ja.use[ind[j,1]]),unlist(ja.use[ind[j,2]])))
    ja_met[j]=inte/uni
}

evo.test=data.frame(matrix(NA,nrow=940,ncol = 4))

for (j in 1:940) {
  a=t.test(dat.met.bhqc.use[j,1:5],dat.met.bhqc.use[j,6:23])
  temp[j,1]=a$estimate[2]
  temp[j,2]=a$statistic
  temp[j,3]=a$p.value
}
temp[,4]=p.adjust(temp[,3])

regulated=read.csv(file = "Report_regulated_BHC_2.csv")

dat.logfc.anov=matrix(NA,nrow = 940,ncol = 18)
for (i in 1:18) {
  dat.logfc.anov[,i]=log2(dat.met.bhqc.use[,(i+5)]/apply(dat.met.bhqc.use[,1:5],1,mean))
}
colnames(dat.logfc.anov)=colnames(dat.met.bhqc.use)[6:23]

p.val=c()
trt=gl(6,3)
for (i in 1:940) {
  p.val[i]=summary(aov(dat.logfc.anov[i,]~trt))[[1]][1,5]
}
p.adj=p.adjust(p.val)

####permutation test####
dat.met.bhqc.use.shuffle.pre=dat.met.bhqc.use[,1:23]
dat.met.bhqc.use.shuffle=dat.met.bhqc.use.shuffle.pre[,c(1:5,sample(6:23,18,replace = F))]
met.test.permute=list()
for (i in 1:6) {
  temp_permute=matrix(NA,nrow = 940,ncol = 4)
  for (j in 1:940) {
    a=t.test(dat.met.bhqc.use.shuffle[j,1:5],dat.met.bhqc.use.shuffle[j,c(ind[i]:(ind[i]+2))])
    temp_permute[j,1]=a$estimate[2]
    temp_permute[j,2]=a$statistic
    temp_permute[j,3]=a$p.value
  }
  temp_permute[,4]=p.adjust(temp_permute[,3])
  met.test.permute[[i]]=temp_permute
}


head(dat.met.bhqc.use)

dat.met.bhqc.use.shuffle.pre=dat.met.bhqc.use[,1:23]
ind_evo=seq(6,23,3)

met.test.permute.emp.p=list()
for (z in 1:6) {
  met.test.permute.evo=list()
  for (i in 1:1000) {
    print(i)
    temp_permute=matrix(NA,nrow = 940,ncol = 3)
    ind_all=c(1:5,ind_evo[z]:(ind_evo[z]+2))
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
for (i in 1:15) {
  inte=length(intersect(unlist(ja.use[ind[i,2]]),unlist(ja.use[ind[i,2]])))
  uni=length(union(unlist(ja.use[ind[i,1]]),unlist(ja.use[ind[i,2]])))
  ja_metabolite_pa[i]=inte/uni
}





####making table####
met_ID_anno=read.csv(file = "metabolome_ID_matching_complete.csv",stringsAsFactors = F,header = T)

for (i in 1:6) {
 tmp=data.frame(matrix(ncol=0,nrow=940))
 tmp[,1:4]=met_ID_anno
 tmp$Avg=apply(dat.met.bhqc.use[,ind[i]:(ind[i]+2)],1,mean)
 tmp$Diff=apply(dat.met.bhqc.use[,ind[i]:(ind[i]+2)],1,mean)-apply(dat.met.bhqc.use[,1:5],1,mean)
 tmp$t_tatistics=apply(dat.met.bhqc.use[,c(1:5,ind[i]:(ind[i]+2))],1,function(x) t.test(x[1:5],x[6:8])$statistic)
 tmp$p.adj=permute.p.adj[[i]]
 write.table(tmp,paste0("/Users/weiyun/Dropbox (PopGen)/backup/Wei-Yun/project4/",i,".txt"),
              col.names = T,row.names = F,quote = F,sep = "\t")
 print(c(sum(tmp$Diff>0&tmp$p.adj<0.05),sum(tmp$Diff<0&tmp$p.adj<0.05)))
}


