
# stolen from methylkit::assocComp
# only important change is the exact=F argument for the wilcox.test

mBase <- mkit_merged
annot <- sample_annotation


scale=TRUE
center=TRUE
mat=percMethylation(mBase) # get matrix
pr=prcomp(mat,scale.=scale,center=center) # get PCA
vars=100*pr$sdev**2/sum(pr$sdev**2) # calc variation explained

# get association p-values using different tests
res=list()

# loop over each Principal Component
for(i in seq(1,8)){
  res[i] <- wilcox.test(split(pr$rotation[,i],annot)[[1]], split(pr$rotation[,i],annot)[[2]], exact=F)$p.value
}

