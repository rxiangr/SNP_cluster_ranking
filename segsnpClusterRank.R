#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<7) {
  stop("7 argument must be supplied, \nUsage: Rscript --vanilla segsnpClusterRank.R <../yourpath/filename_of_multi_trait_snp_effect.gz> <../yourpath/filename_ranking_of_snps.gz> <../yourpath/filename_of_plink_ld_file_of_segment_snps.ld> <ld_cutoff_for_clustering> <proportion_of_top_snps_to_keep> <../yourpath/prefix_of_output> <Number_of_cores_for_analysis>", call.=FALSE)
}

library(data.table)
library(WGCNA)
library('igraph')

#---define input file path/names
snpefffn <- args[1] #(e.g., test.snpeff.txt.gz)
snprankfn <- args[2] #(e.g., test.snprank.txt.gz)
snpldfn <- args[3] #(e.g., test.ld)
#---define variables
ldcut <- as.numeric(args[4]) #(e.g., 0.95)
rkperc <- as.numeric(args[5]) #(e.g., 0.5)
outputpref <- args[6] #(e.g., test)
Ncore <- as.numeric(args[7]) #(e.g., 10)
#---potential variable to edit
NClustCut <- 50 #(number of clusters to cut, currently set to defualt of 50)
if (Ncore>1) {allowWGCNAThreads(Ncore)}

#---read in data
all.t <- fread(snpefffn,header=T)
all.t3 <- setDF(all.t)[,-1]
rownames(all.t3) <- all.t[,1]
snprank <- fread(snprankfn,header=T)
#---we need to 3rd, 6th and the 7th column from the ld file
ld <- fread(snpldfn,header=T)[,c(3,6,7)]
colnames(ld) <- c('Var1','Var2','ld')
cat(paste0('Data reading finished at ',Sys.time()),sep='\n')
setkey(ld,Var1,Var2)
#---calculate correlation matrix between effects of targed SNPs
cormat <- WGCNA::cor(t(all.t3[rownames(all.t3) %in% unique(c(ld[,unique(Var1)],ld[,unique(Var2)])),]),use='pairwise.complete.obs')
cat(paste0('Correlation matrix of snp effects calculated at ',Sys.time()),sep='\n')
#---reshape the correlation matrix (to merge with the plink LD table)
if (nrow(cormat)^2<(2^31-1)) {
cat(paste0('nrow(cormat)^2<(2^31-1) counted at ',Sys.time()),sep='\n')
cormat.t <- setDT(reshape2::melt(cormat))
colnames(cormat.t)[3] <- 'cor'
setkey(cormat.t,Var1,Var2)
ld[cormat.t,p.cor:=i.cor]
ld[,pr:=ld*p.cor]
} else {
cat(paste0('nrow(cormat)^2>(2^31-1) counted at ',Sys.time()),sep='\n')
corlist <- list()
for (i in seq(1,nrow(ld))){
var1<- unlist(ld[i,1])
var2<- unlist(ld[i,2])
corlist[[i]] <- cormat[which(rownames(cormat)==var1),which(colnames(cormat)==var2)]
}
ld[,p.cor:=cbind(corlist)]
ld[, p.cor:=as.numeric(p.cor)]
ld[,pr:=ld*p.cor]
}
cat(paste0('LD and snp effect correlation combined at ',Sys.time()),sep='\n')

#----1st part: if SNPs have too much LD (e.g., ld2>0.95), they are then ranked and selected directly
dat <- ld[ld^2>=ldcut]
dat1 <- merge(snprank,dat[,1:5],by.x='SNP',by.y='Var1')
colnames(dat1)[1] <- 'Var1'
colnames(dat1)[2] <- 'Var1.rank'
dat2 <- merge(snprank,dat1[,1:6],by.x='SNP',by.y='Var2')
colnames(dat2)[1] <- 'Var2'
colnames(dat2)[2] <- 'Var2.rank'
dat2.1 <- dat2
dat2.1[,diff:=Var2.rank-Var1.rank]
dat2.1[,selected:=ifelse(diff<0,Var2,ifelse(diff>0,Var1,Var1))]
resk1 <- setDT(unique(merge(setDF(dat2.1[,9]),snprank,by.x='selected',by.y='SNP')))
resk1[,rk:= frank(snpRank,ties.method='dense')/length(snpRank)]
resk1.1 <- resk1[rk<=rkperc][,-3]
pld <- unique(dat2.1[diff==0][,c(3,4,1,2,5:7)])
colnames(pld)[1] <- 'selected'
pld1 <- pld[unlist(pld[,1]) %in% unlist(resk1.1[,1])]
flist <- unique((c(unlist(pld1[,1]),unlist(pld1[,3]))))
cat(paste0('Part 1 of clustering finished at ',Sys.time()),sep='\n')
#2nd part: the rest of the SNPs (less LD than the set cutoff) participate in the clustering analysis, and then rank and select
fgl <- graph_from_data_frame(setDF(ld[ld^2<ldcut])[,c(1,2,5)],directed=F,vertices=NULL)
fgl1<- simplify(fgl,remove.multiple=T,remove.loops=T, edge.attr.comb = igraph_opt('edge.attr.comb'))
glclus <- cluster_walktrap(fgl1,merges=T,modularity=T,membership=T,weights=(E(fgl1)$pr))
lhc <- as.hclust(glclus)
if (length(lhc$merge[,1])+1>=NClustCut){
k <- data.frame(cutree(lhc, k = NClustCut))
colnames(k)[1] <- 'k'
data2.1 <- merge(k,snprank,by.x='row.names',by.y='SNP')
data2.2 <- data2.1[!(unlist(data2.1[,1]) %in% flist),]
colnames(data2.2)[1] <- 'SNP'
setDT(data2.2)
data2.2[,rk:=frank(snpRank,ties.method='dense')/length(snpRank),by=k]
resk2 <- data2.2[rk<=rkperc][,c(1,grep('snpRank',colnames(data2.2))),with=F]
colnames(resk2)[1] <- 'selected'
resk2[,selected:=as.character(selected)]
reskall <- unique(rbind(resk1.1,resk2))
} else {
k <- data.frame(cutree(lhc, k = length(lhc$merge[,1])+1))
colnames(k)[1] <- 'k'
data2.1 <- merge(k,snprank,by.x='row.names',by.y='SNP')
data2.2 <- data2.1[!(unlist(data2.1[,1]) %in% flist),]
colnames(data2.2)[1] <- 'SNP'
setDT(data2.2)
data2.2[,rk:=frank(snpRank,ties.method='dense')/length(snpRank),by=k]
resk2 <- data2.2[rk<=rkperc][,c(1,grep('snpRank',colnames(data2.2))),with=F]
colnames(resk2)[1] <- 'selected'
resk2[,selected:=as.character(selected)]
reskall <- unique(rbind(resk1.1,resk2))
}
cat(paste0('Part 2 of clustering finished at ',Sys.time()),sep='\n')

#three files tooutput:
write.table(reskall,gzfile(paste0(outputpref,'.selectsnp.txt.gz')),row.names=F,quote=F,sep=' ')
cat(paste0('Results of selected top cluster SNPs saved to ',paste0(outputpref,'.selectsnp.txt.gz')),sep='\n')
write.table(pld1,gzfile(paste0(outputpref,'.perfldsnp.txt.gz')),row.names=F,quote=F,sep=' ')
cat(paste0('Results SNPs with perfect LD saved to ',paste0(outputpref,'.perfldsnp.txt.gz')),sep='\n')
write.table(data2.2,gzfile(paste0(outputpref,'.clustsnp.txt.gz')),row.names=F,quote=F,sep=' ')
cat(paste0('Detailed results of SNP clusters saved to ',paste0(outputpref,'.clustsnp.txt.gz')),sep='\n')
