#!/usr/bin/env Rscript
args=commandArgs()

#################################################
# ATAC-seq and CUT&Tag-seq DBR finding workflow #
#################################################
library("edgeR")
library("NOISeq")

# do TMM normalization
counttable=read.table("edgeR_input_table.txt",header=F)
mobData_data=counttable[,-1]	
tmm_table=tmm(mobData_data, long = 1000, lc = 0, k = 0)
tmm_table =round(tmm_table,digits=0)
new_TMM_count=data.frame(mobData_data[,1],tmm_table)

# prepare for edgeR
group_edgeR=factor(c(rep(0,as.numeric(args[8])),rep(1,as.numeric(args[9]))),label=c(args[6],args[7]))
# create DGE List
wtest_raw=DGEList(counts=new_TMM_count[,-1], group=group_edgeR, genes=new_TMM_count[,1])
# Filter the peaks with CPM > $cpm_value
keep <- rowSums(cpm(wtest_raw)>as.numeric(args[10])) >=sum(as.numeric(args[8])+as.numeric(args[9]))  #delete the the peaks with 0 count in all samples
wtest <- wtest_raw[keep,]
# Normalizing the data (default method is TMM)
# wtest <- calcNormFactors(wtest) 
# Estimate dispersion
wtest<-estimateCommonDisp(wtest, verbose=TRUE)
wtest<-estimateTagwiseDisp(wtest)
# calculate the CPM
wcpm=cpm(wtest$pseudo.counts)
# Differential Analysis (from the tutorial of edgeR)
wet=exactTest(wtest)
wntop=topTags(wet,n=dim(wtest)[1])
wstop=wntop$table[order(rownames(wntop$table)),]
head(wstop)
wsexp=wcpm[order(rownames(wcpm)),]
wnftab=cbind(wstop,wsexp)
wnnftab=wnftab[order(wnftab[,4]),]
wftab=wnnftab
colnames(wftab)=c("genes","logFC","logCPM","PValue","FDR",rep(args[6], as.numeric(args[8])), rep(args[7], as.numeric(args[9])))
# select peaks with log(Fold-Change) > $log2FC and FDR < $qvalue
wftab_selected=wftab[abs(wftab[,2])>as.numeric(args[12]) & wftab[,5]<as.numeric(args[11]),]
#write out the data
write.table(wftab, file = 'edgeR_DBR_full_list.txt', sep = '\t', quote = FALSE, row.names = FALSE)
write.table(wftab_selected, file = 'edgeR_DBR_significant_only.txt', sep = '\t', quote = FALSE, row.names = FALSE)







