library(dplyr)

# Read arguments
# options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

# Read file for normal tissue (tophat+htseq)
norm_s=read.csv(args[1],sep='\t')

# Read file with tumor sample (TCGA htseq workflow)
tumor_s=read.csv(args[2],sep='\t')

# Merge two tables
merged_s=merge(norm_s,tumor_s,by=1,all=F)
# print(merged_s)

# Write result into two tables
cols1=c(1,2)
# cols2=c(1,3)

write.table(select(merged_s,1,2),args[3],sep='\t',col.names=F,row.names=F)
# write.table(merged_s[,cols2],args[4],sep='\t',col.names=F,row.names=F)