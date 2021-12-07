read_NTRK_data=function(dir) {
	# Saving files to the table
	files=grep('counts',list.files(paste(dir,
	                                        'htseq_counts_without_trv_filter',
	                                        sep='')),value=T)
	ntrk_table=data.frame(sampleName=files,fileName=files)

	# Reading sample_sheet file
	sample_sheet=read.csv(paste(dir,'gdc_sample_sheet.2020-12-01.tsv',sep=''),sep='\t')
	# Matching file names with adding sample IDs
	ntrk_table2=merge(ntrk_table,sample_sheet,by.x='fileName',by.y='File.Name',all=T)

	# Reading clinics file
	clinics_file=read.csv(paste(dir,'clinical.cart.2020-12-01/clinical.tsv',sep=''),sep='\t')
	# Matching sample IDs with adding clinical info
	ntrk_table3=merge(ntrk_table2,clinics_file,by.x='Case.ID',by.y='case_submitter_id',all=T)

	return(ntrk_table3)
}

library('dplyr')
library("DESeq2")
library(ggplot2)
library('tidyr')


inputDir='/media/andrey/WD_4Tb/projects/NTRK_fusions_RSF_2020_05_26/TCGA/HTSeq_data/'
ntrk1_table=read_NTRK_data(paste(dir,'NTRK1/',sep=''))
ntrk2_table=read_NTRK_data(paste(dir,'NTRK2/',sep=''))
ntrk3_table=read_NTRK_data(paste(dir,'NTRK3/',sep=''))

# Merge all thyroid samples with NTRK-fusions
ntrk_table=merge(ntrk1_table,ntrk2_table,all=T)
ntrk_table=merge(ntrk_table,ntrk3_table,all=T)
print(paste('NTRK table size after joining:',
            length(ntrk_table$fileName),sep=' '))
write.csv(ntrk_table,
          './DESeq2/thyroid_vs_ntrk/ntrk_table.csv')

#Extract only thyroid samples, leave only columns upto tissue of origin
ntrk_table_thyroid=unique(ntrk_table[(!is.na(ntrk_table$tissue_or_organ_of_origin) & ntrk_table$tissue_or_organ_of_origin=='Thyroid gland'),1:128])
print(paste('NTRK table size after removing samples with treatment and extracting thyroid:',
            length(ntrk_table_thyroid$fileName),sep=' '))
# ntrk_table_thyroid=ntrk_table

# Read file with NTRK-fusions
ntrk_fusions=read.csv('/media/andrey/WD_4Tb/projects/NTRK_fusions_RSF_2020_05_26/TCGA/NTRK_all_fusion_samples.csv',sep=';')
# Match with thyroid data
ntrk_table_thyroid=merge(ntrk_table_thyroid,ntrk_fusions,by.x='case_id',by.y='SampleID_long',all.x=T,all.y=F)
print(paste('NTRK table size after adding NTRK-fusions info:',
            length(ntrk_table_thyroid$fileName),sep=' '))
write.csv(ntrk_table_thyroid,
          paste(inputDir,'DESeq2/thyroid_vs_ntrk/ntrk_table_thyroid.csv'))

# Read thyroid cancer data
# Saving files to the table
dir=paste(inputDir,'thyroid_tumors/')
files=grep('counts.gz',list.files(paste(dir,'htseq_counts_without_trv',sep='')),value=T)
thyroid_table=data.frame(sampleName=files,fileName=files,
                         Fusions='no',FusedGene='no')
print(paste('Initial thyroid table size:',
            length(thyroid_table$fileName),sep=' '))
# Remove files that are in the NTRK table
thyroid_table=thyroid_table[!(thyroid_table$fileName %in% ntrk_table_thyroid$fileName),]
print(paste('Thyroid table size after removing files with NTRK:',
            length(thyroid_table$fileName),sep=' '))

# Reading sample_sheet file
sample_sheet=read.csv(paste(dir,'thyroid_gdc_sample_sheet.2020-12-01.tsv',sep=''),sep='\t')
# Matching file names with adding sample IDs
thyroid_table=merge(thyroid_table,sample_sheet,by.x='fileName',by.y='File.Name',all=T)

# Reading clinics file
clinics_file=read.csv(paste(dir,'clinical.cart.2020-12-01/clinical.tsv',sep=''),sep='\t')
# Matching sample IDs with adding clinical info
thyroid_table=merge(thyroid_table,clinics_file,by.x='Case.ID',by.y='case_submitter_id',all=T)

# Extract only thyroid samples, leave only columns upto tissue of origin and leave only unique rows
thyroid_table=unique(thyroid_table[(!is.na(thyroid_table$tissue_or_organ_of_origin) & thyroid_table$tissue_or_organ_of_origin=='Thyroid gland'),1:128])
print(paste('Thyroid table size after removing samples duplicated:',
            length(thyroid_table$fileName),sep=' '))
# Replace NA to 'no' for Fusions
thyroid_table$Fusions[is.na(thyroid_table$Fusions)]='no'

# Read normal thyroid gland
dir='/media/andrey/WD_4Tb/projects/NTRK_fusions_RSF_2020_05_26/normal_tissues/reads_tophat_htseq_filter/thyroid/'
files=grep('filter',list.files(dir),value=T)
norm_thyroid_table=data.frame(Sample.ID=files,fileName=files,
                              Fusions='normal',FusedGene='normal',Case.ID=files)
print(paste('Initial normal thyroid table size:',
            length(norm_thyroid_table$fileName),sep=' '))
write.csv(norm_thyroid_table,
          paste(inputDir,'DESeq2/thyroid_vs_ntrk/norm_thyroid_table.csv'))

# Merge NTRK-samples with normal thyroid gland
ntrk_with_normal=merge(ntrk_table_thyroid,norm_thyroid_table,
                       by=intersect(colnames(ntrk_table_thyroid),colnames(norm_thyroid_table)),all=T)
ntrk_with_thyroid_with_normal=merge(thyroid_table,ntrk_with_normal,
                          by=intersect(colnames(thyroid_table),colnames(ntrk_with_normal)),all=T)

# Write to table
# Remove repeated samples
ntrk_with_thyroid_with_normal=distinct(ntrk_with_thyroid_with_normal,ntrk_with_thyroid_with_normal$Case.ID,.keep_all=T)
# Change sample name if it doesn't have it
ntrk_with_thyroid_with_normal$sampleName=ntrk_with_thyroid_with_normal$fileName
# Replace NA to 'no' for Fusions, FusedGene and Protein.Change
ntrk_with_thyroid_with_normal$Fusions[is.na(ntrk_with_thyroid_with_normal$Fusions)]='normal'
ntrk_with_thyroid_with_normal$FusedGene[is.na(ntrk_with_thyroid_with_normal$FusedGene)]='normal'
ntrk_with_thyroid_with_normal$Protein.Change[is.na(ntrk_with_thyroid_with_normal$Protein.Change)]='normal'
# Select necessary columns
ntrk_with_thyroid_with_normal_fusions=ntrk_with_thyroid_with_normal[,c('Sample.ID','fileName','Fusions')]
ntrk_with_thyroid_with_normal_fused_gene=ntrk_with_thyroid_with_normal[,c('Sample.ID','fileName','FusedGene')]
# Change column name
colnames(ntrk_with_thyroid_with_normal_fusions)[1]='sampleName'
colnames(ntrk_with_thyroid_with_normal_fused_gene)[1]='sampleName'
# Write to table
write.csv(ntrk_with_thyroid_with_normal_fusions,
          paste(inputDir,'DESeq2/thyroid_vs_ntrk/ntrk_vs_normal.csv'))
# Analyze data with DESeq2
ntrk_with_thyroid_with_normal_dds_fused=DESeqDataSetFromHTSeqCount(sampleTable=ntrk_with_thyroid_with_normal_fused_gene,
                                                                   directory="/media/andrey/WD_4Tb/projects/NTRK_fusions_RSF_2020_05_26/TCGA/HTSeq_data/thyroid_with_ntrk_filter",
                                                                   design=~FusedGene)
ntrk_with_thyroid_with_normal_dds_fused$FusedGene=relevel(ntrk_with_thyroid_with_normal_dds_fused$FusedGene,ref='normal')
ntrk_with_thyroid_with_normal_dds_fused=DESeq(ntrk_with_thyroid_with_normal_dds_fused,fitType='parametric')
ntrk_with_thyroid_with_normal_dds_fused_res_ntrk1=results(ntrk_with_thyroid_with_normal_dds_fused,contrast=c('FusedGene','NTRK1','normal'))
ntrk_with_thyroid_with_normal_dds_fused_res_ntrk1=ntrk_with_thyroid_with_normal_dds_fused_res_ntrk1[(!is.na(ntrk_with_thyroid_with_normal_dds_fused_res_ntrk1$log2FoldChange) &
                                                                                         !is.na(ntrk_with_thyroid_with_normal_dds_fused_res_ntrk1$pvalue) &
                                                                                         !is.na(ntrk_with_thyroid_with_normal_dds_fused_res_ntrk1$padj) &
                                                                                           ntrk_with_thyroid_with_normal_dds_fused_res_ntrk1$pvalue!='' &
                                                                                           ntrk_with_thyroid_with_normal_dds_fused_res_ntrk1$padj!=''),]
write.csv(as.data.frame(ntrk_with_thyroid_with_normal_dds_fused_res_ntrk1),
          paste(paste(inputDir,'DESeq2/thyroid_vs_ntrk/ntrk1_vs_tumor_and_normal_thyroid_dds_res_ntrk.csv'),
                sep=''))
ntrk_with_thyroid_with_normal_dds_fused_res_ntrk3=results(ntrk_with_thyroid_with_normal_dds_fused,contrast=c('FusedGene','NTRK3','normal'))
ntrk_with_thyroid_with_normal_dds_fused_res_ntrk3=ntrk_with_thyroid_with_normal_dds_fused_res_ntrk3[(!is.na(ntrk_with_thyroid_with_normal_dds_fused_res_ntrk3$log2FoldChange) &
                                                                                                       !is.na(ntrk_with_thyroid_with_normal_dds_fused_res_ntrk3$pvalue) &
                                                                                                       !is.na(ntrk_with_thyroid_with_normal_dds_fused_res_ntrk3$padj) &
                                                                                                       ntrk_with_thyroid_with_normal_dds_fused_res_ntrk3$pvalue!='' &
                                                                                                       ntrk_with_thyroid_with_normal_dds_fused_res_ntrk3$padj!=''),]
write.csv(as.data.frame(ntrk_with_thyroid_with_normal_dds_fused_res_ntrk3),
          paste(paste(inputDir,'DESeq2/thyroid_vs_ntrk/ntrk3_vs_tumor_and_normal_thyroid_dds_res_ntrk.csv'),
                sep=''))
ntrk_with_thyroid_with_normal_dds_fused_res_no=results(ntrk_with_thyroid_with_normal_dds_fused,contrast=c('FusedGene','no','normal'))
ntrk_with_thyroid_with_normal_dds_fused_res_no=ntrk_with_thyroid_with_normal_dds_fused_res_no[(!is.na(ntrk_with_thyroid_with_normal_dds_fused_res_no$log2FoldChange) &
                                                                                                       !is.na(ntrk_with_thyroid_with_normal_dds_fused_res_no$pvalue) &
                                                                                                       !is.na(ntrk_with_thyroid_with_normal_dds_fused_res_no$padj) &
                                                                                                 ntrk_with_thyroid_with_normal_dds_fused_res_no$pvalue!='' &
                                                                                                 ntrk_with_thyroid_with_normal_dds_fused_res_no$padj!=''),]
write.csv(as.data.frame(ntrk_with_thyroid_with_normal_dds_fused_res_no),
          paste(paste(inputDir,'DESeq2/thyroid_vs_ntrk/no_ntrk_vs_tumor_and_normal_thyroid_dds_res_ntrk.csv'),
                sep=''))

#Draw some counts
genes=c("SCN4A","IGF1","PPP2R5B","ACTB","RPH3A","TBP","PASK","HDAC1","FLRT3","DTNA","RPL32","AUTS2","ERBB4","FLRT2","NTRK1","KDR","PRSS1","CEBP","CAMK2A")
trans=c("ENSG00000007314","ENSG00000017427","ENSG00000068971","ENSG00000075624","ENSG00000089169","ENSG00000112592","ENSG00000115687","ENSG00000116478","ENSG00000125848","ENSG00000134769","ENSG00000144713","ENSG00000158321","ENSG00000178568","ENSG00000185070","ENSG00000198400","ENSG00000128052","ENSG00000204983","ENSG00000245848","ENSG00000070808")
tcounts=t(counts(ntrk_with_thyroid_with_normal_dds_fused[trans,],normalized=T,replaced=F)) %>%
  merge(colData(ntrk_with_thyroid_with_normal_dds_fused),.,by='row.names') %>%
  gather(gene,expression,(ncol(.)-length(trans)+1):ncol(.))
for (i in seq_along(genes)) {
  tcounts$gene[tcounts$gene==trans[i]]=genes[i]
}
png(paste(paste(inputDir,'DESeq2/thyroid_vs_ntrk/ntrk_vs_thyroid_vs_normal_'),
          paste(genes,collapse='_'),
          '_counts.png',
          sep=''),
    width=1000,height=1000)
ggplot(tcounts, aes(FusedGene, expression, fill=FusedGene)) +
  geom_boxplot() +
  scale_fill_manual(values=c('coral2','skyblue4','gold3','#AF58BA')) +
  facet_wrap(~gene, scales="free_y") +
  labs(y="Expression",
       fill="FusedGene") +
  theme(axis.text.x = element_text(angle=45,
                                   vjust=0.5,
                                   face="bold",
                                   size=10))
dev.off()