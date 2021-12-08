# ERBB4_NTRK_fusions
This repository contains all scripts used for the analysis of differentially expressed genes in thyroid cancer samples with NTRK fusions.

NTRK_all_fusion_samples.thyroid_cancer.csv - table with list of TCGA sample IDs with NTRK fusions
deseq2_htseq_counts.R - R script that was used for the analysis of differentially expressed genes 
filterHtseqCounts.sh - bash script that was used for the intersection of transcript list for the tumor and normal tissue samples
intersect_HTSeq_data.R - R script that was used for the intersection of transcript list for the tumor and normal tissue samples
removeTranscriptVersions.sh - bash script that removes trascript version from the htseq-counts files
