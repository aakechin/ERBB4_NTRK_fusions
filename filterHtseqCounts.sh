dir='thyroid_cancer'
mkdir $dir/htseq_counts_without_trv_filter/
for i in $dir/htseq_counts_without_trv/*.gz

	do
	
	f=$(basename -- $i)
	f=${f%.*}

	Rscript intersect_HTSeq_data.R $i normal_tissues/reads_tophat_htseq/thyroid/1_12_accepted_hits.bam.htseq.counts.gz $dir/htseq_counts_without_trv_filter/${f}

	done
