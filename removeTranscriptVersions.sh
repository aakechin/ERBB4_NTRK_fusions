dir='thyroid_cancer'
mkdir $dir/htseq_counts_without_trv/
for i in $dir/htseq_counts/*

	do

	f=$(basename -- $i)

	zcat $i | sed 's|\.[0-9]\+||g' | gzip > $dir/htseq_counts_without_trv/$f

	done
