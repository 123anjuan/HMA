#Step 1: Transform single cell Atac-seq peaks bed files in from hg38 to hg19
for i in {1..22}
	do
	/liftOver \
     ./differtest_peak_$i.txt \
     hg38ToHg19.over.chain.gz \
     differtest$i.hg38_to_hg19.txt \
     differtest$i.unlifted.txt
    done
