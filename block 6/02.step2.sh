#Step 2: Creating an annot file------------------------source activate ldsc------------------
for i in {1..22}
  do
    for j in {1..22}
      do 
        python make_annot.py \
        --bed-file differtest${i}.hg38_to_hg19.txt \
        --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${j}.bim \
        --annot-file differtest${i}.${j}.annot.gz
      done
  done
