#Step 3: Computing LD scores with an annot file.
for i in {1..22}
  do
    for j in {1..22}
      do 
       python ldsc.py \
       --l2 \
       --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${j} \
       --ld-wind-cm 1 \
       --annot differtest${i}.${j}.annot.gz \
       --thin-annot \
       --print-snps hapmap3_snps/hm.${j}.snp \
       --out differtest${i}.${j}
      done
  done
