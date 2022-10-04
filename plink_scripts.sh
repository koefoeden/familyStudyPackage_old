# Make relatedness matrix
plink2 --make-rel square --vcf all_chroms.vcf.gz --out relatedness_matrix

# PCA of genetic components
plink2 --pca \
--pfile fam_study \
--maf 0.05 \
--require-info ER2 \
--out fam_study_only_genotyped

