# convert vcf files to bed file on Nidhogg
ls *dose.vcf.gz -1 | \
xargs -P 22 -I {} \
bash -c "plink2 --vcf {} --make-bed --out ../plink_files/{}"

# combine chromosomes
for i in {1..22}; do echo chr$i.dose.vcf.gz >> mergelist.txt; done
plink --merge-list mergelist.txt --make-bed --out all_chromsomes

# flip missnps
cat mergelist.txt | \
xargs -P 22 -I {} \
bash -c "plink --bfile {} --flip all_chromsomes-merge.missnp --make-bed --out {}.flipped"

# Combine chromosomes now excluding biallelic
for i in {1..22}; do echo chr$i.dose.vcf.gz >> mergelist.txt; done
plink --merge-list mergelist.txt --make-bed --out all_chrom_biallelic --exclude all_chromsomes-merge.missnp

# Construct standardized genetic relatedness matrix
./gemma -bfile /emc/cbmr/projects/tqb695/family_study/plink_files/all_chrom_biallelic -gk 2 -o relatedness_matrix


# Simple linear model, performs Wald, likelihood ratio and score test
./gemma -bfile [prefix] -lm 4 -o [prefix]



# Univariate Linear Mixed Models (needs relatedness matrix k)
# performs Wald, likelihood ratio and score test
./gemma -bfile [prefix] -k [filename] -lmm 4 -o [prefix]

# Multivariate linear mixed models
./gemma -bfile [prefix] -k [filename] -lmm 4 \
-n [num1] [num2] [num3] -o [prefix]

#**** UNUSED
# convert vcf files to bed with max 2 alleles
ls *dose.vcf.gz -1 | \
xargs -P 22 -I {} \
bash -c "plink2 --vcf {} --make-bed --max-alleles 2 --out ../plink_files/biallelic/{}"
