
#!/bin/bash

wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20220402.zip
unzip plink_linux_x86_64_20220402.zip

dir=$1
dx cd "/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release"

dx download ukb23158_cX_b0_v1*

./plink --bfile ukb23158_cX_b0_v1 --snps 23:154536002:C:T,23:154534419:G:A --recode A --output-missing-genotype . --out var_plink_chrX_genotype
mv var_plink_chrX_genotype.raw var_plink_chrX_genotype.txt

dx cd $dir
dx upload -r var_plink_chrX_genotype.txt
