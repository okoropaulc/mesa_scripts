#prepare data for imputation with predixcan EN
library(data.table)

#check the dosage files header to be sure it in the format
#chromosome rsid position allele1 allele2 MAF id1 ..... idn.
caupheno1 <- fread(file = "Z:/data/mesa_models/mesa_pheno/thrombotic/cau_dosage_chr1_all_sample.txt", header=T, nrows=3)

#use the bash script below to gzip the dosage files. -c is to keep the original files
# for i in {1..22}
# do
# echo $i
# gzip -c cau_dosage_chr${i}_all_sample.txt > cau_dosage_chr${i}_all_sample.txt.gz
# done

#--dosages /home/pokoro/data/mesa_models/mesa_pheno/thrombotic/predixcan_imputation/#cau_dosage_chr1_all_sample.txt.gz
#--samples /home/pokoro/data/mesa_models/mesa_pheno/thrombotic/predixcan_imputation/cau_rankplt5_samples.txt
#--dosages_prefix cau
#--weights /home/pokoro/mesa_models/all/ALL_unfiltered_cpos_hg38.db
#--output_prefix home/pokoro/data/mesa_models/mesa_pheno/thrombotic/predixcan_imputation/en_rankplt5
#--assoc
#--linear
#--pheno /hom/pokoro/data/mesa_models/mesa_pheno/thrombotic/predixcan_imputation/pheno_rankplt5.txt

#make the sample file. the sample id in the pheno must be in the same order they appear in the dosage
thrombomodulin <- fread(file = "Z:/data/mesa_models/mesa_pheno/thrombotic/MESA_thrombotic_phenotypes.txt", header=T)

#also remake the pheno file to have cols IID, pheno only
library(tidyverse)
library(dplyr)

#drop NA
#keep ID and rankplt5
thrombomodulin <- thrombomodulin[,c(1,5)]
thrombomodulin <- drop_na(thrombomodulin) #drop NA
thrombomodulin$sidno <- as.character(thrombomodulin$sidno)
colnames(thrombomodulin)[2] <- "phenotype" #rank_plt5
colnames(thrombomodulin)[1] <- "IID" #rank_plt5
thrombomodulin$IID <- as.character(thrombomodulin$IID)

dos_id_col <- as.data.frame(colnames(caupheno1)[6:length(caupheno1)])
names(dos_id_col) <- "IID"
dos_id_col$IID <- as.character(dos_id_col$IID)
merge_id <- inner_join(dos_id_col, thrombomodulin, by = c("IID" ="IID"))

#make the sample list to have two columns FID and IID
dos_id_col <- cbind(dos_id_col, dos_id_col)
names(dos_id_col) <- c("FID", "IID")

#predixcan does not want the sample file to have columns
fwrite(dos_id_col, file="Z:/data/mesa_models/mesa_pheno/thrombotic/predixcan/cau_rankplt5_samples.txt", row.name=F, quote=F, sep="\t", col.names=F)


#make the pheno file in the format required for predixcan FID IID phenotype
pheno_file <- cbind(dos_id_col, merge_id)
pheno_file <- pheno_file[,c(1,2,4)] #keep only FID, IID, phenotype cols
fwrite(pheno_file, file="Z:/data/mesa_models/mesa_pheno/thrombotic/predixcan/pheno_rankplt5.txt", row.name=F, quote=F, sep="\t")
