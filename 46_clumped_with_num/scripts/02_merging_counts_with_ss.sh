# define varibales
WD_SCRIPT=/lustre1/project/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/46_clumped_with_num/scripts
WD=/lustre1/project/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/36_clump_ss
REGENIE=clumping_regenie
TRACTOR=clumping_tractor
OUT=/lustre1/project/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/46_clumped_with_num/output

# merge 
Rscript $WD_SCRIPT/03_merging.R $WD/$REGENIE/ham_clumped_ss.txt $OUT/$REGENIE/ham_temp.txt $WD/$REGENIE/ham_clumped_ss_counts.txt
Rscript $WD_SCRIPT/03_merging.R $WD/$REGENIE/pvl_clumped_ss.txt $OUT/$REGENIE/pvl_temp.txt $WD/$REGENIE/pvl_clumped_ss_counts.txt

Rscript $WD_SCRIPT/03_merging.R $WD/$TRACTOR/afr_ham_clumped_ss.txt $OUT/$TRACTOR/afr_ham_temp.txt $WD/$TRACTOR/afr_ham_clumped_ss_counts.txt
Rscript $WD_SCRIPT/03_merging.R $WD/$TRACTOR/afr_pvl_clumped_ss.txt $OUT/$TRACTOR/afr_pvl_temp.txt $WD/$TRACTOR/afr_pvl_clumped_ss_counts.txt

Rscript $WD_SCRIPT/03_merging.R $WD/$TRACTOR/amr_ham_clumped_ss.txt $OUT/$TRACTOR/amr_ham_temp.txt $WD/$TRACTOR/amr_ham_clumped_ss_counts.txt
Rscript $WD_SCRIPT/03_merging.R $WD/$TRACTOR/amr_pvl_clumped_ss.txt $OUT/$TRACTOR/amr_pvl_temp.txt $WD/$TRACTOR/amr_pvl_clumped_ss_counts.txt

Rscript $WD_SCRIPT/03_merging.R $WD/$TRACTOR/asn_ham_clumped_ss.txt $OUT/$TRACTOR/asn_ham_temp.txt $WD/$TRACTOR/asn_ham_clumped_ss_counts.txt
Rscript $WD_SCRIPT/03_merging.R $WD/$TRACTOR/asn_pvl_clumped_ss.txt $OUT/$TRACTOR/asn_pvl_temp.txt $WD/$TRACTOR/asn_pvl_clumped_ss_counts.txt

Rscript $WD_SCRIPT/03_merging.R $WD/$TRACTOR/eur_ham_clumped_ss.txt $OUT/$TRACTOR/eur_ham_temp.txt $WD/$TRACTOR/eur_ham_clumped_ss_counts.txt
Rscript $WD_SCRIPT/03_merging.R $WD/$TRACTOR/eur_pvl_clumped_ss.txt $OUT/$TRACTOR/eur_pvl_temp.txt $WD/$TRACTOR/eur_pvl_clumped_ss_counts.txt