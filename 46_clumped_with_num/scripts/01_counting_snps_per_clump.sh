# define varibales
WD=/lustre1/project/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/36_clump_ss
REGENIE=clumping_regenie
TRACTOR=clumping_tractor
OUT=/lustre1/project/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/46_clumped_with_num/output

### create the tmp files to save the counts 
#REGNIE
awk 'NR>1 {print $3, $6+1}' $WD/$REGENIE/ham_ss_clumped.clumped | head -n -2 > $OUT/$REGENIE/ham_temp.txt
awk 'NR>1 {print $3, $6+1}' $WD/$REGENIE/pvl_ss_clumped.clumped | head -n -2 > $OUT/$REGENIE/pvl_temp.txt

#TRACTOR
awk 'NR>1 {print $3, $6+1}' $WD/$TRACTOR/afr_ham_clump.clumped | head -n -2 > $OUT/$TRACTOR/afr_ham_temp.txt
awk 'NR>1 {print $3, $6+1}' $WD/$TRACTOR/afr_pvl_clump.clumped | head -n -2 > $OUT/$TRACTOR/afr_pvl_temp.txt

awk 'NR>1 {print $3, $6+1}' $WD/$TRACTOR/amr_ham_clump.clumped | head -n -2 > $OUT/$TRACTOR/amr_ham_temp.txt
awk 'NR>1 {print $3, $6+1}' $WD/$TRACTOR/amr_pvl_clump.clumped | head -n -2 > $OUT/$TRACTOR/amr_pvl_temp.txt

awk 'NR>1 {print $3, $6+1}' $WD/$TRACTOR/asn_ham_clump.clumped | head -n -2 > $OUT/$TRACTOR/asn_ham_temp.txt
awk 'NR>1 {print $3, $6+1}' $WD/$TRACTOR/asn_pvl_clump.clumped | head -n -2 > $OUT/$TRACTOR/asn_pvl_temp.txt

awk 'NR>1 {print $3, $6+1}' $WD/$TRACTOR/eur_ham_clump.clumped | head -n -2 > $OUT/$TRACTOR/eur_ham_temp.txt
awk 'NR>1 {print $3, $6+1}' $WD/$TRACTOR/eur_pvl_clump.clumped | head -n -2 > $OUT/$TRACTOR/eur_pvl_temp.txt