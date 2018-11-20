source activate snp_calling
source /master/nplatt/schisto_hybridization/scripts/set_env.sh

cd $RESULTS_DIR

mkdir sprime_vs_pcadmix

cd sprime_vs_pcadmix

#get pcadmix and sprime data
cp ../sprime/2018-06-12/niger_tz-out_sprime_2018-06-12.score .
ln -s niger_tz-out_sprime_2018-06-12.score sprime.score

cp ../pcadmix/control-bov-admixed_maf00.ia.csv .
ln -s control-bov-admixed_maf00.ia.csv pcadmix_chr.score

cp ../pcadmix/control-bov-admixed_maf00.vit.csv .
ln -s control-bov-admixed_maf00.vit.csv pcadmix_window.score

cp ../pcadmix/control-bov-admixed_maf00.markers.csv .
ln -s pcadmix_markers.list

#pcadmix sizes per window
cat control-bov-admixed_maf00.markers.csv | sed 's/,/ /g' | awk '{print $1, $2, $NF}' >pcadmix_window.start_stop 

python find_pcadmix_window_coords.py >window_coords.txt

#used excel to calculate percentages and MBs in file
pcadmix_mb.txt




