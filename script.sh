
NTH=$(nproc --all)
CPU_TYPE=$(uname -n)

while read -r PATIENT; do # loop over patients from samples_config.yaml
	case $PATIENT in 
		"patients"*)
		echo skip
		;;
		*)
		echo current patient $PATIENT
		echo -e "patients:\n    $PATIENT" > sample_tmp.yaml

		for i in {1..3} # loop over single patient
		do

		#single thread run
		snakemake --use-conda --config n_sim=$i cpu_type=$CPU_TYPE sample="sample_tmp.yaml"
		python ~/biopipe/script_benchmark.py
		rm -fr ~/biopipe/results/
		rm -fr ~/biopipe/benchmarks/

		#multi-thread run
		for n_th in $(seq 2 2 $NTH) # loop over threads
		do
		snakemake --use-conda --config n_sim=$i map_thrs=$n_th RT_thrs=$n_th BQSR1_thrs=$n_th BQSR2_thrs=$n_th BQSR4_thrs=$n_th HC_thrs=$n_th cpu_type=$CPU_TYPE sample="sample_tmp.yaml"
		python ~/biopipe/script_benchmark.py
		rm -fr ~/biopipe/results/
		rm -fr ~/biopipe/benchmarks/

		#How to let Snakefile run index_bwa rule
		#Removing rule's products

		#mkdir ~/data_ref/trash
		#mv ~/data_ref/human_g1k_v37.* ~/data_ref/trash/
		#mv ~/data_ref/trash/human_g1k_v37.fasta ~/data_ref/
		#rm -fr ~/data_ref/trash/

		done

		
		;;
	esac
done < samples_config.yaml

rm samples_tmp.yaml

