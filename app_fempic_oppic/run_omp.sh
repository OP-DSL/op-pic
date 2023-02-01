
folderName="Log_omp_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo $folderName

mkdir $folderName

echo "OMP TEST START"

export OMP_NUM_THREADS=2
export OMP_PROC_BIND=close

for thr in {1,2,4,8,16,24,32,40,48}; do

	export OMP_NUM_THREADS=${thr}

	for mult in {1,2,4,8,16,32,50}; do
		echo $mult $thr "START" 
		bin/fempic_openmp OPP_ALLOC_MULT=${mult} > $folderName/log_thr${thr}_mult${mult}.log;
		echo "END" 

	done
done

echo "TEST END"

echo $folderName