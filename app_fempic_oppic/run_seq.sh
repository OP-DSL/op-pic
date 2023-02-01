
folderName="Log_seq_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo $folderName

mkdir $folderName

echo "SEQ TEST START"

thr=1

for mult in {1,2,4,8,16,32,50}; do
	echo $mult $thr "START" 
	bin/fempic_genseq OPP_ALLOC_MULT=${mult} > $folderName/log_thr${thr}_mult${mult}.log;
	echo "END" 

done

echo "TEST END"

echo $folderName