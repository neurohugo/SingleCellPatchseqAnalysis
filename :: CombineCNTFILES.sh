// Combine CNT FILES

ls /projects/ps-zhenglab/Hugo/SequencingData/HugosData/Analysis/RNB4/*cnt >listCounts.txt
x=$(wc -l listCounts.txt|awk '{print $1}')
Fullnames=$(echo 'Combine_'$x)
echo '#!/bin/csh' >$Fullnames.pbs
echo '#PBS -q hotel'>>$Fullnames.pbs
echo '#PBS -N' $Fullnames>>$Fullnames.pbs
echo '#PBS -l nodes=6:ppn=1'>>$Fullnames.pbs
echo '#PBS -l walltime=2:50:00'>>$Fullnames.pbs
echo '#PBS -o' $Fullnames'.out'>>$Fullnames.pbs
echo '#PBS -e' $Fullnames'.err'>>$Fullnames.pbs
echo '#PBS -V'>>$Fullnames.pbs
echo '#PBS -M jmk002@ucsd.edu'>>$Fullnames.pbs
echo '#PBS -m abe'>>$Fullnames.pbs
echo '#PBS -A jmk002'>>$Fullnames.pbs
echo 'cd /oasis/tscc/scratch/jmk002'>>$Fullnames.pbs
echo 'mpirun -v -machinefile $PBS_NODEFILE -np 20./jmk002mpi.out'>>$Fullnames.pbs
echo "ls -t -v /projects/ps-zhenglab/Hugo/SequencingData/HugosData/Analysis/RNB4/*cnt | tr '\n' ' '>/home/jmk002/PBSfiles/htseq/RNB4/listCounts.txt;">>$Fullnames.pbs
echo "awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}'" $(cat /home/jmk002/PBSfiles/htseq/RNB4/listCounts.txt) ">/home/jmk002/PBSfiles/htseq/RNB4/Combined.txt">>$Fullnames.pbs

ls -t -v /projects/ps-zhenglab/Hugo/SequencingData/HugosData/Analysis/RNB4/*cnt | tr '\n' ' '>/home/jmk002/PBSfiles/htseq/RNB4/listCounts.txt
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}'  $(cat /home/jmk002/PBSfiles/htseq/RNB4/listCounts.txt) >/home/jmk002/PBSfiles/htseq/RNB4/Combined.txt


x=$(wc -l listCounts.txt|awk '{print $1}')
for i in $(seq 1 1 $x)
do
Thisname=$(head -$i listCounts.txt |tail -1 )
Samplename=$(head -$i listCounts.txt|tail -1 |awk '{split($0,a,"/"); print a[9]}'|awk '{split($0,a,"."); print a[1]}' )
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $Thisname $Thisname | awk '{print $2}' > $Samplename.tmps
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $Thisname $Thisname | awk '{print $1}' > 00names.tmps
done
ls -t -v *tmps | tr '\n' ' '>Listtmps.txt
paste $(cat Listtmps.txt) >CombinedRNBCombined.txt


