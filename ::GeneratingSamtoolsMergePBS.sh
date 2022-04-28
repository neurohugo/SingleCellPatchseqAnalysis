//Generating Samtools Merge PBS

x=$(wc -l listRNB3Bam.txt|awk '{print $1}')

for i in $(seq 1 2 $x)
do

Fullnames=$(head -$i listRNB3Bam.txt |tail -1 |awk '{split($0,a,"/"); print a[10]}'|awk '{split($0,a,"."); print a[1]}' |awk '{split($0,a,"_"); print a[1] "_" a[2]}')
LaneOne=$(head -$i listRNB3Bam.txt |tail -1)
LaneTwo=$(head -$(($i+1)) listRNB3Bam.txt |tail -1)
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
echo 'samtools merge /projects/ps-zhenglab/Hugo/SequencingData/HugosData/Analysis/RNB2/'$Fullnames'_Merge.bam' $LaneOne $LaneTwo >>$Fullnames.pbs
echo 'samtools sort /projects/ps-zhenglab/Hugo/SequencingData/HugosData/Analysis/RNB2/'$Fullnames'_Merge_sorted.bam' '/projects/ps-zhenglab/Hugo/SequencingData/HugosData/Analysis/RNB2/'$Fullnames'_Merge.bam' >>$Fullnames.pbs
done

x=$(wc -l listRNB4Bam.txt|awk '{print $1}')

for i in $(seq 1 1 $x)
do
Fullnames=$(head -$i listRNB4Bam.txt |tail -1 |awk '{split($0,a,"/"); print a[9]}'|awk '{split($0,a,"."); print a[1]}' |awk '{split($0,a,"_"); print a[1] "_" a[2]}')
Unsorted=$(head -$i listRNB4Bam.txt |tail -1)
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
echo 'samtools sort -o /projects/ps-zhenglab/Hugo/SequencingData/HugosData/Analysis/RNB4/'$Fullnames'_sorted.bam' $Unsorted >>$Fullnames.pbs
done

