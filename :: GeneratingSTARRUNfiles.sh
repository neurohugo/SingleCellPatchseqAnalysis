// Generating STAR RUN files

// check
head -1 listRNB5.txt |awk '{split($0,a,"/"); print a[12]}'|awk '{split($0,a,"."); print a[1]}'
// run
x=$(wc -l listRNB5.txt|awk '{print $1}')

for i in $(seq 1 2 $x)
do
Fullnames=$(head -$i listRNB5.txt |tail -1 |awk '{split($0,a,"/"); print a[12]}'|awk '{split($0,a,"."); print a[1]}')
ReadOne=$(head -$i listRNB5.txt |tail -1)
ReadTwo=$(head -$(($i+1)) listRNB5.txt |tail -1)
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
echo 'STAR --genomeDir /projects/ps-zhenglab/Hugo/SequencingData/REF/MouseGenome/INDEX \'>>$Fullnames.pbs
echo '--runThreadN 6 \'>>$Fullnames.pbs
echo '--readFilesIn' $ReadOne $ReadTwo '\'>>$Fullnames.pbs
echo '--outFileNamePrefix /projects/ps-zhenglab/Hugo/SequencingData/HugosData/Analysis/RNB5/'$Fullnames '\'>>$Fullnames.pbs
echo '--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 \'>>$Fullnames.pbs
echo '--outSAMtype BAM SortedByCoordinate \'>>$Fullnames.pbs
echo '--outSAMunmapped Within \'>>$Fullnames.pbs
echo '--outSAMattributes Standard \'>>$Fullnames.pbs
echo '--readFilesCommand zcat'>>$Fullnames.pbs
done
