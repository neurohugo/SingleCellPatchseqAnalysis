//Generating htseqCountsfrom merged

ls /projects/ps-zhenglab/Hugo/SequencingData/HugosData/Analysis/RNB4/*.bam >listMergeBam.txt
x=$(wc -l listMergeBam.txt|awk '{print $1}')
for i in $(seq 1 1 $x)
do
Fullnames=$(head -$i listMergeBam.txt |tail -1 |awk '{split($0,a,"/"); print a[9]}'|awk '{split($0,a,"."); print a[1]}')
ThisRead=$(head -$i listMergeBam.txt |tail -1)
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
echo 'module load pysam'>>$Fullnames.pbs
echo 'htseq-count -f bam' $ThisRead '/projects/ps-zhenglab/Hugo/SequencingData/REF/MouseGenome/GTF/gencode.vM20.annotation.gtf>>/projects/ps-zhenglab/Hugo/SequencingData/HugosData/Analysis/RNB5/'$Fullnames.cnt >>$Fullnames.pbs
done