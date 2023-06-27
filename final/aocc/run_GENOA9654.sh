#!/usr/bin/bash

#SBATCH --nodes=1
#SBATCH --partition=9654_RHEL8.3_Ashish
#SBATCH --exclusive
#SBATCH --job-name=ldc-aocc
##SBATCH --time=1:00:00
#SBATCH --error=aocc9654-%j.err
#SBATCH --output=aocc9654-%j.out
#SBATCH --ntasks-per-node=192
#SBATCH --cpus-per-task=1

source /mnt/share/tools/aocc/AOCC_4.0/aocc-compiler-rel-4.0-3206-145/setenv_AOCC.sh

size=40000
iter=10

echo "====================================================================================" | tee -a run_GENOA9654.log
printf "\n \n" | tee -a run_GENOA9654.log

ncpus=$(lscpu | grep "^CPU(s):" | awk '{print $2}')

export OMP_NUM_THREADS=$ncpus
export OMP_PROC_BIND=close

echo " AOCC AVX512: " | tee -a run_GENOA9654.log
echo "        MAX THREADS=$ncpus" | tee -a run_GENOA9654.log
echo "        OMP_PROC_BIND=close" | tee -a run_GENOA9654.log
printf "\n"
echo "        ./bin/ldc_aoccZNVER4AVX512 $size $iter" | tee -a run_GENOA9654.log

echo "  Time        |     Elem      |      Thrds   " | tee -a run_GENOA9654.log
echo "---------------------------------------------" | tee -a run_GENOA9654.log
padl=15

OUTPUT=`./bin/ldc_aoccZNVER4AVX512 $size $iter`
out=($(echo $OUTPUT | grep -oE "[0-9]+\.[0-9]+"))
printf "% *s" $(((padl-${#out[0]})/2)) ${out[-1]} | tee -a run_GENOA9654.log
printf "% *s" $((padl+(padl-${#out[-1]})/2)) $size | tee -a run_GENOA9654.log
printf "% *s" $(((padl-${#ncpus})/2)) $ncpus | tee -a run_GENOA9654.log
printf "\n" | tee -a run_GENOA9654.log
printf "\n" | tee -a run_GENOA9654.log

echo " AOCC AVX2: " | tee -a run_GENOA9654.log
echo "        MAX THREADS=$ncpus" | tee -a run_GENOA9654.log
echo "        OMP_PROC_BIND=close" | tee -a run_GENOA9654.log
printf "\n"
echo "        ./bin/ldc_aoccZNVER4AVX2 $size $iter" | tee -a run_GENOA9654.log

echo "  Time        |     Elem      |      Thrds   " | tee -a run_GENOA9654.log
echo "---------------------------------------------" | tee -a run_GENOA9654.log
padl=15

OUTPUT=`./bin/ldc_aoccZNVER4AVX2 $size $iter`
out=($(echo $OUTPUT | grep -oE "[0-9]+\.[0-9]+"))
printf "% *s" $(((padl-${#out[0]})/2)) ${out[-1]} | tee -a run_GENOA9654.log
printf "% *s" $((padl+(padl-${#out[-1]})/2)) $size | tee -a run_GENOA9654.log
printf "% *s" $(((padl-${#ncpus})/2)) $ncpus | tee -a run_GENOA9654.log
printf "\n" | tee -a run_GENOA9654.log

