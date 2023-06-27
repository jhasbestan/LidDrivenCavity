#!/usr/bin/bash

#SBATCH --nodes=1
#SBATCH --partition=9654_RHEL8.3_Ashish
#SBATCH --exclusive
#SBATCH --job-name=ldc-intel
##SBATCH --time=1:00:00
#SBATCH --error=intel9654-%j.err
#SBATCH --output=intel9654-%j.out
#SBATCH --ntasks-per-node=192
#SBATCH --cpus-per-task=1

module use /mnt/share/intel/oneapi/2022.2/modulefiles
module load icc/latest

size=40000
iter=10

echo "====================================================================================" | tee -a run_GENOA9654.log
printf "\n \n" | tee -a run_GENOA9654.log

ncpus=$(lscpu | grep "^CPU(s):" | awk '{print $2}')

export OMP_NUM_THREADS=$ncpus
export KMP_AFFINITY=granularity=fine,compact

echo " INTEL AVX512: " | tee -a run_GENOA9654.log
echo "        MAX THREADS=$ncpus" | tee -a run_GENOA9654.log
echo "        KMP_AFFINITY=granularity=fine,compact" | tee -a run_GENOA9654.log
printf "\n"
echo "        ./bin/incomp_avx512 $size $iter" | tee -a run_GENOA9654.log

echo "  Time        |     Elem      |      Thrds   " | tee -a run_GENOA9654.log
echo "---------------------------------------------" | tee -a run_GENOA9654.log
padl=15

OUTPUT=`./bin/incomp_avx512 $size $iter`
out=($(echo $OUTPUT | grep -oE "[0-9]+\.[0-9]+"))
printf "% *s" $(((padl-${#out[0]})/2)) ${out[-1]} | tee -a run_GENOA9654.log
printf "% *s" $((padl+(padl-${#out[1]})/2)) $size | tee -a run_GENOA9654.log
printf "% *s" $(((padl-${#ncpus})/2)) $ncpus | tee -a run_GENOA9654.log
printf "\n" | tee -a run_GENOA9654.log
printf "\n" | tee -a run_GENOA9654.log

echo " INTEL AVX2: " | tee -a run_GENOA9654.log
echo "        MAX THREADS=$ncpus" | tee -a run_GENOA9654.log
echo "        KMP_AFFINITY=granularity=fine,compact" | tee -a run_GENOA9654.log
printf "\n"
echo "        ./bin/incomp_avx2 $size $iter" | tee -a run_GENOA9654.log

echo "  Time        |     Elem        |     Thrds  " | tee -a run_GENOA9654.log
echo "---------------------------------------------" | tee -a run_GENOA9654.log
padl=15

OUTPUT=`./bin/incomp_avx2 $size $iter`
out=($(echo $OUTPUT | grep -oE "[0-9]+\.[0-9]+"))
printf "% *s" $(((padl-${#out[0]})/2)) ${out[-1]} | tee -a run_GENOA9654.log
printf "% *s" $((padl+(padl-${#out[1]})/2)) $size | tee -a run_GENOA9654.log
printf "% *s" $(((padl-${#ncpus})/2)) $ncpus | tee -a run_GENOA9654.log
printf "\n" | tee -a run_GENOA9654.log

