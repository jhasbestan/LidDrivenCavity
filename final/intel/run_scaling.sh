#!/usr/bin/bash

#SBATCH --nodes=1
#SBATCH --partition=9654_RHEL8.3_Ashish
#SBATCH --exclusive
#SBATCH --job-name=ldc-intel
##SBATCH --time=1:00:00
#SBATCH --error=scaling9654-%j.err
#SBATCH --output=scaling9654-%j.out
#SBATCH --ntasks-per-node=192
#SBATCH --cpus-per-task=1

module use /mnt/share/intel/oneapi/2022.2/modulefiles
module load icc/latest

size=40000
iter=10

echo "====================================================================================" | tee -a run_scaling.log
printf "\n \n" | tee -a run_scaling.log

ncpus=$(lscpu | grep "^CPU(s):" | awk '{print $2}')
export KMP_AFFINITY=granularity=fine,compact

echo " Intel AVX512: " | tee -a run_scaling.log
echo "        MAX THREADS=$ncpus" | tee -a run_scaling.log
echo "        KMP_AFFINITY=granularity=fine,compact" | tee -a run_scaling.log
printf "\n"
echo "        ./bin/incomp_avx512 $size $iter" | tee -a run_scaling.log

echo "  Time        |     Elem        |     Thrds  "
echo "---------------------------------------------"
padl=15

for THRDS in 192 96 48 24 12 
do
        export OMP_NUM_THREADS=$THRDS
	OUTPUT=`./bin/incomp_avx512 $size $iter`
	out=($(echo $OUTPUT | grep -oE "[0-9]+\.[0-9]+"))
#	printf "% *f" $(((padl-${#out[0]})/2)) ${out[0]};
        printf "% *s" $(((padl-${#out[0]})/2)) ${out[-1]};
	printf "% *s" $((padl+(padl-${#out[-1]})/2)) $size;
	printf "% *s" $(((padl-${#THRDS})/2)) $THRDS;
	printf "\n"
done

