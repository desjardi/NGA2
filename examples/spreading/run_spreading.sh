#!/bin/bash
CA=($(seq 10 10 170))
OH=(1e-3 3e-3 1e-2 3e-2 1e-1 3e-1 1e+0 3e+0 1e+1 3e+1 1e+2)
#for j in ${!CA[@]}; do
#   for i in ${!OH[@]}; do
#      echo "Running a simulation with Oh=${OH[$i]} and CA=${CA[$j]}"
#      mpiexec -n 4 ./nga.dp.gnu.opt.mpi.exe -i input  --Oh=${OH[$i]} --CA=${CA[$j]} >> result_${CA[$j]}.txt
#   done
#done

# Just one test run with Oh=1e-3
for j in ${!CA[@]}; do
   echo "Running a simulation with Oh=1e-3 and CA=${CA[$j]}"
   mpiexec -n 4 ./nga.dp.gnu.opt.mpi.exe -i input --Oh=1e-3 --CA=${CA[$j]} >> result.txt
done
