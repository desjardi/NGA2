#!/bin/bash
one_over_two_pi_squared=0.025330295910584444
g=1
rho_l=2
rho_g=1
Lx=1
klc_vals=($(seq 0 0.1 1.2))
for i in ${!klc_vals[@]}; do
   echo "Running a simulation with k*lc=${klc_vals[$i]}"
   sigma=$(echo "$g * ($rho_l - $rho_g) * ${klc_vals[$i]} * ${klc_vals[$i]} * $Lx * $Lx * $one_over_two_pi_squared" | bc -l)
   echo "sigma is =$sigma"
   mpiexec -n 4 ./nga.dp.gnu.opt.mpi.exe -i input --"Surface tension coefficient"=$sigma >> result.txt
done
