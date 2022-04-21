import numpy as np
import os
import re

# Input range of cases to run
minklc=0.0
maxklc=1.0
nklc=11
klc=np.linspace(minklc,maxklc,num=nklc)

# Process the input file to obtain simulation parameters
with open('input') as fobj:
    for line in fobj:
        line_data = re.split(':',line)
        if line_data[0].rstrip()=='Gravity':
            gx=float(line_data[1].split()[0])
            gy=float(line_data[1].split()[1])
            gz=float(line_data[1].split()[2])
        if line_data[0].rstrip()=='Liquid density':
            rho_l=float(line_data[1])
        if line_data[0].rstrip()=='Gas density':
            rho_g=float(line_data[1])
        if line_data[0].rstrip()=='Liquid dynamic viscosity':
            mu_l=float(line_data[1])
        if line_data[0].rstrip()=='Gas dynamic viscosity':
            mu_g=float(line_data[1])
        if line_data[0].rstrip()=='Lx':
            Lx=float(line_data[1])

# Data will be stored in 'result.txt' file
os.system('cp result.txt oldresult.txt')
os.system('rm result.txt')

# Loop over cases, generate sigma, and run
for myklc in klc:
    sigma=abs(gy)*(rho_l-rho_g)*myklc**2*Lx**2/(2*np.pi)**2
    print('Running simulation for klc=',myklc,' and sigma=',sigma)
    command=f'mpiexec -n 4 ./nga.dp.gnu.opt.mpi.exe -i input --"Surface tension coefficient={sigma}" >> result.txt'
    os.system(command)
