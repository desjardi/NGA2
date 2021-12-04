import numpy as np
from matplotlib import pyplot as plt
import re

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

# Create non-dimensional numbers
At=(rho_l-rho_g)/(rho_l+rho_g)
iGa=(2*np.pi)**3/(np.abs(gy)*Lx**3)*((mu_l+mu_g)/(rho_l+rho_g))**2

# Prepare plot
plt.title('Multiphase Rayleigh-Taylor Instability')
plt.xlabel('Normalized wavenumber')
plt.ylabel('Normalized growth rate')

# Plot reference inviscid and viscous data
ref_klc=np.linspace(0,1,num=1000)
inviscid_ngr=np.sqrt(At*(ref_klc-ref_klc**3))
viscous_ngr =np.sqrt(At*(ref_klc-ref_klc**3)+iGa*ref_klc)-np.sqrt(iGa*ref_klc)
plt.plot(ref_klc,inviscid_ngr,'-',lw=2,color='blue',label='Inviscid')
plt.plot(ref_klc,viscous_ngr,'-',lw=2,color='green',label='Viscous')
plt.ylim(bottom=0,top=1.2*max(inviscid_ngr))

# Plot results from result.txt file
result=np.fromfile('result.txt',dtype=float,sep=' ')
result=np.reshape(result,(-1,4))
lc =result[:,0]
tau=result[:,1]
klc=result[:,2]
ngr=result[:,3]
plt.plot(klc,np.clip(ngr,0,None),'+',ms=10,mew=2,color='red',label='nga2')

# Plot
plt.legend()
plt.show()
