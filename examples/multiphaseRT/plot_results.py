import numpy as np
from scipy.optimize import fsolve
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
Pl=rho_l/(rho_l+rho_g)
Pg=rho_g/(rho_l+rho_g)
At=Pl-Pg
iGal=(2*np.pi)**3/(np.abs(gy)*Lx**3)*((     mu_l)/(     +rho_l))**2
iGag=(2*np.pi)**3/(np.abs(gy)*Lx**3)*((mu_g     )/(rho_g      ))**2
iGam=(2*np.pi)**3/(np.abs(gy)*Lx**3)*((mu_g-mu_l)/(rho_g+rho_l))**2
iGap=(2*np.pi)**3/(np.abs(gy)*Lx**3)*((mu_g+mu_l)/(rho_g+rho_l))**2

# Prepare plot
plt.title('Multiphase Rayleigh-Taylor Instability')
plt.xlabel('Normalized wavenumber')
plt.ylabel('Normalized growth rate')

# Plot reference inviscid and viscous data
ref_klc=np.linspace(0,1,num=1000)
inviscid_ngr=np.sqrt(At*(ref_klc-ref_klc**3))
viscous_ngr =np.sqrt(At*(ref_klc-ref_klc**3)+iGap*ref_klc)-np.sqrt(iGap*ref_klc)
plt.plot(ref_klc,inviscid_ngr,'-',lw=2,color='blue',label='Inviscid exact')
plt.plot(ref_klc, viscous_ngr,'-',lw=2,color='green',label='Viscous approx.')
plt.ylim(bottom=0,top=1.2*max(inviscid_ngr))

# Define dispersion relation for two-phase viscous RT
def viscous_dispersion(x,klc,At,Pl,Pg,iGal,iGag,iGam):
    if klc <= 0.0:
        return (0.0)
    else:
        return ((At*klc*(1-klc**2)-x**2)*(Pl*np.sqrt(1+x/np.sqrt(iGag*klc))+Pg*np.sqrt(1+x/np.sqrt(iGal*klc))-1)
        -4*Pl*Pg*x**2+4*np.sqrt(iGam*klc)*x*(Pl*np.sqrt(1+x/np.sqrt(iGag*klc))-Pg*np.sqrt(1+x/np.sqrt(iGal*klc))-At)
        +4*iGam*klc*(np.sqrt(1+x/np.sqrt(iGag*klc))-1)*(np.sqrt(1+x/np.sqrt(iGal*klc))-1))

# Plot exact viscous RT solution
viscous2_ngr=0*ref_klc
for ik, myklc in np.ndenumerate(ref_klc):
    viscous2_ngr[ik]=fsolve(func=viscous_dispersion,x0=viscous_ngr[ik],args=(myklc,At,Pl,Pg,iGal,iGag,iGam))
plt.plot(ref_klc,viscous2_ngr,'-',lw=2,color='purple',label='Viscous exact')

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
