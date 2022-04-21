import struct
import numpy as np

# This reads nga's config.grid and config.geom files
config=open('config.grid',mode='rb').read()
(name,coord,xper,yper,zper,nx,ny,nz)=struct.unpack_from('64siiiiiii',config,offset=0)
x=struct.unpack_from('d'*(nx+1),config,offset=92)
y=struct.unpack_from('d'*(ny+1),config,offset=92+8*(nx+1))
z=struct.unpack_from('d'*(nz+1),config,offset=92+8*(nx+1)+8*(ny+1))
geom=open('config.geom',mode='rb').read()
(n1,n2,n3,nval,nvar)=struct.unpack_from('iiiii',geom,offset=0)
varname=struct.unpack_from('8s',geom,offset=20)
var=struct.unpack_from('d'*n1*n2*n3,geom,offset=28)
wall=np.array(var)
wall=np.reshape(wall,(nx,ny,nz))

# This reads nga's data file
data=open('config.geom',mode='rb').read()
(n1,n2,n3,nval,nvar)=struct.unpack_from('iiiii',data,offset=0)
valname=struct.unpack_from('8s'*nval,data,offset=20)
val=struct.unpack_from('d'*nval,data,offset=20+nval*8)
varname=struct.unpack_from('8s'*nvar,data,offset=20+nval*(8+8))
var=struct.unpack_from('d'*nvar*n1*n2*n3,data,offset=20+nval*(8+8)+nvar*8)
