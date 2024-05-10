#!/usr/bin/env python3

# Generates the plic_net_class.f90 file for a given Pytorch model

import torch
import numpy as np
import struct
import os

torch.set_default_dtype(torch.float64)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=np.inf)
model = torch.jit.load('./model.pt')

file = open("plic_net.f90", "w")
print("!> PLIC-Net File\n!> Provides the architecture for the neural network and the weights/biases\n!> Use gen_plic_net.py in NGA2/tools/scripts/plic_net to generate this file for a given Pytorch model",file=file)
print("module plic_net",file=file)
print("\tuse precision, only: WP\n\timplicit none\n",file=file)

count = 0
for param in model.parameters():
    count = count + 1
    name = ""
    if count%2 != 0:
        name = "lay" + str(int(count/2)+1) + "_weight"
        weight = "reshape(["
        for i in range(param.size()[0]-1):
            weight = weight + np.array2string(param[i].detach().numpy(), separator=', ') + "&\n\t,"
        weight = weight + np.array2string(param[i+1].detach().numpy(), separator=', ') + "],shape(" + name + "))"
        size1 = param.transpose(0,1).size()[0]
        size2 = param.transpose(0,1).size()[1]
        print("\treal(WP), dimension(" + str(size1) + "," + str(size2) + "), parameter :: " + name + "=" + weight,file=file)
    else:
        size1 = param.size()[0]
        name = "lay" + str(int((count-1)/2)+1) + "_bias"
        print("\treal(WP), dimension(" + str(size1) + "), parameter :: " + name + "=" + np.array2string(param.detach().numpy(), separator=', '),file=file)

print("\n\tcontains\n\n\tsubroutine get_normal(moments,normal)\n\t\timplicit none", file=file)
print("\t\treal(WP), dimension(:), intent(in) :: moments", file=file)
print("\t\treal(WP), dimension(:), intent(out) :: normal", file=file)
print("\t\treal(WP), dimension(:), allocatable :: tmparr", file=file)
print("\t\ttmparr=max(0.0_WP,matmul(moments,lay1_weight)+lay1_bias)",file=file)
for i in range(int((count-1)/2)-1):
    name1 = "lay" + str(i+2) + "_weight"
    name2 = "lay" + str(i+2) + "_bias"
    print("\t\ttmparr=max(0.0_WP,matmul(tmparr," + name1 + ")+" + name2 + ")",file=file)
name1 = "lay" + str(i+3) + "_weight"
name2 = "lay" + str(i+3) + "_bias"
print("\t\tnormal=matmul(tmparr," + name1 + ")+" + name2,file=file)
print("\t\tdeallocate(tmparr)", file=file)
print("\tend subroutine", file=file)

print("\n\tsubroutine reflect_moments(moments,center,direction)\n\t\timplicit none",file=file)
print("\t\treal(WP), dimension(0:188), intent(inout) :: moments\n\t\treal(WP), dimension(0:2), intent(in) :: center\n\t\tinteger, intent(out) :: direction", file=file)
print("\n\t\tdirection=0\n\t\tif (center(0).lt.0.and.center(1).ge.0.and.center(2).ge.0) then\n\t\t\tdirection = 1\n\t\t\tcall reflect_moments_x(moments)",file=file)
print("\t\telse if (center(0).ge.0.and.center(1).lt.0.and.center(2).ge.0) then\n\t\t\tdirection = 2\n\t\t\tcall reflect_moments_y(moments)",file=file)
print("\t\telse if (center(0).ge.0.and.center(1).ge.0.and.center(2).lt.0) then\n\t\t\tdirection = 3\n\t\t\tcall reflect_moments_z(moments)",file=file)
print("\t\telse if (center(0).lt.0.and.center(1).lt.0.and.center(2).ge.0) then\n\t\t\tdirection = 4\n\t\t\tcall reflect_moments_x(moments)\n\t\t\tcall reflect_moments_y(moments)",file=file)
print("\t\telse if (center(0).lt.0.and.center(1).ge.0.and.center(2).lt.0) then\n\t\t\tdirection = 5\n\t\t\tcall reflect_moments_x(moments)\n\t\t\tcall reflect_moments_z(moments)",file=file)
print("\t\telse if (center(0).ge.0.and.center(1).lt.0.and.center(2).lt.0) then\n\t\t\tdirection = 6\n\t\t\tcall reflect_moments_y(moments)\n\t\t\tcall reflect_moments_z(moments)",file=file)
print("\t\telse if (center(0).lt.0.and.center(1).lt.0.and.center(2).lt.0) then\n\t\t\tdirection = 7\n\t\t\tcall reflect_moments_x(moments)\n\t\t\tcall reflect_moments_y(moments)\n\t\t\tcall reflect_moments_z(moments)",file=file)
print("\t\tend if\n\tend subroutine", file=file)

print("\n\tsubroutine reflect_moments_x(moments)\n\t\timplicit none",file=file)
print("\t\treal(WP), dimension(0:188), intent(inout) :: moments", file=file)
print("\t\tinteger :: i,j,k,n\n\t\treal(WP) :: temp", file=file)
print("\t\tdo k=0,2\n\t\t\tdo j=0,2\n\t\t\t\tdo i=0,2\n\t\t\t\t\tif (i.eq.0) then\n\t\t\t\t\tdo n=0,6",file=file)
print("\t\t\t\t\t\tif (n.eq.1.or.n.eq.4) then\n\t\t\t\t\t\t\ttemp = moments(7*(i*9+j*3+k)+n)\n\t\t\t\t\t\t\tmoments(7*(i*9+j*3+k)+n) = -moments(7*(2*9+j*3+k)+n)\n\t\t\t\t\t\t\tmoments(7*(2*9+j*3+k)+n) = -temp",file=file)
print("\t\t\t\t\t\telse\n\t\t\t\t\t\t\ttemp = moments(7*(i*9+j*3+k)+n)\n\t\t\t\t\t\t\tmoments(7*(i*9+j*3+k)+n) = moments(7*(2*9+j*3+k)+n)\n\t\t\t\t\t\t\tmoments(7*(2*9+j*3+k)+n) = temp\n\t\t\t\t\t\tend if\n\t\t\t\t\tend do",file=file)
print("\t\t\t\t\telse if (i.eq.1) then\n\t\t\t\t\t\tmoments(7*(i*9+j*3+k)+1) = -moments(7*(i*9+j*3+k)+1)",file=file)
print("\t\t\t\t\t\tmoments(7*(i*9+j*3+k)+4) = -moments(7*(i*9+j*3+k)+4)\n\t\t\t\t\tend if\n\t\t\t\tend do\n\t\t\tend do\n\t\tend do\n\tend subroutine", file=file)

print("\n\tsubroutine reflect_moments_y(moments)\n\t\timplicit none",file=file)
print("\t\treal(WP), dimension(0:188), intent(inout) :: moments", file=file)
print("\t\tinteger :: i,j,k,n\n\t\treal(WP) :: temp", file=file)
print("\t\tdo k=0,2\n\t\t\tdo j=0,2\n\t\t\t\tdo i=0,2\n\t\t\t\t\tif (j.eq.0) then\n\t\t\t\t\tdo n=0,6",file=file)
print("\t\t\t\t\t\tif (n.eq.2.or.n.eq.5) then\n\t\t\t\t\t\t\ttemp = moments(7*(i*9+j*3+k)+n)\n\t\t\t\t\t\t\tmoments(7*(i*9+j*3+k)+n) = -moments(7*(i*9+2*3+k)+n)\n\t\t\t\t\t\t\tmoments(7*(i*9+2*3+k)+n) = -temp",file=file)
print("\t\t\t\t\t\telse\n\t\t\t\t\t\t\ttemp = moments(7*(i*9+j*3+k)+n)\n\t\t\t\t\t\t\tmoments(7*(i*9+j*3+k)+n) = moments(7*(i*9+2*3+k)+n)\n\t\t\t\t\t\t\tmoments(7*(i*9+2*3+k)+n) = temp\n\t\t\t\t\t\tend if\n\t\t\t\t\tend do",file=file)
print("\t\t\t\t\telse if (j.eq.1) then\n\t\t\t\t\t\tmoments(7*(i*9+j*3+k)+2) = -moments(7*(i*9+j*3+k)+2)",file=file)
print("\t\t\t\t\t\tmoments(7*(i*9+j*3+k)+5) = -moments(7*(i*9+j*3+k)+5)\n\t\t\t\t\tend if\n\t\t\t\tend do\n\t\t\tend do\n\t\tend do\n\tend subroutine", file=file)

print("\n\tsubroutine reflect_moments_z(moments)\n\t\timplicit none",file=file)
print("\t\treal(WP), dimension(0:188), intent(inout) :: moments", file=file)
print("\t\tinteger :: i,j,k,n\n\t\treal(WP) :: temp", file=file)
print("\t\tdo k=0,2\n\t\t\tdo j=0,2\n\t\t\t\tdo i=0,2\n\t\t\t\t\tif (k.eq.0) then\n\t\t\t\t\tdo n=0,6",file=file)
print("\t\t\t\t\t\tif (n.eq.3.or.n.eq.6) then\n\t\t\t\t\t\t\ttemp = moments(7*(i*9+j*3+k)+n)\n\t\t\t\t\t\t\tmoments(7*(i*9+j*3+k)+n) = -moments(7*(i*9+j*3+2)+n)\n\t\t\t\t\t\t\tmoments(7*(i*9+j*3+2)+n) = -temp",file=file)
print("\t\t\t\t\t\telse\n\t\t\t\t\t\t\ttemp = moments(7*(i*9+j*3+k)+n)\n\t\t\t\t\t\t\tmoments(7*(i*9+j*3+k)+n) = moments(7*(i*9+j*3+2)+n)\n\t\t\t\t\t\t\tmoments(7*(i*9+j*3+2)+n) = temp\n\t\t\t\t\t\tend if\n\t\t\t\t\tend do",file=file)
print("\t\t\t\t\telse if (k.eq.1) then\n\t\t\t\t\t\tmoments(7*(i*9+j*3+k)+3) = -moments(7*(i*9+j*3+k)+3)",file=file)
print("\t\t\t\t\t\tmoments(7*(i*9+j*3+k)+6) = -moments(7*(i*9+j*3+k)+6)\n\t\t\t\t\tend if\n\t\t\t\tend do\n\t\t\tend do\n\t\tend do\n\tend subroutine", file=file)

print("end module plic_net", file=file)