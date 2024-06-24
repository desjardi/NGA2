#!/usr/bin/env python3

# Generates the plicnet.f90 file for a given Pytorch model

import torch
import numpy as np
import struct
import os

torch.set_default_dtype(torch.float64)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=np.inf)
model = torch.jit.load('./model.pt')

file = open("plicnet.f90", "w")
print("!> PLIC-Net File\n!> Provides the architecture for the neural network and the weights/biases\n!> Use generate_plicnet.py in NGA2/tools/scripts/plicnet to generate this file for a given Pytorch model",file=file)
print("module plicnet",file=file)
print("   use precision, only: WP\n   implicit none",file=file)

count = 0
hidden = 0
for param in model.parameters():
    count = count + 1
    name = ""
    if count == 3:
        hidden = param.transpose(0,1).size()[0] #Assumes each hidden layer has same number of neurons
    if count%2 != 0:
        name = "lay" + str(int(count/2)+1) + "_weight"
        weight = "reshape(["
        for i in range(param.size()[0]-1):
            weight = weight + np.array2string(param[i].detach().numpy(), separator=', ') + "&\n   ,"
        weight = weight + np.array2string(param[i+1].detach().numpy(), separator=', ') + "],shape(" + name + "))"
        size1 = param.transpose(0,1).size()[0]
        size2 = param.transpose(0,1).size()[1]
        print("   real(WP), dimension(" + str(size1) + "," + str(size2) + "), parameter :: " + name + "=" + weight,file=file)
    else:
        size1 = param.size()[0]
        name = "lay" + str(int((count-1)/2)+1) + "_bias"
        print("   real(WP), dimension(" + str(size1) + "), parameter :: " + name + "=" + np.array2string(param.detach().numpy(), separator=', '),file=file)

print("   contains\n   subroutine get_normal(moments,normal)\n      implicit none", file=file)
print("      real(WP), dimension(:), intent(in) :: moments !< Needs to be of size 189", file=file)
print("      real(WP), dimension(:), intent(out) :: normal !< Needs to be of size 3", file=file)
print("      real(WP), dimension(" + str(hidden) + ") :: tmparr", file=file)
print("      tmparr=max(0.0_WP,matmul(moments,lay1_weight)+lay1_bias)",file=file)
for i in range(int((count-1)/2)-1):
    name1 = "lay" + str(i+2) + "_weight"
    name2 = "lay" + str(i+2) + "_bias"
    print("      tmparr=max(0.0_WP,matmul(tmparr," + name1 + ")+" + name2 + ")",file=file)
name1 = "lay" + str(i+3) + "_weight"
name2 = "lay" + str(i+3) + "_bias"
print("      normal=matmul(tmparr," + name1 + ")+" + name2,file=file)
print("   end subroutine", file=file)

print("   subroutine reflect_moments(moments,center,direction)\n      implicit none",file=file)
print("      real(WP), dimension(0:), intent(inout) :: moments !< Needs to be of size (0:188)\n      real(WP), dimension(0:), intent(in) :: center     !< Needs to be of size (0:2)\n      integer, intent(out) :: direction", file=file)
print("      direction=0\n      if (center(0).lt.0.and.center(1).ge.0.and.center(2).ge.0) then\n         direction = 1\n         call reflect_moments_x(moments)",file=file)
print("      else if (center(0).ge.0.and.center(1).lt.0.and.center(2).ge.0) then\n         direction = 2\n         call reflect_moments_y(moments)",file=file)
print("      else if (center(0).ge.0.and.center(1).ge.0.and.center(2).lt.0) then\n         direction = 3\n         call reflect_moments_z(moments)",file=file)
print("      else if (center(0).lt.0.and.center(1).lt.0.and.center(2).ge.0) then\n         direction = 4\n         call reflect_moments_x(moments)\n         call reflect_moments_y(moments)",file=file)
print("      else if (center(0).lt.0.and.center(1).ge.0.and.center(2).lt.0) then\n         direction = 5\n         call reflect_moments_x(moments)\n         call reflect_moments_z(moments)",file=file)
print("      else if (center(0).ge.0.and.center(1).lt.0.and.center(2).lt.0) then\n         direction = 6\n         call reflect_moments_y(moments)\n         call reflect_moments_z(moments)",file=file)
print("      else if (center(0).lt.0.and.center(1).lt.0.and.center(2).lt.0) then\n         direction = 7\n         call reflect_moments_x(moments)\n         call reflect_moments_y(moments)\n         call reflect_moments_z(moments)",file=file)
print("      end if\n   end subroutine", file=file)

print("   subroutine reflect_moments_x(moments)\n      implicit none",file=file)
print("      real(WP), dimension(0:), intent(inout) :: moments !< Needs to be of size (0:188)", file=file)
print("      integer :: i,j,k,n\n      real(WP) :: temp", file=file)
print("      do k=0,2\n         do j=0,2\n            do i=0,2\n               if (i.eq.0) then\n               do n=0,6",file=file)
print("                  if (n.eq.1.or.n.eq.4) then\n                     temp = moments(7*(i*9+j*3+k)+n)\n                     moments(7*(i*9+j*3+k)+n) = -moments(7*(2*9+j*3+k)+n)\n                     moments(7*(2*9+j*3+k)+n) = -temp",file=file)
print("                  else\n                     temp = moments(7*(i*9+j*3+k)+n)\n                     moments(7*(i*9+j*3+k)+n) = moments(7*(2*9+j*3+k)+n)\n                     moments(7*(2*9+j*3+k)+n) = temp\n                  end if\n               end do",file=file)
print("               else if (i.eq.1) then\n                  moments(7*(i*9+j*3+k)+1) = -moments(7*(i*9+j*3+k)+1)",file=file)
print("                  moments(7*(i*9+j*3+k)+4) = -moments(7*(i*9+j*3+k)+4)\n               end if\n            end do\n         end do\n      end do\n   end subroutine", file=file)

print("   subroutine reflect_moments_y(moments)\n      implicit none",file=file)
print("      real(WP), dimension(0:), intent(inout) :: moments !< Needs to be of size (0:188)", file=file)
print("      integer :: i,j,k,n\n      real(WP) :: temp", file=file)
print("      do k=0,2\n         do j=0,2\n            do i=0,2\n               if (j.eq.0) then\n               do n=0,6",file=file)
print("                  if (n.eq.2.or.n.eq.5) then\n                     temp = moments(7*(i*9+j*3+k)+n)\n                     moments(7*(i*9+j*3+k)+n) = -moments(7*(i*9+2*3+k)+n)\n                     moments(7*(i*9+2*3+k)+n) = -temp",file=file)
print("                  else\n                     temp = moments(7*(i*9+j*3+k)+n)\n                     moments(7*(i*9+j*3+k)+n) = moments(7*(i*9+2*3+k)+n)\n                     moments(7*(i*9+2*3+k)+n) = temp\n                  end if\n               end do",file=file)
print("               else if (j.eq.1) then\n                  moments(7*(i*9+j*3+k)+2) = -moments(7*(i*9+j*3+k)+2)",file=file)
print("                  moments(7*(i*9+j*3+k)+5) = -moments(7*(i*9+j*3+k)+5)\n               end if\n            end do\n         end do\n      end do\n   end subroutine", file=file)

print("   subroutine reflect_moments_z(moments)\n      implicit none",file=file)
print("      real(WP), dimension(0:), intent(inout) :: moments !< Needs to be of size (0:188)", file=file)
print("      integer :: i,j,k,n\n      real(WP) :: temp", file=file)
print("      do k=0,2\n         do j=0,2\n            do i=0,2\n               if (k.eq.0) then\n               do n=0,6",file=file)
print("                  if (n.eq.3.or.n.eq.6) then\n                     temp = moments(7*(i*9+j*3+k)+n)\n                     moments(7*(i*9+j*3+k)+n) = -moments(7*(i*9+j*3+2)+n)\n                     moments(7*(i*9+j*3+2)+n) = -temp",file=file)
print("                  else\n                     temp = moments(7*(i*9+j*3+k)+n)\n                     moments(7*(i*9+j*3+k)+n) = moments(7*(i*9+j*3+2)+n)\n                     moments(7*(i*9+j*3+2)+n) = temp\n                  end if\n               end do",file=file)
print("               else if (k.eq.1) then\n                  moments(7*(i*9+j*3+k)+3) = -moments(7*(i*9+j*3+k)+3)",file=file)
print("                  moments(7*(i*9+j*3+k)+6) = -moments(7*(i*9+j*3+k)+6)\n               end if\n            end do\n         end do\n      end do\n   end subroutine", file=file)

print("end module plicnet", file=file)