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

def format_array_content(content, max_len=1900):
    content = content.strip()
    if len(content) <= max_len:
        return content
    lines = []
    while len(content) > 0:
        if len(content) <= max_len:
            lines.append(content)
            break
        break_pos = content.rfind(',', 0, max_len)
        if break_pos == -1:
            break_pos = max_len
        lines.append(content[:break_pos+1] + "&")
        content = "   " + content[break_pos+1:].lstrip()
    return "\n".join(lines)

file = open("plicnet.f90", "w")
print("!> PLIC-Net File\n!> Provides the architecture for the neural network and the weights/biases\n!> Use generate_plicnet.py in NGA2/tools/scripts/plicnet to generate this file for a given Pytorch model",file=file)
print("module plicnet",file=file)
print("   use precision, only: WP\n   implicit none",file=file)
print("   integer :: idx", file=file)

param_info = []
count = 0
hidden = 0

for param in model.parameters():
    count = count + 1
    name = ""
    if count == 3:
        hidden = param.transpose(0,1).size()[0] #Assumes each hidden layer has same number of neurons
    
    if count%2 != 0:
        name = "lay" + str(int(count/2)+1) + "_weight"
        size1 = param.transpose(0,1).size()[0]
        size2 = param.transpose(0,1).size()[1]
        print(f"   real(WP), dimension({size1},{size2}), save :: {name}", file=file)
        param_info.append({'name': name, 'type': 'weight', 'data': param.detach().numpy()})
    else:
        name = "lay" + str(int((count-1)/2)+1) + "_bias"
        size1 = param.size()[0]
        print(f"   real(WP), dimension({size1}), save :: {name}", file=file)
        param_info.append({'name': name, 'type': 'bias', 'data': param.detach().numpy()})

bias_size = 1000
for info in param_info:
    name = info['name']
    data = info['data']
    if info['type'] == 'weight':
        num_rows = data.shape[1]
        num_cols = data.shape[0]
        for j in range(num_cols):
            col_data = data[j, :]
            values_str = np.array2string(col_data, separator=', ')[1:-1].strip()
            formatted_values = format_array_content(values_str)
            var_list_str = f"({name}(idx, {j+1}), idx=1, {num_rows})"
            print(f"   DATA {var_list_str} /&\n   {formatted_values}/", file=file)
    elif info['type'] == 'bias':
        total_size = data.shape[0]
        for i in range(0, total_size, bias_size):
            start_index = i
            end_index = min(i + bias_size, total_size)
            chunk_data = data[start_index:end_index]
            values_str = np.array2string(chunk_data, separator=', ')[1:-1].strip()
            formatted_values = format_array_content(values_str)
            var_list_str = f"({name}(idx), idx={start_index+1}, {end_index})"
            print(f"   DATA {var_list_str} /&\n   {formatted_values}/", file=file)

print("\n   contains\n   subroutine get_normal(moments,normal)\n      implicit none", file=file)
print("      real(WP), dimension(:), intent(in) :: moments !< Needs to be of size 189", file=file)
print("      real(WP), dimension(:), intent(out) :: normal !< Needs to be of size 3", file=file)
print("      real(WP), dimension(" + str(hidden) + ") :: tmparr", file=file)
print("      tmparr=max(0.0_WP,matmul(moments,lay1_weight)+lay1_bias)",file=file)
i=-1
for i in range(int((count-1)/2)-1):
    name1 = "lay" + str(i+2) + "_weight"
    name2 = "lay" + str(i+2) + "_bias"
    print("      tmparr=max(0.0_WP,matmul(tmparr," + name1 + ")+" + name2 + ")",file=file)
name1 = "lay" + str(i+3) + "_weight"
name2 = "lay" + str(i+3) + "_bias"
print("      normal=matmul(tmparr," + name1 + ")+" + name2,file=file)
print("   end subroutine", file=file)

reflect_subroutines = """   subroutine reflect_moments(moments,center,direction,direction2)
      implicit none
      real(WP), dimension(0:), intent(inout) :: moments !< Needs to be of size (0:188)
      real(WP), dimension(0:), intent(in) :: center     !< Needs to be of size (0:2)
      integer, intent(out) :: direction, direction2
      real(WP), dimension(0:2) :: new_center
      real(WP) :: temp
      direction=0
      direction2=0
      new_center = center
      if (abs(new_center(0)).le.1e-12) new_center(0)=0
      if (abs(new_center(1)).le.1e-12) new_center(1)=0
      if (abs(new_center(2)).le.1e-12) new_center(2)=0
      if (new_center(0).lt.0.and.new_center(1).ge.0.and.new_center(2).ge.0) then
         direction=1
         call reflect_moments_x(moments)
         new_center(0) = -new_center(0)
      else if (new_center(0).ge.0.and.new_center(1).lt.0.and.new_center(2).ge.0) then
         direction=2
         call reflect_moments_y(moments)
         new_center(1) = -new_center(1)
      else if (new_center(0).ge.0.and.new_center(1).ge.0.and.new_center(2).lt.0) then
         direction=3
         call reflect_moments_z(moments)
         new_center(2) = -new_center(2)
      else if (new_center(0).lt.0.and.new_center(1).lt.0.and.new_center(2).ge.0) then
         direction=4
         call reflect_moments_x(moments)
         call reflect_moments_y(moments)
         new_center(0) = -new_center(0)
         new_center(1) = -new_center(1)
      else if (new_center(0).lt.0.and.new_center(1).ge.0.and.new_center(2).lt.0) then
         direction=5
         call reflect_moments_x(moments)
         call reflect_moments_z(moments)
         new_center(0) = -new_center(0)
         new_center(2) = -new_center(2)
      else if (new_center(0).ge.0.and.new_center(1).lt.0.and.new_center(2).lt.0) then
         direction=6
         call reflect_moments_y(moments)
         call reflect_moments_z(moments)
         new_center(1) = -new_center(1)
         new_center(2) = -new_center(2)
      else if (new_center(0).lt.0.and.new_center(1).lt.0.and.new_center(2).lt.0) then
         direction=7
         call reflect_moments_x(moments)
         call reflect_moments_y(moments)
         call reflect_moments_z(moments)
         new_center(0) = -new_center(0)
         new_center(1) = -new_center(1)
         new_center(2) = -new_center(2)
      end if

      if (abs(new_center(0)-new_center(1)).le.1e-12.and.(new_center(0)-new_center(2)).gt.1e-12) then
         direction2=0
      else if (abs(new_center(1)-new_center(2)).le.1e-12.and.(new_center(0)-new_center(1)).gt.1e-12) then
         direction2=0
      else if (abs(new_center(0)-new_center(1)).le.1e-12.and.(new_center(2)-new_center(0)).gt.1e-12) then
         direction2=3
         call reflect_moments_xz(moments)
         temp = new_center(0)
         new_center(0) = new_center(2)
         new_center(2) = temp
      else if (abs(new_center(0)-new_center(2)).le.1e-12.and.(new_center(1)-new_center(0)).gt.1e-12) then
         direction2=1
         call reflect_moments_xy(moments)
         temp = new_center(0)
         new_center(0) = new_center(1)
         new_center(1) = temp
      else if (abs(new_center(0)-new_center(2)).le.1e-12.and.(new_center(0)-new_center(1)).gt.1e-12) then
         direction2=2
         call reflect_moments_yz(moments)
         temp = new_center(1)
         new_center(1) = new_center(2)
         new_center(2) = temp
      else if (abs(new_center(1)-new_center(2)).le.1e-12.and.(new_center(1)-new_center(0)).gt.1e-12) then
         direction2=3
         call reflect_moments_xz(moments)
         temp = new_center(0)
         new_center(0) = new_center(2)
         new_center(2) = temp
      else if (new_center(1).gt.new_center(0).and.new_center(0).ge.new_center(2)) then
         direction2=1
         call reflect_moments_xy(moments)
         temp = new_center(0)
         new_center(0) = new_center(1)
         new_center(1) = temp
      else if (new_center(2).gt.new_center(1).and.new_center(0).ge.new_center(2)) then
         direction2=2
         call reflect_moments_yz(moments)
         temp = new_center(1)
         new_center(1) = new_center(2)
         new_center(2) = temp
      else if (new_center(2).gt.new_center(1).and.new_center(1).ge.new_center(0)) then
         direction2=3
         call reflect_moments_xz(moments)
         temp = new_center(0)
         new_center(0) = new_center(2)
         new_center(2) = temp
      else if (new_center(1).gt.new_center(0)) then
         direction2=4
         call reflect_moments_xy(moments)
         call reflect_moments_yz(moments)
         temp = new_center(0)
         new_center(0) = new_center(1)
         new_center(1) = temp
         temp = new_center(1)
         new_center(1) = new_center(2)
         new_center(2) = temp
      else if (new_center(2).gt.new_center(1)) then
         direction2=5
         call reflect_moments_xy(moments)
         call reflect_moments_xz(moments)
         temp = new_center(0)
         new_center(0) = new_center(1)
         new_center(1) = temp
         temp = new_center(0)
         new_center(0) = new_center(2)
         new_center(2) = temp
      end if
   end subroutine reflect_moments
   subroutine reflect_moments_x(moments)
      implicit none
      real(WP), dimension(0:), intent(inout) :: moments !< Needs to be of size (0:188)
      integer :: i,j,k,n
      real(WP) :: temp
      do k=0,2
         do j=0,2
            do i=0,2
               if (i.eq.0) then
                  do n=0,6
                     if (n.eq.1.or.n.eq.4) then
                        temp=moments(7*(i*9+j*3+k)+n)
                        moments(7*(i*9+j*3+k)+n)=-moments(7*(2*9+j*3+k)+n)
                        moments(7*(2*9+j*3+k)+n)=-temp
                     else
                        temp=moments(7*(i*9+j*3+k)+n)
                        moments(7*(i*9+j*3+k)+n)=+moments(7*(2*9+j*3+k)+n)
                        moments(7*(2*9+j*3+k)+n)=+temp
                     end if
                  end do
               else if (i.eq.1) then
                  moments(7*(i*9+j*3+k)+1)=-moments(7*(i*9+j*3+k)+1)
                  moments(7*(i*9+j*3+k)+4)=-moments(7*(i*9+j*3+k)+4)
               end if
            end do
         end do
      end do
   end subroutine reflect_moments_x
   subroutine reflect_moments_y(moments)
      implicit none
      real(WP), dimension(0:), intent(inout) :: moments !< Needs to be of size (0:188)
      integer :: i,j,k,n
      real(WP) :: temp
      do k=0,2
         do j=0,2
            do i=0,2
               if (j.eq.0) then
                  do n=0,6
                     if (n.eq.2.or.n.eq.5) then
                        temp=moments(7*(i*9+j*3+k)+n)
                        moments(7*(i*9+j*3+k)+n)=-moments(7*(i*9+2*3+k)+n)
                        moments(7*(i*9+2*3+k)+n)=-temp
                     else
                        temp = moments(7*(i*9+j*3+k)+n)
                        moments(7*(i*9+j*3+k)+n)=+moments(7*(i*9+2*3+k)+n)
                        moments(7*(i*9+2*3+k)+n)=+temp
                     end if
                  end do
               else if (j.eq.1) then
                  moments(7*(i*9+j*3+k)+2)=-moments(7*(i*9+j*3+k)+2)
                  moments(7*(i*9+j*3+k)+5)=-moments(7*(i*9+j*3+k)+5)
               end if
            end do
         end do
      end do
   end subroutine reflect_moments_y
   subroutine reflect_moments_z(moments)
      implicit none
      real(WP), dimension(0:), intent(inout) :: moments !< Needs to be of size (0:188)
      integer :: i,j,k,n
      real(WP) :: temp
      do k=0,2
         do j=0,2
            do i=0,2
               if (k.eq.0) then
                  do n=0,6
                     if (n.eq.3.or.n.eq.6) then
                        temp=moments(7*(i*9+j*3+k)+n)
                        moments(7*(i*9+j*3+k)+n)=-moments(7*(i*9+j*3+2)+n)
                        moments(7*(i*9+j*3+2)+n)=-temp
                     else
                        temp=moments(7*(i*9+j*3+k)+n)
                        moments(7*(i*9+j*3+k)+n)=+moments(7*(i*9+j*3+2)+n)
                        moments(7*(i*9+j*3+2)+n)=+temp
                     end if
                  end do
               else if (k.eq.1) then
                  moments(7*(i*9+j*3+k)+3)=-moments(7*(i*9+j*3+k)+3)
                  moments(7*(i*9+j*3+k)+6)=-moments(7*(i*9+j*3+k)+6)
               end if
            end do
         end do
      end do
   end subroutine reflect_moments_z
   subroutine reflect_moments_xy(moments)
      implicit none
      real(WP), dimension(0:), intent(inout) :: moments !< Needs to be of size (0:188)
      integer :: i,j,k,n
      real(WP) :: temp
      do k=0,2
         do j=0,2
            do i=0,2
               if (i.eq.j) then
                  do n=0,6
                     if (n.eq.1.or.n.eq.4) then
                        temp=moments(7*(i*9+j*3+k)+n)
                        moments(7*(i*9+j*3+k)+n)=moments(7*(i*9+j*3+k)+n+1)
                        moments(7*(i*9+j*3+k)+n+1)=temp
                     end if
                  end do
               else if (i.gt.j) then
                  do n=0,6
                     if (n.eq.1.or.n.eq.4) then
                        temp=moments(7*(i*9+j*3+k)+n)
                        moments(7*(i*9+j*3+k)+n)=moments(7*(j*9+i*3+k)+n+1)
                        moments(7*(j*9+i*3+k)+n+1)=temp
                        temp = moments(7*(j*9+i*3+k)+n)
                        moments(7*(j*9+i*3+k)+n)=moments(7*(i*9+j*3+k)+n+1)
                        moments(7*(i*9+j*3+k)+n+1)=temp
                     else if (n.eq.0.or.n.eq.3.or.n.eq.6) then
                        temp=moments(7*(i*9+j*3+k)+n)
                        moments(7*(i*9+j*3+k)+n)=moments(7*(j*9+i*3+k)+n)
                        moments(7*(j*9+i*3+k)+n)=temp
                     end if
                  end do
               end if
            end do
         end do
      end do
   end subroutine reflect_moments_xy
   subroutine reflect_moments_yz(moments)
      implicit none
      real(WP), dimension(0:), intent(inout) :: moments !< Needs to be of size (0:188)
      integer :: i,j,k,n
      real(WP) :: temp
      do k=0,2
         do j=0,2
            do i=0,2
               if (j.eq.k) then
                  do n=0,6
                     if (n.eq.2.or.n.eq.5) then
                        temp=moments(7*(i*9+j*3+k)+n)
                        moments(7*(i*9+j*3+k)+n)=moments(7*(i*9+j*3+k)+n+1)
                        moments(7*(i*9+j*3+k)+n+1)=temp
                     end if
                  end do
               else if (j.gt.k) then
                  do n=0,6
                     if (n.eq.2.or.n.eq.5) then
                        temp=moments(7*(i*9+j*3+k)+n)
                        moments(7*(i*9+j*3+k)+n)=moments(7*(i*9+k*3+j)+n+1)
                        moments(7*(i*9+k*3+j)+n+1)=temp
                        temp = moments(7*(i*9+k*3+j)+n)
                        moments(7*(i*9+k*3+j)+n)=moments(7*(i*9+j*3+k)+n+1)
                        moments(7*(i*9+j*3+k)+n+1)=temp
                     else if (n.eq.0.or.n.eq.1.or.n.eq.4) then
                        temp=moments(7*(i*9+j*3+k)+n)
                        moments(7*(i*9+j*3+k)+n)=moments(7*(i*9+k*3+j)+n)
                        moments(7*(i*9+k*3+j)+n)=temp
                     end if
                  end do
               end if
            end do
         end do
      end do
   end subroutine reflect_moments_yz
   subroutine reflect_moments_xz(moments)
      implicit none
      real(WP), dimension(0:), intent(inout) :: moments !< Needs to be of size (0:188)
      integer :: i,j,k,n
      real(WP) :: temp
      do k=0,2
         do j=0,2
            do i=0,2
               if (i.eq.k) then
                  do n=0,6
                     if (n.eq.1.or.n.eq.4) then
                        temp=moments(7*(i*9+j*3+k)+n)
                        moments(7*(i*9+j*3+k)+n)=moments(7*(i*9+j*3+k)+n+2)
                        moments(7*(i*9+j*3+k)+n+2)=temp
                     end if
                  end do
               else if (i.gt.k) then
                  do n=0,6
                     if (n.eq.1.or.n.eq.4) then
                        temp=moments(7*(i*9+j*3+k)+n)
                        moments(7*(i*9+j*3+k)+n)=moments(7*(k*9+j*3+i)+n+2)
                        moments(7*(k*9+j*3+i)+n+2)=temp
                        temp = moments(7*(k*9+j*3+i)+n)
                        moments(7*(k*9+j*3+i)+n)=moments(7*(i*9+j*3+k)+n+2)
                        moments(7*(i*9+j*3+k)+n+2)=temp
                     else if (n.eq.0.or.n.eq.2.or.n.eq.5) then
                        temp=moments(7*(i*9+j*3+k)+n)
                        moments(7*(i*9+j*3+k)+n)=moments(7*(k*9+j*3+i)+n)
                        moments(7*(k*9+j*3+i)+n)=temp
                     end if
                  end do
               end if
            end do
         end do
      end do
   end subroutine reflect_moments_xz"""

print(reflect_subroutines, file=file)
print("end module plicnet", file=file)

file.close()
