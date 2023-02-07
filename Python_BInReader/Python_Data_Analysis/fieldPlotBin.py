# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 15:37:10 2022

@author: irp1y14
"""

#import numpy as np

#z = np.fromfile('ux-2400.bin', dtype=float, count=256, sep='', offset=0)


import struct
import numpy as np
import matplotlib.pyplot as plt

with open("frac-4.000000.bin","rb") as file:
    # read first number, which gives the size N
    N = int(struct.unpack('f'*1, file.read(4*1))[0])
    size = (N+1)*(N+1)
    print('UX BIN init size:\t', N, '\tCELLS\t', size)
    

with open("frac-4.000000.bin","rb") as file:
    numbers = struct.unpack('f'*size, file.read(4*size))

A = np.array(numbers)
A = np.reshape(A, (N+1,N+1))
x = A[1:,0]
y = A[0,1:]
f = A[1:,1:]
f = np.rot90(f, k=1)
#f=np.flipud(f)
plt.clf()
plt.contourf(x,y,f, antialiased=False)
#plt.imshow(f)
#plt.imshow(f, vmin=0, vmax=1)
plt.colorbar()
plt.show()

vol_arr=[]
count_x=0
for x_i in x[1:]:
    count_y=0
    for y_i in y[1:]:
        dr=y_i-y[count_y]
        dx=x_i-x[count_x]
        vf_i=2*np.pi*dr*dx*y_i*f[count_x,count_y]
        #print(dr,dx,vf_i,f[count_x,count_y])
        vol_arr.append(vf_i)
        count_y+=1
        #print('countx:',str(count_x),'\t','county:',str(count_y),'\t','x_i',str(x_i),'\t','y_i',str(y_i),'\t','dr',str(dr),'\t','dx',str(dx),'\t','vf_i',str(vf_i),file=f)
        #print(fra_arr[count_x,count_y])
    count_x+=1
vol=np.sum(vol_arr)
print(vol)
