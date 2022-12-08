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

with open("ux-2400.bin","rb") as file:
    # read first number, which gives the size N
    N = int(struct.unpack('f'*1, file.read(4*1))[0])
    
print(N)
size = (N+1)*(N+1)
print(size)

with open("ux-2400.bin","rb") as file:
    numbers = struct.unpack('f'*size, file.read(4*size))

A = np.array(numbers)
A = np.reshape(A, (N+1,N+1))
x = A[1:,0]
y = A[0,1:]
f = A[1:,1:]
f = np.rot90(f, k=3)


plt.contourf(x,y,f, antialiased=False)
#plt.imshow(f)
#plt.imshow(f, vmin=0, vmax=1)
plt.show()
print (x,y)
