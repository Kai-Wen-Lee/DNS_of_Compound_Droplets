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

with open("frac-12.400000.bin","rb") as file:
    # read first number, which gives the size N
    N = int(struct.unpack('f'*1, file.read(4*1))[0])
    size = (N+1)*(N+1)
    print('UX BIN init size:\t', N, '\tCELLS\t', size)

'''
with open("frac-12.400000.bin","rb") as file:
    numbers = np.fromfile("frac-0.100000.bin", dtype=np.float32)
'''

# Create a memory-mapped array for the binary file
numbers = np.memmap("frac-12.400000.bin", dtype=np.float32, mode="r")

# Print the first 10 values
print(numbers[:10])


#A = np.array(numbers)
A = np.reshape(numbers, (N+1,N+1))
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
