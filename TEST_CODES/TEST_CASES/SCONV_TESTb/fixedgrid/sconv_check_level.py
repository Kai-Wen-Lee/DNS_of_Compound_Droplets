import struct
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
import re

lvl_1 = "level_1"
lvl_2 = "level_2"
lvl_3 = "level_3"
lvl_4 = "level_4"
lvl_5 = "level_5"
lvl_6 = "level_6"
lvl_7 = "level_7"
lvl_8 = "level_8"
lvl_9 = "level_9"
lvl_10 = "level_10"
dir_arr=[lvl_1,lvl_2,lvl_3,lvl_4,lvl_5,lvl_6,lvl_7,lvl_8,lvl_9,lvl_10]

numbers_1 = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers_1.split(value)
    parts[1::2] = map(float, parts[1::2])
    return parts
    
def file_open (directory):
    output_val=[]
    
    filesindir=sorted(glob.iglob(f'{directory}/*'),key=numericalSort)
    for filename in filesindir:
        
        with open(filename,"rb") as file:
            N = int(struct.unpack('f'*1, file.read(4*1))[0])
            size = (N+1)*(N+1)
        with open(filename,"rb") as file:
            numbers = struct.unpack('f'*size, file.read(4*size))

        output_val.append(numbers)

    return output_val,N,size

def resolve_xyf (frac_arr,N):
    A = np.array(frac_arr)
    A = np.reshape(A, (N+1,N+1))
    x = A[1:,0]
    y = A[0,1:]
    f = A[1:,1:]
    f = np.rot90(f, k=1)
    f=np.flipud(f)
    return x,y,f

def compute_vol (x,y,f):
    vol_arr=[]
    count_x=0
    for x_i in x[1:]:
        count_y=0
        for y_i in y[1:]:
            dr=y_i-y[count_y]
            dx=x_i-x[count_x]
            vf_i=2*np.pi*dr*dx*y_i*f[count_x,count_y]
            vol_arr.append(vf_i)
            count_y+=1
        count_x+=1
    vol=np.sum(vol_arr)

    return vol

if __name__ == "__main__":
    plt.clf()
    U_l=1.0/4.356435644
    D_i=1.0
    D_rel=0.732484076
    D_o=D_i/D_rel
    x_t=np.linspace(0.0,15.0,num=1000)
    x_p=np.linspace(0.0,150.0,num=1000)
    y_t=[]
    
    for xi in x_t:
        y_t.append(np.pi*((D_o**2)-(D_i**2))*U_l*xi)
    for dir_i in dir_arr:
        f_all,N,size=file_open(dir_i)
        vol_vf=[]

        for fi in f_all:
            x,y,f=resolve_xyf(fi,N)
            vol_vf.append(compute_vol(x,y,f))
        plt.plot(vol_vf)

    plt.plot(x_p,y_t)
    dir_arr.append('analytical')

    plt.legend(dir_arr)
    plt.show()
