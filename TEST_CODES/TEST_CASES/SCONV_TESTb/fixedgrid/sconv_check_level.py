import struct
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

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
def file_open (directory):
    output_val=[]

    filesindir=sorted(glob.iglob(f'{directory}/*'),key=os.path.basename)
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
    f = np.rot90(f, k=3)
    return x,y,f

if __name__ == "__main__":
    plt.clf()
    for dir_i in dir_arr:
        f_all,N,size=file_open(dir_i)
        vol_vf=[]

        for fi in f_all:
            x,y,f=resolve_xyf(fi,N)
            vol_vf.append(np.sum(f))
        plt.plot(vol_vf)
            
    plt.legend(dir_arr)
    plt.show()
