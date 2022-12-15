import struct
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
dir_512 = "data_512"
dir_256 = "data_256"
dir_128 = "data_128"
dir_64 = "data_64"
dir_32 = "data_32"
dir_16 = "data_16"

dir_arr=[dir_16,dir_32,dir_64,dir_128,dir_256,dir_512]
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
    count=0
    
    #x=np.linspace(0.0,9.99,num=100)
    for dir_i in dir_arr:
        
        f_all,N,size=file_open(dir_i)
        vol_vf=[]
        count+=1
        print ((2.0**count)*8.0)

        for fi in f_all:
            x,y,f=resolve_xyf(fi,N)
            vol_vf.append(np.sum(f)*(20/(2.0**count)*8.0)*(20/(2.0**count)*8.0))
        plt.plot(vol_vf)
            

    plt.legend(dir_arr)
    plt.show()
