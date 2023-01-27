import struct
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
directory = "data"

for filename in sorted(glob.iglob(f'{directory}/*'),key=os.path.basename):
    print(filename)
    with open(filename,"rb") as file:
    # read first number, which gives the size N
        N = int(struct.unpack('f'*1, file.read(4*1))[0])
        size = (N+1)*(N+1)
        print('UX BIN init size:\t', N, '\tCELLS\t', size)
    

    with open(filename,"rb") as file:
        numbers = struct.unpack('f'*size, file.read(16*size))
        print (numbers)

    A = np.array(numbers)
    A = np.reshape(A, (N+1,N+1))
    x = A[1:,0]
    y = A[0,1:]
    f = A[1:,1:]
    f = np.rot90(f, k=3)
    f = np.flip(f, 1)

 

    plt.clf()

    plt.contourf(x,y,f, antialiased=False)
    plt.colorbar()
    plt.title(filename)
    plt.show()


