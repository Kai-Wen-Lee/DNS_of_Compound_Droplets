import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
import re
import vtk

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
#dir_arr=[lvl_8,lvl_9,lvl_10]
numbers_1 = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers_1.split(value)
    parts[1::2] = map(float, parts[1::2])
    return parts
    
def file_open (directory):
    vol=[]
    filesindir=sorted(glob.iglob(f'{directory}/*'),key=numericalSort)
    for filename in filesindir:
        # Read the source file.
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(filename)
        reader.Update()  # Needed because of GetScalarRange
        output = reader.GetOutput()
        volume=0.0
        for n in range(output.GetNumberOfCells()):
            f_tuple = output.GetCellData().GetArray("f").GetTuple(n)
            cell = output.GetCell(n)
            cell_type = cell.GetCellType()
            area = vtk.vtkTriangle.TriangleArea(cell.GetPoints().GetPoint(0),cell.GetPoints().GetPoint(1),cell.GetPoints().GetPoint(2))
            + vtk.vtkTriangle.TriangleArea(cell.GetPoints().GetPoint(0),cell.GetPoints().GetPoint(2),cell.GetPoints().GetPoint(3))
            #print(area)
            volume += f_tuple[0] * area * 2.0 * np.pi
        #print(filename,'\t',volume)
        vol.append(volume)
    return vol

def plot_err(e,t,dir_arr):
    plt.clf()    
    # Plot the magnitude of the difference
    for i in range(len(e)):
        diff =np.abs(np.array(e[i])-np.array(t))
        plt.plot(diff)
    plt.legend(dir_arr)
    plt.show()
    
def plot_cum_err(e,t,dir_arr):
    plt.clf()    
    # Plot the magnitude of the difference
   
    for i in range(len(e)):
        diff=np.abs(np.array(e[i])-np.array(t))
        cum=0
        cum_list=[]
        for j in range(len(diff)):
             cum+=diff[j]
             cum_list.append(cum)
        plt.plot(cum_list)
    plt.legend(dir_arr)
    plt.show()
    
if __name__ == "__main__":
    plt.clf()
    U_l=1.0/4.356435644
    D_i=1.0
    D_rel=0.732484076
    D_o=D_i/D_rel
    x_t=np.linspace(0.0,15.0,num=151)
    x_p=np.linspace(0.0,150.0,num=151)
    y_t=[]
    v_all=[]
    for xi in x_t:
        y_t.append(np.pi*((D_o**2)-(D_i**2))*U_l*xi)
    for dir_i in dir_arr:
        v=file_open(dir_i)
        plt.plot(v)
        v_all.append(v)
        print(len(v))
        print(dir_i)
    plt.plot(x_p,y_t)
    dir_arr.append('analytical')
    plt.legend(dir_arr)
    plt.show()
    plot_err(v_all,y_t,dir_arr)
    plot_cum_err(v_all,y_t,dir_arr)
