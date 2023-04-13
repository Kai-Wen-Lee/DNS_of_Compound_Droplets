import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import glob
import os
import sys
import re
import vtk

lvl_1 = "1"
lvl_2 = "2"
lvl_3 = "3"
lvl_4 = "4"
lvl_5 = "5"
lvl_6 = "6"
lvl_7 = "7"
lvl_8 = "8"
lvl_9 = "9"
lvl_10 = "10"
dir_arr=[lvl_1,lvl_2,lvl_3,lvl_4,lvl_5,lvl_6,lvl_7,lvl_8,lvl_9,lvl_10]
#dir_arr=[lvl_1,lvl_2,lvl_3,lvl_4,lvl_5,lvl_6,lvl_7,lvl_8,lvl_9]
numbers_1 = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers_1.split(value)
    parts[1::2] = map(float, parts[1::2])
    return parts
    
def get_input (directory):
    vol=[]
    cell_area=[]
    cell_num=0
    filesindir=sorted(glob.iglob(f'{directory}/*'),key=numericalSort)
    for filename in filesindir:
        # Read the source file.
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(filename)
        reader.Update()  # Needed because of GetScalarRange
        output = reader.GetOutput()
        volume=0.0
        cell_size_filter = vtk.vtkCellSizeFilter()
        cell_size_filter.SetInputData(output)
        cell_size_filter.Update()
        output_cellsize = cell_size_filter.GetOutput()
        cell_sizes = output_cellsize.GetCellData().GetArray('Area')
        
        cell_centers = vtk.vtkCellCenters()
        cell_centers.SetInputData(output)
        cell_centers.Update()
        
        for n in range(output.GetNumberOfCells()):
            f_tuple = output.GetCellData().GetArray("f").GetTuple(n)
            cell = output.GetCell(n)
            area=cell_sizes.GetValue(n)
            centroid = cell_centers.GetOutput().GetPoint(n)
            volume += f_tuple[0] * area * 2.0 * np.pi * centroid[1]
            cell_area.append(area)
        vol.append(volume)
        cell_num+=output.GetNumberOfCells()
    return vol,cell_num, np.array(cell_area)

def numfmt(x, pos): # your custom formatter function: divide by 10.0
    s = '{}'.format(x / 10.0)
    return s
yfmt = tkr.FuncFormatter(numfmt)

def cell_no(v_dict):
    cat_arr = list(v_dict.keys())
    val_arr = np.array(list(v_dict.values()))
    plt.rcParams["figure.figsize"] = [8,6]
    plt.rcParams['figure.dpi'] = 160
    plt.bar(cat_arr,val_arr)
    plt.ylabel("Number of Cells")
    plt.show()


def cell_size_dist(input_dict):
    cat_arr = list(input_dict.keys())
    c_size_t = (20/256)**2
    val_arr = np.array(list(input_dict.values()), dtype=object)/c_size_t

    for i3 in range(len(cat_arr)):
        bins = 10**(np.linspace(np.log10(0.001), np.log10(100000), 100))
        plt.yscale('log')
        plt.xscale('log')
        plt.hist(val_arr[i3], bins=bins)
        plt.title(str(cat_arr[i3]))
        plt.xlabel("Relative Cell area")
        plt.ylabel("Number of Cells")
        plt.show()
def vol_analytical():
    U_l=1.0/3.07425133
    D_i=1.0
    D_rel=0.533333333
    D_o=D_i/D_rel
    
    x_t=np.linspace(0.0,30.0,num=600)
    x_t2=np.linspace(0.0,300.0,num=600)
    y_t=[]
    for xi in x_t:
        y_t1=(np.pi*((D_o**2)-(D_i**2))*U_l*xi)+np.pi*((D_o-D_i)/2)*(D_i+(D_o-D_i)/2)
        y_t.append(y_t1)
    return x_t2,y_t

def vol_vs_time(v_dict):
    plt.clf()
    cat_arr = list(v_dict.keys())
    val_arr = list(v_dict.values())
    #x_t,y_t=vol_analytical()
    for i in range(len(cat_arr)):
        plt.plot(val_arr[i], label='Case %i'%i+1)

    #plt.plot(x_t,y_t, label='Analytical Case')
    plt.legend()
    plt.title("Volume in computational domain against time")
    plt.xlabel("Dimensionless time, T")
    plt.ylabel("Dimensionless volume, V")
    plt.gca().xaxis.set_major_formatter(yfmt)
    plt.show()


def vol_err(v_dict):
    plt.clf()
    cat_arr = list(v_dict.keys())
    val_arr = list(v_dict.values())
    x_t,y_t=vol_analytical()  
    # Plot the magnitude of the difference
    for i in range(len(val_arr)):
        diff =(np.abs(np.array(val_arr[i])-np.array(y_t))/np.array(y_t))*100
        plt.plot(diff)
    plt.legend(cat_arr)
    plt.title("Percentage error per iteration against time")
    plt.xlabel("Dimensionless time, T")
    plt.ylabel("Percentage error, %")
    
    plt.gca().xaxis.set_major_formatter(yfmt)
    plt.show()

def vol_cum_err(v_dict):
    plt.clf()
    cat_arr = list(v_dict.keys())
    val_arr = list(v_dict.values())
    x_t,y_t=vol_analytical()
    
    for i in range(len(val_arr)):
        diff=np.abs(np.array(val_arr[i])-np.array(y_t))
        cum=0
        cum_list=[]
        for j in range(len(diff)):
             cum+=diff[j]
             cum_list.append(cum)
        plt.plot(cum_list)
    plt.legend(cat_arr)
    plt.title("Cumulative error (compared to analytical case) against time")
    plt.xlabel("Dimensionless time, T")
    plt.ylabel("Cumulative error")
    
    plt.gca().xaxis.set_major_formatter(yfmt)
    plt.show()
        
if __name__ == "__main__":
    plt.rcParams["figure.figsize"] = [8,6]
    plt.rcParams['figure.dpi'] = 160
    v_arr=[]
    c_num_arr=[]
    c_area_arr=[]
    for dir_i in dir_arr:
        v,c_num,a=get_input(dir_i)
        v_arr.append(v)
        c_num_arr.append(c_num)
        c_area_arr.append(a)
        print(dir_i)
    c_num_dict = dict(zip(dir_arr, c_num_arr))
    c_area_dict= dict(zip(dir_arr, c_area_arr))
    v_dict=dict(zip(dir_arr, v_arr))
    
    print(c_num_dict)
    #cell_no(c_num_dict)
    #cell_size_dist(c_area_dict)
    vol_vs_time(v_dict)
    #vol_err(v_dict)
    #vol_cum_err(v_dict)
