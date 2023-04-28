import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import vtk

filename ='/home/kaiwen/Documents/Data/TypeI_We85_v02/3/frac-45.300000.bin'

# Read the source file.
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(filename)
reader.Update()  # Needed because of GetScalarRange
output = reader.GetOutput()

cell_centers = vtk.vtkCellCenters()
cell_centers.SetInputData(output)
cell_centers.Update()

print(output.GetCellData().GetArray("Vect-u.x"))

# Get the cell center points and create a new point data array for the velocity vectors
center_points = cell_centers.GetOutput().GetPoints()
cell_velocities = vtk.vtkDoubleArray()
cell_velocities.SetName("cell_velocity")
cell_velocities.SetNumberOfComponents(3)
cell_velocities.SetNumberOfTuples(output.GetNumberOfCells())

coord=[]
u_arr=[]
for n in range(output.GetNumberOfCells()):
    f_tuple = output.GetCellData().GetArray("f").GetTuple(n)
    u_tuple = output.GetCellData().GetArray("Vect-u.x").GetTuple(n)
    coord.append(center_points.GetPoint(n))
    u_arr.append(u_tuple[0])  # use the x-component of velocity as the scalar value for contour plot

# Create the contour plot using Matplotlib
fig, ax = plt.subplots()
contour_plot = ax.tricontourf(np.array(coord)[:,0], np.array(coord)[:,1], u_arr)
fig.colorbar(contour_plot)
plt.show()
