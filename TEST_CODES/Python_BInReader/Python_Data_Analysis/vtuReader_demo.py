import vtk
import matplotlib.pyplot as plt
import numpy as np
# The source file
file_name = "frac-8.000000.bin"

# Read the source file.
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(file_name)
reader.Update()  # Needed because of GetScalarRange
output = reader.GetOutput()

print("Number of points:", output.GetNumberOfPoints())
print("Number of cells:", output.GetNumberOfCells())
f=[]
x=[]
y=[]

for n in range(output.GetNumberOfCells()):
    f_tuple = output.GetCellData().GetArray("f").GetTuple(n)
    f.append(f_tuple[0])

for m in range(output.GetNumberOfPoints()):
    cell_coord_tuple = output.GetPoints().GetData().GetTuple(m)
    x.append(cell_coord_tuple[0])
    y.append(cell_coord_tuple[1])

count=0

vol_frac = vtk.vtkContourFilter()
vol_frac.SetInputConnection(reader.GetOutputPort())
#vol_frac.SetValue(0, 0.5)  # Set the isosurface value to 0.5
vol_frac.Update()

# Extract the axisymmetric surface
axisym = vtk.vtkExtractEdges()
axisym.SetInputConnection(vol_frac.GetOutputPort())
axisym.Update()

# Calculate the axisymmetric volume
geom_filter = vtk.vtkGeometryFilter()
geom_filter.SetInputConnection(axisym.GetOutputPort())
geom_filter.Update()

vol = vtk.vtkMassProperties()
vol.SetInputData(geom_filter.GetOutput())
vol.Update()

axisym_volume = vol.GetVolume()

print("The axisymmetric volume is:", axisym_volume)
print(vol_frac)
X,Y=np.meshgrid(x,y)

'''
print(len(x),len(y))
with open('output.txt', 'w') as f:
    for p in range(len(x)):
        f.write(str(x[p]))
        f.write('\t')
        f.write(str(y[p]))

X,Y=np.meshgrid(x,y)
Z=[]
count=0
for j in Y:
    Z_inner=[]
    for i in X:
        count+=1
        Z_inner.append(count)
    Z.append(Z_inner)
        

plt.clf()
plt.pcolormesh(X,Y,Z, antialiased=False)
plt.colorbar()
plt.show()
'''
