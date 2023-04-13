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
volume=0.0


for n in range(output.GetNumberOfCells()):
    f_tuple = output.GetCellData().GetArray("f").GetTuple(n)
    f.append(f_tuple[0])
    cell = output.GetCell(n)
    cell_type = cell.GetCellType()
    
    # Calculate the area of the cell as the product of its width and heigh
    area = vtk.vtkTriangle.TriangleArea(cell.GetPoints().GetPoint(0),cell.GetPoints().GetPoint(1),cell.GetPoints().GetPoint(2))
    + vtk.vtkTriangle.TriangleArea(cell.GetPoints().GetPoint(0),cell.GetPoints().GetPoint(2),cell.GetPoints().GetPoint(3))
    volume += f_tuple[0] * area * 2.0 * np.pi
    

for m in range(output.GetNumberOfPoints()):
    cell_coord_tuple = output.GetPoints().GetData().GetTuple(m)
    x.append(cell_coord_tuple[0])
    y.append(cell_coord_tuple[1])
    

print("Total volume:", volume)
