#!/usr/bin/env python

#/// \ingroup examples
#/// \author Lutz Hofmann <lutz.hofmann@iwr.uni-heidelberg.de>
#/// \date August 2019.
#///
#/// \brief Minimalist python TTK example pipeline, including:
#///  -# The computation of a persistence curve
#///  -# The computation of a persistence diagram
#///  -# The selection of the most persistent pairs of the diagram
#///  -# The pre-simplification of the data according to this selection
#///  -# The computation of the Morse-Smale complex on this simplified data
#///  -# The storage of the output of this pipeline to disk.
#///
#/// This example reproduces the Figure 1 of the TTK companion paper:
#/// "The Topology ToolKit", J. Tierny, G. Favelier, J. Levine, C. Gueunet, M.
#/// Michaux., IEEE Transactions on Visualization and Computer Graphics, Proc.
#/// of IEEE VIS 2017.
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
import sys

from vtk import (
    vtkDataObject,
    vtkTableWriter,
    vtkThreshold,
    vtkXMLPolyDataWriter,
    vtkXMLUnstructuredGridReader,
    vtkXMLUnstructuredGridWriter,
    vtkDataSetReader,
)

from topologytoolkit import (
    ttkMorseSmaleComplex,
    ttkPersistenceCurve,
    ttkPersistenceDiagram,
    ttkTopologicalSimplification,
)

# if len(sys.argv) == 2:
inputFilePath = "D:\\Research Stuff\\VTK_Python_API\\PyhtonAPI_Examples\\inputData.vtu"
# else:
    # print("Missing mandatory argument: Path to input VTU file")
    # sys.exit()

# 1. loading the input data
reader = vtkXMLUnstructuredGridReader()
reader.SetFileName(inputFilePath)

# 2. computing the persistence curve
curve = ttkPersistenceCurve()
curve.SetInputConnection(reader.GetOutputPort())
curve.SetInputArrayToProcess(0, 0, 0, 0, "data")
curve.SetDebugLevel(3)

# 3. computing the persitence diagram
diagram = ttkPersistenceDiagram()
diagram.SetInputConnection(reader.GetOutputPort())
diagram.SetInputArrayToProcess(0, 0, 0, 0, "data")
diagram.SetDebugLevel(3)

# 4. selecting the critical point pairs
criticalPairs = vtkThreshold()
criticalPairs.SetInputConnection(diagram.GetOutputPort())
criticalPairs.SetInputArrayToProcess(
    0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_CELLS, "PairIdentifier")
criticalPairs.ThresholdBetween(-0.1, 999999)

# 5. selecting the most persistent pairs
persistentPairs = vtkThreshold()
persistentPairs.SetInputConnection(criticalPairs.GetOutputPort())
persistentPairs.SetInputArrayToProcess(
    0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_CELLS, "Persistence")
persistentPairs.ThresholdBetween(0.05, 999999)

# 6. simplifying the input data to remove non-persistent pairs
topologicalSimplification = ttkTopologicalSimplification()
topologicalSimplification.SetInputConnection(0, reader.GetOutputPort())
topologicalSimplification.SetInputArrayToProcess(0, 0, 0, 0, "data")
topologicalSimplification.SetInputConnection(1, persistentPairs.GetOutputPort())
topologicalSimplification.SetDebugLevel(3)

# 7. computing the Morse-Smale complex
morseSmaleComplex = ttkMorseSmaleComplex()
morseSmaleComplex.SetInputConnection(topologicalSimplification.GetOutputPort())
morseSmaleComplex.SetInputArrayToProcess(0, 0, 0, 0, "data")
morseSmaleComplex.SetDebugLevel(3)

#% try to obtain the number of arrays in the intput data set before saving 
number_of_arrays = reader.GetOutputDataObject(0).GetPointData().GetNumberOfArrays()
print("The number of arrays in the input vtu file (before saving outputs of MS complex): {}".format(number_of_arrays))

#%% 8. saving the output data
curveWriter = vtkTableWriter()
curveWriter.SetInputConnection(curve.GetOutputPort())
curveWriter.SetFileName("curve.vtk")
curveWriter.Write()

segWriter = vtkXMLUnstructuredGridWriter()
segWriter.SetInputConnection(morseSmaleComplex.GetOutputPort(1))
segWriter.SetFileName("separatrices.vtu")
segWriter.Write()

# segWriter = vtkXMLUnstructuredGridWriter()
# segWriter.SetInputConnection(morseSmaleComplex.GetOutputPort(0))
# segWriter.SetFileName("critical_points.vtu")
# segWriter.Write()

segWriter = vtkXMLUnstructuredGridWriter()
segWriter.SetInputConnection(morseSmaleComplex.GetOutputPort(3))
segWriter.SetFileName("segmentation.vtu")
segWriter.Write()

#% try to obtain the number of arrays in the intput data set after saving 
number_of_arrays = reader.GetOutputDataObject(0).GetPointData().GetNumberOfArrays()
print("The number of arrays in the input vtu file (after saving outputs of MS complex): {}".format(number_of_arrays))

#%% Input

#conver input data into numpy array

input_Data = {}
number_of_arrays = reader.GetOutputDataObject(0).GetPointData().GetNumberOfArrays()

for i in range(number_of_arrays):
    input_Data[reader.GetOutputDataObject(0).GetPointData().GetArrayName(i)]=vtk_to_numpy(reader.GetOutputDataObject(0).GetPointData().GetArray(i))

data = vtk_to_numpy(reader.GetOutputDataObject(0).GetPoints().GetData())
Xmin,Xmax,Ymin,Ymax,Zmin,Zmax = reader.GetOutputDataObject(0).GetPoints().GetBounds()

# plot the surface
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.plot_wireframe(data[:,0], data[:,1], data[:,2], color='black')


#%% Outputs
from vtk.util.numpy_support import vtk_to_numpy

# convert critical Points to numpy arrays

critical_points_output = {}
number_of_arrays = morseSmaleComplex.GetOutputDataObject(0).GetPointData().GetNumberOfArrays()

for i in range(number_of_arrays):
    critical_points_output[morseSmaleComplex.GetOutputDataObject(0).GetPointData().GetArrayName(i)]=vtk_to_numpy(morseSmaleComplex.GetOutputDataObject(0).GetPointData().GetArray(i))

# convert 1-D separatrices to numpy arrays
D1_sep_output_PD = {}  # pointdata
D1_sep_output_CD = {}  # celldata

number_of_arrays_PD = morseSmaleComplex.GetOutputDataObject(1).GetPointData().GetNumberOfArrays()
number_of_arrays_CD = morseSmaleComplex.GetOutputDataObject(1).GetCellData().GetNumberOfArrays()

for i in range(number_of_arrays_PD):
    D1_sep_output_PD[morseSmaleComplex.GetOutputDataObject(1).GetPointData().GetArrayName(i)]=vtk_to_numpy(morseSmaleComplex.GetOutputDataObject(1).GetPointData().GetArray(i))

for i in range(number_of_arrays_CD):
    D1_sep_output_CD[morseSmaleComplex.GetOutputDataObject(1).GetCellData().GetArrayName(i)]=vtk_to_numpy(morseSmaleComplex.GetOutputDataObject(1).GetCellData().GetArray(i))

# convert 2-D separatrices to numpy arrays
D2_sep_output_PD = {}  # pointdata
D2_sep_output_CD = {}  # celldata

number_of_arrays_PD = morseSmaleComplex.GetOutputDataObject(2).GetPointData().GetNumberOfArrays()
number_of_arrays_CD = morseSmaleComplex.GetOutputDataObject(2).GetCellData().GetNumberOfArrays()

for i in range(number_of_arrays_PD):
    D2_sep_output_PD[morseSmaleComplex.GetOutputDataObject(2).GetPointData().GetArrayName(i)]=vtk_to_numpy(morseSmaleComplex.GetOutputDataObject(2).GetPointData().GetArray(i))

for i in range(number_of_arrays_CD):
    D2_sep_output_CD[morseSmaleComplex.GetOutputDataObject(2).GetCellData().GetArrayName(i)]=vtk_to_numpy(morseSmaleComplex.GetOutputDataObject(2).GetCellData().GetArray(i))

# convert segmentation data into numpy arrays
Seg_output_PD = {}  # pointdata
Seg_output_CD = {}  # celldata

number_of_arrays_PD = morseSmaleComplex.GetOutputDataObject(3).GetPointData().GetNumberOfArrays()
number_of_arrays_CD = morseSmaleComplex.GetOutputDataObject(3).GetCellData().GetNumberOfArrays()

for i in range(number_of_arrays_PD):
    Seg_output_PD[morseSmaleComplex.GetOutputDataObject(3).GetPointData().GetArrayName(i)]=vtk_to_numpy(morseSmaleComplex.GetOutputDataObject(3).GetPointData().GetArray(i))

for i in range(number_of_arrays_CD):
    Seg_output_CD[morseSmaleComplex.GetOutputDataObject(3).GetCellData().GetArrayName(i)]=vtk_to_numpy(morseSmaleComplex.GetOutputDataObject(3).GetCellData().GetArray(i))

Scalars_segmentation = vtk_to_numpy(morseSmaleComplex.GetOutputDataObject(3).GetPointData().GetScalars())


#get the diagram output
number_of_arrays_PD = diagram.GetOutputDataObject(0).GetPointData().GetNumberOfArrays()
number_of_arrays_CD = diagram.GetOutputDataObject(0).GetCellData().GetNumberOfArrays()
Diag_output_PD = {}  # pointdata
Diag_output_CD = {}  # celldata

for i in range(number_of_arrays_PD):
    Diag_output_PD[diagram.GetOutputDataObject(0).GetPointData().GetArrayName(i)]=vtk_to_numpy(diagram.GetOutputDataObject(0).GetPointData().GetArray(i))

for i in range(number_of_arrays_CD):
    Diag_output_CD[diagram.GetOutputDataObject(0).GetCellData().GetArrayName(i)]=vtk_to_numpy(diagram.GetOutputDataObject(0).GetCellData().GetArray(i))


