from pyevtk.hl import pointsToVTK
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
import os

from vtk import (
    vtkDataObject,
    vtkTableWriter,
    vtkThreshold,
    vtkXMLPolyDataWriter,
    vtkUnstructuredGridReader,
    vtkXMLUnstructuredGridReader,
    vtkXMLUnstructuredGridWriter,
    vtkXMLStructuredGridWriter,
    vtkDataSetReader,
    vtkUnstructuredGridWriter,
    vtkXMLStructuredGridReader
)

from topologytoolkit import (
    ttkMorseSmaleComplex,
    ttkPersistenceCurve,
    ttkPersistenceDiagram,
    ttkTopologicalSimplification,
)

# Example 1
npoints = 100
x = np.random.rand(npoints)
y = np.random.rand(npoints)
z = np.random.rand(npoints)
pressure = np.random.rand(npoints)
temp = np.random.rand(npoints)
pointsToVTK("./rnd_points", x, y, z, data={"temp": temp, "pressure": pressure})

# Example 2
x = np.arange(1.0, 10.0, 0.1)
y = np.arange(1.0, 10.0, 0.1)
z = np.arange(1.0, 10.0, 0.1)
pointsToVTK("./line_points", x, y, z, data={"elev": z})

#%%


# if len(sys.argv) == 2:
    
inputFilePath = os.path.join(os.path.dirname(__file__),"rnd_points.vtu")
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

# 8. saving the output data
# curveWriter = vtkTableWriter()
# curveWriter.SetInputConnection(curve.GetOutputPort())
# curveWriter.SetFileName("curve.vtk")
# curveWriter.Write()

# sepWriter = vtkXMLPolyDataWriter()
# sepWriter.SetInputConnection(morseSmaleComplex.GetOutputPort(1))
# sepWriter.SetFileName("separatrices.vtp")
# sepWriter.Write()

segWriter = vtkXMLUnstructuredGridWriter()
segWriter.SetInputConnection(morseSmaleComplex.GetOutputPort(3))
segWriter.SetFileName("segmentation.vtu")
segWriter.Write()




