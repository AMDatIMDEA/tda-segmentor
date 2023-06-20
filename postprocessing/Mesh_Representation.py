import vtkmodules.all as vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np

Material = 'FAU'

segmentation = vtk.vtkDataSetReader()
segmentation.SetFileName("../Results/"+Material+"/"+Material+"_Segmentation.vtk")
segmentation.Update()

accessible = vtk.vtkThreshold()
accessible.SetInputConnection(segmentation.GetOutputPort(0))
accessible.SetInputArrayToProcess(0,0,0,0,"This is distance grid")
accessible.ThresholdBetween(-9e9,-1.2)
accessible.Update()

accessible_array = accessible.GetOutput().GetPointData().GetArray(2)

segmentIDs_list = np.unique(vtk_to_numpy(accessible_array))

(xmax, ymax, zmax) = segmentation.GetOutput().GetDimensions()

append = vtk.vtkAppendFilter()

for x in [0, xmax-1]:
    for y in [0, ymax-1]:
        for z in [0, zmax-1]:
            aTransform = vtk.vtkTransform()
            aTransform.Translate(x, y, z)
            
            transform = vtk.vtkTransformFilter()
            transform.SetInputConnection(accessible.GetOutputPort(0))
            transform.SetTransform(aTransform)
            transform.Update()
            
            append.AddInputConnection(transform.GetOutputPort(0))
            append.MergePointsOn()
            append.SetTolerance(0.0)
            append.Update()

triangulation = vtk.vtkDataSetTriangleFilter()
triangulation.SetInputConnection(append.GetOutputPort())
triangulation.Update()

final = vtk.vtkAppendFilter()


for segmentID in segmentIDs_list:
    
    currentSegment = vtk.vtkThreshold()
    currentSegment.SetInputConnection(triangulation.GetOutputPort(0))
    currentSegment.SetInputArrayToProcess(0,0,0,0,"DescendingManifold")
    currentSegment.ThresholdBetween(segmentID, segmentID)
    currentSegment.Update()
    
    surface2 = vtk.vtkDataSetSurfaceFilter()
    surface2.SetInputConnection(currentSegment.GetOutputPort(0))
    surface2.Update()
    
    region = vtk.vtkConnectivityFilter()
    region.SetExtractionModeToLargestRegion()
    region.SetInputConnection(surface2.GetOutputPort(0))
    region.Update()
    
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName('../Data/'+Material+'_surface'+str(segmentID)+'.vtk')
    writer.SetInputData(region.GetOutput())
    writer.Write()