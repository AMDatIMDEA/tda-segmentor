import vtkmodules.all as vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np

segmentation = vtk.vtkDataSetReader()
segmentation.SetFileName("../Results/DDR/DDR_Segmentation.vtk")
segmentation.Update()

segmentation_grid = segmentation.GetOutput()
segmentation_array = segmentation_grid.GetPointData().GetArray(2)

density = vtk.vtkGaussianCubeReader2()
density.SetFileName("../Data/DDR_CH4_density.cube")
density.Update()

criticalPoints = vtk.vtkPolyDataReader()
criticalPoints.SetFileName("../Results/DDR/DDR_CriticalPoints.vtk")
criticalPoints.Update()

writer = vtk.vtkPolyDataWriter()
writer.SetFileName('../Data/DDR_graph.vtk')

density_grid = density.GetGridOutput()
density_array = density_grid.GetPointData().GetArray(0)

density_grid.SetSpacing(0.15, 0.15, 0.15)
segmentation_grid.SetSpacing(0.15, 0.15, 0.15)

accessible_criticals = vtk.vtkThresholdPoints()
accessible_criticals.SetInputConnection(criticalPoints.GetOutputPort(0))

accessible_criticals.SetInputArrayToProcess(0,0,0,0,"This is distance grid")

accessible_criticals.ThresholdBetween(-9e9,-1.2)

accessible_criticals.Update()

accessible_saddles = vtk.vtkThresholdPoints()
accessible_saddles.SetInputConnection(accessible_criticals.GetOutputPort(0))
accessible_saddles.SetInputArrayToProcess(0,0,0,0,"CellDimension")
accessible_saddles.ThresholdBetween(1.0,1.0)
accessible_saddles.Update()

accessible_min = vtk.vtkThresholdPoints()
accessible_min.SetInputConnection(accessible_criticals.GetOutputPort(0))
accessible_min.SetInputArrayToProcess(0,0,0,0,"CellDimension")
accessible_min.ThresholdBetween(0.0,0.0)
accessible_min.Update()

saddlesDataSet = vtk.vtkDataSet.SafeDownCast(accessible_saddles.GetOutputDataObject(0))

minDataSet = vtk.vtkDataSet.SafeDownCast(accessible_min.GetOutputDataObject(0))

n_saddles = saddlesDataSet.GetNumberOfPoints()
saddles_segments = []

pointLocator = vtk.vtkPointLocator()
pointLocator.SetDataSet(segmentation_grid)
pointLocator.BuildLocator()

for i in range(n_saddles):
    coords = (saddlesDataSet.GetPoint(i))
    closest_points = vtk.vtkIdList()
    pointLocator.FindPointsWithinRadius(0.3,np.array(coords)*0.15, closest_points)
    close_segments = set()
    for j in range(closest_points.GetNumberOfIds()):
        close_segments.add(segmentation_array.GetValue(closest_points.GetId(j)))
    saddles_segments.append(close_segments)
    
def get_relative_coords(saddleID, minID, saddlesDataSet, minDataSet, Resolution):

    min_coords = np.array(minDataSet.GetPoint(minID)) * Resolution
    saddle_coords = np.array(saddlesDataSet.GetPoint(saddleID)) * Resolution

    Lx = segmentation_grid.GetDimensions()[0]*Resolution
    Ly = segmentation_grid.GetDimensions()[1]*Resolution
    Lz = segmentation_grid.GetDimensions()[2]*Resolution

    X_rel_candidates = np.array([min_coords[0] - saddle_coords[0],
                            min_coords[0] - saddle_coords[0] + Lx,
                            min_coords[0] - saddle_coords[0] - Lx])

    X_rel = [i for i in X_rel_candidates if abs(i) == min(abs(X_rel_candidates))]

    Y_rel_candidates = np.array([min_coords[1] - saddle_coords[1],
                            min_coords[1] - saddle_coords[1] + Ly,
                            min_coords[1] - saddle_coords[1] - Ly])

    Y_rel = [i for i in Y_rel_candidates if abs(i) == min(abs(Y_rel_candidates))]

    Z_rel_candidates = np.array([min_coords[2] - saddle_coords[2],
                            min_coords[2] - saddle_coords[2] + Lz,
                            min_coords[2] - saddle_coords[2] - Lz])

    Z_rel = [i for i in Z_rel_candidates if abs(i) == min(abs(Z_rel_candidates))]
    
    return (X_rel[0], Y_rel[0], Z_rel[0])

n_mins = minDataSet.GetNumberOfPoints()
mins_segments = []

pointLocator = vtk.vtkPointLocator()
pointLocator.SetDataSet(segmentation_grid)
pointLocator.BuildLocator()

node_points = vtk.vtkPoints()
edges = vtk.vtkCellArray()
segmentID = vtk.vtkIntArray()
segmentID.SetName("SegmentID")
c = 0

for i in range(n_mins):
    
    coords = (minDataSet.GetPoint(i))
    node_points.InsertPoint(c, coords[0], coords[1], coords[2])
    
    current_node = c
    c += 1 
    
    closest_points = vtk.vtkIdList()
    pointLocator.FindPointsWithinRadius(0.3,np.array(coords)*0.15, closest_points)
    close_segments = set()
    for j in range(closest_points.GetNumberOfIds()):
        close_segments.add(segmentation_array.GetValue(closest_points.GetId(j)))
    print("segmentID: ", list(close_segments)[0])
    segmentID.InsertNextValue(list(close_segments)[0])
    print("minID: ", i)
    saddles_in_this_segment = []
    print("saddles (ID, relative coordinates, distance from min): ")
    for k in range(len(saddles_segments)):
        if list(close_segments)[0] in saddles_segments[k]:
            saddles_in_this_segment.append(k)
            rel_coords = get_relative_coords(k, i, saddlesDataSet, minDataSet, 1)
            print("    ",k, " , ", rel_coords, " , ", np.linalg.norm(np.array(rel_coords)))
            
            abs_coords =  np.array(coords) - np.array(rel_coords)
            
            node_points.InsertPoint(c, abs_coords[0], abs_coords[1], abs_coords[2])
            
            edges.InsertNextCell(2)
            edges.InsertCellPoint(current_node)
            edges.InsertCellPoint(c)
            
            segmentID.InsertNextValue(list(close_segments)[0])
            
            c+=1
    
    print("---")
    
g = vtk.vtkPolyData()
g.SetPoints(node_points)
g.SetLines(edges)
g.GetPointData().AddArray(segmentID)

writer.SetInputData(g)
writer.Write()
print("Results saved to " + writer.GetFileName())