/*********************************************************************

poretda - A segmentation tool for porous structures using the topology
          toolkit (https://topology-tool-kit.github.io/)

Authors: Jorge Zorrilla Prieto, Aditya Vasudevan, Maciek Haranczyk 
IMDEA Materiales Institute

contact : maciej.haranczyk@imdea.org
          aditya.vasudevan@imdea.org
**********************************************************************/


//Visualization Toolkit(VTK) libraries
//-------------------------------------------------------------------------------------------------
#include <vtkNew.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkFieldData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkThreshold.h>
#include <vtkPointSet.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkPolyDataWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkDataSetWriter.h>
#include <vtkResampleToImage.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkFeatureEdges.h>
#include <vtkDataSetReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkThresholdPoints.h>
#include <vtkGaussianCubeReader2.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkAppendFilter.h>
#include <vtkAppendDataSets.h>
#include <vtkConnectivityFilter.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkArrayCalculator.h>
#include <vtkPointLocator.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiBlockDataGroupFilter.h>
#include <vtkSTLReader.h>
#include <vtkMultiBlockDataSet.h>
//#include <vtkGroupDataSetsFilter.h>
#include "vtkLogger.h"
#include "vtkInformation.h"
#include "vtkPartitionedDataSetCollection.h"
//-------------------------------------------------------------------------------------------------



//Topology Toolkit(TTK) libraries 
//-------------------------------------------------------------------------------------------------
#include <ttkMorseSmaleComplex.h>
#include <ttkArrayPreconditioning.h>
#include <ttkPeriodicGrid.h>
#include <ttkPersistenceDiagram.h>
#include <ttkTopologicalSimplification.h>
#include <ttkExtract.h>
#include <ttkFTMTree.h>
#include <ttkEigenField.h>
#include <ttkBottleneckDistance.h>
#include <ttkPersistenceDiagramClustering.h>
//-------------------------------------------------------------------------------------------------



//C++ Standard libraries included
//-------------------------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>
#include <cmath>
//#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <chrono>
#include <filesystem>
#include <dirent.h>
#include <sstream>
//-------------------------------------------------------------------------------------------------


//Namespaces used
//-------------------------------------------------------------------------------------------------
using namespace std;
using namespace std::chrono;
//-------------------------------------------------------------------------------------------------


#ifndef headers_h
#define headers_h


#endif /* headers_h */

