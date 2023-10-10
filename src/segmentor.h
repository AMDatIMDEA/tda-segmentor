/*********************************************************************

TDA-Segmentor    A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:         Aditya Vasudevan (adityavv.iitkgp@gmail.com)J
                 Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
                 Maciek Haranczyk (maciej.haranczyk@imdea.org)
                 IMDEA Materiales Institute
 
**********************************************************************/

#ifndef segmentor_h
#define segmentor_h

#include <stdio.h>
#include "headers.h"
#include "logger.h"
#include "distanceGrid.h"
#include "PEgrid.h"
#include "grid.h"

/*
 The main analysis class that contains the functions of TTK.
 
 - Grid stores an instance of the grid, which can be either
   distanceGrid or PEgrid.
 - theMSC stores the results of the segmentation as a 3D grid, however this is in
   normalized coordinates of a unit cube. The segmentation results in the
   actual coordinates of the nanoporous materials is stored in grid->segmentation
 - thePersistenceDiagram and thePersistenceCurve are variables that stores the diagram
   and curve respectively which can be reused if necessary for other functions.
 - If the name of input file is FAU.cube, then BaseFileName = FAU, extensionName = .cube,
   Directory = segmentor + Basefile + .results.
 
 */

class segmentor
{
    
public:
        
    segmentor(const parameters &p);
    ~segmentor();
    
    grid*                                              readInputFile(const parameters &p, bool writeGridFile = true);
    vtkSmartPointer<vtkImageData>                      readFromCubeFile();
    void                                               readDistanceGrid(vtkSmartPointer<vtkImageData> imageData);
    void                                               readPEgrid(vtkSmartPointer<vtkImageData> imageData);
    vtkSmartPointer<ttkTriangulationManager>           generatePeriodicGrid(vtkSmartPointer<ttkScalarFieldSmoother> smoothgrid,bool periodicConditions,
                                                                            bool useAllCores);
    vtkSmartPointer<ttkScalarFieldSmoother>            smoothInputGrid(vtkSmartPointer<vtkImageData> grid, int niterations, bool useAllCores); 
    void                                               computePersistenceDiagram(vtkSmartPointer<ttkTriangulationManager> grid, bool useAllCores);
    
    void                                               MSC(vtkSmartPointer<ttkTriangulationManager> grid,double persistenceThreshold,
                                                           double saddlesaddleIncrement, bool writeOutputs, bool useAllCores);
    void                                               persistencecurve(vtkSmartPointer<ttkTriangulationManager> grid, bool useAllCores);
    // ----
    
    //auto                                               ftmtree(vtkSmartPointer<ttkTriangulationManager> grid, double persistenceThreshold, bool useAllCores);
    //void                                               accessibleVoidGraph(vtkSmartPointer<ttkFTMTree> ftmTree, double moleculeRadius, bool useAllCores);
    //void                                               accessibleSolidGraph(vtkSmartPointer<ttkFTMTree> ftmTree, bool useAllCores);

    

    
private:
    
    vtkIdType                                          getNumberOfDescendingManifolds();
    vtkIdType                                          getNumberOfAscendingManifolds();

    void                                               getGridResolutionFromCubeFile(double gridResolution[3][3]);
    void                                               getArrayName(std::string &arrayname);
    void                                               getArrayNameFromCubeFile(std::string &arrayname);
    void                                               getArrayNameFromVTIFile(std::string &arrayname);
    
    vtkSmartPointer<ttkMorseSmaleComplex>              theMSC;
    vtkSmartPointer<ttkPersistenceDiagram>             thePersistenceDiagram;
    vtkSmartPointer<ttkPersistenceCurve>               thePersistenceCurve;
    
    // Variables
    
    string                                              BaseFileName, fileName, extensionName;
    string                                              arrayName, Directory;
    grid*                                               Grid;
    bool                                                isPersistenceDiagramComputed;
    bool                                                writeFractionalGrid; 

};


#endif 

