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

class segmentor
{
private:
    

    
    
public:
    
    //Cleaned functions
    
    segmentor(const parameters &p);
    ~segmentor();
    
    vtkIdType                                          getNumberOfDescendingManifolds(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex);
    vtkIdType                                          getNumberOfAscendingManifolds(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex);

    void                                               getGridResolutionFromCubeFile(double gridResolution[3][3]);
    void                                               getArrayName(std::string &arrayname);
    void                                               getArrayNameFromCubeFile(std::string &arrayname);
    void                                               getArrayNameFromVTIFile(std::string &arrayname);
    grid*                                              readInputFile(const parameters &p, bool writeGridFile = true);
    
    
    // --- The main TTK functions
    vtkSmartPointer<ttkTriangulationManager>           generatePeriodicGrid(vtkSmartPointer<vtkImageData> grid,bool periodicConditions,
                                                                            bool useAllCores);
    void                                               computePersistenceDiagram(vtkSmartPointer<ttkTriangulationManager> grid, bool useAllCores);
    
    void                                               MSC(vtkSmartPointer<ttkTriangulationManager> grid,double persistenceThreshold,
                                                           double saddlesaddleIncrement, bool writeOutputs, bool useAllCores);
    void                                               persistencecurve(vtkSmartPointer<ttkTriangulationManager> grid, bool useAllCores);
    // ----
    
    //auto                                               ftmtree(vtkSmartPointer<ttkTriangulationManager> grid, double persistenceThreshold, bool useAllCores);
    //void                                               accessibleVoidGraph(vtkSmartPointer<ttkFTMTree> ftmTree, double moleculeRadius, bool useAllCores);
    //void                                               accessibleSolidGraph(vtkSmartPointer<ttkFTMTree> ftmTree, bool useAllCores);

    
    // Variables
    
    string                                               BaseFileName, fileName, extensionName;
    string                                               arrayName, Directory;
    grid*                                                Grid;
    bool                                                 isPersistenceDiagramComputed;
    
    
    
    vtkSmartPointer<ttkMorseSmaleComplex>                theMSC;
    vtkSmartPointer<ttkPersistenceDiagram>               thePersistenceDiagram;
    vtkSmartPointer<ttkPersistenceCurve>                 thePersistenceCurve;

};


#endif 

