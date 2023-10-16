/*********************************************************************

TDA-Segmentor    A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:         Aditya Vasudevan (adityavv.iitkgp@gmail.com)
                 Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
                 Maciek Haranczyk (maciej.haranczyk@imdea.org)
          
                 IMDEA Materiales Institute
 
**********************************************************************/
#ifndef grid_hpp
#define grid_hpp

#include <stdio.h>
#include "headers.h"
#include "parameters.h"
#include "segmenteddata.h"

#define DEBUG 1


double    determinant(double matrix[3][3]);
void      abcToxyz (const double abc[3], double xyz[3], const double unitCellVectors[3][3]);
void      xyzToabc (const double xyz[3], double abc[3], const double invUnitCellVectors[3][3]);
/*
 
 This is the geometry class that stores the grid:
 
 - nx, ny, nz are the number of grid points
 - gridResolution[3][3] is the vector along the grid cell. For an orthorhombic
   lattice, this is purely diagonal, but for a general triclinic lattice, it has
   non-zero off diagonal terms
 - unitCellVectors[3][3] the points in the same direction of the gridResolution,
   but has a magnitude of the box lengths of the grid.
 - invUnitCellVectors[3][3] is the map that maps the cube in fractional coordinates to
   the general triclinic unit cell.
 
 - cubicGrid is a vtkImageData grid that has the scalar values but defined on a unit cube.
 - originalGrid is the grid in the actual coordinates of the nanoporous material.
 - gridPointsXYZ are the coordinates of the grid points of the originalGrid.
 
 - segmentation stores the segmentation data on actual coordinates.
 - criticalPoints stores the critical points data on actual coordinates.
 
 - criticalPointsData stores the same data as criticalPoints, but without the VTK wrapper, in native C++
 - segmentationData stores the segmentation data, but again without the VTK wrapper. 
 */



class grid {
    
public:
    
                                              grid(const parameters &p);
                                              
    
    virtual                                   ~grid();
    
    
    void                                      convertToSuperCell();
    void                                      defineUnitCellVectors();

    
    
    virtual void                              voidSegmentation() = 0;
    virtual void                              solidSegmentation() = 0;
    virtual void                              accessibleVoidSpace(double moleculeRadius, bool useAllCores) = 0;
    virtual void                              accessibleVoidGraph(double moleculeRadius, bool useAllCores) = 0;
    virtual void                              accessibleSolidGraph(double moleculeRadius, bool useAllCores) = 0;
    
    
    
    size_t                                    nx, ny, nz;
    double                                    gridResolution[3][3];
    double                                    unitCellVectors[3][3]; 
    double                                    invUnitCellVectors[3][3]; // This matrix maps a general triclinic unit cell to a cube of length unity.
    string                                    baseFileName, fileName, extensionName, Directory, arrayName;
    
    vtkSmartPointer<vtkImageData>             cubicGrid; // this is the general lattice mapped to a cube of length unity.
    vtkSmartPointer<vtkStructuredGrid>        originalGrid; // this is the grid, and in the most general case can also be triclinic.
    vtkSmartPointer<vtkPoints>                gridPointsXYZ; // grid points of the originalGrid.

    // Data from segmentation that can be used for post processing //
    vtkSmartPointer<vtkStructuredGrid>        segmentation;
    vtkSmartPointer<vtkPolyData>              criticalPoints;
    
    // Data from segmentation stored in native C++ format which can be used for post processing too //
    
    std::vector<criticalPoint>                criticalPointsData;
    segmentedData*                            segmentationData;

    protected:
    
    void                                      getSaddleConnectivity(vtkDataSet* saddlesDataSet, vtkDataSet* structureDataSet, 
                                                                    std::string manifoldName, double cellSize, 
                                                                    vector<set<int>> & saddlesConnectivity);
    void                                      writeGraphinNT2format(vtkSmartPointer<vtkUnstructuredGrid> graph);
    vtkSmartPointer<vtkUnstructuredGrid>      saveGraphForVisualization(vtkSmartPointer<vtkUnstructuredGrid> graph);
    
    
};



#endif /* grid_hpp */
