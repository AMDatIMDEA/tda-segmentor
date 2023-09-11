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

#define DEBUG 1


double    determinant(double matrix[3][3]);
void      abcToxyz (double abc[3], double xyz[3], double unitCellVectors[3][3]);


class grid {
    
public:
    
                                              grid(const parameters &p);
                                              
    
    virtual                                   ~grid();
    
    
    void                                      generateSuperCell( );
    void                                      defineUnitCellVectors();

    
    
    virtual void                              voidSegmentation() = 0;
    virtual void                              solidSegmentation() = 0;
    virtual void                              accessibleVoidSpace(double moleculeRadius, bool useAllCores) = 0;
    
    
    
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
    
    
};



#endif /* grid_hpp */
