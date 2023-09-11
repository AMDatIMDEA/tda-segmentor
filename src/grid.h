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

/*
 
 This is the geometry class that stores the grid:
 
 - nx, ny, nz are the number of grid points
 - gridResolution[3][3] is the vector along the grid cell. For an orthorhombic
   lattice, this is purely diagonal, but for a general triclinic lattice, it has
   non-zero off diagonal terms
 - unitCellVectors[3][3] the points in the same direction of the gridResolution,
   but has a magnitude of the box lengths of the grid.
 - invUnitCellVectors[3][3] is the map that maps the normalized unit cube to
   the general triclinic unit cell.
 
 - cubicGrid is a vtkImageData grid that has the scalar values but defined on a unit cube
 - originalGrid is the grid in the actual coordinates of the nanoporous material.
 - gridPointsXYZ are the coordinates of the grid points of the originalGrid.
 
 - segmentation stores the segmentation data on actual coordinates.
 - criticalPoints stores the critical points data on actual coordinates. 
 
 */



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
