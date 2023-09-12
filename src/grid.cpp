/*********************************************************************

TDA-Segmentor    A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:         Aditya Vasudevan (adityavv.iitkgp@gmail.com)
          Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
          Maciek Haranczyk (maciej.haranczyk@imdea.org)
          
          IMDEA Materiales Institute
 
**********************************************************************/

#include "grid.h"
#include "logger.h"


grid::grid (const parameters& p):
nx(),ny(),nz(),
baseFileName(),
extensionName(),
fileName(),
arrayName(),
cubicGrid(vtkSmartPointer<vtkImageData>::New()),
originalGrid(vtkSmartPointer<vtkStructuredGrid>::New()),
gridPointsXYZ(vtkSmartPointer<vtkPoints>::New()),
segmentation(vtkSmartPointer<vtkStructuredGrid>::New()),
criticalPoints(vtkSmartPointer<vtkPolyData>::New()),
criticalPointsData()
{
    
    fileName = p.inputfilename;
    baseFileName = p.basefilename;
    extensionName = p.extensionname;
    Directory = p.Directory;
    arrayName = p.arrayName;
    
    for (size_t i = 0; i < 3; i++){
        for (size_t j = 0; j < 3; j++){
            gridResolution[i][j] = 0.0;
            unitCellVectors[i][j] = 0.0;
            invUnitCellVectors[i][j] = 0.0;
        }
    }
    
    segmentationData = new segmentedData();

};




grid::~grid()
{
    delete segmentationData;
}




void grid::defineUnitCellVectors(){
    
    // Unit vector matrix of the unit cell
    unitCellVectors[0][0] = gridResolution[0][0]*(nx-1); unitCellVectors[1][0] = gridResolution[1][0]*(nx-1); unitCellVectors[2][0] = gridResolution[2][0]*(nx-1);
    unitCellVectors[0][1] = gridResolution[0][1]*(ny-1); unitCellVectors[1][1] = gridResolution[1][1]*(ny-1);; unitCellVectors[2][1] = gridResolution[2][1]*(ny-1);
    unitCellVectors[0][2] = gridResolution[0][2]*(nz-1); unitCellVectors[1][2] = gridResolution[1][2]*(nz-1); unitCellVectors[2][2] = gridResolution[2][2]*(nz-1);
    

    double D =  determinant(unitCellVectors);
    
    if (abs(D) < 1e-10){
        logger::mainlog << "Value of determinant is really small!" << endl;
        logger::errlog << "Value of determinant is really small!" << endl;
    }
    
    double invDet = 1/D;
    // Inverse unit cell matrix from abc to xyz
    invUnitCellVectors[0][0] = invDet*   (unitCellVectors[2][2]*unitCellVectors[1][1]-unitCellVectors[2][1]*unitCellVectors[1][2]);
    invUnitCellVectors[0][1] = invDet*-1*(unitCellVectors[2][2]*unitCellVectors[0][1]-unitCellVectors[2][1]*unitCellVectors[0][2]);
    invUnitCellVectors[0][2] = invDet*   (unitCellVectors[1][2]*unitCellVectors[0][1]-unitCellVectors[1][1]*unitCellVectors[0][2]);
    invUnitCellVectors[1][0] = invDet*-1*(unitCellVectors[2][2]*unitCellVectors[1][0]-unitCellVectors[2][0]*unitCellVectors[1][2]);
    invUnitCellVectors[1][1] = invDet*   (unitCellVectors[2][2]*unitCellVectors[0][0]-unitCellVectors[2][0]*unitCellVectors[0][2]);
    invUnitCellVectors[1][2] = invDet*-1*(unitCellVectors[1][2]*unitCellVectors[0][0]-unitCellVectors[1][0]*unitCellVectors[0][2]);
    invUnitCellVectors[2][0] = invDet*   (unitCellVectors[2][1]*unitCellVectors[1][0]-unitCellVectors[2][0]*unitCellVectors[1][1]);
    invUnitCellVectors[2][1] = invDet*-1*(unitCellVectors[2][1]*unitCellVectors[0][0]-unitCellVectors[2][0]*unitCellVectors[0][1]);
    invUnitCellVectors[2][2] = invDet*   (unitCellVectors[1][1]*unitCellVectors[0][0]-unitCellVectors[1][0]*unitCellVectors[0][1]);

    
}




void grid::generateSuperCell(){
    
    ttk::Timer superCellTimer;
    logger::mainlog << "\n\nGenerating Super Cell:         " << "\n";

    vtkSmartPointer<vtkImageAppend> appendX = vtkSmartPointer<vtkImageAppend>::New();
    appendX->AddInputData(cubicGrid);
    appendX->SetAppendAxis(0);
    appendX->AddInputData(cubicGrid);
    appendX->Update();
    
    vtkSmartPointer<vtkImageAppend> appendXY = vtkSmartPointer<vtkImageAppend>::New();
    appendXY->AddInputData(appendX->GetOutput());
    appendXY->SetAppendAxis(1);
    appendXY->AddInputData(appendX->GetOutput());
    appendXY->Update();
    
    vtkSmartPointer<vtkImageAppend> appendXYZ = vtkSmartPointer<vtkImageAppend>::New();
    appendXYZ->AddInputData(appendXY->GetOutput());
    appendXYZ->SetAppendAxis(2);
    appendXYZ->AddInputData(appendXY->GetOutput());
    appendXYZ->Update();

    vtkSmartPointer<vtkImageData> appendedImage = appendXYZ->GetOutput();
    vtkIdType cellDimsOriginal[3];
    cubicGrid->GetDimensions(cellDimsOriginal);
    logger::mainlog << "Number of points in the original grid  : (" << cellDimsOriginal[0]
                                                                   << " X " << cellDimsOriginal[1]
                                                                   << " X "<< cellDimsOriginal[2] << ")" << endl;
    vtkIdType cellDimsSuperCell[3];
    appendedImage->GetDimensions(cellDimsSuperCell);
    logger::mainlog << "Number of points in the super cell grid : (" << cellDimsSuperCell[0]
                                                                   << " X " << cellDimsSuperCell[1]
                                                                   << " X "<< cellDimsSuperCell[2] << ")" << endl;

    double superCellCreationTime = superCellTimer.getElapsedTime();
    superCellTimer.reStart();
    logger::mainlog << "Time taken for Super Cell creation            : " << superCellCreationTime << endl;
    
    cubicGrid->Initialize();
    cubicGrid = appendedImage;
    // We reinitialize the grid points, nx, ny, nz that now belongs to the super cell
    nx = cellDimsSuperCell[0]; ny = cellDimsSuperCell[1]; nz = cellDimsSuperCell[2];
    gridPointsXYZ->Initialize();
    // Store all the locations of the grid points for a general triclinic lattice.
    double x = 0.0, y = 0.0, z = 0.0;
    for (unsigned int k = 0; k < cellDimsSuperCell[2]; k++){
        for (unsigned int j = 0; j < cellDimsSuperCell[1]; j++){
            for (unsigned int i = 0; i < cellDimsSuperCell[0]; i++){
                
              x = gridResolution[0][0]*i + gridResolution[0][1]*j + gridResolution[0][2] * k;
              y = gridResolution[1][0]*i + gridResolution[1][1]*j + gridResolution[1][2] * k;
              z = gridResolution[2][0]*i + gridResolution[2][1]*j + gridResolution[2][2] * k;
              gridPointsXYZ->InsertNextPoint(x, y, z);
          }
        }
      }

    originalGrid->Initialize();
    vtkNew<vtkDoubleArray> pointValues;
    pointValues->SetNumberOfComponents(1);
    pointValues->SetNumberOfTuples(nx*ny*nz);
    
    for (size_t i = 0; i < (nx*ny*nz); ++i)
    {
      pointValues->SetValue(i, cubicGrid->GetPointData()->GetArray(arrayName.c_str())->GetVariantValue(i).ToDouble());
    }

    originalGrid->SetDimensions(static_cast<int>(nx), static_cast<int>(ny),
                                  static_cast<int>(nz));
    originalGrid->SetPoints(gridPointsXYZ);
    originalGrid->GetPointData()->SetScalars(pointValues);

   vtkNew<vtkStructuredGridWriter> strucGridWriter;
   strucGridWriter->SetInputData(originalGrid);
   strucGridWriter->SetFileName((Directory+"/"+baseFileName+"_superCellgrid.vtk").c_str());
   strucGridWriter->Write();
    
    double superCellWriteTime = superCellTimer.getElapsedTime();
    logger::mainlog << "Time taken to write Super Cell creation       : " << superCellWriteTime << endl;
    
    
    double totalTime = superCellCreationTime + superCellWriteTime;
    logger::mainlog << "Total Time in the supercell module            : " << totalTime << endl;

    
}



double determinant(double matrix[3][3]){
    
    double determinant =  matrix[0][0]*(matrix[2][2]*matrix[1][1] - matrix[2][1]*matrix[1][2])
                        - matrix[1][0]*(matrix[2][2]*matrix[0][1] - matrix[2][1]*matrix[0][2])
                        + matrix[2][0]*(matrix[1][2]*matrix[0][1] - matrix[1][1]*matrix[0][2]);
    
    return determinant;
    
}



void abcToxyz (double coordABC[3], double coordXYZ[3], double unitCellVectors[3][3]){
    
    coordXYZ[0] = coordABC[0]*unitCellVectors[0][0]+coordABC[1]*unitCellVectors[0][1]+coordABC[2]*unitCellVectors[0][2];
    coordXYZ[1] = coordABC[1]*unitCellVectors[1][1]+coordABC[2]*unitCellVectors[1][2];
    coordXYZ[2] = coordABC[2]*unitCellVectors[2][2];
    
    
}
