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




void grid::convertToSuperCell(){
    
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
    // We reinitialize the grid points, nx, ny, nz that now belongs to the super cell, and reset the spacing
    nx = cellDimsSuperCell[0]; ny = cellDimsSuperCell[1]; nz = cellDimsSuperCell[2];
    cubicGrid->SetSpacing(1.0/(nx-1), 1.0/(ny-1), 1.0/(nz-1));
    // We also need to reinitialize the unit cell vectors that should now map the super cell
    // to the unit cube
    this->defineUnitCellVectors();
    
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



void abcToxyz (const double coordABC[3], double coordXYZ[3], const double unitCellVectors[3][3]){
    
    coordXYZ[0] = coordABC[0]*unitCellVectors[0][0]+coordABC[1]*unitCellVectors[0][1]+coordABC[2]*unitCellVectors[0][2];
    coordXYZ[1] = coordABC[1]*unitCellVectors[1][1]+coordABC[2]*unitCellVectors[1][2];
    coordXYZ[2] = coordABC[2]*unitCellVectors[2][2];
    
}



void xyzToabc(const double coordXYZ[3], double coordABC[3], const double invUnitCellVectors[3][3]){
    
    coordABC[0] = coordXYZ[0]*invUnitCellVectors[0][0]+coordXYZ[1]*invUnitCellVectors[0][1]+coordXYZ[2]*invUnitCellVectors[0][2];
    coordABC[1] = coordXYZ[1]*invUnitCellVectors[1][1]+coordXYZ[2]*invUnitCellVectors[1][2];
    coordABC[2] = coordXYZ[2]*invUnitCellVectors[2][2];
    
}



/**
 * @brief Typically saddle should lie on the boudary of two or more segments, and this routine 
 * finds all the segments which the saddle is connected to. 
 * 
 * **** Updates saddlesConnectivity variable. ****
 * 
 * Connectivtity is found out by drawing a sphere of certain radius and looking for grid points
 * within that radius, and then get its corresponding segment ID. 
 * 
 * @param saddlesDataSet           - the vtkDataSet of the saddles
 * @param structureDataSet         - the vtkDataSet of the structure, that can correspond to either void,
 *                                   solid, or accessible void
 * @param manifoldName             - Based on the nature of the analysis, this can be either 
 *                                   an ascending or descending manifold. 
 * @param cellSize                 - Size of the grid cell. 
 * @param saddlesConnectivity      - Stores all the segments the saddle is connected to.  
 * 
*/
void grid::getSaddleConnectivity(vtkDataSet* saddlesDataSet, vtkDataSet* structureDataSet, 
                            std::string manifoldName, double cellSize, 
                            vector<set<int>> & saddlesConnectivity)
{
    
    // Loop one by one on each saddle, and look for the grid points 
    // within a radius, find its segment ID, and store them in saddlesConnectivity for each saddle. 
    for (size_t k = 0; k < saddlesDataSet->GetNumberOfPoints(); k++) 
    {

        double currentSaddleCoords[3];                    
        saddlesDataSet->GetPoint(k, currentSaddleCoords); 

        vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
        pointLocator->SetDataSet(structureDataSet);
        pointLocator->BuildLocator();
        vtkSmartPointer<vtkIdList> closestPoints = vtkSmartPointer<vtkIdList>::New(); 
        double searchRadius = sqrt(3.0) * cellSize * 2; 
        pointLocator->FindPointsWithinRadius(searchRadius, currentSaddleCoords, closestPoints);

        std::set<int> closestRegionsToSaddle; 
        for (size_t kk = 0; kk < closestPoints->GetNumberOfIds(); kk++)
        {
            auto currentClosestRegion = structureDataSet->GetPointData()->GetAbstractArray(manifoldName.c_str())->GetVariantValue(closestPoints->GetId(kk)).ToInt();
            closestRegionsToSaddle.insert(currentClosestRegion);
        }

        saddlesConnectivity.push_back(closestRegionsToSaddle);

    }

    if (DEBUG)
    {
        for (size_t i = 0; i < saddlesDataSet->GetNumberOfPoints(); i++)
        {
            if (saddlesConnectivity[i].size() == 1 ){
                logger::mainlog << "Saddle " << i << " is inside segment " << *saddlesConnectivity[i].begin() << " and is not connected to any segment. " << endl; 
            } else {
                logger::mainlog << "Saddle " << i << " is connected to segments : ";
                for (auto it : saddlesConnectivity[i]){
                    logger::mainlog << it << ", ";
                }
                logger::mainlog << "\n";

            }

        }

    }
}




/**
 * @brief The graph generated by graph modules is periodic. When such a graph is plotted, edges across the box
 *        are plotted, which make it look ugly. Such periodic edges are identified and plotted by connecting it to
 *        the neighboring box instead. 
 * @param graph          : the is the periodic graph generated by the accessibleVoidGraph/accessibleSolidGraph modules.
*/
vtkSmartPointer<vtkUnstructuredGrid> grid::saveGraphForVisualization(vtkSmartPointer<vtkUnstructuredGrid> graph){

    //Create a new graph just for visualization, this shows the nodes outside the periodic box.
    vtkSmartPointer<vtkUnstructuredGrid> vizGraph = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkNew<vtkPoints> allNodes;
    vizGraph->SetPoints(allNodes);
    vtkCellIterator *it = graph->NewCellIterator();
    for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextCell())
     {
         if (it->GetCellType() == VTK_LINE)
         {
             vtkIdList *pointIds = it->GetPointIds();
             
             double p1[3], p2[3];
             graph->GetPoint(pointIds->GetId(0),p1); // coordinates of birth point
             graph->GetPoint(pointIds->GetId(1),p2); // coordinates of death point
             
             double pabc1[3], pabc2[3];
             xyzToabc(p1, pabc1, invUnitCellVectors);
             xyzToabc(p2, pabc2, invUnitCellVectors);
             
             vtkIdType id1 = allNodes->InsertNextPoint(p1);
             vtkIdType id2 = allNodes->InsertNextPoint(p2);
             
             
             int periodicity[3] = {0,0,0};
             double dp[3];for (size_t i = 0; i < 3; i++){
                 dp[i] = pabc2[i] - pabc1[i];
             }
             
             for (size_t i = 0; i < 3 ; i++){
                 if ( abs(dp[i]) > 0.5 )
                 {
                     if (dp[i] > 0.0) periodicity[i] = -1;
                     else if (dp[i] <  0.0) periodicity[i] = 1;
                 }
             }
             
             bool periodicityFlag = false;
             if (periodicity[0] != 0 || periodicity[1] != 0 || periodicity[2] != 0) periodicityFlag = true;
             
             double periodicNode1[3], periodicNode2[3];
             if (periodicityFlag){
                 
                 for(size_t i = 0; i < 3; i++){
                     
                     periodicNode2[i] = pabc2[i] + periodicity[i];
                     periodicNode1[i] = pabc1[i] - periodicity[i];
                 }
                 
                 double periodicNodeXYZ1[3], periodicNodeXYZ2[3];
                 abcToxyz(periodicNode1, periodicNodeXYZ1,unitCellVectors);
                 abcToxyz(periodicNode2, periodicNodeXYZ2,unitCellVectors);
                 vtkIdType id3 = allNodes->InsertNextPoint(periodicNodeXYZ1);
                 vtkIdType id4 = allNodes->InsertNextPoint(periodicNodeXYZ2);
                 vtkSmartPointer<vtkLine> imageLine1 = vtkSmartPointer<vtkLine>::New();
                 imageLine1->GetPointIds()->SetId(0,id1);
                 imageLine1->GetPointIds()->SetId(1,id4);
                 vtkSmartPointer<vtkLine> imageLine2 = vtkSmartPointer<vtkLine>::New();
                 imageLine2->GetPointIds()->SetId(0,id2);
                 imageLine2->GetPointIds()->SetId(1,id3);
                 
                 vizGraph->InsertNextCell(imageLine1->GetCellType(), imageLine1->GetPointIds());
                 vizGraph->InsertNextCell(imageLine2->GetCellType(), imageLine2->GetPointIds());
                 
             } else {
                 
                 vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                 line->GetPointIds()->SetId(0,id1);
                 line->GetPointIds()->SetId(1,id2);
                 
                 vizGraph->InsertNextCell(line->GetCellType(), line->GetPointIds());
             }
             
             
         }
         else {
             logger::mainlog << " Error in accessible Void Graph module: graph is not VTK_LINE" << endl;
             logger::errlog << " Error in accessible Void Graph module: graph is not VTK_LINE" << endl;
         }
         
     }
    it->Delete();

    return vizGraph;
}




/**
 * @brief Saves the graph in .nt2 format which is the following. 
 *        First all the nodes of the graph are stored with an ID, and its X, Y, Z, coordinates respectively.
 *        Then the edges are stored with the first two columns indicating the nodes they connect; first the birth and then the death.  
 *        Then the next three columns indicates periodicity in X, Y, Z, directions respectively. 
 *        If the value is 0, then it is not periodic along this axis. 
 *        If the value is 1, then the death connects the birth along the positive direction of the axis. 
 *        If the value is -1, then the death connects the birth along the negative direction of the axis. 
*/
void grid::writeGraphinNT2format(vtkSmartPointer<vtkUnstructuredGrid> graph){

    // save the graph in .nt2 format
    ofstream graphFile;
    std::string graphFileName = Directory + "/" + baseFileName + "-voidGraph" + ".nt2";
    graphFile.open((graphFileName).c_str());
    assert(graphFile.is_open());
    
    // We first write all the nodes
    graphFile << "Nodes: " << "\n";
    for (size_t i = 0; i < graph->GetNumberOfPoints(); i++){
        
        double coord[3];
        graph->GetPoint(i,coord);
        graphFile << i << " " << coord[0] << " " << coord[1] << " " << coord[2] << "\n";
    }
    
    
    vtkIdType cellDims[3];
    double spacing[3];
    cubicGrid->GetDimensions(cellDims);
    cubicGrid->GetSpacing(spacing);
    double boxLength[3];
    for (size_t i = 0; i < 3; i++){
        boxLength[i] = spacing[i] * (cellDims[i]-1);
    }

    
    // Next we store all the edges
    graphFile << "Edges: " << "\n";
    vtkCellIterator* it = graph->NewCellIterator();
    for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextCell())
     {
         if (it->GetCellType() == VTK_LINE)
         {
             vtkIdList *pointIds = it->GetPointIds();
             
             double p1[3], p2[3];
             graph->GetPoint(pointIds->GetId(0),p1); // coordinates of birth point
             graph->GetPoint(pointIds->GetId(1),p2); // coordinates of death point
             
             double pabc1[3], pabc2[3];
             xyzToabc(p1, pabc1, invUnitCellVectors);
             xyzToabc(p2, pabc2, invUnitCellVectors);

             int periodicity[3] = {0,0,0};
             double dp[3];for (size_t i = 0; i < 3; i++){
                 dp[i] = pabc2[i] - pabc1[i];
             }
             
             for (size_t i = 0; i < 3 ; i++){
                 if ( abs(dp[i]) > 0.5 * boxLength[i] )
                 {
                     if (dp[i] > 0.0) periodicity[i] = -1;
                     else if (dp[i] <  0.0) periodicity[i] = 1;
                 }
             }
             
             graphFile << pointIds->GetId(0) << " -> " << pointIds->GetId(1) << " "
                           << periodicity[0] << " " << periodicity[1] << " " << periodicity[2] << "\n";
         }
         else {
             logger::mainlog << " Error in accessible Void Graph module: graph is not VTK_LINE" << endl;
             logger::errlog << " Error in accessible Void Graph module: graph is not VTK_LINE" << endl;
         }
         
     }
    it->Delete();
    
    graphFile.close();
    logger::mainlog << "Graph is stored in the file " <<  graphFileName << endl;

}