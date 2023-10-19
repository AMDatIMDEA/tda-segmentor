/*********************************************************************

TDA-Segmentor    A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:         Aditya Vasudevan (adityavv.iitkgp@gmail.com)
                 Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
                 Maciek Haranczyk (maciej.haranczyk@imdea.org)
          
                 IMDEA Materiales Institute
 
**********************************************************************/
#include "PEgrid.h"
#include "logger.h"


PEgrid::PEgrid(const parameters& p):
grid(p)
{
    arrayName = "Potential Energy";
    
};




PEgrid::~PEgrid(){
    
    
}




/**
 * @brief This routine writes the segment information corresponding to the void space 
 *        for an energy grid. 
 * For every segment, this routine also identifies the number of segments it is connected to. 
 * 
 * A file baseFileName + _Void_Segments.csv is written that has many columns:
 * regionID            - ID of the segment
 * x,y,z               - coordinates of the grid point
 * Scalar              - Value of the distance function at the grid point
 * RegionMaxValue      - Maximum value of the distance function in the segment
 * isMaximum, isSaddle - Notes if the grid point corresponds to a maxima or saddle
 * numberOfPoints      - number of grid points in the segment
 * Volume              - Volume of the segment 
 * numberOfConnections - number of segments the segmentID is connected to. 
 * 
 * Note : if DEBUG is set to 1 in grid.h, then the connectivity is printed in the log file. 
*/
void PEgrid::voidSegmentation(){
    
    logger::mainlog << "\n\nSegmentor: Void Segmentation Module" << "\n" << flush;
    
    ttk::Timer VoidSegmentationTimer;
    //Writer of the .csv results file
    ofstream misDatos;
    misDatos.open((Directory + "/" + baseFileName + "_Void_Segments.csv").c_str());
    assert(misDatos.is_open());
    misDatos << "regionID,x,y,z,Scalar,RegionMaxValue,isMaximum,isSaddle,numberOfPoints,Volume,numberOfConnections" << "\n";
    
    ofstream materialInfo;
    materialInfo.open( (Directory+ "/"+ baseFileName + "-materialInfo.csv").c_str());
    assert(materialInfo.is_open());
    materialInfo << "regionID,AverageMinimumEnergy,NumberOfPoints,NumberOfConexions,MaxMinEnergy,MinMinEnergy" << "\n";
    
    //--------------------------------------------------------------------------------------------------------
    double maxMinima=-9e10;
    double minMinima=0;
    double acumulatedEnergy = 0.0;
    double averageMinimaEnergy = 0.0;

    //Find the minimum critical points
    vtkSmartPointer<vtkThresholdPoints> minima = vtkSmartPointer<vtkThresholdPoints>::New();
    minima->SetInputData(criticalPoints);
    minima->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"CellDimension");
    minima->ThresholdBetween(0,0);
    minima->Update();
    //Find the minima on the void structure
    vtkSmartPointer<vtkThresholdPoints> minimaVoid = vtkSmartPointer<vtkThresholdPoints>::New();
    minimaVoid->SetInputConnection(minima->GetOutputPort(0));
    minimaVoid->SetInputArrayToProcess(0,0,0,0,"Potential Energy");
    minimaVoid->ThresholdBetween(-9e10,0.0);
    minimaVoid->Update();

    auto minimaDataSet = vtkDataSet::SafeDownCast(minimaVoid->GetOutputDataObject(0))->GetPointData()->GetAbstractArray("Potential Energy");
        
    for (size_t ii = 0; ii < minimaDataSet->GetNumberOfValues(); ii++)
    {
        
        double currentEnergy = minimaDataSet->GetVariantValue(ii).ToDouble();
        if (currentEnergy < minMinima){
            minMinima = currentEnergy;
        }
        if (currentEnergy > maxMinima)
        {
            maxMinima = currentEnergy;
        }
        acumulatedEnergy += currentEnergy;
            
    }

    averageMinimaEnergy = acumulatedEnergy/minimaDataSet->GetNumberOfValues();
        
    // Get the max bounding box of the voxel grid.     
    double cellDimensions[6];
    segmentation->GetCellBounds(0,cellDimensions);
    double cellSizeX = cellDimensions[1] - cellDimensions[0];
    double cellSizeY = cellDimensions[3] - cellDimensions[2];
    double cellSizeZ = cellDimensions[5] - cellDimensions[4];

    double cellSize = 0.0;
    if (cellSizeX > cellSizeY) cellSize = cellSizeX;
    else cellSize = cellSizeY;
    if (cellSizeZ > cellSize) cellSize = cellSizeZ;

    //Triangulate the segmentation to improve precision
    vtkSmartPointer<vtkDataSetTriangleFilter> triangulation = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
    triangulation->SetInputData(segmentation);
    triangulation->Update();

    // Each voxel can be fit with 6 tetrahedrons; so the unitCellVolume is 1/6th of the voxel volume. 
    double unitCellVolume = determinant(gridResolution)/6.0;

    //Segmentation corresponding to the void structure
    vtkSmartPointer<vtkThreshold> voidSegmentation = vtkSmartPointer<vtkThreshold>::New();
    voidSegmentation->SetInputConnection(triangulation->GetOutputPort());
    voidSegmentation->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    voidSegmentation->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    voidSegmentation->SetLowerThreshold(-9e10);
    voidSegmentation->SetUpperThreshold(0.0);
    voidSegmentation->Update();

    auto currentVoidDataSet = vtkDataSet::SafeDownCast(voidSegmentation->GetOutputDataObject(0));
    // Set of ID list of the manifolds
    vtkSmartPointer<vtkAbstractArray> descendingManifoldArray = currentVoidDataSet->GetPointData()->GetAbstractArray("DescendingManifold");
    
    std::set<int> descendingManifoldIDList;
    for (size_t i = 0; i < descendingManifoldArray->GetNumberOfValues(); i++ ){
        
        descendingManifoldIDList.insert(descendingManifoldArray->GetVariantValue(i).ToInt());
        
    }
    
    //Find the 1-saddle critical points
    vtkSmartPointer<vtkThresholdPoints> saddles = vtkSmartPointer<vtkThresholdPoints>::New();
    saddles->SetInputData(criticalPoints);
    saddles->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"CellDimension");
    saddles->ThresholdBetween(1,1);
    saddles->Update();
    
    //Find the 2-saddles on the void structure
    vtkSmartPointer<vtkThresholdPoints> positiveSaddles = vtkSmartPointer<vtkThresholdPoints>::New();
    positiveSaddles->SetInputConnection(saddles->GetOutputPort(0));
    positiveSaddles->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    positiveSaddles->ThresholdBetween(-9e10,0.0);
    positiveSaddles->Update();

    //DataSet of the saddles of the Ascending Segmentation of the void structure
    auto saddlesDataSet = vtkDataSet::SafeDownCast(positiveSaddles->GetOutputDataObject(0));
    logger::mainlog << "Number of accessible saddles: " << saddlesDataSet->GetNumberOfPoints() << endl;
    
    //2d vector to store the saddles id and the regions connected to them
    vector<set<int>> saddlesConnectivity;
    std::string manifoldName("DescendingManifold");
    getSaddleConnectivity(saddlesDataSet, currentVoidDataSet, manifoldName, cellSize, saddlesConnectivity);
    
    for (auto i : descendingManifoldIDList) 
    {
        int currentRegion = i;
        
        //Current Region of the Descending Segmentation
        vtkSmartPointer<vtkThreshold> sectionID = vtkSmartPointer<vtkThreshold>::New();
        sectionID->SetInputConnection(voidSegmentation->GetOutputPort());
        sectionID->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"DescendingManifold");
        sectionID->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
        sectionID->SetLowerThreshold(currentRegion);
        sectionID->SetUpperThreshold(currentRegion);
        sectionID->Update();

        //DataSet of the specific region of the Descending Segmentation
        auto sectionIDDataset = vtkDataSet::SafeDownCast(sectionID->GetOutputDataObject(0));

        //---------------------------------------------------------------------------
        int numberOfConnections = 0; //Number of connections of the current region
        vector<int> regionsSaddlesID; //ID of the points that work as Saddle in the region
        std::set <int> connectedSegments;
        if(sectionIDDataset->GetNumberOfPoints() > 0) //If not an empty region
        {
            for (size_t k = 0; k < saddlesDataSet->GetNumberOfPoints(); k++) 
            {
                
                auto it = find(saddlesConnectivity[k].begin(), saddlesConnectivity[k].end(), currentRegion);
 
                // If element was found and is not an isolated saddle
                if ((it != saddlesConnectivity[k].end()) && (saddlesConnectivity[k].size() > 1) )
                {
                    double currentSaddleCoords[3]; //Coordinates of the current saddle
                    saddlesDataSet->GetPoint(k,currentSaddleCoords);
                    int closestRegionPoint = sectionIDDataset->FindPoint(currentSaddleCoords); //ID of the closest point to the saddle from the region
                    regionsSaddlesID.push_back(closestRegionPoint);

                    ++numberOfConnections;

                    for (auto it2 : saddlesConnectivity[k]){
                        if (it2 != currentRegion){
                            connectedSegments.insert(it2);
                        }

                    }

                }

            }
        }

        //Array corresponding to the scalar values of the Region
        auto scalarValues = sectionIDDataset->GetPointData()->GetArray(arrayName.c_str());
        int segmentNumberOfCells = sectionIDDataset->GetNumberOfCells();
        double segmentVolume = segmentNumberOfCells * unitCellVolume;

        double minimumValue = 0.0;
        int minID; 
        for (size_t j = 0; j < sectionIDDataset->GetNumberOfPoints(); j++) 
        {
            double currentValue = scalarValues->GetVariantValue(j).ToDouble(); 
            if (currentValue <= minimumValue) 
            {
                minimumValue = currentValue;
                minID = j;
            }
        }
        bool isMinima;
        for (size_t j = 0; j < sectionIDDataset->GetNumberOfPoints(); j++) //For each of the points of the segment
        {
            bool isSaddle = false;
            
            for (size_t jj = 0; jj < regionsSaddlesID.size(); jj++)
            {
                if (j == regionsSaddlesID[jj])
                {
                    isSaddle = true;
                }
                
            }
            double pointCoords[3];
            sectionIDDataset->GetPoint(j,pointCoords);
            if(j == minID)
            {
                isMinima = 1;
            }
            else
            {
                isMinima = 0;
            }

            misDatos << currentRegion <<","<< pointCoords[0]<<","<<pointCoords[1]<<","<<pointCoords[2]<<","<<scalarValues->GetVariantValue(j).ToDouble()
                                      << "," << minimumValue<<","<< isMinima << "," << isSaddle <<","<< sectionIDDataset->GetNumberOfPoints()<< "," 
                                      << segmentVolume << "," << numberOfConnections<<"\n";

        }
        
        if (DEBUG)
        {
            logger::mainlog << "Number of connections for segment " << currentRegion << " is " << numberOfConnections  << " and number neighbor segments identified : " << connectedSegments.size();
            
            if (numberOfConnections > 0 ){
                logger::mainlog << " These are ";
                for (auto in : connectedSegments ){logger::mainlog << in << ", "; }
            }
            
            if (numberOfConnections != connectedSegments.size()){
                logger::mainlog << " It is possible that the two segments are connected more than once. ";
            }
            logger::mainlog << "\n";
        }
        
            materialInfo << currentRegion << "," << averageMinimaEnergy << "," << sectionIDDataset->GetNumberOfPoints()<< "," << numberOfConnections << "," << maxMinima << "," << minMinima << "\n";

        
    }
    
    misDatos.close();
    materialInfo.close();
    
    double elapsedTime = VoidSegmentationTimer.getElapsedTime();
    logger::mainlog << "Time elapsed in the void segmentation module: " << elapsedTime << "(s)\n" << flush;
    
};




void PEgrid::solidSegmentation(){
    
    logger::mainlog << "This method is not implemented for the potential energy grids!!" << endl;
    logger::errlog << "This method is not implemented for the potential energy grids!!" << endl;
    
    
};




void PEgrid::accessibleVoidSpace(double moleculeRadius, bool useAllCores){
    
    logger::mainlog << "This method is not implemented for the potential energy grids!!" << endl;
    logger::errlog << "This method is not implemented for the potential energy grids!!" << endl;
    
    
};




/**
 * @brief After the morse smale segmentation, this routine generates a graph
 *        representation of the accessible void space. 
 * @param moleculeRadius : Radius of the probe atom. 
 * @param useAllcores    : Use all the cores for the TTK modules. 
*/
void PEgrid::accessibleVoidGraph(double moleculeRadius, bool useAllCores){
    
    logger::mainlog << "\n\nSegmentor: Accessible Void Graph Module" << "\n" << flush;
    
    ttk::Timer voidGraphTimer;
    
    // Get the max bounding box of the voxel grid.     
    double cellDimensions[6];
    segmentation->GetCellBounds(0,cellDimensions);
    double cellSizeX = cellDimensions[1] - cellDimensions[0];
    double cellSizeY = cellDimensions[3] - cellDimensions[2];
    double cellSizeZ = cellDimensions[5] - cellDimensions[4];

    double cellSize = 0.0;
    if (cellSizeX > cellSizeY) cellSize = cellSizeX;
    else cellSize = cellSizeY;
    if (cellSizeZ > cellSize) cellSize = cellSizeZ;
    
    //Triangulate the segmentation to improve precision
    vtkSmartPointer<vtkDataSetTriangleFilter> triangulation = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
    triangulation->SetInputData(segmentation);
    triangulation->Update();

    //Segmentation corresponding to the void structure
    vtkSmartPointer<vtkThreshold> voidSegmentation = vtkSmartPointer<vtkThreshold>::New();
    voidSegmentation->SetInputConnection(triangulation->GetOutputPort());
    voidSegmentation->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    voidSegmentation->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    voidSegmentation->SetLowerThreshold(-9e9);
    voidSegmentation->SetUpperThreshold(0.0);
    voidSegmentation->Update();
    
    //Same structure segmentation but with a Field Data added
    
    auto currentVoidDataSet = vtkDataSet::SafeDownCast(voidSegmentation->GetOutputDataObject(0));

    vtkSmartPointer<vtkAbstractArray> descendingManifoldArray = currentVoidDataSet->GetPointData()->GetAbstractArray("DescendingManifold");
    
    std::set<int> descendingManifoldIDList;
    for (size_t i = 0; i < descendingManifoldArray->GetNumberOfValues(); i++ ){
        
        descendingManifoldIDList.insert(descendingManifoldArray->GetVariantValue(i).ToInt());
        
    }

    //Find the 1-saddle critical points
    vtkSmartPointer<vtkThresholdPoints> saddles = vtkSmartPointer<vtkThresholdPoints>::New();
    saddles->SetInputData(criticalPoints);
    saddles->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"CellDimension");
    saddles->ThresholdBetween(1.0,1.0);
    saddles->Update();
    //Find the 2-saddles on the void structure
    vtkSmartPointer<vtkThresholdPoints> accessibleSaddles = vtkSmartPointer<vtkThresholdPoints>::New();
    accessibleSaddles->SetInputConnection(saddles->GetOutputPort(0));
    accessibleSaddles->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    accessibleSaddles->ThresholdBetween(-9e9, 0.0);
    accessibleSaddles->Update();

    //DataSet of the accessible saddles to the molecule
    auto saddlesDataSet = vtkDataSet::SafeDownCast(accessibleSaddles->GetOutputDataObject(0));
    logger::mainlog << "Number of accessible saddles: " << saddlesDataSet->GetNumberOfPoints() << endl;
    
    // Next we find the minima critical points, and those that are accessible
    vtkSmartPointer<vtkThresholdPoints> minimas = vtkSmartPointer<vtkThresholdPoints>::New();
    minimas->SetInputData(criticalPoints);
    minimas->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"CellDimension");
    minimas->ThresholdBetween(0.0,0.0);
    minimas->Update();
    //Find the minima on the solid structure
    vtkSmartPointer<vtkThresholdPoints> accessibleMinimas = vtkSmartPointer<vtkThresholdPoints>::New();
    accessibleMinimas->SetInputConnection(minimas->GetOutputPort(0));
    accessibleMinimas->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    accessibleMinimas->ThresholdBetween(-9e9, 0.0);
    accessibleMinimas->Update();
    //DataSet of the accessible maximas to the molecule
    auto minimaDataSet = vtkDataSet::SafeDownCast(accessibleMinimas->GetOutputDataObject(0));
    logger::mainlog << "Number of accessible minimas: " << minimaDataSet->GetNumberOfPoints() << endl;

    //2d vector to store the saddles id and the regions connected to them
    vector<set<int>> saddlesConnectivity;
    std::string manifoldName("DescendingManifold");
    getSaddleConnectivity(saddlesDataSet, currentVoidDataSet, manifoldName, cellSize, saddlesConnectivity);
    
    // Next we find the segment that each minima belongs to:
    std::map<int,int> segmentIDofMinima;
    for (size_t k = 0; k < minimaDataSet->GetNumberOfPoints(); k++){
        
        double currentSaddleCoords[3]; //Coordinates of the current saddle
        minimaDataSet->GetPoint(k,currentSaddleCoords); //Save its coordinates
        vtkIdType closestPoint = currentVoidDataSet->FindPoint(currentSaddleCoords[0],currentSaddleCoords[1],currentSaddleCoords[2]);
        int currentClosestRegion = currentVoidDataSet->GetPointData()->GetAbstractArray("DescendingManifold")->GetVariantValue(closestPoint).ToInt();
        
        segmentIDofMinima.insert((std::pair<int,int> (k,currentClosestRegion)));
        
    }
    
    for (auto im : segmentIDofMinima){
        if (DEBUG) logger::mainlog << "For minima ID : " << im.first << ", segment ID is " << im.second << endl;
    }
    
    // Create the critical points as a point data
    vtkIdType nsaddles = saddlesDataSet->GetNumberOfPoints();
    vtkIdType nminimas = minimaDataSet->GetNumberOfPoints();
    
    vtkNew<vtkPoints> graphNodes;
    graphNodes->SetNumberOfPoints(nsaddles + nminimas);
    
    // We first set the saddle points
    for (size_t i = 0; i < nsaddles; i++){
        double coordinates[3];
        saddlesDataSet->GetPoint(i, coordinates);
        graphNodes->SetPoint(i, coordinates);
    }
    
    for (size_t i = 0; i < nminimas; i++){
        double coordinates[3];
        minimaDataSet->GetPoint(i, coordinates);
        int ipoint = nsaddles + i;
        graphNodes->SetPoint(ipoint, coordinates);
    }

    
    vtkNew<vtkUnstructuredGrid> graph;
    graph->SetPoints(graphNodes);
    
    for (size_t i = 0; i < nsaddles; i++){
        
        for (auto it : saddlesConnectivity[i]) {

            int connectedSegmentID = it;
            
            for (auto ip : segmentIDofMinima){
                if (ip.second == connectedSegmentID){
                    
                    int maximaID = ip.first;
                    // we join saddle ID and maxima ID
                    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                    int ipoint = nsaddles + maximaID;
                    line->GetPointIds()->SetId(0,i);
                    line->GetPointIds()->SetId(1,ipoint);
                    if (DEBUG) logger::mainlog << "Connecting saddle << " << i << " and maxima " << ipoint << endl;
                    graph->InsertNextCell(line->GetCellType(), line->GetPointIds());
                }
            }
        }
    }

    /* this generates a new graph that shows the edges that connect periodically by joining it with the neighoring box */
    vtkSmartPointer<vtkUnstructuredGrid> vizGraph = saveGraphForVisualization(graph);
    
    vtkSmartPointer<vtkUnstructuredGridWriter> vizGraphWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    vizGraphWriter->SetInputData(vizGraph);
    vizGraphWriter->SetFileName((Directory+"/" + baseFileName+"_viz_graph.vtk").c_str());
    vizGraphWriter->Write();

    /* save in a .nt2 format that can be used to post process*/
    writeGraphinNT2format(graph);
    
    double elapsedTime = voidGraphTimer.getElapsedTime();
    logger::mainlog << "Time elapsed in the accessible void graph module: " << elapsedTime << "(s)\n" << flush;
    
}




void PEgrid::accessibleSolidGraph(double moleculeRadius, bool useAllCores){

    logger::mainlog << "This method is not implemented for the potential energy grids!!" << endl;
    logger::errlog << "This method is not implemented for the potential energy grids!!" << endl;
    
}
