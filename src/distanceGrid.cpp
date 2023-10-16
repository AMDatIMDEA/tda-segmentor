/*********************************************************************

TDA-Segmentor    A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:         Aditya Vasudevan (adityavv.iitkgp@gmail.com)
                 Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
                 Maciek Haranczyk (maciej.haranczyk@imdea.org)
          
                 IMDEA Materiales Institute
 
**********************************************************************/

#include "distanceGrid.h"
#include "logger.h"

distanceGrid::distanceGrid(const parameters& p):
grid(p)
{
    arrayName = "This is distance grid";
};




distanceGrid::~distanceGrid(){
    
}



/**
 * @brief This routine writes the segment information for the voids (distance function > 0.0)
 * For every segment, this routine also identifies the number of segments it is connected to. 
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
void distanceGrid::voidSegmentation(){
    
    logger::mainlog << "\n\nSegmentor: Void Segmentation Module" << "\n" << flush;
    
    ttk::Timer VoidSegmentationTimer;
    //Writer of the .csv results file
    ofstream misDatos;
    misDatos.open((Directory + "/" + baseFileName + "_Void_Segments.csv").c_str());
    assert(misDatos.is_open());
    misDatos << "regionID,x,y,z,Scalar,RegionMaxValue,isMaximum,isSaddle,numberOfPoints,Volume,numberOfConnections" << "\n";
    

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
    voidSegmentation->SetInputArrayToProcess(0,0,0,0,"This is distance grid");
    voidSegmentation->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    voidSegmentation->SetLowerThreshold(0.0);
    voidSegmentation->SetUpperThreshold(9e9);
    voidSegmentation->Update();
    
    auto currentVoidDataSet = vtkDataSet::SafeDownCast(voidSegmentation->GetOutputDataObject(0));

    vtkSmartPointer<vtkAbstractArray> ascendingManifoldArray = currentVoidDataSet->GetPointData()->GetAbstractArray("AscendingManifold");
    
    std::set<int> ascendingManifoldIDList;
    for (size_t i = 0; i < ascendingManifoldArray->GetNumberOfValues(); i++ ){
        
        ascendingManifoldIDList.insert(ascendingManifoldArray->GetVariantValue(i).ToInt());
        
    }
    
    //Find the 2-saddle critical points
    vtkSmartPointer<vtkThresholdPoints> saddles = vtkSmartPointer<vtkThresholdPoints>::New();
    saddles->SetInputData(criticalPoints);
    saddles->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"CellDimension");
    saddles->ThresholdBetween(2,2);
    saddles->Update();
    
    //Find the 2-saddles on the void structure
    vtkSmartPointer<vtkThresholdPoints> positiveSaddles = vtkSmartPointer<vtkThresholdPoints>::New();
    positiveSaddles->SetInputConnection(saddles->GetOutputPort(0));
    positiveSaddles->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    positiveSaddles->ThresholdBetween(0.0,9e9);
    positiveSaddles->Update();

    //DataSet of the saddles of the Ascending Segmentation of the void structure
    auto saddlesDataSet = vtkDataSet::SafeDownCast(positiveSaddles->GetOutputDataObject(0));
    logger::mainlog << "Number of accessible saddles: " << saddlesDataSet->GetNumberOfPoints() << endl;
    
    //2d vector to store the saddles id and the regions connected to them
    vector<set<int>> saddlesConnectivity;
    std::string manifoldName("AscendingManifold");
    getSaddleConnectivity(saddlesDataSet, currentVoidDataSet, manifoldName, cellSize, saddlesConnectivity);
    
    for (auto i : ascendingManifoldIDList) 
    {
        int currentRegion = i;
        
        vtkSmartPointer<vtkThreshold> sectionID = vtkSmartPointer<vtkThreshold>::New();
        sectionID->SetInputConnection(voidSegmentation->GetOutputPort());
        sectionID->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"AscendingManifold");
        sectionID->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
        sectionID->SetLowerThreshold(currentRegion);
        sectionID->SetUpperThreshold(currentRegion);
        sectionID->Update();

        auto sectionIDDataset = vtkDataSet::SafeDownCast(sectionID->GetOutputDataObject(0));

        int numberOfConnections = 0; 
        vector<int> regionsSaddlesID; 
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

        // Get the distance array in the segment    
        auto scalarValues = sectionIDDataset->GetPointData()->GetArray(arrayName.c_str());
        int segmentNumberOfCells = sectionIDDataset->GetNumberOfCells();
        double segmentVolume = segmentNumberOfCells * unitCellVolume; 

        // Get the maximumValue  and its index of the distance array in the segment. 
        double maximumValue = 0.0;
        int maxID;                                                         
        for (size_t j = 0; j < sectionIDDataset->GetNumberOfPoints(); j++) 
        {
            double currentValue = scalarValues->GetVariantValue(j).ToDouble(); 
            if (currentValue >= maximumValue)                                  
            {
                maximumValue = currentValue;
                maxID = j;
            }
        }

        bool isMaxima;
        for (size_t j = 0; j < sectionIDDataset->GetNumberOfPoints(); j++) // For each of the points of the segment
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
            sectionIDDataset->GetPoint(j, pointCoords);
            if (j == maxID)
            {
                isMaxima = 1;
            }
            else
            {
                isMaxima = 0;
            }

            misDatos << currentRegion << "," << pointCoords[0] << "," << pointCoords[1] << "," << pointCoords[2] << "," << scalarValues->GetVariantValue(j).ToDouble() << "," << maximumValue << "," << isMaxima << "," << isSaddle << "," << sectionIDDataset->GetNumberOfPoints() << "," << segmentVolume << "," << numberOfConnections << "\n";
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
    
        
    }
    
    misDatos.close();
    
    double elapsedTime = VoidSegmentationTimer.getElapsedTime();
    logger::mainlog << "Time elapsed in the void segmentation module: " << elapsedTime << "(s)\n" << flush;
    

    
    
    
};




/**
 * @brief This routine writes the segment information for the solids (distance function < 0.0)
 * For every segment, this routine also identifies the number of segments it is connected to. 
 * A file baseFileName + _Solid_Segments.csv is written that has many columns:
 * regionID            - ID of the segment
 * x,y,z               - coordinates of the grid point
 * Scalar              - Value of the distance function at the grid point
 * RegionMaxValue      - Maximum value of the distance function in the segment
 * isMinima, isSaddle - Notes if the grid point corresponds to a minima or saddle
 * numberOfPoints      - number of grid points in the segment
 * Volume              - Volume of the segment 
 * numberOfConnections - number of segments the segmentID is connected to. 
 * 
 * Note : if DEBUG is set to 1 in grid.h, then the connectivity is printed in the log file. 
*/
void distanceGrid::solidSegmentation(){
    
    logger::mainlog << "\n\nSegmentor: Solid Segmentation Module" << "\n" << flush;
    
    ttk::Timer SolidSegmentationTimer;
    //Writer of the .csv results file
    ofstream misDatos;
    misDatos.open((Directory+"/"+ baseFileName +"_Solid_Segments.csv").c_str());
    assert(misDatos.is_open());
    misDatos << "regionID,x,y,z,Scalar,RegionMinValue,isMinima,isSaddle,numberOfPoints,Volume,numberOfConnections" << "\n";

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
    vtkSmartPointer<vtkThreshold> solidSegmentation = vtkSmartPointer<vtkThreshold>::New();
    solidSegmentation->SetInputConnection(triangulation->GetOutputPort());
    solidSegmentation->SetInputArrayToProcess(0,0,0,0,"This is distance grid");
    solidSegmentation->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    solidSegmentation->SetLowerThreshold(-9e9);
    solidSegmentation->SetUpperThreshold(-1e-10);
    solidSegmentation->Update();

    auto currentSolidDataSet = vtkDataSet::SafeDownCast(solidSegmentation->GetOutputDataObject(0));
    
    // Set of ID list of the manifolds
    vtkSmartPointer<vtkAbstractArray> descendingManifoldArray = currentSolidDataSet->GetPointData()->GetAbstractArray("DescendingManifold");
    
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
    //Find the 1-saddles on the solid structure
    vtkSmartPointer<vtkThresholdPoints> negativeSaddles = vtkSmartPointer<vtkThresholdPoints>::New();
    negativeSaddles->SetInputConnection(saddles->GetOutputPort(0));
    negativeSaddles->SetInputArrayToProcess(0,0,0,0,"This is distance grid");
    negativeSaddles->ThresholdBetween(-9e10,-1e-10);
    negativeSaddles->Update();

    //DataSet of the saddles of the Descending Segmentation of the void structure
    auto saddlesDataSet = vtkDataSet::SafeDownCast(negativeSaddles->GetOutputDataObject(0));
    logger::mainlog << "Number of accessible saddles: " << saddlesDataSet->GetNumberOfPoints() << endl;

    //  A vector of sets to store the segments connected by each saddle. 
    vector<set<int>> saddlesConnectivity;
    std::string manifoldName("DescendingManifold");
    getSaddleConnectivity(saddlesDataSet, currentSolidDataSet, manifoldName, cellSize, saddlesConnectivity);


    for (auto i : descendingManifoldIDList) 
    {
        int currentRegion = i;
        
        //Current Region of the Descending Segmentation
        vtkSmartPointer<vtkThreshold> sectionID = vtkSmartPointer<vtkThreshold>::New();
        sectionID->SetInputConnection(solidSegmentation->GetOutputPort());
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
            for (size_t k = 0; k < saddlesDataSet->GetNumberOfPoints(); k++) //For each of the saddles
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
        
        int segmentNumberOfCells = sectionIDDataset->GetNumberOfCells();
        double segmentVolume = segmentNumberOfCells * unitCellVolume; 

        //Array corresponding to the scalar values of the Region
        auto scalarValues = sectionIDDataset->GetPointData()->GetArray("This is distance grid");
        // Find minimum value in the segment.
        double minimumValue = 10e9;
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

        
    }
    
    misDatos.close();
    
    double elapsedTime = SolidSegmentationTimer.getElapsedTime();
    logger::mainlog << "Time elapsed in the solid segmentation module: " << elapsedTime << "(s)\n" << flush;

    
    
};




/**
 * @brief This routine writes the segment information for the voids that are accessible to a probe
 *        of certain radius (distance function  > probe-radius)
 * @param This routine takes the radius of the probe atom as a parameter.  
 * For every segment, this routine also identifies the number of segments it is connected to. 
 * A file baseFileName + -avs-radius.csv is written that has many columns:
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
void distanceGrid::accessibleVoidSpace(double moleculeRadius, bool useAllCores){
    
    logger::mainlog << "\n\nSegmentor: Accessible Void Space Module" << "\n" << flush;
    
    ttk::Timer VoidSpaceTimer;
    //Writer of the .csv results file
    ofstream segmentResults;
    segmentResults.open((Directory + "/" + baseFileName + "-avs-" + to_string(moleculeRadius) + ".csv").c_str());
    assert(segmentResults.is_open());
    segmentResults << "regionID,x,y,z,Scalar,RegionMaxValue,isMaximum,isSaddle,numberOfPoints,Volume,numberOfConnections" << "\n";

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


    //Segmentation corresponding to the void structure accessible to the void space
    vtkSmartPointer<vtkThreshold> voidSegmentation = vtkSmartPointer<vtkThreshold>::New();
    voidSegmentation->SetInputConnection(triangulation->GetOutputPort());
    voidSegmentation->SetInputArrayToProcess(0,0,0,0,"This is distance grid");
    voidSegmentation->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    voidSegmentation->SetLowerThreshold(1.0*moleculeRadius);
    voidSegmentation->SetUpperThreshold(9e9);
    voidSegmentation->Update();

    //Same structure segmentation but with a Field Data added
    vtkSmartPointer<ttkExtract> accessibleSpace = vtkSmartPointer<ttkExtract>::New();
    accessibleSpace->SetUseAllCores(true);
    accessibleSpace->SetInputConnection(voidSegmentation->GetOutputPort());
    accessibleSpace->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"AscendingManifold");
    accessibleSpace->SetExtractionMode(3); //Array Values
    accessibleSpace->SetExtractUniqueValues(true);
    accessibleSpace->Update();

    auto accessibleSpaceDataSet = vtkDataSet::SafeDownCast(accessibleSpace->GetOutputDataObject(0));
    auto segmentsID = accessibleSpaceDataSet->GetFieldData()->GetAbstractArray("UniqueAscendingManifold");


    //Find the 2-saddle critical points
    vtkSmartPointer<vtkThresholdPoints> saddles = vtkSmartPointer<vtkThresholdPoints>::New();
    saddles->SetInputData(criticalPoints);
    saddles->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"CellDimension");
    saddles->ThresholdBetween(2.0,2.0);
    saddles->Update();
    //Find the 2-saddles on the void structure
    vtkSmartPointer<vtkThresholdPoints> accessibleSaddles = vtkSmartPointer<vtkThresholdPoints>::New();
    accessibleSaddles->SetInputConnection(saddles->GetOutputPort(0));
    accessibleSaddles->SetInputArrayToProcess(0,0,0,0,"This is distance grid");
    accessibleSaddles->ThresholdBetween(1.0*moleculeRadius, 9e9);
    accessibleSaddles->Update();

    //DataSet of the accessible saddles to the molecule
    auto saddlesDataSet = vtkDataSet::SafeDownCast(accessibleSaddles->GetOutputDataObject(0));
    logger::mainlog << "Number of accessible saddles: " << saddlesDataSet->GetNumberOfPoints() << endl;

        //2d vector to store the saddles id and the regions connected to them
    vector<set<int>> saddlesConnectivity;
    std::string manifoldName("AscendingManifold");
    getSaddleConnectivity(saddlesDataSet, accessibleSpaceDataSet, manifoldName, cellSize, saddlesConnectivity);

    for (size_t i = 0; i < segmentsID->GetNumberOfValues(); i++) //For each of the void segments
    {
        int currentRegion = segmentsID->GetVariantValue(i).ToInt();

        //Current Region of the Ascending Segmentation
        vtkSmartPointer<vtkThreshold> segment = vtkSmartPointer<vtkThreshold>::New();
        segment->SetInputConnection(voidSegmentation->GetOutputPort());
        segment->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"AscendingManifold");
        segment->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
        segment->SetLowerThreshold(currentRegion);
        segment->SetUpperThreshold(currentRegion);
        segment->Update();

        //DataSet of the specific region of the Descending Segmentation
        auto segmentDataset = vtkDataSet::SafeDownCast(segment->GetOutputDataObject(0));
        
        //---------------------------------------------------------------------------
        int numberOfConnections = 0; //Number of connections of the current region
        vector<int> regionsSaddlesID; 
        std::set <int> connectedSegments;
        if(segmentDataset->GetNumberOfPoints() > 0) //If not an empty region
        {
        
            for (size_t k = 0; k < saddlesDataSet->GetNumberOfPoints(); k++) //For each of the saddles
            {
                //Check if the segment is connected to a saddle
                auto it = find(saddlesConnectivity[k].begin(), saddlesConnectivity[k].end(), currentRegion);
                
                // If element was found and is not an isolated saddle
                if ((it != saddlesConnectivity[k].end()) && (saddlesConnectivity[k].size() > 1) )
                {
                    double currentSaddleCoords[3]; //Coordinates of the current saddle
                    saddlesDataSet->GetPoint(k,currentSaddleCoords);
                    int closestRegionPoint = segmentDataset->FindPoint(currentSaddleCoords); //ID of the closest point to the saddle from the region
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
        auto scalarValues = segmentDataset->GetPointData()->GetArray("This is distance grid");

        int segmentNumberOfCells = segmentDataset->GetNumberOfCells();
        double segmentsVolume = segmentNumberOfCells * unitCellVolume;


        // Get the maximumValue  and its index of the distance array in the segment. 
        double maximumValue = 0.0;
        int maxID;                                                         
        for (size_t j = 0; j < segmentDataset->GetNumberOfPoints(); j++) 
        {
            double currentValue = scalarValues->GetVariantValue(j).ToDouble(); 
            if (currentValue >= maximumValue)                                  
            {
                maximumValue = currentValue;
                maxID = j;
            }
        }

        bool isMaxima;
        for (size_t j = 0; j < segmentDataset->GetNumberOfPoints(); j++) // For each of the points of the segment
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
            segmentDataset->GetPoint(j, pointCoords);
            if (j == maxID)
            {
                isMaxima = 1;
            }
            else
            {
                isMaxima = 0;
            }

            segmentResults << currentRegion << "," << pointCoords[0] << "," << pointCoords[1] << "," << pointCoords[2] << "," << scalarValues->GetVariantValue(j).ToDouble() << "," << maximumValue << "," << isMaxima << "," 
                           << isSaddle << "," << segmentDataset->GetNumberOfPoints() << "," << segmentsVolume << "," << numberOfConnections << "\n";
        }
        
        // Print the segment connectivity for debugging
        if (DEBUG)
        {
            logger::mainlog << "Number of connections for segment " << currentRegion << " is " << numberOfConnections  << " and number neighbor segments identified : " << connectedSegments.size();
            
            if (numberOfConnections > 0 ){
                logger::mainlog << " These are ";
                for (auto in : connectedSegments ){logger::mainlog << in << ", "; }
            }
            
            if (numberOfConnections != connectedSegments.size()){
                logger::mainlog << "- It is possible that they are connected at multiple points ";
            }
            logger::mainlog << "\n";
        }
        
    }
    
    segmentResults.close();
    
    double elapsedTime = VoidSpaceTimer.getElapsedTime();
    logger::mainlog << "Time elapsed in the accesible void space module: " << elapsedTime << "(s) \n" << flush;
    
}




/**
 * @brief After the morse smale segmentation, this routine generates a graph
 *        representation of the accessible void space. 
 * @param moleculeRadius : Radius of the probe atom. 
 * @param useAllcores    : Use all the cores for the TTK modules. 
*/
void distanceGrid::accessibleVoidGraph(double moleculeRadius, bool useAllCores){
    
    logger::mainlog << "\n\nSegmentor: Accessible Void Graph Module" << "\n" << flush;
    
    ttk::Timer VoidGraphTimer;
    
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

    //Segmentation corresponding to the void structure accessible to the void space
    vtkSmartPointer<vtkThreshold> voidSegmentation = vtkSmartPointer<vtkThreshold>::New();
    voidSegmentation->SetInputConnection(triangulation->GetOutputPort());
    voidSegmentation->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    voidSegmentation->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    voidSegmentation->SetLowerThreshold(1.0*moleculeRadius);
    voidSegmentation->SetUpperThreshold(9e9);
    voidSegmentation->Update();
    
    //Same structure segmentation but with a Field Data added
    
    auto currentVoidDataSet = vtkDataSet::SafeDownCast(voidSegmentation->GetOutputDataObject(0));

    vtkSmartPointer<vtkAbstractArray> ascendingManifoldArray = currentVoidDataSet->GetPointData()->GetAbstractArray("AscendingManifold");
    
    std::set<int> ascendingManifoldIDList;
    for (size_t i = 0; i < ascendingManifoldArray->GetNumberOfValues(); i++ ){
        
        ascendingManifoldIDList.insert(ascendingManifoldArray->GetVariantValue(i).ToInt());
        
    }

    //Find the 2-saddle critical points
    vtkSmartPointer<vtkThresholdPoints> saddles = vtkSmartPointer<vtkThresholdPoints>::New();
    saddles->SetInputData(criticalPoints);
    saddles->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"CellDimension");
    saddles->ThresholdBetween(2.0,2.0);
    saddles->Update();
    //Find the 2-saddles on the void structure
    vtkSmartPointer<vtkThresholdPoints> accessibleSaddles = vtkSmartPointer<vtkThresholdPoints>::New();
    accessibleSaddles->SetInputConnection(saddles->GetOutputPort(0));
    accessibleSaddles->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    accessibleSaddles->ThresholdBetween(1.0*moleculeRadius, 9e9);
    accessibleSaddles->Update();

    //DataSet of the accessible saddles to the molecule
    auto saddlesDataSet = vtkDataSet::SafeDownCast(accessibleSaddles->GetOutputDataObject(0));
    logger::mainlog << "Number of accessible saddles: " << saddlesDataSet->GetNumberOfPoints() << endl;
    
    // Next we find the maxima critical points, and those that are accessible
    vtkSmartPointer<vtkThresholdPoints> maximas = vtkSmartPointer<vtkThresholdPoints>::New();
    maximas->SetInputData(criticalPoints);
    maximas->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"CellDimension");
    maximas->ThresholdBetween(3.0,3.0);
    maximas->Update();
    //Find the 2-saddles on the void structure
    vtkSmartPointer<vtkThresholdPoints> accessibleMaximas = vtkSmartPointer<vtkThresholdPoints>::New();
    accessibleMaximas->SetInputConnection(maximas->GetOutputPort(0));
    accessibleMaximas->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    accessibleMaximas->ThresholdBetween(1.0*moleculeRadius, 9e9);
    accessibleMaximas->Update();
    //DataSet of the accessible maximas to the molecule
    auto maximaDataSet = vtkDataSet::SafeDownCast(accessibleMaximas->GetOutputDataObject(0));
    logger::mainlog << "Number of accessible maximas: " << maximaDataSet->GetNumberOfPoints() << endl;

    //2d vector to store the saddles id and the regions connected to them
    vector<set<int>> saddlesConnectivity;
    std::string manifoldName("AscendingManifold");
    getSaddleConnectivity(saddlesDataSet, currentVoidDataSet, manifoldName, cellSize, saddlesConnectivity);

    // Next we find the segment that each maxima belongs to:
    std::map<int,int> segmentIDofMaxima;
    for (size_t k = 0; k < maximaDataSet->GetNumberOfPoints(); k++){
        
        double currentSaddleCoords[3]; //Coordinates of the current saddle
        maximaDataSet->GetPoint(k,currentSaddleCoords); //Save its coordinates
        vtkIdType closestPoint = currentVoidDataSet->FindPoint(currentSaddleCoords[0],currentSaddleCoords[1],currentSaddleCoords[2]);
        int currentClosestRegion = currentVoidDataSet->GetPointData()->GetAbstractArray("AscendingManifold")->GetVariantValue(closestPoint).ToInt();
        
        segmentIDofMaxima.insert((std::pair<int,int> (k,currentClosestRegion)));
        
    }
    
    for (auto im : segmentIDofMaxima){
        if (DEBUG) logger::mainlog << "For maxima ID : " << im.first << ", segment ID is " << im.second << endl;
    }
    
    // Create the critical points as a point data
    vtkIdType nsaddles = saddlesDataSet->GetNumberOfPoints();
    vtkIdType nmaximas = maximaDataSet->GetNumberOfPoints();
    
    vtkNew<vtkPoints> graphNodes;
    graphNodes->SetNumberOfPoints(nsaddles + nmaximas);
    
    // We first set the saddle points
    for (size_t i = 0; i < nsaddles; i++){
        double coordinates[3];
        saddlesDataSet->GetPoint(i, coordinates);
        graphNodes->SetPoint(i, coordinates);
    }
    
    for (size_t i = 0; i < nmaximas; i++){
        double coordinates[3];
        maximaDataSet->GetPoint(i, coordinates);
        int ipoint = nsaddles + i;
        graphNodes->SetPoint(ipoint, coordinates);
    }

    
    vtkNew<vtkUnstructuredGrid> graph;
    graph->SetPoints(graphNodes);
    
    for (size_t i = 0; i < nsaddles; i++){
        
        for (auto it : saddlesConnectivity[i]) {

            int connectedSegmentID = it;
            
            for (auto ip : segmentIDofMaxima){
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

    double elapsedTime = VoidGraphTimer.getElapsedTime();
    logger::mainlog << "Time elapsed in the accesible void graph module: " << elapsedTime << "(s) \n" << flush;

}




/**
 * @brief After the morse smale segmentation, this routine generates a graph
 *        representation of the accessible solid space. 
 * @param moleculeRadius : Radius of the probe atom. 
 *                        (this does not make sense for the solid, but still can be optionally 
 *                         entered as a small value like 1.0 or 2.0 to clean out some noise.)
 * @param useAllcores    : Use all the cores for the TTK modules. 
*/
void distanceGrid::accessibleSolidGraph(double moleculeRadius, bool useAllCores){
    
    /* We note that a saddle lies on the boundary of two or more segments and the minima
     lies inside the segment. So we ientify the saddle that connects the two segments
     and simply join it with the minima of that segment */
    
    logger::mainlog << "\n\nSegmentor: Accessible Solid Graph Module" << "\n" << flush;
    
    ttk::Timer solidGraphTimer;
    
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

    //Segmentation corresponding to the solid structure
    vtkSmartPointer<vtkThreshold> solidSegmentation = vtkSmartPointer<vtkThreshold>::New();
    solidSegmentation->SetInputConnection(triangulation->GetOutputPort());
    solidSegmentation->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    solidSegmentation->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    solidSegmentation->SetLowerThreshold(-9e9);
    solidSegmentation->SetUpperThreshold(-moleculeRadius);
    solidSegmentation->Update();
    
    //Same structure segmentation but with a Field Data added
    
    auto currentSolidDataSet = vtkDataSet::SafeDownCast(solidSegmentation->GetOutputDataObject(0));

    vtkSmartPointer<vtkAbstractArray> descendingManifoldArray = currentSolidDataSet->GetPointData()->GetAbstractArray("DescendingManifold");
    
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
    accessibleSaddles->ThresholdBetween(-9e9, -moleculeRadius);
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
    accessibleMinimas->ThresholdBetween(-9e9, -moleculeRadius);
    accessibleMinimas->Update();
    //DataSet of the accessible maximas to the molecule
    auto minimaDataSet = vtkDataSet::SafeDownCast(accessibleMinimas->GetOutputDataObject(0));
    logger::mainlog << "Number of accessible minimas: " << minimaDataSet->GetNumberOfPoints() << endl;

    //2d vector to store the saddles id and the regions connected to them
    vector<set<int>> saddlesConnectivity;
    std::string manifoldName("DescendingManifold");
    getSaddleConnectivity(saddlesDataSet, currentSolidDataSet, manifoldName, cellSize, saddlesConnectivity);

    // Next we find the segment that each minima belongs to:
    std::map<int,int> segmentIDofMinima;
    for (size_t k = 0; k < minimaDataSet->GetNumberOfPoints(); k++){
        
        double currentSaddleCoords[3]; //Coordinates of the current saddle
        minimaDataSet->GetPoint(k,currentSaddleCoords); //Save its coordinates
        vtkIdType closestPoint = currentSolidDataSet->FindPoint(currentSaddleCoords[0],currentSaddleCoords[1],currentSaddleCoords[2]);
        int currentClosestRegion = currentSolidDataSet->GetPointData()->GetAbstractArray("DescendingManifold")->GetVariantValue(closestPoint).ToInt();
        
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
                    
                    int minimaID = ip.first;
                    // we join saddle ID and maxima ID
                    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                    int ipoint = nsaddles + minimaID;
                    line->GetPointIds()->SetId(0,i);
                    line->GetPointIds()->SetId(1,ipoint);
                    if (DEBUG) logger::mainlog << "Connecting saddle << " << i << " and minima " << ipoint << endl;
                    graph->InsertNextCell(line->GetCellType(), line->GetPointIds());
                }
            }
        }
    }
    
    vtkSmartPointer<vtkUnstructuredGrid> vizGraph = saveGraphForVisualization(graph);
    
    vtkSmartPointer<vtkUnstructuredGridWriter> vizGraphWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    vizGraphWriter->SetInputData(vizGraph);
    vizGraphWriter->SetFileName((Directory+"/" + baseFileName+"_viz_graph.vtk").c_str());
    vizGraphWriter->Write();

    /* save in a .nt2 format that can be used to post process*/
    writeGraphinNT2format(graph);

    double elapsedTime = solidGraphTimer.getElapsedTime();
    logger::mainlog << "Time elapsed in the accesible void graph module: " << elapsedTime << "(s) \n" << flush;
}
