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




void distanceGrid::voidSegmentation(){
    
    logger::mainlog << "\n\nSegmentor: Void Segmentation Module" << "\n" << flush;
    
    ttk::Timer VoidSegmentationTimer;
    //Writer of the .csv results file
    ofstream misDatos;
    misDatos.open((Directory + "/" + baseFileName + "_Void_Segments.csv").c_str());
    assert(misDatos.is_open());
    misDatos << "regionID,x,y,z,Scalar,RegionMaxValue,isMaximum,isSaddle,numberOfPoints,numberOfConnections,xScaled,yScaled,zScaled" << "\n";
    
    
    //Compute cell dimensions of the input file
    //---------------------------------------------------------------------------------------------
    double cellDimensions[6];
    segmentation->GetCellBounds(0,cellDimensions);
    //Cell size of the current dataset
    double cellSize = cellDimensions[1] - cellDimensions[0];
    //---------------------------------------------------------------------------------------------


    //Segmentation corresponding to the void structure
    vtkSmartPointer<vtkThresholdPoints> voidSegmentation = vtkSmartPointer<vtkThresholdPoints>::New();
    voidSegmentation->SetInputData(segmentation);
    voidSegmentation->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    voidSegmentation->ThresholdBetween(0.0,9e9);
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
    vector<vector<int>> saddlesConnectivity;
    //Set default values to -1.0
    size_t nmaxneighbors = 20;
    saddlesConnectivity.resize(saddlesDataSet->GetNumberOfPoints(),vector<int>(nmaxneighbors,-1.0));
    std::map<int,int> regionsWithSaddleInside;
    for (size_t k = 0; k < saddlesDataSet->GetNumberOfPoints(); k++) //For each of the saddles
    {
                
        double currentSaddleCoords[3]; //Coordinates of the current saddle
        saddlesDataSet->GetPoint(k,currentSaddleCoords);//Save the current saddle coordinates
        
        //Check that this saddle is not noise inside the region
        vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
        pointLocator->SetDataSet(currentVoidDataSet);
        pointLocator->BuildLocator();
        vtkSmartPointer<vtkIdList> closestPoints = vtkSmartPointer<vtkIdList>::New(); //IDs of the closest points to the saddle in the void structure
        //Find the in the void structure the closest points to the saddle inside a sphere of radius equal to the cell size
        pointLocator->FindPointsWithinRadius(sqrt(2.0) * cellSize,currentSaddleCoords,closestPoints);

        vector<int> closestRegionsToSaddle; //Closest Regions ID to the saddle
        //logger::mainlog << "Current Saddle ID: " << k << endl;
        for (size_t kk = 0; kk < closestPoints->GetNumberOfIds(); kk++)
        {
            auto currentClosestRegion = currentVoidDataSet->GetPointData()->GetAbstractArray("AscendingManifold")->GetVariantValue(closestPoints->GetId(kk)).ToInt();
            //logger::mainlog << currentClosestRegion << endl;
            closestRegionsToSaddle.push_back(currentClosestRegion);
        }
        sort(closestRegionsToSaddle.begin(), closestRegionsToSaddle.end()); //Order the values of the segmentation
        vector<int>::iterator it;
        it = unique(closestRegionsToSaddle.begin(), closestRegionsToSaddle.end());  //Delete repeated values
        closestRegionsToSaddle.resize(distance(closestRegionsToSaddle.begin(),it)); //Resize with the unique values
        if (closestRegionsToSaddle.size() > 1) //If the number of connected regions to this saddle is greater than 1
        {
            int contador = 0;
            for (size_t mm = 0; mm < closestRegionsToSaddle.size(); mm++)
            {
                saddlesConnectivity[k][contador] = closestRegionsToSaddle[mm];
                ++contador;
            }
            
        }
        
        if (closestRegionsToSaddle.size() == 1)
        {
            regionsWithSaddleInside.insert(std::pair<int,int> (k,closestRegionsToSaddle[0]));
        }

    }
    
    if (DEBUG)
    {
        for (size_t i = 0; i < saddlesDataSet->GetNumberOfPoints(); i++){
            logger::mainlog << "Saddle " << i << " is connected to segments : " ;
            for (size_t j = 0; j < nmaxneighbors; j++){
                if (saddlesConnectivity[i][j] != -1){
                    logger::mainlog << saddlesConnectivity[i][j] << ", ";
                }
            }
            
            // if saddle is totally inside a particular segment, output that region
            for (auto ip : regionsWithSaddleInside){
                if (ip.first == i) logger::mainlog << ip.second;
            }
            logger::mainlog << "\n" << flush;
        }
    }
    
    for (auto i : ascendingManifoldIDList) //For each of the void segments
    {
        int currentRegion = i;
        
        //Current Region of the Ascending Segmentation
        vtkSmartPointer<vtkThresholdPoints> sectionID = vtkSmartPointer<vtkThresholdPoints>::New();
        sectionID->SetInputConnection(voidSegmentation->GetOutputPort());
        sectionID->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"AscendingManifold");
        sectionID->ThresholdBetween(currentRegion,currentRegion);
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
 
                // If element was found
                if (it != saddlesConnectivity[k].end())
                {
                    double currentSaddleCoords[3]; //Coordinates of the current saddle
                    saddlesDataSet->GetPoint(k,currentSaddleCoords);
                    int closestRegionPoint = sectionIDDataset->FindPoint(currentSaddleCoords); //ID of the closest point to the saddle from the region
                    regionsSaddlesID.push_back(closestRegionPoint);

                    ++numberOfConnections;
                    for (size_t c = 0; c < nmaxneighbors; c++){
                        if ((saddlesConnectivity[k][c] != currentRegion) && (saddlesConnectivity[k][c] != -1))
                            connectedSegments.insert(saddlesConnectivity[k][c]);
                    }

                }

            }
                
        }
        
        //Check the number of Connections of each segment
        if (numberOfConnections == 0)
        {
            bool found = false;
            for (auto it : regionsWithSaddleInside){
                if (it.second == currentRegion) found = true;
            }
            if (found)
            {
                ++numberOfConnections;
                connectedSegments.insert(currentRegion);
            }
            
        }
            

        //Array corresponding to the scalar values of the Region
        auto scalarValues = sectionIDDataset->GetPointData()->GetArray(arrayName.c_str());

        double maximumValue = 0.0;
        
        int maxID; //ID of the maximum point of the region
        for (size_t j = 0; j < sectionIDDataset->GetNumberOfPoints(); j++) //For each of the points of the segment
        {
            double currentValue = scalarValues->GetVariantValue(j).ToDouble(); //Current scalar value of the point
            if (currentValue >= maximumValue) //Check if this point distance is greater than maximim
            {
                maximumValue = currentValue;
                maxID = j;
            }
        }
        bool isMaxima;
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
            if(j == maxID)
            {
                isMaxima = 1;
            }
            else
            {
                isMaxima = 0;
            }

            misDatos << currentRegion <<","<< pointCoords[0]<<","<<pointCoords[1]<<","<<pointCoords[2]<<","<<scalarValues->GetVariantValue(j).ToDouble()<< "," << maximumValue<<","<< isMaxima << "," << isSaddle <<","<< sectionIDDataset->GetNumberOfPoints()<< "," << numberOfConnections<<","<< pointCoords[0]<<","<<pointCoords[1]<<","<< pointCoords[2]<<"\n";

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




void distanceGrid::solidSegmentation(){
    
    logger::mainlog << "\n\nSegmentor: Solid Segmentation Module" << "\n" << flush;
    
    ttk::Timer SolidSegmentationTimer;
    //Writer of the .csv results file
    ofstream misDatos;
    misDatos.open((Directory+"/"+ baseFileName +"_Solid.csv").c_str());
    assert(misDatos.is_open());
    misDatos << "regionID,x,y,z,Scalar,RegionMinValue,isMinima,isSaddle,numberOfPoints,numberOfConnections,xScaled,yScaled,zScaled" << "\n";

    
    //Compute cell dimensions
    double dimensionesCelda[6];
    segmentation->GetCellBounds(0,dimensionesCelda);
    //Cell size of the current dataset
    double cellSize = dimensionesCelda[1] - dimensionesCelda[0];

    //Segmentation corresponding to the solid structure
    vtkSmartPointer<vtkThresholdPoints> solidSegmentation = vtkSmartPointer<vtkThresholdPoints>::New();
    solidSegmentation->SetInputData(segmentation);
    solidSegmentation->SetInputArrayToProcess(0,0,0,0,"This is distance grid");
    solidSegmentation->ThresholdBetween(-9e10,-1e-10);
    solidSegmentation->Update();

    auto currentSolidDataSet = vtkDataSet::SafeDownCast(solidSegmentation->GetOutputDataObject(0));
    
    // Set of ID list of the manifolds
    vtkSmartPointer<vtkDataSet> solidSegmentationDataSet = vtkDataSet::SafeDownCast(solidSegmentation->GetOutputDataObject(0));
    vtkSmartPointer<vtkAbstractArray> descendingManifoldArray = solidSegmentationDataSet->GetPointData()->GetAbstractArray("DescendingManifold");
    
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
    //Find the 1-saddles on the void structure
    vtkSmartPointer<vtkThresholdPoints> negativeSaddles = vtkSmartPointer<vtkThresholdPoints>::New();
    negativeSaddles->SetInputConnection(saddles->GetOutputPort(0));
    negativeSaddles->SetInputArrayToProcess(0,0,0,0,"This is distance grid");
    negativeSaddles->ThresholdBetween(-9e10,-1e-10);
    negativeSaddles->Update();

    //DataSet of the saddles of the Descending Segmentation of the void structure
    auto saddlesDataSet = vtkDataSet::SafeDownCast(negativeSaddles->GetOutputDataObject(0));
    logger::mainlog << "Number of accessible saddles: " << saddlesDataSet->GetNumberOfPoints() << endl;

    vector<vector<int>> saddlesConnectivity;
    size_t nmaxneighbors = 20;
    saddlesConnectivity.resize(saddlesDataSet->GetNumberOfPoints(),vector<int>(nmaxneighbors,-1.0));
    std::map<int,int> regionsWithSaddleInside;
    for (size_t k = 0; k < saddlesDataSet->GetNumberOfPoints(); k++) //For each of the saddles
    {
                
        double currentSaddleCoords[3]; //Coordinates of the current saddle
        saddlesDataSet->GetPoint(k,currentSaddleCoords);
        
        //Check that this saddle is not noise inside the region
        vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
        pointLocator->SetDataSet(currentSolidDataSet);
        pointLocator->BuildLocator();
        vtkSmartPointer<vtkIdList> closestPoints = vtkSmartPointer<vtkIdList>::New(); //IDs of the closest points to the saddle in the void structure
        //Find the in the void structure the closest points to the saddle inside a sphere of radius
        pointLocator->FindPointsWithinRadius(sqrt(2.0) * cellSize,currentSaddleCoords,closestPoints);

        vector<int> closestRegionsToSaddle; //Closest Regions ID to the saddle
        //logger::mainlog << "Current Saddle ID: " << k << endl;
        for (size_t kk = 0; kk < closestPoints->GetNumberOfIds(); kk++)
        {
            auto currentClosestRegion = currentSolidDataSet->GetPointData()->GetAbstractArray("DescendingManifold")->GetVariantValue(closestPoints->GetId(kk)).ToInt();
            //logger::mainlog << currentClosestRegion << endl;
            closestRegionsToSaddle.push_back(currentClosestRegion);
        }
        sort(closestRegionsToSaddle.begin(), closestRegionsToSaddle.end()); //Order the values of the segmentation
        vector<int>::iterator it;
        it = unique(closestRegionsToSaddle.begin(), closestRegionsToSaddle.end());  //Delete repeated values
        closestRegionsToSaddle.resize(distance(closestRegionsToSaddle.begin(),it)); //Resize with the unique values
        if (closestRegionsToSaddle.size() > 1) //If the number of connected regions to this saddle is greater than 1
        {
            int contador = 0;
            for (size_t mm = 0; mm < closestRegionsToSaddle.size(); mm++)
            {
                saddlesConnectivity[k][contador] = closestRegionsToSaddle[mm];
                ++contador;
            }
            
        }
        
        if (closestRegionsToSaddle.size() == 1)
        {
            regionsWithSaddleInside.insert(std::pair<int,int> (k,closestRegionsToSaddle[0]));
        }
        

    }

    
    
    if (DEBUG)
    {
        for (size_t i = 0; i < saddlesDataSet->GetNumberOfPoints(); i++){
            logger::mainlog << "Saddle " << i << " is connected to segments : " ;
            for (size_t j = 0; j < nmaxneighbors; j++){
                if (saddlesConnectivity[i][j] != -1){
                    logger::mainlog << saddlesConnectivity[i][j] << ", ";
                }
            }
            
            // if saddle is totally inside a particular segment, output that region
            for (auto ip : regionsWithSaddleInside){
                if (ip.first == i) logger::mainlog << ip.second;
            }
            logger::mainlog << "\n" << flush;
        }
    }
    

    

    for (auto i : descendingManifoldIDList) //For each of the void segments
    {
        int currentRegion = i;
        
        //Current Region of the Descending Segmentation
        vtkSmartPointer<vtkThresholdPoints> sectionID = vtkSmartPointer<vtkThresholdPoints>::New();
        sectionID->SetInputConnection(solidSegmentation->GetOutputPort());
        sectionID->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"DescendingManifold");
        sectionID->ThresholdBetween(currentRegion,currentRegion);
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
 
                // If element was found
                if (it != saddlesConnectivity[k].end())
                {
                    double currentSaddleCoords[3]; //Coordinates of the current saddle
                    saddlesDataSet->GetPoint(k,currentSaddleCoords);
                    int closestRegionPoint = sectionIDDataset->FindPoint(currentSaddleCoords); //ID of the closest point to the saddle from the region
                    regionsSaddlesID.push_back(closestRegionPoint);

                    ++numberOfConnections;
                    for (size_t c = 0; c < nmaxneighbors; c++){
                        if ((saddlesConnectivity[k][c] != currentRegion) && (saddlesConnectivity[k][c] != -1))
                            connectedSegments.insert(saddlesConnectivity[k][c]);
                    }

                }

            }

                
                
        }
        
        //Check the number of Connections of each segment
        if (numberOfConnections == 0)
        {
            bool found = false;
            for (auto it : regionsWithSaddleInside){
                if (it.second == currentRegion) found = true;
            }
            if (found)
            {
                ++numberOfConnections;
                connectedSegments.insert(currentRegion);
            }
            
        }

        //Array corresponding to the scalar values of the Region
        auto scalarValues = sectionIDDataset->GetPointData()->GetArray("This is distance grid");

        double minimumValue = 10e2;
        
        int minID; //ID of the maximum point of the region
        for (size_t j = 0; j < sectionIDDataset->GetNumberOfPoints(); j++) //For each of the points of the segment
        {
            bool isMaximum = false;
            double currentValue = scalarValues->GetVariantValue(j).ToDouble(); //Current scalar value of the point
            if (currentValue <= minimumValue) //Check if this point distance is smaller than minimum
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

            misDatos << currentRegion <<","<< pointCoords[0]<<","<<pointCoords[1]<<","<<pointCoords[2]<<","<<scalarValues->GetVariantValue(j).ToDouble()<< "," << minimumValue<<","<< isMinima << "," << isSaddle <<","<< sectionIDDataset->GetNumberOfPoints()<< "," << numberOfConnections<<","<< pointCoords[0]<<","<< pointCoords[1]<<","<< pointCoords[2]<<"\n";

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




void distanceGrid::accessibleVoidSpace(double moleculeRadius, bool useAllCores){
    
    logger::mainlog << "\n\nSegmentor: Accessible Void Space Module" << "\n" << flush;
    
    ttk::Timer VoidSpaceTimer;
    //Writer of the .csv results file
    ofstream segmentResults;
    segmentResults.open((Directory + "/" + baseFileName + "-accessible-void-space-mRad-" + to_string(moleculeRadius) + ".csv").c_str());
    assert(segmentResults.is_open());
    segmentResults << "regionID,Scalar,Volume,NumberOfConexions" << "\n";

    //Compute cell dimensions of the input file
    //---------------------------------------------------------------------------------------------
    double cellDimensions[6];
    segmentation->GetCellBounds(0,cellDimensions);
    //Cell size of the current dataset
    double cellSize = cellDimensions[1] - cellDimensions[0];
    //---------------------------------------------------------------------------------------------

    //Volume of each tetrahedron. As we know the volume of an unit cubic cell and each
    //cubic cell is made of 6 tetrahedrons. We set their volume to be a sixth part of the total
    
    double unitCellVolume = determinant(gridResolution)/6.0;

    //Triangulate the segmentation to improve precision
    vtkSmartPointer<vtkDataSetTriangleFilter> triangulation = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
    triangulation->SetInputData(segmentation);
    triangulation->Update();

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

    vector<vector<int>> saddlesConnectivity;
    size_t nmaxneighbors = 20;
    saddlesConnectivity.resize(saddlesDataSet->GetNumberOfPoints(),vector<int>(4,-1.0));
    std::map<int,int> regionsWithSaddleInside;
    for (size_t k = 0; k < saddlesDataSet->GetNumberOfPoints(); k++) //For each of the saddles
    {
        //logger::mainlog << "Current Saddle ID:" << endl;
        double currentSaddleCoords[3]; //Coordinates of the current saddle
        saddlesDataSet->GetPoint(k,currentSaddleCoords); //Save its coordinates
        
        //Check that this saddle is not noise inside the region
        vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
        pointLocator->SetDataSet(accessibleSpaceDataSet);
        pointLocator->BuildLocator();
        vtkSmartPointer<vtkIdList> closestPoints = vtkSmartPointer<vtkIdList>::New(); //IDs of the closest points to the saddle in the accessible void structure
        //Find  in the void structure the closest points to the saddle inside a sphere of radius
        pointLocator->FindPointsWithinRadius(sqrt(2.0)*cellSize,currentSaddleCoords,closestPoints);

       
        //Find the closest segments to each of the saddles that work as connectors between segments
        vector<int> closestRegionsToSaddle; //Closest Regions ID to the saddle
        //logger::mainlog << "Current Saddle ID: " << k << endl;
        for (size_t kk = 0; kk < closestPoints->GetNumberOfIds(); kk++)
        {
            auto currentClosestRegion = accessibleSpaceDataSet->GetPointData()->GetAbstractArray("AscendingManifold")->GetVariantValue(closestPoints->GetId(kk)).ToInt();
            //logger::mainlog << currentClosestRegion << endl;
            closestRegionsToSaddle.push_back(currentClosestRegion);
        }
        sort(closestRegionsToSaddle.begin(), closestRegionsToSaddle.end()); //Order the values of the connected segments
        vector<int>::iterator it;
        it = unique(closestRegionsToSaddle.begin(), closestRegionsToSaddle.end());  //Delete repeated values
        closestRegionsToSaddle.resize(distance(closestRegionsToSaddle.begin(),it)); //Resize with the unique values
        if (closestRegionsToSaddle.size() > 1) //If the number of connected regions to this saddle is greater than 1
        {
            //logger::mainlog << "YES" <<endl;
            int contador = 0;
            for (size_t mm = 0; mm < closestRegionsToSaddle.size(); mm++)
            {
                //logger::mainlog << closestRegionsToSaddle[mm] << endl;

                saddlesConnectivity[k][contador] = closestRegionsToSaddle[mm];
                ++contador;
            }
            
        }
        if (closestRegionsToSaddle.size() == 1)
        {
            regionsWithSaddleInside.insert(std::pair<int,int> (k,closestRegionsToSaddle[0]));
        }
        

    }
    
    // Print the saddleconnectivity for debugging
    if (DEBUG)
    {
        for (size_t i = 0; i < saddlesDataSet->GetNumberOfPoints(); i++){
            logger::mainlog << "Saddle " << i << " is connected to segments : " ;
            for (size_t j = 0; j < nmaxneighbors; j++){
                if (saddlesConnectivity[i][j] != -1){
                    logger::mainlog << saddlesConnectivity[i][j] << ", ";
                }
            }
            // if saddle is totally inside a particular segment, output that region
            for (auto ip : regionsWithSaddleInside){
                if (ip.first == i) logger::mainlog << ip.second;
            }
            logger::mainlog << "\n" << flush;
        }
    }

    for (size_t i = 0; i < segmentsID->GetNumberOfValues(); i++) //For each of the void segments
    {
        int currentRegion = segmentsID->GetVariantValue(i).ToInt();
        //logger::mainlog << "Current Region: " <<  currentRegion << endl;
        
        //Current Region of the Descending Segmentation
        vtkSmartPointer<vtkThreshold> segment = vtkSmartPointer<vtkThreshold>::New();
        segment->SetInputConnection(voidSegmentation->GetOutputPort());
        segment->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"AscendingManifold");
        segment->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
        segment->SetLowerThreshold(currentRegion);
        segment->SetUpperThreshold(currentRegion);
        segment->Update();

        //DataSet of the specific region of the Descending Segmentation
        auto segmentDataset = vtkDataSet::SafeDownCast(segment->GetOutputDataObject(0));
        
        
        int segmentNumberOfCells = segmentDataset->GetNumberOfCells();
        //logger::mainlog << segmentNumberOfCells << endl;
        //---------------------------------------------------------------------------
        int numberOfConnections = 0; //Number of connections of the current region
        std::set <int> connectedSegments;
        if(segmentDataset->GetNumberOfPoints() > 0) //If not an empty region
        {
        
            for (size_t k = 0; k < saddlesDataSet->GetNumberOfPoints(); k++) //For each of the saddles
            {
                //Check if the segment is connected to a saddle
                auto it = find(saddlesConnectivity[k].begin(), saddlesConnectivity[k].end(), currentRegion);
 
                // If element was found
                if (it != saddlesConnectivity[k].end())
                {
                    ++numberOfConnections;
                    for (size_t c = 0; c < nmaxneighbors; c++){
                        if ((saddlesConnectivity[k][c] != currentRegion) && (saddlesConnectivity[k][c] != -1))
                            connectedSegments.insert(saddlesConnectivity[k][c]);
                    }
                }

            }
        }
        //logger::mainlog << numberOfConnections << endl;

        //Check the number of Connections of each segment
        if (numberOfConnections == 0)
        {
            bool found = false;
            for (auto it : regionsWithSaddleInside){
                if (it.second == currentRegion) found = true;
            }
            if (found)
            {
                ++numberOfConnections;
                connectedSegments.insert(currentRegion);
            }
            
        }

        //logger::mainlog << numberOfConnections << endl;

        
        //Array corresponding to the scalar values of the Region
        auto scalarValues = segmentDataset->GetPointData()->GetArray("This is distance grid");

        //Write the output file
        for (size_t j = 0; j < segmentDataset->GetNumberOfPoints(); j++) //For each of the points of the segment
        {
            double segmentsVolume = segmentNumberOfCells * unitCellVolume;
            segmentResults << currentRegion<<","<<scalarValues->GetVariantValue(j).ToDouble()<< "," << segmentsVolume << "," << numberOfConnections <<"\n";

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




void distanceGrid::accessibleVoidGraph(double moleculeRadius, bool useAllCores){
    
    /* We note that a saddle lies on the boundary of two or more segments and the maxima
     lies inside the segment. So we identify the saddle that connects the two segments
     and simply join it with the maxima of that segment */
    
    logger::mainlog << "\n\nSegmentor: Accessible Void Graph Module" << "\n" << flush;
    
    ttk::Timer VoidGraphTimer;
    
    //Compute cell dimensions of the input file
    //---------------------------------------------------------------------------------------------
    double cellDimensions[6];
    segmentation->GetCellBounds(0,cellDimensions);
    //Cell size of the current dataset
    double cellSize = cellDimensions[1] - cellDimensions[0];
    //---------------------------------------------------------------------------------------------
    
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

    vector<vector<int>> saddlesConnectivity;
    size_t nmaxneighbors = 20;
    saddlesConnectivity.resize(saddlesDataSet->GetNumberOfPoints(),vector<int>(nmaxneighbors,-1.0));
    std::map<int,int> regionsWithSaddleInside;
    for (size_t k = 0; k < saddlesDataSet->GetNumberOfPoints(); k++) //For each of the saddles
    {
        //logger::mainlog << "Current Saddle ID:" << endl;
        double currentSaddleCoords[3]; //Coordinates of the current saddle
        saddlesDataSet->GetPoint(k,currentSaddleCoords); //Save its coordinates
        
        //Check that this saddle is not noise inside the region
        vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
        pointLocator->SetDataSet(currentVoidDataSet);
        pointLocator->BuildLocator();
        vtkSmartPointer<vtkIdList> closestPoints = vtkSmartPointer<vtkIdList>::New(); //IDs of the closest points to the saddle in the accessible void structure
        //Find  in the void structure the closest points to the saddle inside a sphere of radius
        pointLocator->FindPointsWithinRadius(sqrt(2.0)*cellSize,currentSaddleCoords,closestPoints);

       
        //Find the closest segments to each of the saddles that work as connectors between segments
        vector<int> closestRegionsToSaddle; //Closest Regions ID to the saddle
        logger::mainlog << "Current Saddle ID: " << k << endl;
        for (size_t kk = 0; kk < closestPoints->GetNumberOfIds(); kk++)
        {
            auto currentClosestRegion = currentVoidDataSet->GetPointData()->GetAbstractArray("AscendingManifold")->GetVariantValue(closestPoints->GetId(kk)).ToInt();
            logger::mainlog << currentClosestRegion << endl;
            closestRegionsToSaddle.push_back(currentClosestRegion);
        }
        sort(closestRegionsToSaddle.begin(), closestRegionsToSaddle.end()); //Order the values of the connected segments
        vector<int>::iterator it;
        it = unique(closestRegionsToSaddle.begin(), closestRegionsToSaddle.end());  //Delete repeated values
        closestRegionsToSaddle.resize(distance(closestRegionsToSaddle.begin(),it)); //Resize with the unique values
        if (closestRegionsToSaddle.size() > 1) //If the number of connected regions to this saddle is greater than 1
        {
            logger::mainlog << "YES" <<endl;
            int contador = 0;
            for (size_t mm = 0; mm < closestRegionsToSaddle.size(); mm++)
            {
                logger::mainlog << closestRegionsToSaddle[mm] << endl;

                saddlesConnectivity[k][contador] = closestRegionsToSaddle[mm];
                ++contador;
            }
            
        }
        if (closestRegionsToSaddle.size() == 1)
        {
            regionsWithSaddleInside.insert(std::pair<int,int> (k,closestRegionsToSaddle[0]));
        }
        

    }
    
    for (auto ip: regionsWithSaddleInside) {
        
        logger::mainlog << "Saddle ID for regions with saddle inside = " << ip.first << " and segment ID is " << ip.second << endl;
    }
    
    // Print the saddleconnectivity for debugging
    if (DEBUG)
    {
        for (size_t i = 0; i < saddlesDataSet->GetNumberOfPoints(); i++){
            logger::mainlog << "Saddle " << i << " is connected to segments : " ;
            for (size_t j = 0; j < nmaxneighbors; j++){
                if (saddlesConnectivity[i][j] != -1){
                    logger::mainlog << saddlesConnectivity[i][j] << ", ";
                }
            }
            // if saddle is totally inside a particular segment, output that region
            for (auto ip : regionsWithSaddleInside){
                if (ip.first == i) logger::mainlog << ip.second;
            }
            logger::mainlog << "\n" << flush;
        }
    }

    // Next we find the segment that each maxima belongs to:
    std::map<int,int> segmentIDofMaxima;
    for (size_t k = 0; k < maximaDataSet->GetNumberOfPoints(); k++){
        
        double currentSaddleCoords[3]; //Coordinates of the current saddle
        maximaDataSet->GetPoint(k,currentSaddleCoords); //Save its coordinates
        vtkIdType closestPoint = currentVoidDataSet->FindPoint(currentSaddleCoords[0],currentSaddleCoords[1],currentSaddleCoords[2]);
        int currentClosestRegion = currentVoidDataSet->GetPointData()->GetAbstractArray("AscendingManifold")->GetVariantValue(closestPoint).ToInt();
        
        segmentIDofMaxima.insert((std::pair<int,int> (k,currentClosestRegion)));
        logger::mainlog << " For maxima ID : " << k << ", closestPoint is " << closestPoint << " and segment ID is " << currentClosestRegion << endl;
        
    }
    
    for (auto im : segmentIDofMaxima){
        logger::mainlog << "For maxima ID : " << im.first << ", segment ID is " << im.second << endl;
    }
    
    for (auto i : ascendingManifoldIDList) //For each of the void segments
    {
        int currentRegion = i;
        //Current Region of the Ascending Segmentation
        vtkSmartPointer<vtkThresholdPoints> sectionID = vtkSmartPointer<vtkThresholdPoints>::New();
        sectionID->SetInputConnection(voidSegmentation->GetOutputPort());
        sectionID->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"AscendingManifold");
        sectionID->ThresholdBetween(currentRegion,currentRegion);
        sectionID->Update();

        //DataSet of the specific region of the Descending Segmentation
        auto sectionIDDataset = vtkDataSet::SafeDownCast(sectionID->GetOutputDataObject(0));
        int numberOfPoints = sectionIDDataset->GetNumberOfPoints();
        
        logger::mainlog << "For ascending manifold ID: " << i << " number of points is " << numberOfPoints << endl;
        
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
        
        for (size_t j =0; j < nmaxneighbors; j++){
            if (saddlesConnectivity[i][j] != -1){
                
                int connectedSegmentId = saddlesConnectivity[i][j];
                
                for (auto ip : segmentIDofMaxima){
                    if (ip.second == connectedSegmentId){
                        
                        int maximaID = ip.first;
                        // we join saddle ID and maxima ID
                        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                        int ipoint = nsaddles + maximaID;
                        line->GetPointIds()->SetId(0,i);
                        line->GetPointIds()->SetId(1,ipoint);
                        logger::mainlog << "Connecting saddle << " << i << " and maxima " << ipoint << endl;
                        graph->InsertNextCell(line->GetCellType(), line->GetPointIds());
                    }
                }
            }
        }
        
        
        for (auto ip : regionsWithSaddleInside){
            if (ip.first == i){
                int connectedSegmentId = ip.second;
                
                for (auto im : segmentIDofMaxima){
                    if (im.second == connectedSegmentId){
                        
                        int maximaID = im.first;
                        // we join saddle ID and maxima ID
                        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                        line->GetPointIds()->SetId(0,i);
                        line->GetPointIds()->SetId(1,i+maximaID);
                        graph->InsertNextCell(line->GetCellType(), line->GetPointIds());
                    }
                }
            }
            
        }
        
    }
    
    vtkSmartPointer<vtkUnstructuredGridWriter> graphWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    graphWriter->SetInputData(graph);
    graphWriter->SetFileName((Directory+"/" + baseFileName+"_graph.vtk").c_str());
    graphWriter->Write();
    
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
                 logger::mainlog << "dp[" << i << "] = " << dp[i] << endl;
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
    
    vtkSmartPointer<vtkUnstructuredGridWriter> vizGraphWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    vizGraphWriter->SetInputData(vizGraph);
    vizGraphWriter->SetFileName((Directory+"/" + baseFileName+"_viz_graph.vtk").c_str());
    vizGraphWriter->Write();


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
    it = graph->NewCellIterator();
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




void distanceGrid::accessibleSolidGraph(double moleculeRadius, bool useAllCores){
    
    /* We note that a saddle lies on the boundary of two or more segments and the minima
     lies inside the segment. So we ientify the saddle that connects the two segments
     and simply join it with the minima of that segment */
    
    logger::mainlog << "\n\nSegmentor: Accessible Solid Graph Module" << "\n" << flush;
    
    ttk::Timer solidGraphTimer;
    
    //Compute cell dimensions of the input file
    //---------------------------------------------------------------------------------------------
    double cellDimensions[6];
    segmentation->GetCellBounds(0,cellDimensions);
    //Cell size of the current dataset
    double cellSize = cellDimensions[1] - cellDimensions[0];
    logger::mainlog << "Cell size is : " << cellSize << endl;
    //---------------------------------------------------------------------------------------------
    
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

    vector<vector<int>> saddlesConnectivity;
    size_t nmaxneighbors = 20;
    saddlesConnectivity.resize(saddlesDataSet->GetNumberOfPoints(),vector<int>(nmaxneighbors,-1.0));
    std::map<int,int> regionsWithSaddleInside;
    for (size_t k = 0; k < saddlesDataSet->GetNumberOfPoints(); k++) //For each of the saddles
    {
        //logger::mainlog << "Current Saddle ID:" << endl;
        double currentSaddleCoords[3]; //Coordinates of the current saddle
        saddlesDataSet->GetPoint(k,currentSaddleCoords); //Save its coordinates
        //Check that this saddle is not noise inside the region
        vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
        pointLocator->SetDataSet(currentSolidDataSet);
        pointLocator->BuildLocator();
        vtkSmartPointer<vtkIdList> closestPoints = vtkSmartPointer<vtkIdList>::New(); //IDs of the closest points to the saddle in the accessible void structure
        //Find  in the void structure the closest points to the saddle inside a sphere of radius
        pointLocator->FindPointsWithinRadius(sqrt(2.0)*cellSize,currentSaddleCoords,closestPoints);

       
        //Find the closest segments to each of the saddles that work as connectors between segments
        vector<int> closestRegionsToSaddle; //Closest Regions ID to the saddle
        for (size_t kk = 0; kk < closestPoints->GetNumberOfIds(); kk++)
        {
            auto currentClosestRegion = currentSolidDataSet->GetPointData()->GetAbstractArray("DescendingManifold")->GetVariantValue(closestPoints->GetId(kk)).ToInt();
            closestRegionsToSaddle.push_back(currentClosestRegion);
        }
        sort(closestRegionsToSaddle.begin(), closestRegionsToSaddle.end()); //Order the values of the connected segments
        vector<int>::iterator it;
        it = unique(closestRegionsToSaddle.begin(), closestRegionsToSaddle.end());  //Delete repeated values
        closestRegionsToSaddle.resize(distance(closestRegionsToSaddle.begin(),it)); //Resize with the unique values
        if (closestRegionsToSaddle.size() > 1) //If the number of connected regions to this saddle is greater than 1
        {
            int contador = 0;
            for (size_t mm = 0; mm < closestRegionsToSaddle.size(); mm++)
            {
                logger::mainlog << closestRegionsToSaddle[mm] << endl;

                saddlesConnectivity[k][contador] = closestRegionsToSaddle[mm];
                ++contador;
            }
            
        }
        if (closestRegionsToSaddle.size() == 1)
        {
            regionsWithSaddleInside.insert(std::pair<int,int> (k,closestRegionsToSaddle[0]));
        }
        

    }
    
    for (auto ip: regionsWithSaddleInside) {
        
        logger::mainlog << "Saddle ID for regions with saddle inside = " << ip.first << " and segment ID is " << ip.second << endl;
    }
    
    // Print the saddleconnectivity for debugging
    if (DEBUG)
    {
        for (size_t i = 0; i < saddlesDataSet->GetNumberOfPoints(); i++){
            logger::mainlog << "Saddle " << i << " is connected to segments : " ;
            for (size_t j = 0; j < nmaxneighbors; j++){
                if (saddlesConnectivity[i][j] != -1){
                    logger::mainlog << saddlesConnectivity[i][j] << ", ";
                }
            }
            // if saddle is totally inside a particular segment, output that region
            for (auto ip : regionsWithSaddleInside){
                if (ip.first == i) logger::mainlog << ip.second;
            }
            logger::mainlog << "\n" << flush;
        }
    }

    // Next we find the segment that each minima belongs to:
    std::map<int,int> segmentIDofMinima;
    for (size_t k = 0; k < minimaDataSet->GetNumberOfPoints(); k++){
        
        double currentSaddleCoords[3]; //Coordinates of the current saddle
        minimaDataSet->GetPoint(k,currentSaddleCoords); //Save its coordinates
        vtkIdType closestPoint = currentSolidDataSet->FindPoint(currentSaddleCoords[0],currentSaddleCoords[1],currentSaddleCoords[2]);
        int currentClosestRegion = currentSolidDataSet->GetPointData()->GetAbstractArray("DescendingManifold")->GetVariantValue(closestPoint).ToInt();
        
        segmentIDofMinima.insert((std::pair<int,int> (k,currentClosestRegion)));
        logger::mainlog << " For minima ID : " << k << ", closestPoint is " << closestPoint << " and segment ID is " << currentClosestRegion << endl;
        
    }
    
    for (auto im : segmentIDofMinima){
        logger::mainlog << "For minima ID : " << im.first << ", segment ID is " << im.second << endl;
    }
    
    for (auto i : descendingManifoldIDList) //For each of the void segments
    {
        int currentRegion = i;
        //Current Region of the Ascending Segmentation
        vtkSmartPointer<vtkThresholdPoints> sectionID = vtkSmartPointer<vtkThresholdPoints>::New();
        sectionID->SetInputConnection(solidSegmentation->GetOutputPort());
        sectionID->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"DescendingManifold");
        sectionID->ThresholdBetween(currentRegion,currentRegion);
        sectionID->Update();

        //DataSet of the specific region of the Descending Segmentation
        auto sectionIDDataset = vtkDataSet::SafeDownCast(sectionID->GetOutputDataObject(0));
        int numberOfPoints = sectionIDDataset->GetNumberOfPoints();
        
        logger::mainlog << "For descending manifold ID: " << i << " number of points is " << numberOfPoints << endl;
        
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
        
        for (size_t j =0; j < nmaxneighbors; j++){
            if (saddlesConnectivity[i][j] != -1){
                
                int connectedSegmentId = saddlesConnectivity[i][j];
                
                for (auto ip : segmentIDofMinima){
                    if (ip.second == connectedSegmentId){
                        
                        int minimaID = ip.first;
                        // we join saddle ID and maxima ID
                        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                        int ipoint = nsaddles + minimaID;
                        line->GetPointIds()->SetId(0,i);
                        line->GetPointIds()->SetId(1,ipoint);
                        logger::mainlog << "Connecting saddle << " << i << " and minima " << (ipoint-nsaddles) << endl;
                        graph->InsertNextCell(line->GetCellType(), line->GetPointIds());
                    }
                }
            }
        }
        
        
        for (auto ip : regionsWithSaddleInside){
            if (ip.first == i){
                int connectedSegmentId = ip.second;
                
                for (auto im : segmentIDofMinima){
                    if (im.second == connectedSegmentId){
                        
                        int minimaID = im.first;
                        // we join saddle ID and maxima ID
                        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                        line->GetPointIds()->SetId(0,i);
                        line->GetPointIds()->SetId(1,i+minimaID);
                        graph->InsertNextCell(line->GetCellType(), line->GetPointIds());
                    }
                }
            }
            
        }
        
    }
    
    vtkSmartPointer<vtkUnstructuredGridWriter> graphWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    graphWriter->SetInputData(graph);
    graphWriter->SetFileName((Directory+"/" + baseFileName+"_graph.vtk").c_str());
    graphWriter->Write();
    
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
    
    vtkSmartPointer<vtkUnstructuredGridWriter> vizGraphWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    vizGraphWriter->SetInputData(vizGraph);
    vizGraphWriter->SetFileName((Directory+"/" + baseFileName+"_viz_graph.vtk").c_str());
    vizGraphWriter->Write();


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
    it = graph->NewCellIterator();
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
             logger::mainlog << " Error in accessible Solid Graph module: graph is not VTK_LINE" << endl;
             logger::errlog << " Error in accessible Solid Graph module: graph is not VTK_LINE" << endl;
         }
         
     }
    it->Delete();
    
    graphFile.close();
    logger::mainlog << "Graph is stored in the file " <<  graphFileName << endl;
    
}
