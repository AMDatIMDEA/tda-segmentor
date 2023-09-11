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
    saddlesConnectivity.resize(saddlesDataSet->GetNumberOfPoints(),vector<int>(4,-1.0));
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
            for (size_t j = 0; j < 4; j++){
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
                    for (size_t c = 0; c < 4; c++){
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
    saddlesConnectivity.resize(saddlesDataSet->GetNumberOfPoints(),vector<int>(4,-1.0));
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
            for (size_t j = 0; j < 4; j++){
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
                    for (size_t c = 0; c < 4; c++){
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
            for (size_t j = 0; j < 4; j++){
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
                    for (size_t c = 0; c < 4; c++){
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
