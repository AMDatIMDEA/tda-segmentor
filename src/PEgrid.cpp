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




void PEgrid::voidSegmentation(){
    
    logger::mainlog << "\n\nSegmentor: Void Segmentation Module" << "\n" << flush;
    
    ttk::Timer VoidSegmentationTimer;
    //Writer of the .csv results file
    ofstream misDatos;
    misDatos.open((Directory + "/" + baseFileName + "_Void_Segments.csv").c_str());
    assert(misDatos.is_open());
    misDatos << "regionID,x,y,z,Scalar,RegionMaxValue,isMaximum,isSaddle,numberOfPoints,numberOfConnections,xScaled,yScaled,zScaled" << "\n";
    
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
    voidSegmentation->Update();
    
    auto currentVoidDataSet = vtkDataSet::SafeDownCast(voidSegmentation->GetOutputDataObject(0));
    // Set of ID list of the manifolds
    vtkSmartPointer<vtkAbstractArray> ascendingManifoldArray = currentVoidDataSet->GetPointData()->GetAbstractArray("DescendingManifold");
    
    std::set<int> ascendingManifoldIDList;
    for (size_t i = 0; i < ascendingManifoldArray->GetNumberOfValues(); i++ ){
        
        ascendingManifoldIDList.insert(ascendingManifoldArray->GetVariantValue(i).ToInt());
        
    }
    
    //Find the 2-saddle critical points
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
            auto currentClosestRegion = currentVoidDataSet->GetPointData()->GetAbstractArray("DescendingManifold")->GetVariantValue(closestPoints->GetId(kk)).ToInt();
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




void PEgrid::accessibleVoidGraph(double moleculeRadius, bool useAllCores){
    
    
    
}




void PEgrid::accessibleSolidGraph(double moleculeRadius, bool useAllCores){

    
    
}
