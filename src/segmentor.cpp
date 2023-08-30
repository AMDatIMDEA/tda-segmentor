/*********************************************************************

TDA-Segmentor    A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:         Aditya Vasudevan (adityavv.iitkgp@gmail.com)
                 Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
                 Maciek Haranczyk (maciej.haranczyk@imdea.org)

                 IMDEA Materiales Institute
 
**********************************************************************/


segmentor::segmentor(const parameters &p){

    fileName = p.inputfilename;
    BaseFileName = p.basefilename;
    extensionName = p.extensionname;
    arrayName = p.arrayName; 
    // All the results will be saved in cwd/segmentor-BaseFileName.results/
    std::stringstream str;
    str << "segmentor-" << BaseFileName << ".results";
    Directory = str.str();

    std::ifstream check(Directory.c_str());
    if (!check)
    {
        std::stringstream com;
        com << "mkdir " << Directory;
        if( system( com.str().c_str() ) != 0)
            std::cout << "\n Unable to make directory";
    }
    else
    {
        std::stringstream com;
        com << "rm -f " << Directory << "/*";
        if( system( com.str().c_str() ) != 0)
            std::cout << "\n Unable to remove directory";
    }
    
    for (size_t i = 0; i < 3; i++){
        for (size_t j = 0; j < 3; j++){
            GridResolution[i][j] = 0.0;
            ucVectors[i][j] = 0.0;
            invUCVectors[i][j] = 0.0;
        }
    }

}

/**
 * @brief segmentor Destructor class
 *
 */
segmentor::~segmentor()
{
    logger::mainlog << "segmentor: Closing class" << "\n";

}



/**
 * @brief Compute a supercell based on an input unit cell
 * CAUTION: The results of this computation could have big storage size
 *
 * @param grid Input grid from a reader function of this class
 * @return auto
 */
auto segmentor::superCell(vtkSmartPointer<vtkImageData> grid){
    
    ttk::Timer superCellTimer;
    logger::mainlog << "\n\nSegmentor: Super Cell module:         " << "\n";

    vtkSmartPointer<vtkImageAppend> appendX = vtkSmartPointer<vtkImageAppend>::New();
    appendX->AddInputData(grid);
    appendX->SetAppendAxis(0);
    appendX->AddInputData(grid);
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
    grid->GetDimensions(cellDimsOriginal);
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
    
    // We reinitialize the grid points, nz, ny, nz that now belongs to the super cell
    nx = cellDimsSuperCell[0]; ny = cellDimsSuperCell[1]; nz = cellDimsSuperCell[2];
    gridPointsXYZ->Initialize();
    // Store all the locations of the grid points for a general triclinic lattice.
    double x = 0.0, y = 0.0, z = 0.0;
    for (unsigned int k = 0; k < cellDimsSuperCell[2]; k++){
        for (unsigned int j = 0; j < cellDimsSuperCell[1]; j++){
            for (unsigned int i = 0; i < cellDimsSuperCell[0]; i++){
                
              x = GridResolution[0][0]*i + GridResolution[0][1]*j + GridResolution[0][2] * k;
              y = GridResolution[1][0]*i + GridResolution[1][1]*j + GridResolution[1][2] * k;
              z = GridResolution[2][0]*i + GridResolution[2][1]*j + GridResolution[2][2] * k;
              gridPointsXYZ->InsertNextPoint(x, y, z);
          }
        }
      }

    
    vtkSmartPointer<vtkXMLImageDataWriter> imageWriter = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    imageWriter->SetInputData(appendedImage);
    imageWriter->SetFileName((Directory+"/"+BaseFileName+"_superCellgrid.vti").c_str());
    imageWriter->Write();
    
    double superCellWriteTime = superCellTimer.getElapsedTime();
    logger::mainlog << "Time taken to write Super Cell creation       : " << superCellWriteTime << endl;
    
    
    double totalTime = superCellCreationTime + superCellWriteTime;
    logger::mainlog << "Total Time in the supercell module            : " << totalTime << endl;

    return appendedImage;
}


auto segmentor::segmentsShapes(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex,int numberOfEigenFunctions, bool writeSegments,string scalar, bool useAllCores)
{
    logger::mainlog << "Segment Shapes Module" << "\n";
    logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";


    //Compute the outer box bounds of the material
    auto initialData = vtkDataSet::SafeDownCast(morseSmaleComplex->GetOutputDataObject(3));
    double materialBounds[6]; //Material outer box bounds
    initialData->GetBounds(materialBounds);

    int xlim = materialBounds[1];
    int ylim = materialBounds[3];
    int zlim = materialBounds[5];

    logger::mainlog << xlim << "," << ylim << "," << zlim << endl;

    //Segmentation corresponding to the void structure
    vtkSmartPointer<vtkThreshold> voidStructure = vtkSmartPointer<vtkThreshold>::New();
    voidStructure->SetInputConnection(morseSmaleComplex->GetOutputPort(3));
    voidStructure->SetAllScalars(1);
    voidStructure->SetInputArrayToProcess(0,0,0,0,"This is distance grid");
    voidStructure->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    voidStructure->SetLowerThreshold(-99999999);
    voidStructure->SetUpperThreshold(0.0);
    voidStructure->Update();

    //Same structure segmentation but with a Field Data added
    vtkSmartPointer<ttkExtract> voidStructureID = vtkSmartPointer<ttkExtract>::New();
    //descendingManifoldIDList->SetDebugLevel(1);
    voidStructureID->SetUseAllCores(useAllCores);
    voidStructureID->SetInputConnection(voidStructure->GetOutputPort());
    voidStructureID->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"DescendingManifold");
    voidStructureID->SetExtractionMode(3); //Array Values
    voidStructureID->SetExtractUniqueValues(true);
    voidStructureID->Update();

    auto voidDataSet = vtkDataSet::SafeDownCast(voidStructureID->GetOutputDataObject(0)); //Data Set of the void structure
    auto idList = voidDataSet->GetFieldData()->GetAbstractArray("UniqueDescendingManifold"); //Array of the unique segment IDs

    
    for (size_t i = 0; i < idList->GetNumberOfValues(); i++)
    {
        int segmentID = idList->GetVariantValue(i).ToInt();
        logger::mainlog << "Current Segment:" << segmentID << endl;
        
        //Current Region of the Descending Segmentation
        vtkSmartPointer<vtkThreshold> segment = vtkSmartPointer<vtkThreshold>::New();
        segment->SetInputConnection(voidStructure->GetOutputPort());
        segment->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"DescendingManifold");
        segment->SetAllScalars(1);
        // segment->SetLowerThreshold(segmentID);
        // segment->SetUpperThreshold(segmentID);
        // segment->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
        segment->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
        segment->SetLowerThreshold(segmentID);
        segment->SetUpperThreshold(segmentID);
        segment->Update();
        //DataSet of the Segment
        auto segmentData = vtkDataSet::SafeDownCast(segment->GetOutputDataObject(0));

        //Check for isolated regions in the segment
        vtkSmartPointer<vtkConnectivityFilter> isolatedRegions = vtkSmartPointer<vtkConnectivityFilter>::New();
        isolatedRegions->SetInputConnection(segment->GetOutputPort());
        isolatedRegions->SetExtractionModeToAllRegions();
        isolatedRegions->SetRegionIdAssignmentMode(0);
        isolatedRegions->ColorRegionsOn();
        isolatedRegions->Update();


        //Number of isolated regions in the segment
        int numberOfIsolated = isolatedRegions->GetNumberOfExtractedRegions();

        bool signal0 = false; //Segment has multiple isolated regions

        if (numberOfIsolated > 1) //If there is more than one isolated region
        {
            signal0 = true;
        }

        if (signal0 == false) //Segment with only one isolated region
        {
            if (segmentData->GetNumberOfPoints() > 200)
            {
                
                logger::mainlog << "Segment: " << segmentID << " is completed\n";
                //Triangulate the segment
                vtkSmartPointer<vtkDataSetTriangleFilter> triangulation = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
                triangulation->SetInputConnection(segment->GetOutputPort());
                triangulation->Update();

                vtkSmartPointer<ttkEigenField> eigenField = vtkSmartPointer<ttkEigenField>::New();
                eigenField->SetInputConnection(triangulation->GetOutputPort());
                eigenField->SetUseAllCores(useAllCores);
                eigenField->SetEigenNumber(numberOfEigenFunctions);
                eigenField->SetOutputFieldName("EigenFunctions");
                eigenField->Update();

                //Extract the last EigenFunction
                vtkSmartPointer<vtkArrayCalculator> extraction = vtkSmartPointer<vtkArrayCalculator>::New();
                extraction->SetInputConnection(eigenField->GetOutputPort());
                extraction->SetAttributeTypeToPointData();
                extraction->AddScalarArrayName("EigenFunctions",numberOfEigenFunctions-1);
                extraction->SetFunction("EigenFunctions");
                extraction->SetResultArrayName("maxEigenFunction");
                extraction->Update();

                if (writeSegments)
                {
                    //logger::mainlog << "Writing segment" << acceptedRegions[i] << endl;
                    
                    vtkSmartPointer<vtkUnstructuredGridWriter> regionWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
                    regionWriter->SetInputConnection(extraction->GetOutputPort());
                    regionWriter->SetFileName(("../Results/MaterialsRegions/" + BaseFileName + "_EigenField_" + to_string(segmentID) + ".vtk").c_str());
                    regionWriter->Write();
                }
            }
            
        }

        else
        {
            logger::mainlog << "Segment: " << segmentID << " is splitted in pieces\n";
            
            
            vtkSmartPointer<vtkAppendFilter> append = vtkSmartPointer<vtkAppendFilter>::New(); //Filter used to append several fragments of the region
            for (size_t j = 0; j < numberOfIsolated; j++)
            {
                vtkSmartPointer<vtkThreshold> piece = vtkSmartPointer<vtkThreshold>::New();
                piece->SetInputConnection(isolatedRegions->GetOutputPort());
                piece->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"RegionId");
                piece->SetAllScalars(1);

                piece->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
                piece->SetLowerThreshold(j);
                piece->SetUpperThreshold(j);
                piece->Update();

                auto pieceData = vtkDataSet::SafeDownCast(piece->GetOutputDataObject(0));

                double pieceBounds[6];
                pieceData->GetBounds(pieceBounds);

                double xMin = pieceBounds[0];
                double xMax = pieceBounds[1];
                double yMin = pieceBounds[2];
                double yMax = pieceBounds[3];
                double zMin = pieceBounds[4];
                double zMax = pieceBounds[5];

                int xTransform = 0;
                int yTransform = 0;
                int zTransform = 0;

                
                if (int(xMax) == xlim && int(xMin) != 0)
                {
                    xTransform -= xlim;
                }
                
                if (int(yMax) == ylim && int(yMin) != 0)
                {
                    yTransform -= ylim;
                }
                
                if (int(zMax) == zlim && int(zMin) != 0)
                {
                    zTransform -= zlim;
                }
                
                

                //Translate isolated fragments to a common place
                vtkSmartPointer<vtkTransform> aTransform = vtkSmartPointer<vtkTransform>::New();
                aTransform->Translate(xTransform,yTransform,zTransform);

                vtkSmartPointer<vtkTransformFilter> transform = vtkSmartPointer<vtkTransformFilter>::New();
                transform->SetInputConnection(piece->GetOutputPort());
                transform->SetTransform(aTransform);
                transform->Update();

                append->AddInputConnection(transform->GetOutputPort());
                append->MergePointsOn();
                append->SetTolerance(0.0);
                append->Update();



            }

            if (writeSegments)
            {
                //logger::mainlog << "Writing segment" << acceptedRegions[i] << endl;
                
                vtkSmartPointer<vtkUnstructuredGridWriter> regionWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
                regionWriter->SetInputConnection(append->GetOutputPort());
                regionWriter->SetFileName(("../Results/MaterialsRegions/" + BaseFileName + "_EigenField_" + to_string(segmentID) + ".vtk").c_str());
                regionWriter->Write();
            }
            

        }

    }
    



}

auto segmentor::segmentsShapes2(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex,int numberOfEigenFunctions, bool writeSegments,string scalar, bool useAllCores)
{
    logger::mainlog << "Segment Shapes Module" << "\n";
    logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";


    //Compute the outer box bounds of the material
    auto initialData = vtkDataSet::SafeDownCast(morseSmaleComplex->GetOutputDataObject(3));
    double materialBounds[6]; //Material outer box bounds
    initialData->GetBounds(materialBounds);

    int xlim = materialBounds[1];
    int ylim = materialBounds[3];
    int zlim = materialBounds[5];

    logger::mainlog << xlim << "," << ylim << "," << zlim << endl;

    //Segmentation corresponding to the void structure
    vtkSmartPointer<vtkThreshold> voidStructure = vtkSmartPointer<vtkThreshold>::New();
    voidStructure->SetInputConnection(morseSmaleComplex->GetOutputPort(3));
    voidStructure->SetAllScalars(1);
    voidStructure->SetInputArrayToProcess(0,0,0,0,"This is distance grid");
    voidStructure->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    voidStructure->SetLowerThreshold(-9999999999);
    voidStructure->SetUpperThreshold(0.0);
    voidStructure->Update();

    //Same structure segmentation but with a Field Data added
    vtkSmartPointer<ttkExtract> voidStructureID = vtkSmartPointer<ttkExtract>::New();
    //descendingManifoldIDList->SetDebugLevel(1);
    voidStructureID->SetUseAllCores(useAllCores);
    voidStructureID->SetInputConnection(voidStructure->GetOutputPort());
    voidStructureID->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"DescendingManifold");
    voidStructureID->SetExtractionMode(3); //Array Values
    voidStructureID->SetExtractUniqueValues(true);
    voidStructureID->Update();

    auto voidDataSet = vtkDataSet::SafeDownCast(voidStructureID->GetOutputDataObject(0)); //Data Set of the void structure
    auto idList = voidDataSet->GetFieldData()->GetAbstractArray("UniqueDescendingManifold"); //Array of the unique segment IDs

    
    for (size_t i = 0; i < idList->GetNumberOfValues(); i++)
    {
        int segmentID = idList->GetVariantValue(i).ToInt();
        logger::mainlog << "Current Segment:" << segmentID << endl;
        
        //Current Region of the Descending Segmentation
        vtkSmartPointer<vtkThreshold> segment = vtkSmartPointer<vtkThreshold>::New();
        segment->SetInputConnection(voidStructure->GetOutputPort());
        segment->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"DescendingManifold");
        segment->SetAllScalars(1);
        // segment->SetLowerThreshold(segmentID);
        // segment->SetUpperThreshold(segmentID);
        // segment->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
        segment->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
        segment->SetLowerThreshold(segmentID);
        segment->SetUpperThreshold(segmentID);
        segment->Update();
        //DataSet of the Segment
        auto segmentData = vtkDataSet::SafeDownCast(segment->GetOutputDataObject(0));

        //Check for isolated regions in the segment
        vtkSmartPointer<vtkConnectivityFilter> isolatedRegions = vtkSmartPointer<vtkConnectivityFilter>::New();
        isolatedRegions->SetInputConnection(segment->GetOutputPort());
        isolatedRegions->SetExtractionModeToAllRegions();
        isolatedRegions->SetRegionIdAssignmentMode(0);
        isolatedRegions->ColorRegionsOn();
        isolatedRegions->Update();


        //Number of isolated regions in the segment
        int numberOfIsolated = isolatedRegions->GetNumberOfExtractedRegions();

        bool signal0 = false; //Segment has multiple isolated regions
        bool signal1 = false;

        if (numberOfIsolated > 1) //If there is more than one isolated region
        {
            signal0 = true;
        }

        if (signal0 == false) //Segment with only one isolated region
        {
            if (segmentData->GetNumberOfPoints() > 200)
            {
                
                logger::mainlog << "Segment: " << segmentID << " is completed\n";
                //Triangulate the segment
                vtkSmartPointer<vtkDataSetTriangleFilter> triangulation = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
                triangulation->SetInputConnection(segment->GetOutputPort());
                triangulation->Update();

                vtkSmartPointer<ttkEigenField> eigenField = vtkSmartPointer<ttkEigenField>::New();
                eigenField->SetInputConnection(triangulation->GetOutputPort());
                eigenField->SetUseAllCores(useAllCores);
                eigenField->SetEigenNumber(numberOfEigenFunctions);
                eigenField->SetOutputFieldName("EigenFunctions");
                eigenField->Update();

                //Extract the last EigenFunction
                vtkSmartPointer<vtkArrayCalculator> extraction = vtkSmartPointer<vtkArrayCalculator>::New();
                extraction->SetInputConnection(eigenField->GetOutputPort());
                extraction->SetAttributeTypeToPointData();
                extraction->AddScalarArrayName("EigenFunctions",numberOfEigenFunctions-1);
                extraction->SetFunction("EigenFunctions");
                extraction->SetResultArrayName("maxEigenFunction");
                extraction->Update();

                if (writeSegments)
                {
                    //logger::mainlog << "Writing segment" << acceptedRegions[i] << endl;
                    
                    vtkSmartPointer<vtkUnstructuredGridWriter> regionWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
                    regionWriter->SetInputConnection(extraction->GetOutputPort());
                    regionWriter->SetFileName(("../Results/MaterialsRegions/" + BaseFileName + "_EigenField_" + to_string(segmentID) + ".vtk").c_str());
                    regionWriter->Write();
                }
            }
            
        }

        else
        {
            logger::mainlog << "Segment: " << segmentID << " is splitted in pieces\n";
            //Check how the pieces are splitted
            double segmentBounds[6];
            segmentData->GetBounds(segmentBounds);

            int segmentXLim = segmentBounds[1];
            int segmentYLim = segmentBounds[3];
            int segmentZLim = segmentBounds[5];

            logger::mainlog << "Segment Limits: " << segmentXLim << "," << segmentYLim << "," << segmentZLim << "\n";

            bool signalX = false;
            bool signalY = false;
            bool signalZ = false;

            //Possibilities to check when computing the super cell
            vector<int> xVec{ 0, 0, 0 };
            vector<int> yVec{ 0, 0, 0 };
            vector<int> zVec{ 0, 0, 0 };

            if (segmentXLim == xlim)
            {
                vector<int> xVec{ +xlim, -xlim, 0 };
                signalX = true;
            }

            if (segmentYLim == ylim)
            {
                
                vector<int> yVec{ +ylim, -ylim, 0 };
                signalY = true;
            }

            if (segmentZLim == zlim)
            {
                
                vector<int> zVec{ +zlim, -zlim, 0 };
                signalZ = true;
            }


            //Filter used to merge the new instances
            //vtkSmartPointer<vtkAppendFilter> append = vtkSmartPointer<vtkAppendFilter>::New();
            
            vtkSmartPointer<vtkMultiBlockDataGroupFilter> group = vtkSmartPointer<vtkMultiBlockDataGroupFilter>::New();

            //vtkSmartPointer<vtkGroupDataSetsFilter> group = vtkSmartPointer<vtkGroupDataSetsFilter>::New();
            
            if ((signalX==true) && (signalY==true) && (signalZ==false))
            {
                logger::mainlog << "SignalX: " << signalX << " SignalY: " << signalY << " SignalZ: " << signalZ << endl;
                
                for(auto x : xVec)
                {
                    for(auto y : yVec)
                    {
                        
                        int xTransform = 0; //Initial Coords
                        int yTransform = 0; //Initial Coords
                        int zTransform = 0; //Initial Coords

                        xTransform += x;
                        yTransform += y;

                        vtkSmartPointer<vtkTransform> aTransform = vtkSmartPointer<vtkTransform>::New();
                        aTransform->Translate(xTransform,yTransform,zTransform);

                        vtkSmartPointer<vtkTransformFilter> transform = vtkSmartPointer<vtkTransformFilter>::New();
                        transform->SetInputConnection(segment->GetOutputPort());
                        transform->SetTransform(aTransform);
                        transform->Update();

                        //group->AddInputConnection(transform->GetOutputPort());
                        // append->AddInputConnection(transform->GetOutputPort());
                        // append->MergePointsOn();
                        // append->SetTolerance(0.0);
                        // append->Update();
                        group->AddInputConnection(transform->GetOutputPort());


                    }
                }
                
            }
            if (signalX == true && signalZ == true && (signalY == false))
            {
                logger::mainlog << "SignalX: " << signalX << " SignalY: " << signalY << " SignalZ: " << signalZ << endl;
                
                for(auto x : xVec)
                {
                    for(auto z : zVec)
                    {
                        
                        int xTransform = 0; //Initial Coords
                        int yTransform = 0; //Initial Coords
                        int zTransform = 0; //Initial Coords

                        xTransform += x;
                        zTransform += z;

                        vtkSmartPointer<vtkTransform> aTransform = vtkSmartPointer<vtkTransform>::New();
                        aTransform->Translate(xTransform,yTransform,zTransform);

                        vtkSmartPointer<vtkTransformFilter> transform = vtkSmartPointer<vtkTransformFilter>::New();
                        transform->SetInputConnection(voidStructure->GetOutputPort());
                        transform->SetTransform(aTransform);
                        transform->Update();

                        //group->AddInputConnection(transform->GetOutputPort());
                        // append->AddInputConnection(transform->GetOutputPort());
                        // append->MergePointsOn();
                        // append->SetTolerance(0.0);
                        // append->Update();
                        group->AddInputConnection(transform->GetOutputPort());


                        
                    }
                }
                
            }
            if (signalY == true && signalZ == true && (signalX == false))
            {
                logger::mainlog << "SignalX: " << signalX << " SignalY: " << signalY << " SignalZ: " << signalZ << endl;
                
                for(auto z : zVec)
                {
                    for(auto y : yVec)
                    {
                        
                        int xTransform = 0; //Initial Coords
                        int yTransform = 0; //Initial Coords
                        int zTransform = 0; //Initial Coords

                        zTransform += z;
                        yTransform += y;

                        vtkSmartPointer<vtkTransform> aTransform = vtkSmartPointer<vtkTransform>::New();
                        aTransform->Translate(xTransform,yTransform,zTransform);

                        vtkSmartPointer<vtkTransformFilter> transform = vtkSmartPointer<vtkTransformFilter>::New();
                        transform->SetInputConnection(voidStructure->GetOutputPort());
                        transform->SetTransform(aTransform);
                        transform->Update();

                        //group->AddInputConnection(transform->GetOutputPort());
                        // append->AddInputConnection(transform->GetOutputPort());
                        // append->MergePointsOn();
                        // append->SetTolerance(0.0);
                        // append->Update();
                        group->AddInputConnection(transform->GetOutputPort());

                        
                    }
                }
                
            }
            if (signalX == true && signalY == true && signalZ == true)
            {
                logger::mainlog << "SignalX: " << signalX << " SignalY: " << signalY << " SignalZ: " << signalZ << endl;
                
                for(auto x : xVec)
                {
                    for(auto y : yVec)
                    {
                        for(auto z : zVec)
                        {
                            
                            
                            
                            int xTransform = 0; //Initial Coords
                            int yTransform = 0; //Initial Coords
                            int zTransform = 0; //Initial Coords

                            xTransform += x;
                            yTransform += y;
                            zTransform += z;

                            vtkSmartPointer<vtkTransform> aTransform = vtkSmartPointer<vtkTransform>::New();
                            aTransform->Translate(xTransform,yTransform,zTransform);

                            vtkSmartPointer<vtkTransformFilter> transform = vtkSmartPointer<vtkTransformFilter>::New();
                            transform->SetInputConnection(voidStructure->GetOutputPort());
                            transform->SetTransform(aTransform);
                            transform->Update();

                            //group->AddInputConnection(transform->GetOutputPort());
                            // append->AddInputConnection(transform->GetOutputPort());
                            // append->MergePointsOn();
                            // append->SetTolerance(0.0);
                            // append->Update();

                            group->AddInputConnection(transform->GetOutputPort());
                            logger::mainlog << group->GetNumberOfInputConnections(0) << endl;


                        }

                        
                    }
                }
                
            }
            group->Update();
            
            
            logger::mainlog <<group->GetNumberOfOutputPorts() << endl;



            logger::mainlog << "Merging points" << endl;
            vtkSmartPointer<vtkAppendFilter> append = vtkSmartPointer<vtkAppendFilter>::New();
            append->SetInputConnection(group->GetOutputPort());
            append->MergePointsOn();
            append->SetTolerance(0.0);
            append->Update();



            

            // logger::mainlog << group->GetNumberOfInputPorts() << endl;
            // logger::mainlog << group->GetNumberOfInputConnections() << endl;


            auto appendDataSet = vtkDataSet::SafeDownCast(append->GetOutputDataObject(0));
            logger::mainlog << appendDataSet->GetNumberOfPoints() << endl;

            
        }

    }
    
}




auto segmentor::segmentSelection(string inputFile, int numberOfEigenFunctions, bool writeOutputs,bool useAllCores)
{
    //Delete the extension from the filename
    std::string base_filename = inputFile.substr(inputFile.find_last_of("/\\") + 1);
    std::string::size_type const p(base_filename.find_last_of('.'));
    //Input File name
    std::string file_without_extension = base_filename.substr(0, p);
    string fileDirectory = "../Results/" + file_without_extension + "/" + file_without_extension;

    string filePath = fileDirectory  + "_Segmentation.vtk";

    string filePath2 = fileDirectory  + "_CriticalPoints.vtk";
    
    
    
    logger::mainlog << "Reading segmentation" << endl;
    vtkSmartPointer<vtkUnstructuredGridReader> segmentation = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    segmentation->SetFileName(filePath.c_str());
    segmentation->Update();

    auto segmentationDataSet = vtkDataSet::SafeDownCast(segmentation->GetOutputDataObject(0));
    double materialBounds[6];
    segmentationDataSet->GetBounds(materialBounds);
    

    double xlim = materialBounds[1];
    double ylim = materialBounds[3];
    double zlim = materialBounds[5];

    double centerXMax = 0.75 * xlim;
    double centerXMin = 0.25 * xlim;
    double centerYMax = 0.75 * ylim;
    double centerYMin = 0.25 * ylim;
    double centerZMax = 0.75 * zlim;
    double centerZMin = 0.25 * zlim;

    logger::mainlog << xlim << "," << ylim << "," << zlim << endl;

    logger::mainlog << centerXMin << "," << centerXMax << "|" << centerYMin << "," << centerYMax << "|" << centerZMin << "," << centerZMax << endl;


    //Segmentation corresponding to the solid structure
    vtkSmartPointer<vtkThreshold> voidStructure = vtkSmartPointer<vtkThreshold>::New();
    voidStructure->SetInputConnection(segmentation->GetOutputPort());
    voidStructure->SetAllScalars(1);
    voidStructure->SetInputArrayToProcess(0,0,0,0,"This is distance grid");
    voidStructure->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    voidStructure->SetLowerThreshold(-9e9);
    voidStructure->SetUpperThreshold(0.0);
    voidStructure->Update();

    auto voidStructureDataSet = vtkDataSet::SafeDownCast(voidStructure->GetOutputDataObject(0));

    logger::mainlog << "Reading Critical Points" << endl;

    //Critical points file
    vtkSmartPointer<vtkPolyDataReader> criticalPoints = vtkSmartPointer<vtkPolyDataReader>::New();
    criticalPoints->SetFileName(filePath2.c_str());
    criticalPoints->Update();

    auto criticalPointsDataSet = vtkDataSet::SafeDownCast(criticalPoints->GetOutputDataObject(0));

    vtkSmartPointer<vtkThreshold> minimas = vtkSmartPointer<vtkThreshold>::New();
    minimas->SetInputConnection(criticalPoints->GetOutputPort());
    minimas->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"CellDimension");
    minimas->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    minimas->SetLowerThreshold(-9e9);
    minimas->SetUpperThreshold(0.0);
    minimas->Update();

    vtkSmartPointer<vtkThreshold> minimasVoid = vtkSmartPointer<vtkThreshold>::New();
    minimasVoid->SetInputConnection(minimas->GetOutputPort());
    minimasVoid->SetInputArrayToProcess(0,0,0,0,"This is distance grid");
    minimasVoid->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    minimasVoid->SetLowerThreshold(-9e9);
    minimasVoid->SetUpperThreshold(0.0);
    minimasVoid->Update();

    auto minimasDataSet = vtkDataSet::SafeDownCast(minimasVoid->GetOutputDataObject(0));

    logger::mainlog << "Looking for the segments of interest" << endl;
    for (size_t i = 0; i < minimasDataSet->GetNumberOfPoints(); i++)
    {
        double pointCoords[3];
        minimasDataSet->GetPoint(i,pointCoords);


        if ((pointCoords[0] > centerXMin) && (pointCoords[0] < centerXMax) && (pointCoords[1] > centerYMin) && (pointCoords[1] < centerYMax) && (pointCoords[2] > centerZMin) && (pointCoords[2] < centerZMax))
        {
            logger::mainlog << pointCoords[0] << "," << pointCoords[1] << "," << pointCoords[2] << endl;

            
            int closestPoint = voidStructureDataSet->FindPoint(pointCoords);
            int segmentID = voidStructureDataSet->GetPointData()->GetAbstractArray("DescendingManifold")->GetVariantValue(closestPoint).ToInt();

            vtkSmartPointer<vtkThreshold> segment = vtkSmartPointer<vtkThreshold>::New();
            segment->SetInputConnection(voidStructure->GetOutputPort());
            segment->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"DescendingManifold");
            segment->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
            segment->SetLowerThreshold(segmentID);
            segment->SetUpperThreshold(segmentID);
            segment->Update();

            // //Find the outer surface of the void space
            // vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
            // surfaceFilter->SetInputConnection(segment->GetOutputPort());
            // surfaceFilter->SetPieceInvariant(true);
            // surfaceFilter->Update();

            vtkSmartPointer<vtkDataSetTriangleFilter> triangulation = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
            triangulation->SetInputConnection(segment->GetOutputPort());
            triangulation->Update();

            vtkSmartPointer<ttkEigenField> eigen = vtkSmartPointer<ttkEigenField>::New();
            eigen->SetInputConnection(triangulation->GetOutputPort());
            eigen->SetUseAllCores(useAllCores);
            eigen->SetEigenNumber(numberOfEigenFunctions);
            eigen->SetOutputFieldName("EigenFunctions");
            eigen->Update();

            //Extract the last EigenFunction
            vtkSmartPointer<vtkArrayCalculator> extraction = vtkSmartPointer<vtkArrayCalculator>::New();
            extraction->SetInputConnection(eigen->GetOutputPort());
            extraction->SetAttributeTypeToPointData();
            extraction->AddScalarArrayName("EigenFunctions",numberOfEigenFunctions-1);
            extraction->SetFunction("EigenFunctions");
            extraction->SetResultArrayName("maxEigenFunction");
            extraction->Update();

            if (writeOutputs)
            {
                vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
                writer->SetInputConnection(extraction->GetOutputPort());
                writer->SetFileName(("../Results/MaterialsRegions/" + file_without_extension + "_Segment_" + to_string(segmentID) +".vtk").c_str());
                writer->Write();
            }

            vtkSmartPointer<ttkPersistenceDiagram> persistenceDiagram = vtkSmartPointer<ttkPersistenceDiagram>::New();
            persistenceDiagram->SetUseAllCores(useAllCores);
            persistenceDiagram->SetInputConnection(extraction->GetOutputPort());
            persistenceDiagram->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS, "maxEigenFunction");
            persistenceDiagram->Update();

            vtkSmartPointer<vtkThreshold> criticalPairs = vtkSmartPointer<vtkThreshold>::New();
            criticalPairs->SetInputConnection(persistenceDiagram->GetOutputPort());
            criticalPairs->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"PairIdentifier");
            criticalPairs->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
            criticalPairs->SetLowerThreshold(-0.1);
            criticalPairs->SetUpperThreshold(9e9);
            criticalPairs->Update();

            ofstream comparativeData; //Stream to write the comparatives of each region
            //comparativeData.open("../Results/DiagramsGrids/" + file_without_extension + ".csv"); //Opening the writer
            comparativeData.open("../Results/PersistenceDiagrams/" + file_without_extension + "_" + to_string(segmentID) + ".csv"); //Opening the writer
            assert(comparativeData.is_open());
            comparativeData << "Birth,Death" << "\n";

            auto persistenceDiagramDataSet = vtkDataSet::SafeDownCast(criticalPairs->GetOutputDataObject(0));

            for (size_t j = 1; j < persistenceDiagramDataSet->GetNumberOfPoints(); j=j+2)
            {
                double currentPointCoords[3];
                persistenceDiagramDataSet->GetPoint(j,currentPointCoords);


                
                
                comparativeData << currentPointCoords[0] << "," << currentPointCoords[1]  << "\n";

                
                
            }
            comparativeData.close();
            
    
        }
        
    }
    




}


/**
 * @brief Get the  grid(vtkImageData) from a Gaussian Cube file(.cube)
 *
 * @param inputFilePath Input File Path of the Nanoporous Material
 * @param gridResolution Grid resolution of the Gaussian Cube file
 * @param writeGridFile Write the results to an output file(OPTIONAL)
 * @return auto Grid from the Gaussian Cube file
 */
auto segmentor::readInputFile(bool writeGridFile)
{
    logger::mainlog << "\n\nSegmentor: Reader Module" << "\n";
    logger::mainlog << "Reading " << fileName << endl;
    ttk::Timer readerTime;
    
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    
    if (extensionName == ".cube"){
        
        vtkSmartPointer<vtkGaussianCubeReader2> cubeReader = vtkSmartPointer<vtkGaussianCubeReader2>::New();
        cubeReader->SetFileName(fileName.c_str());
        cubeReader->Update();
        imageData = cubeReader->GetGridOutput();
        getGridResolutionFromCubeFile();
        vtkIdType cellDims[3];
        imageData->GetDimensions(cellDims);
        imageData->SetSpacing(1.0/(cellDims[0]-1), 1.0/(cellDims[1]-1), 1.0/(cellDims[2]-1));
        if (arrayName.empty()) getArrayNameFromCubeFile(arrayName);
        logger::mainlog << "Array that is going to be used for TDA analysis: " << arrayName << endl;
        
    } else if (extensionName == ".vti"){
        
        vtkSmartPointer<vtkXMLImageDataReader> dataReader = vtkSmartPointer<vtkXMLImageDataReader>::New();
        dataReader->SetFileName((fileName).c_str());
        dataReader->Update();
        imageData = dataReader->GetOutput();
        double imageGridRes[3];
        imageData->GetSpacing(imageGridRes);
        GridResolution[0][0] = imageGridRes[0];
        GridResolution[1][1] = imageGridRes[1];
        GridResolution[2][2] = imageGridRes[2];
        if (arrayName.empty()) {
            char * name = imageData->GetPointData()->GetAbstractArray(0)->GetName();
            arrayName = name;
        }
        logger::mainlog << "Array that is going to be used for TDA analysis: " << arrayName << endl;

    } else{
        
        logger::mainlog << "Extension type of the input file is neither .cube nor .vti" <<endl;
        logger::errlog << "Extension type of the input file neither .cube nor .vti" <<endl;
        exit(0);
    }
    
    vtkIdType cellDims[3];
    imageData->GetDimensions(cellDims);
    nx = cellDims[0]; ny = cellDims[1]; nz = cellDims[2];
    
    defineUnitCellVectors();
    
    vtkIdType numberOfCells = imageData->GetNumberOfCells();
    double unitCellVolume = determinant(GridResolution);
    volume = unitCellVolume * numberOfCells;
    
    logger::mainlog << "Grid Vector (X)             : (" << GridResolution[0][0] << ", " << GridResolution[1][0] << ", " << GridResolution[2][0] << ")\n";
    logger::mainlog << "Grid Vector (Y)             : (" << GridResolution[0][1] << ", " << GridResolution[1][1] << ", " << GridResolution[2][1] << ")\n";
    logger::mainlog << "Grid Vector (Z)             : (" << GridResolution[0][2] << ", " << GridResolution[1][2] << ", " << GridResolution[2][2] << ")\n";
    logger::mainlog << "Number of points in the grid        : (" << cellDims[0] << " X " << cellDims[1] << " X "<< cellDims[2] << ")" << endl;
    logger::mainlog << "Volume of the unit cell             : "  << volume << " (A^o)^3" << endl;
    logger::mainlog << "Unit Cell Vector a :   " << ucVectors[0][0] << "    " << ucVectors[1][0] << "    " << ucVectors[2][0] << "\n";
    logger::mainlog << "Unit Cell Vector b :   " << ucVectors[0][1] << "    " << ucVectors[1][1] << "    " << ucVectors[2][1] << "\n";
    logger::mainlog << "Unit Cell Vector c :   " << ucVectors[0][2] << "    " << ucVectors[1][2] << "    " << ucVectors[2][2] << "\n";

    
    // Store all the locations of the grid points for a general triclinic lattice.
    double x = 0.0, y = 0.0, z = 0.0;
    for (unsigned int k = 0; k < cellDims[2]; k++){
        for (unsigned int j = 0; j < cellDims[1]; j++){
            for (unsigned int i = 0; i < cellDims[0]; i++){
                
              x = GridResolution[0][0]*i + GridResolution[0][1]*j + GridResolution[0][2] * k;
              y = GridResolution[1][0]*i + GridResolution[1][1]*j + GridResolution[1][2] * k;
              z = GridResolution[2][0]*i + GridResolution[2][1]*j + GridResolution[2][2] * k;
              gridPointsXYZ->InsertNextPoint(x, y, z);
          }
        }
      }

    
    if (writeGridFile == true)
    {
        vtkNew<vtkDoubleArray> pointValues;
        pointValues->SetNumberOfComponents(1);
        pointValues->SetNumberOfTuples(nx*ny*nz);
        
        for (size_t i = 0; i < (nx*ny*nz); ++i)
        {
          pointValues->SetValue(i, imageData->GetPointData()->GetArray(arrayName.c_str())->GetVariantValue(i).ToDouble());
        }
        vtkNew<vtkStructuredGrid> structuredGrid;
        structuredGrid->SetDimensions(static_cast<int>(nx), static_cast<int>(ny),
                                      static_cast<int>(nz));
        structuredGrid->SetPoints(gridPointsXYZ);
        structuredGrid->GetPointData()->SetScalars(pointValues);

       vtkNew<vtkStructuredGridWriter> strucGridWriter;
       strucGridWriter->SetInputData(structuredGrid);
       strucGridWriter->SetFileName((Directory+"/"+BaseFileName+"_grid.vtk").c_str());
       strucGridWriter->Write();
        
    }
 
    Grid = imageData;
    double elapsedTime = readerTime.getElapsedTime();
    logger::mainlog << "Time elapsed in the reader module: " << elapsedTime << "(s)" << endl;
    return imageData;
}




void segmentor::getGridResolutionFromCubeFile() {
    
    ifstream inputFile;
    
    inputFile.open(fileName);
    int nx, ny, nz;
    if (inputFile.is_open())
    {
            
        string line;
        //Number of the line that is being readed
        int lineNumber = 0;
        while (getline(inputFile,line))
        {
            if (lineNumber == 3)
            {
                stringstream ss(line);
                ss >> nx;
                ss >> GridResolution[0][0]; ss >> GridResolution[1][0]; ss >> GridResolution[2][0];
            }
                
            if (lineNumber == 4)
            {
                stringstream ss(line);
                ss >> ny;
                ss >> GridResolution[0][1]; ss >> GridResolution[1][1]; ss >> GridResolution[2][1];
            }
                
            if (lineNumber == 5)
            {
                stringstream ss(line);
                ss >> nz;
                ss >> GridResolution[0][2]; ss >> GridResolution[1][2]; ss >> GridResolution[2][2];
            }
                
            lineNumber++;
        }
                    
    }
        
    inputFile.close();
    
}




void segmentor::defineUnitCellVectors(){
    
    // Unit vector matrix of the unit cell
    ucVectors[0][0] = GridResolution[0][0]*(nx-1); ucVectors[1][0] = GridResolution[1][0]*(nx-1); ucVectors[2][0] = GridResolution[2][0]*(nx-1);
    ucVectors[0][1] = GridResolution[0][1]*(ny-1); ucVectors[1][1] = GridResolution[1][1]*(ny-1);; ucVectors[2][1] = GridResolution[2][1]*(ny-1);
    ucVectors[0][2] = GridResolution[0][2]*(nz-1); ucVectors[1][2] = GridResolution[1][2]*(nz-1); ucVectors[2][2] = GridResolution[2][2]*(nz-1);
    

    double D =  determinant(ucVectors);
    
    if (abs(D) < 1e-10){
        logger::mainlog << "Value of determinant is really small!" << endl;
        logger::errlog << "Value of determinant is really small!" << endl;
    }
    
    double invDet = 1/D;
    // Inverse unit cell matrix from abc to xyz
    invUCVectors[0][0] = invDet*   (ucVectors[2][2]*ucVectors[1][1]-ucVectors[2][1]*ucVectors[1][2]);
    invUCVectors[0][1] = invDet*-1*(ucVectors[2][2]*ucVectors[0][1]-ucVectors[2][1]*ucVectors[0][2]);
    invUCVectors[0][2] = invDet*   (ucVectors[1][2]*ucVectors[0][1]-ucVectors[1][1]*ucVectors[0][2]);
    invUCVectors[1][0] = invDet*-1*(ucVectors[2][2]*ucVectors[1][0]-ucVectors[2][0]*ucVectors[1][2]);
    invUCVectors[1][1] = invDet*   (ucVectors[2][2]*ucVectors[0][0]-ucVectors[2][0]*ucVectors[0][2]);
    invUCVectors[1][2] = invDet*-1*(ucVectors[1][2]*ucVectors[0][0]-ucVectors[1][0]*ucVectors[0][2]);
    invUCVectors[2][0] = invDet*   (ucVectors[2][1]*ucVectors[1][0]-ucVectors[2][0]*ucVectors[1][1]);
    invUCVectors[2][1] = invDet*-1*(ucVectors[2][1]*ucVectors[0][0]-ucVectors[2][0]*ucVectors[0][1]);
    invUCVectors[2][2] = invDet*   (ucVectors[1][1]*ucVectors[0][0]-ucVectors[1][0]*ucVectors[0][1]);
    
    
}

double segmentor::determinant(double matrix[3][3]){
    
    double determinant =  matrix[0][0]*(matrix[2][2]*matrix[1][1] - matrix[2][1]*matrix[1][2])
                        - matrix[1][0]*(matrix[2][2]*matrix[0][1] - matrix[2][1]*matrix[0][2])
                        + matrix[2][0]*(matrix[1][2]*matrix[0][1] - matrix[1][1]*matrix[0][2]);
    
    return determinant;
    
}

void segmentor::getArrayNameFromCubeFile(std::string &nameOfArray){
    
    ifstream inputFile;
    inputFile.open(fileName);
    
    if (inputFile.is_open())
    {
            
        string line;
        //Number of the line that is being readed
        int lineNumber = 0;
        while (getline(inputFile,line))
        {
            if (lineNumber == 1)
            {
                nameOfArray = line;
            }
                
            lineNumber++;
        }
                    
    }
        
}




void segmentor::abcToxyz (double coordABC[3], double coordXYZ[3]){
    
    coordXYZ[0] = coordABC[0]*ucVectors[0][0]+coordABC[1]*ucVectors[0][1]+coordABC[2]*ucVectors[0][2];
    coordXYZ[1] = coordABC[1]*ucVectors[1][1]+coordABC[2]*ucVectors[1][2];
    coordXYZ[2] = coordABC[2]*ucVectors[2][2];
    
    
}
/**
 * @brief Set periodic conditions for the input grid and compute a distance field for the points
 * based on their Potential Energy Atom. Set positive distance values for the pore structure and negative
 * values for the solid structure
 * @param grid Input grid
 * @param periodicConditions  Set Periodic Boundary Conditions to true or false
 * @param computeDistanceField  Compute a distance field of the grid if need
 * @param writeFile Write an output file of the grid including a distance field(if previously computed)
 * @return Grid with periodic conditions(if true) and distance field(if done) included
 */
auto segmentor::inputPrecondition2(vtkSmartPointer<vtkImageData> grid, bool periodicConditions, bool computeDistanceField,bool writeFile)
{
    logger::mainlog << "segmentor: InputPrecondition Function 2" << "\n";
    logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
    //Periodic Boundary Conditions
    vtkSmartPointer<ttkTriangulationManager> periodGrid = vtkSmartPointer<ttkTriangulationManager>::New();;
    periodGrid->SetDebugLevel(3);
    periodGrid->SetUseAllCores(true);
    periodGrid->SetInputData(grid);
    periodGrid->SetPeriodicity(periodicConditions);
    periodGrid->Update();
    auto periodGridDataSet = vtkDataSet::SafeDownCast(periodGrid->GetOutputDataObject(0));
    if (computeDistanceField)
    {
        //Find the points corresponding to the void space
        vtkSmartPointer<vtkThreshold> structureData = vtkSmartPointer<vtkThreshold>::New();
        structureData->SetInputConnection(periodGrid->GetOutputPort());
        structureData->SetInputArrayToProcess(0,0,0,0,"potentialEnergyAtom");
        structureData->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
        structureData->SetLowerThreshold(-999999999);
        structureData->SetUpperThreshold(-0.00000001);
        structureData->Update();
        
        
        //Find the outer surface of the void space
        vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
        surfaceFilter->SetInputConnection(structureData->GetOutputPort());
        surfaceFilter->SetPieceInvariant(true);
        surfaceFilter->Update();
        auto surfaceFilterData = vtkDataSet::SafeDownCast(surfaceFilter->GetOutputDataObject(0));
        vtkSmartPointer<vtkFeatureEdges> featureEdges = vtkSmartPointer<vtkFeatureEdges>::New();
        featureEdges->SetInputConnection(surfaceFilter->GetOutputPort());
        featureEdges->SetBoundaryEdges(false);
        featureEdges->SetFeatureEdges(true);
        featureEdges->SetNonManifoldEdges(false);
        featureEdges->Update();
        
        auto featuresDataSet = vtkDataSet::SafeDownCast(featureEdges->GetOutputDataObject(0));

        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        for (size_t i = 0; i < featuresDataSet->GetNumberOfPoints(); i++)
        {
            double coordinates[3];
            featuresDataSet->GetPoint(i,coordinates);
            points->InsertPoint(i,coordinates);
            
        }
        logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
        logger::mainlog << "DISTANCE FIELD CALCULATION" << endl;
        vtkSmartPointer<vtkDoubleArray> distanceArray = vtkSmartPointer<vtkDoubleArray>::New();
        distanceArray->SetName("This is distance grid");
        vtkSmartPointer<vtkPointSet> pointSet = vtkSmartPointer<vtkPointSet>::New();
        pointSet->SetPoints(points);
        
        
        for (size_t i = 0; i < grid->GetNumberOfPoints(); i++)
        {
            double gridCoordinates[3]; //Coords of the current point of the grid
            grid->GetPoint(i,gridCoordinates);
            //double gridMinDistance = 1e6;
            
            auto closestId = pointSet->FindPoint(gridCoordinates); //Find the closest point in the surface
            double seedCoordinates[3]; // Coords of the closest point in the surface
            pointSet->GetPoint(closestId,seedCoordinates);
            double gridDistance = sqrt(pow(gridCoordinates[0]-seedCoordinates[0],2)+pow(gridCoordinates[1]-seedCoordinates[1],2)+pow(gridCoordinates[2]-seedCoordinates[2],2));

            
         
            if (grid->GetPointData()->GetArray("potentialEnergyAtom")->GetVariantValue(i).ToDouble() >= 0)
            {
                distanceArray->InsertValue(i,gridDistance);
            }
            else
            {
                distanceArray->InsertValue(i,-gridDistance);

            }
            
            
            
                

        }

        periodGridDataSet->GetPointData()->AddArray(distanceArray);


    }
    //Write output to file
    if(writeFile)
    {
        vtkSmartPointer<vtkXMLImageDataWriter> gridWriter = vtkSmartPointer<vtkXMLImageDataWriter>::New();
        //gridWriter->SetInputData(Grid);
        gridWriter->SetInputConnection(periodGrid->GetOutputPort());
        gridWriter->SetFileName((Directory+"/"+BaseFileName+"distance.vti").c_str());
        gridWriter->Write();
    }
    logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";

    return periodGrid;
}

/**
 * @brief  Set periodic conditions for the input grid. Set positive distance values for the solid structure and negative
 * values for the pore structure
 *
 * @param grid Input grid
 * @param changeValues Set positive distance values for the solid structure and negative
 * values for the pore structure
 * @param periodicConditions  Periodic Boundary Conditions
 * @return auto
 */
auto segmentor::inputPrecondition(vtkSmartPointer<vtkImageData> grid, bool changeValues,bool periodicConditions, bool useAllCores)
{
    
    
    logger::mainlog << "\n\nSegmentor: InputPrecondition Module" << "\n";
    ttk::Timer periodicTimer;
    //VTK function used to set Periodic Boundary Conditions
    vtkSmartPointer<ttkTriangulationManager> periodGrid = vtkSmartPointer<ttkTriangulationManager>::New();
    periodGrid->SetUseAllCores(useAllCores);
    periodGrid->SetInputData(grid);
    periodGrid->SetPeriodicity(periodicConditions);
    periodGrid->Update();

    double elapsedTime = periodicTimer.getElapsedTime();
    logger::mainlog << "Time elapsed in periodic condition setter module: " << elapsedTime << "\n" << flush;

    return periodGrid;
    
}


void segmentor::oneGridFileCreator(string scalarName, string inputFilePath, double persistencePercentage, bool useAllCores)
{
    logger::mainlog << "Grid Files Creation" << endl;
    logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";


    //Creating Results folder:-----------------------------------------------------------------------------
    int state = system("mkdir -p ../Results/GridFiles"); //Create a directory to save the results
    if(!state)
    {
        logger::mainlog << "Grid Files Folder created" << "\n";

    }

    int state2 = system("mkdir -p ../Results/Done"); //CDirectory to write files when a material is completed

    string directory = "../Data/TobaccoPrueba/" ; //Where the files are stored
    

    std::string base_filename = inputFilePath.substr(inputFilePath.find_last_of("-") + 1);
    std::string::size_type const p(base_filename.find_last_of('.'));
    std::string file_without_extension = base_filename.substr(0, p); //Get the material filename without extension nor directory
    

    ifstream checkFile("../Results/Done/" + file_without_extension); //Check if the file was previously computed somewhere

    if (!checkFile)
    {
        
        string gridName = "../Results/GridFiles/" + file_without_extension + ".grid";
        string filePath = directory + inputFilePath; //To read the files
        ofstream outputFile(gridName);

        vtkGaussianCubeReader2 * cubeReader = vtkGaussianCubeReader2::New();
        cubeReader->SetFileName(filePath.data()); //Set the input file
        cubeReader->Update();

        vtkImageData * imageData = vtkImageData::New();
        imageData = cubeReader->GetGridOutput();

        auto energyDataSet = imageData->GetPointData()->GetAbstractArray(scalarName.c_str());

        for (size_t i = 0; i < energyDataSet->GetNumberOfValues(); i++)
        {
            outputFile << energyDataSet->GetVariantValue(i).ToDouble() * 120.0 << endl;

        }

        outputFile.close();

        string signalFile = "../Results/Done/" + file_without_extension;
        ofstream signalR(signalFile);
        signalR << "Checking \n";
        signalR.close();
        
    }
    else
    {
        logger::mainlog << "Already computed\n";
    }
    

    




}


/**
 * @brief From a energy grid(.cube file) creates a .grid file format that could be read from the Randy Snurr module.
 * It only takes the energy values corresponding to the void space and also makes a noise simplification
 *
 * @param inputFilePath Input file path
 */
void segmentor::gridFileCreator(string scalarName, string inputFilePath, double persistencePercentage, bool useAllCores)
{
    logger::mainlog << "GridFileCreator Module" << "\n";
    logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
    //Creating folder:-----------------------------------------------------------------------------
    int state = system("mkdir -p ../Results/GridFiles"); //Create a directory to save the results
    if(!state)
    {
        logger::mainlog << "Grid Files Folder created" << "\n";

    }

    int state2 = system("mkdir -p ../Results/Done"); //CDirectory to write files when a material is completed

    //Save the names of all the files in the folder in one single file
    system("rm -f ../Results/fileList.txt"); //Delete previous one if existed
    string function1 = "ls " + inputFilePath + " >> ../Results/fileList.txt";
    
    auto crystals = system((function1.c_str())); //Create the file with all the file names

    if (!crystals)
    {
        logger::mainlog << "File with the filenames created" << endl;
    }

    string line;
    ifstream inputFile("../Results/fileList.txt");
    vector<string> inputFiles;
    if (inputFile.is_open())
    {
        while(getline(inputFile,line))
        {
            inputFiles.push_back(line); //Save all the file names in a vector
        }
        inputFile.close();
    }


    #pragma omp parallel for
    for (size_t i = 0; i < inputFiles.size(); i++) //For each of the materials
    {
        logger::mainlog << inputFiles[i] << endl;
        string input = inputFilePath + "/" + inputFiles[i]; //Materials .cube path

        std::string base_filename = inputFiles[i].substr(inputFiles[i].find_last_of("-") + 1);
        std::string::size_type const p(base_filename.find_last_of('.'));
        std::string file_without_extension = base_filename.substr(0, p); //Get the material filename without extension nor directory

        ifstream checkFile("../Results/Done/" + file_without_extension); //Check if the file was previously computed somewhere


        if (!checkFile) //If not previously computed
        {
            string mofGrid = "../Results/GridFiles/" + file_without_extension + ".grid";
            ofstream outputFile(mofGrid);


            vtkSmartPointer<vtkGaussianCubeReader2> cubeReader = vtkSmartPointer<vtkGaussianCubeReader2>::New();
            cubeReader->SetFileName(input.data()); //Set the input file
            cubeReader->Update();
            //Image data output from the Gaussian Cube file
            vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
            imageData = cubeReader->GetGridOutput();

            //Periodic Boundary Conditions
            vtkSmartPointer<ttkTriangulationManager> period = vtkSmartPointer<ttkTriangulationManager>::New();;
            //period->SetDebugLevel(3);
            period->SetUseAllCores(useAllCores);
            period->SetInputData(imageData);
            period->SetPeriodicity(true);
            period->Update();


            //Given grid dataSet
            auto inputGridDataSet = vtkDataSet::SafeDownCast(period->GetOutputDataObject(0));
            //Potential Energy DataSet of the input grid
            auto energyDataSet = inputGridDataSet->GetPointData()->GetAbstractArray(scalarName.c_str());

            double minimumEnergy = 0.0; //Minimum energy value of the material
            for (size_t j = 0; j < energyDataSet->GetNumberOfValues(); j++) //For each of the points in the energy DataSet
            {
                double currentEnergy = energyDataSet->GetVariantValue(j).ToDouble(); //Current energy value tested
                if (currentEnergy < minimumEnergy)
                {
                    minimumEnergy = currentEnergy;
                }
            }
            logger::mainlog << "Minimum Energy Value of the material: " <<  minimumEnergy << endl;
            
            //Persistence Diagram of the data
            vtkSmartPointer<ttkPersistenceDiagram> persistenceDiagram = vtkSmartPointer<ttkPersistenceDiagram>::New();
            //persistenceDiagram->SetDebugLevel(3);
            persistenceDiagram->SetUseAllCores(useAllCores);
            //persistenceDiagram->SetInputData(grid);
            persistenceDiagram->SetInputConnection(period->GetOutputPort());
            persistenceDiagram->SetInputArrayToProcess(0,0,0,0,scalarName.c_str());
            
            //We delete the persistence pairs corresponding to the graph diagonal
            vtkSmartPointer<vtkThreshold> criticalPairs = vtkSmartPointer<vtkThreshold>::New();
            criticalPairs->SetInputConnection(persistenceDiagram->GetOutputPort());
            criticalPairs->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"PairIdentifier");
            criticalPairs->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
            criticalPairs->SetLowerThreshold(-0.1);
            criticalPairs->SetUpperThreshold(9e9);
            criticalPairs->Update();
            //Persistence DataSet
            auto persistenceDataSet = vtkDataSet::SafeDownCast(criticalPairs->GetOutputDataObject(0))->GetCellData()->GetArray("Persistence");
            //Persistence maximum to calculate Thresholds
            double maximumPersistence = 0.0;
            for (size_t j = 0; j < persistenceDataSet->GetNumberOfValues(); j++)
            {
                double currentPersistenceValue = persistenceDataSet->GetVariantValue(j).ToDouble();

                if (currentPersistenceValue <= abs(minimumEnergy)) //Check that the current persistence value is less or equal than the minimum energy
                {
                    if(currentPersistenceValue > maximumPersistence)
                    {
                        maximumPersistence = currentPersistenceValue;
                    }
                }
                
                
            }

            logger::mainlog << "Maximum Persistence of the negative values: " << maximumPersistence << endl;
            
            //Persistence Threshold for simplification
            double minimumPersistence = persistencePercentage * maximumPersistence;
            //Persistence threshold for future simplifications
            vtkSmartPointer<vtkThreshold> persistentPairs = vtkSmartPointer<vtkThreshold>::New();
            persistentPairs->SetInputConnection(criticalPairs->GetOutputPort());
            persistentPairs->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "Persistence");
            // persistentPairs->SetLowerThreshold(minimumPersistence);
            // persistentPairs->SetUpperThreshold(9.0e21);
            // persistentPairs->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
            persistentPairs->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
            persistentPairs->SetLowerThreshold(minimumPersistence);
            persistentPairs->SetUpperThreshold(9.0e21);

            //Topological simplification from the persistence results
            vtkSmartPointer<ttkTopologicalSimplification> topologicalSimplification = vtkSmartPointer<ttkTopologicalSimplification>::New();
            //topologicalSimplification->SetDebugLevel(4);
            topologicalSimplification->SetUseAllCores(useAllCores);
            //topologicalSimplification->SetInputData(grid);
            topologicalSimplification->SetInputConnection(0,period->GetOutputPort());
            topologicalSimplification->SetInputArrayToProcess(0,0,0, 0,scalarName.c_str());
            topologicalSimplification->SetInputConnection(1, persistentPairs->GetOutputPort());
            topologicalSimplification->Update();

            // vtkSmartPointer<vtkThreshold> voidSpace = vtkSmartPointer<vtkThreshold>::New();
            // voidSpace->SetInputConnection(topologicalSimplification->GetOutputPort());
            // voidSpace->SetInputArrayToProcess(0,0,0,0,scalarName.c_str());
            // voidSpace->SetUpperThreshold(0.0);
            // voidSpace->Update();

            auto voidSpaceDataSet = vtkDataSet::SafeDownCast(topologicalSimplification->GetOutputDataObject(0));

            auto energyDataSet2 = voidSpaceDataSet->GetPointData()->GetAbstractArray(scalarName.c_str());

            for (size_t j = 0; j < energyDataSet2->GetNumberOfValues(); j++)
            {
                outputFile << energyDataSet2->GetVariantValue(j).ToDouble() *120.27 << "\n";

            }

            outputFile.close();


            string signalFile = "../Results/Done/" + file_without_extension;
            ofstream signalR(signalFile);
            signalR << "Checking \n";
            signalR.close();
        }
        else
        {
            logger::mainlog << "Already computed \n";
        }
        
    }

    

    logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
    
}



/**
 * @brief Check if a values is contained in a vector and finds its index if true
 *
 * @param v Vector to check
 * @param K Value to check
 * @return auto Index of the value if present
 */
auto  segmentor::getIndex(vector<int> v, int K)
{
    auto it = find(v.begin(), v.end(), K);
 
    // If element was found
    if (it != v.end())
    {
     
        // calculating the index
        // of K
        int index = it - v.begin();
        return index;
    }
    else {
        // If the element is not
        // present in the vector
        return -1;
    }
}

/**
 * @brief Computes a topological simplification based on the persistence of the scalar field
 * attached to the input grid. Besides, it computes the Morse Smale Complex Segmentation. It writes 3 files: 1) Critical points file.
 * 2) Segmentation file. 3) Separatrices file.
 * @param grid  Input grid to analyse
 * @param persistenceThreshold Persistence threshold to discard noisy maxima's and minima's.
 * @param saddlesaddleIncrement Persistence threshold increment(if needed) for the
 *  simplification of the sadde-saddle connectors
 * @param writeOutputs Write MSC results to external files
 * @param useAllCores Use all cores available to speed up computations
 * @return auto Morse Smale Complex complete field information
 */
auto segmentor::MSC(vtkSmartPointer<ttkTriangulationManager> grid,double persistenceThreshold, double saddlesaddleIncrement, bool writeOutputs, bool useAllCores)
{
    logger::mainlog << "\n\nSegmentor: Morse Smale Complex Module" << "\n" << flush;
    
    ttk::Timer MSCTimer;
    //Persistence Diagram of the data
    vtkSmartPointer<ttkPersistenceDiagram> persistenceDiagram = vtkSmartPointer<ttkPersistenceDiagram>::New();
    //persistenceDiagram->SetDebugLevel(3);
    persistenceDiagram->SetUseAllCores(useAllCores);
    persistenceDiagram->SetInputConnection(grid->GetOutputPort());
    persistenceDiagram->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    
    //We delete the persistence pairs corresponding to the graph diagonal
    vtkSmartPointer<vtkThreshold> criticalPairs = vtkSmartPointer<vtkThreshold>::New();
    criticalPairs->SetInputConnection(persistenceDiagram->GetOutputPort());
    criticalPairs->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"PairIdentifier");
    criticalPairs->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    criticalPairs->SetLowerThreshold(-0.1);
    criticalPairs->SetUpperThreshold(9e9);
    criticalPairs->Update();
    
    //Persistence DataSet
    auto persistenceDataSet = vtkDataSet::SafeDownCast(criticalPairs->GetOutputDataObject(0))->GetCellData()->GetArray("Persistence");
    
    double* persistenceRange = persistenceDataSet->GetRange();
    
    double minimumPersistence = persistenceRange[0];
    double maximumPersistence = persistenceRange[1];
    
    // If persistenceThreshold is not provided as input, then 10% of max is automatically taken.
    if (persistenceThreshold == 0.0) {
        persistenceThreshold = 0.1 * maximumPersistence;
    }
    
    logger::mainlog << "Maximum persistence = " << maximumPersistence << endl;
    logger::mainlog << "Persistence threshold = " << persistenceThreshold << endl;
    
    //Persistence threshold for future simplifications
    vtkSmartPointer<vtkThreshold> persistentPairs = vtkSmartPointer<vtkThreshold>::New();
    persistentPairs->SetInputConnection(criticalPairs->GetOutputPort());
    persistentPairs->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "Persistence");
    persistentPairs->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    persistentPairs->SetLowerThreshold(persistenceThreshold);
    persistentPairs->SetUpperThreshold(9e21);

    //Topological simplification from the persistence results
    vtkSmartPointer<ttkTopologicalSimplification> topologicalSimplification = vtkSmartPointer<ttkTopologicalSimplification>::New();
    //topologicalSimplification->SetDebugLevel(4);
    topologicalSimplification->SetUseAllCores(useAllCores);
    //topologicalSimplification->SetInputData(grid);
    topologicalSimplification->SetInputConnection(0,grid->GetOutputPort());
    topologicalSimplification->SetInputArrayToProcess(0,0,0, 0,arrayName.c_str());
    topologicalSimplification->SetInputConnection(1, persistentPairs->GetOutputPort());

    //=============================================================================================
    //=============================================================================================
    //3.3 Morse Smale Complex Computation
    //Morse Smale Complex Computation
    vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex = vtkSmartPointer<ttkMorseSmaleComplex>::New();
    morseSmaleComplex->SetUseAllCores(useAllCores);
    morseSmaleComplex->SetReturnSaddleConnectors(1);
    morseSmaleComplex->SetSaddleConnectorsPersistenceThreshold(saddlesaddleIncrement*persistenceThreshold);
    morseSmaleComplex->SetInputConnection(topologicalSimplification->GetOutputPort());
    morseSmaleComplex->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    morseSmaleComplex->SetComputeSaddleConnectors(false);
    morseSmaleComplex->SetComputeAscendingSeparatrices1(false);
    morseSmaleComplex->SetComputeAscendingSeparatrices2(false);
    morseSmaleComplex->SetComputeDescendingSeparatrices1(false);
    morseSmaleComplex->SetComputeDescendingSeparatrices2(false);
    morseSmaleComplex->Update();
    
    double timeTakenForMSC = MSCTimer.getElapsedTime();
    MSCTimer.reStart();
    logger::mainlog << "Time taken for MSC creation: " << timeTakenForMSC << "(s)" << endl;
    
    vtkIdType numberOfDescendingManifolds = getNumberOfDescendingManifolds(morseSmaleComplex);
    vtkIdType numberOfAscendingManifolds = getNumberOfAscendingManifolds(morseSmaleComplex);
    
    logger::mainlog << "Total number of descending manifolds (typically solid segments) : " << numberOfDescendingManifolds << endl;
    logger::mainlog << "Total number of ascending manifolds (typically void segments) : " << numberOfAscendingManifolds << endl;

    if (writeOutputs)
    {
        // We write the critical points for a triclinic lattice:
        auto segmentationDataSet = vtkDataSet::SafeDownCast(morseSmaleComplex->GetOutputDataObject(3));
        size_t numberOfArrays = segmentationDataSet->GetPointData()->GetNumberOfArrays();
        segmentation->SetDimensions((int) nx, (int) ny,(int) nz);
        segmentation->SetPoints(gridPointsXYZ);
        
        for (size_t i = 0; i < numberOfArrays; i++){
            vtkNew<vtkDoubleArray> pointValues;
            char * name = segmentationDataSet->GetPointData()->GetAbstractArray((int)i)->GetName();
            pointValues->SetName(name);
            pointValues->SetNumberOfComponents(1);
            pointValues->SetNumberOfTuples(nx*ny*nz);
            vtkIdType numberOfPoints = segmentationDataSet->GetNumberOfPoints();
            for (size_t j = 0; j < segmentationDataSet->GetNumberOfPoints(); j++){
                
                pointValues->SetValue(j, segmentationDataSet->GetPointData()->GetArray(name)->GetVariantValue(j).ToDouble());
                
            }
            
            segmentation->GetPointData()->AddArray(pointValues);
        }
        
        vtkNew<vtkStructuredGridWriter> segmentationWriter;
        segmentationWriter->SetInputData(segmentation);
        segmentationWriter->SetFileName((Directory+"/"+BaseFileName+"_Segmentation.vtk").c_str());
        segmentationWriter->Write();
        
        auto criticalPointsDataSet = vtkDataSet::SafeDownCast(morseSmaleComplex->GetOutputDataObject(0));
        numberOfArrays = criticalPointsDataSet->GetPointData()->GetNumberOfArrays();
        
        vtkNew<vtkPoints> cPoints;
        for (size_t i = 0; i < criticalPointsDataSet->GetNumberOfPoints(); i++){
            double coordABC[3];
            criticalPointsDataSet->GetPoint(i,coordABC);
            double xt = coordABC[0]*ucVectors[0][0]+coordABC[1]*ucVectors[0][1]+coordABC[2]*ucVectors[0][2];
            double yt = coordABC[1]*ucVectors[1][1]+coordABC[2]*ucVectors[1][2];
            double zt = coordABC[2]*ucVectors[2][2];
            
            cPoints->InsertNextPoint(xt, yt, zt);
            
        }
        
        criticalPoints->SetPoints(cPoints);
        
        for (size_t i = 0; i < numberOfArrays; i++){
            vtkNew<vtkDoubleArray> pointValues;
            char * name = criticalPointsDataSet->GetPointData()->GetAbstractArray((int)i)->GetName();
            pointValues->SetName(name);
            pointValues->SetNumberOfComponents(1);
            vtkIdType numberOfPoints = criticalPointsDataSet->GetNumberOfPoints();
            pointValues->SetNumberOfTuples(numberOfPoints);
            for (size_t j = 0; j < criticalPointsDataSet->GetNumberOfPoints(); j++){
                
                pointValues->SetValue(j, criticalPointsDataSet->GetPointData()->GetArray(name)->GetVariantValue(j).ToDouble());
                
            }
            
            criticalPoints->GetPointData()->AddArray(pointValues);
        }
        
        vtkSmartPointer<vtkPolyDataWriter> critPointWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
        critPointWriter->SetInputData(criticalPoints);
        critPointWriter->SetFileName((Directory+"/" + BaseFileName+"_CriticalPoints.vtk").c_str());
        critPointWriter->Write();

        // //Saddle connectors
        // vtkNew<vtkThreshold> saddleSeparatrices{};
        // saddleSeparatrices->SetInputConnection(morseSmaleComplex->GetOutputPort(1));
        // saddleSeparatrices->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"SeparatrixType");
        // saddleSeparatrices->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
        // saddleSeparatrices->SetLowerThreshold(1);
        // saddleSeparatrices->SetUpperThreshold(1);
        
        
        // vtkNew<vtkUnstructuredGridWriter> saddleSepWriter{};
        // saddleSepWriter->SetInputConnection(saddleSeparatrices->GetOutputPort());
        // //saddleSepWriter->SetFileName("../results/saddleSep.vtk");
        // saddleSepWriter->SetFileName((directory+"/saddleSep.vtk").c_str());
        // saddleSepWriter->Write();
        
            
        // //Ascending separatrices of the MSC
        // vtkSmartPointer<vtkThreshold> ascendingSeparatrices = vtkSmartPointer<vtkThreshold>::New();
        // ascendingSeparatrices->SetInputConnection(morseSmaleComplex->GetOutputPort(1));
        // ascendingSeparatrices->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"SeparatrixType");
        // ascendingSeparatrices->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
        // ascendingSeparatrices->SetLowerThreshold(2);
        // ascendingSeparatrices->SetUpperThreshold(2);
        
        // //Ascending separatrices file
        // vtkSmartPointer<vtkUnstructuredGridWriter> asc1Writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
        // asc1Writer->SetInputConnection(ascendingSeparatrices->GetOutputPort());
        // asc1Writer->SetFileName((Directory+"/" + BaseFileName+"_Asc1Separatrices.vtk").c_str());
        // asc1Writer->Write();

        // //Descending separatrices of the MSC
        // vtkSmartPointer<vtkThreshold> descendingSeparatrices = vtkSmartPointer<vtkThreshold>::New();
        // descendingSeparatrices->SetInputConnection(morseSmaleComplex->GetOutputPort(1));
        // descendingSeparatrices->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"SeparatrixType");
        // descendingSeparatrices->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
        // descendingSeparatrices->SetLowerThreshold(0);
        // descendingSeparatrices->SetUpperThreshold(0);
        
        // //Ascending separatrices file
        // vtkSmartPointer<vtkUnstructuredGridWriter> desc1Writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
        // desc1Writer->SetInputConnection(descendingSeparatrices->GetOutputPort());
        // desc1Writer->SetFileName((Directory+"/" + BaseFileName+"_Des1Separatrices.vtk").c_str());
        // desc1Writer->Write();
    }
    
    
    double writeTime = MSCTimer.getElapsedTime();
    logger::mainlog << "Time taken to write MSC files: " << writeTime << "(s)" << endl;
    double totalTime = timeTakenForMSC + writeTime;
    logger::mainlog << "Total time elapsed in the Morse Smale Complex module: " << totalTime << "(s)" << endl;
    
    return morseSmaleComplex;
    
}




auto segmentor::ftmtree(vtkSmartPointer<ttkTriangulationManager> grid, double persistenceThreshold, bool useAllCores)
{
    
    logger::mainlog << "\n\nSegmentor: FTM tree module" << "\n" << flush;
    
    ttk::Timer graphTimer;
    //Persistence Diagram of the data
    vtkSmartPointer<ttkPersistenceDiagram> persistenceDiagram = vtkSmartPointer<ttkPersistenceDiagram>::New();
    //persistenceDiagram->SetDebugLevel(3);
    persistenceDiagram->SetUseAllCores(useAllCores);
    persistenceDiagram->SetInputConnection(grid->GetOutputPort());
    persistenceDiagram->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    
    //We delete the persistence pairs corresponding to the graph diagonal
    vtkSmartPointer<vtkThreshold> criticalPairs = vtkSmartPointer<vtkThreshold>::New();
    criticalPairs->SetInputConnection(persistenceDiagram->GetOutputPort());
    criticalPairs->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"PairIdentifier");
    criticalPairs->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    criticalPairs->SetLowerThreshold(-0.1);
    criticalPairs->SetUpperThreshold(9e9);
    criticalPairs->Update();
    
    //Persistence DataSet
    auto persistenceDataSet = vtkDataSet::SafeDownCast(criticalPairs->GetOutputDataObject(0))->GetCellData()->GetArray("Persistence");
    
    double* persistenceRange = persistenceDataSet->GetRange();
    
    double minimumPersistence = persistenceRange[0];
    double maximumPersistence = persistenceRange[1];
    
    // If persistenceThreshold is not provided as input, then 10% of max is automatically taken.
    if (persistenceThreshold == 0.0) {
        persistenceThreshold = 0.1 * maximumPersistence;
    }
    
    logger::mainlog << "Maximum persistence = " << maximumPersistence << endl;
    logger::mainlog << "Persistence threshold = " << persistenceThreshold << endl;
    
    //Persistence threshold for future simplifications
    vtkSmartPointer<vtkThreshold> persistentPairs = vtkSmartPointer<vtkThreshold>::New();
    persistentPairs->SetInputConnection(criticalPairs->GetOutputPort());
    persistentPairs->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "Persistence");
    persistentPairs->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    persistentPairs->SetLowerThreshold(persistenceThreshold);
    persistentPairs->SetUpperThreshold(9e21);

    //Topological simplification from the persistence results
    vtkSmartPointer<ttkTopologicalSimplification> topologicalSimplification = vtkSmartPointer<ttkTopologicalSimplification>::New();
    //topologicalSimplification->SetDebugLevel(4);
    topologicalSimplification->SetUseAllCores(useAllCores);
    //topologicalSimplification->SetInputData(grid);
    topologicalSimplification->SetInputConnection(0,grid->GetOutputPort());
    topologicalSimplification->SetInputArrayToProcess(0,0,0, 0,arrayName.c_str());
    topologicalSimplification->SetInputConnection(1, persistentPairs->GetOutputPort());

    // TTK contour tree calculation
    int treetype = 0;
    if (arrayName == "This is distance grid"){
        treetype = 1;
    } else if (arrayName == "Potential Energy"){
        treetype = 0;
    }
    vtkSmartPointer<ttkFTMTree> ftmTree = vtkSmartPointer<ttkFTMTree>::New();
    //ftmTree->SetDebugLevel(3);
    ftmTree->SetUseAllCores(useAllCores);
    ftmTree->SetInputConnection(topologicalSimplification->GetOutputPort());
    ftmTree->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    ftmTree->SetTreeType(treetype);
    ftmTree->SetWithSegmentation(true);
    ftmTree->Update();

    logger::mainlog << "FTM Tree computed" << endl;


    //Critical points file
    auto nodesDataSet = vtkDataSet::SafeDownCast(ftmTree->GetOutputDataObject(0));
    size_t numberOfArrays = nodesDataSet->GetPointData()->GetNumberOfArrays();
    
    vtkNew<vtkPoints> cPoints;
    for (size_t i = 0; i < nodesDataSet->GetNumberOfPoints(); i++){
        double coordABC[3];
        nodesDataSet->GetPoint(i,coordABC);
        double xt = coordABC[0]*ucVectors[0][0]+coordABC[1]*ucVectors[0][1]+coordABC[2]*ucVectors[0][2];
        double yt = coordABC[1]*ucVectors[1][1]+coordABC[2]*ucVectors[1][2];
        double zt = coordABC[2]*ucVectors[2][2];
        
        cPoints->InsertNextPoint(xt, yt, zt);
        
    }
    
    ftmTreeNodes->SetPoints(cPoints);
    
    for (size_t i = 0; i < numberOfArrays; i++){
        vtkNew<vtkDoubleArray> pointValues;
        char * name = nodesDataSet->GetPointData()->GetAbstractArray((int)i)->GetName();
        pointValues->SetName(name);
        pointValues->SetNumberOfComponents(1);
        vtkIdType numberOfPoints = nodesDataSet->GetNumberOfPoints();
        pointValues->SetNumberOfTuples(numberOfPoints);
        for (size_t j = 0; j < nodesDataSet->GetNumberOfPoints(); j++){
            
            pointValues->SetValue(j, nodesDataSet->GetPointData()->GetArray(name)->GetVariantValue(j).ToDouble());
            
        }
        
        ftmTreeNodes->GetPointData()->AddArray(pointValues);
    }
    
    vtkSmartPointer<vtkPolyDataWriter> nodesWriterPoly = vtkSmartPointer<vtkPolyDataWriter>::New();
    nodesWriterPoly->SetInputData(ftmTreeNodes);
    nodesWriterPoly->SetFileName((Directory+"/" + BaseFileName+"_FTM_nodes_poly.vtk").c_str());
    nodesWriterPoly->Write();

    //arcs file
    auto arcsDataSet = vtkDataSet::SafeDownCast(ftmTree->GetOutputDataObject(1));
    vtkNew<vtkPoints> edgeNodes;
    for (size_t i = 0; i < arcsDataSet->GetNumberOfPoints(); i++){
        double coordABC[3];
        arcsDataSet->GetPoint(i,coordABC);
        double xt = coordABC[0]*ucVectors[0][0]+coordABC[1]*ucVectors[0][1]+coordABC[2]*ucVectors[0][2];
        double yt = coordABC[1]*ucVectors[1][1]+coordABC[2]*ucVectors[1][2];
        double zt = coordABC[2]*ucVectors[2][2];
        
        edgeNodes->InsertNextPoint(xt, yt, zt);
        
    }
    
    ftmTreeEdges->SetPoints(edgeNodes);
    
    vtkCellIterator *it = arcsDataSet->NewCellIterator();
    for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextCell())
     {
         if (it->GetCellType() == VTK_LINE)
         {
             vtkIdList *pointIds = it->GetPointIds();
             vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
             line->GetPointIds()->SetId(0,pointIds->GetId(0));
             line->GetPointIds()->SetId(1,pointIds->GetId(1));
             
             ftmTreeEdges->InsertNextCell(line->GetCellType(), line->GetPointIds());
         }
         
     }
    it->Delete();
    
    numberOfArrays = arcsDataSet->GetPointData()->GetNumberOfArrays();
    for (size_t i = 0; i < numberOfArrays; i++){
        vtkNew<vtkDoubleArray> pointValues;
        char * name = arcsDataSet->GetPointData()->GetAbstractArray((int)i)->GetName();
        pointValues->SetName(name);
        pointValues->SetNumberOfComponents(1);
        vtkIdType numberOfPoints = arcsDataSet->GetNumberOfPoints();
        pointValues->SetNumberOfTuples(numberOfPoints);
        for (size_t j = 0; j < arcsDataSet->GetNumberOfPoints(); j++){
            
            pointValues->SetValue(j, arcsDataSet->GetPointData()->GetArray(name)->GetVariantValue(j).ToDouble());
            
        }
        
        ftmTreeEdges->GetPointData()->AddArray(pointValues);
    }
    
    size_t numberOfCellScalars = arcsDataSet->GetCellData()->GetNumberOfArrays();
    
    for (size_t i = 0; i < numberOfCellScalars; i++){
        vtkNew<vtkDoubleArray> pointValues;
        char * name = arcsDataSet->GetCellData()->GetAbstractArray((int)i)->GetName();
        pointValues->SetName(name);
        pointValues->SetNumberOfComponents(1);
        vtkIdType numberOfCells = arcsDataSet->GetNumberOfCells();
        pointValues->SetNumberOfTuples(numberOfCells);
        for (size_t j = 0; j < arcsDataSet->GetNumberOfCells(); j++){
            
            pointValues->SetValue(j, arcsDataSet->GetCellData()->GetArray(name)->GetVariantValue(j).ToDouble());
            
        }
        
        ftmTreeEdges->GetCellData()->AddArray(pointValues);
    }
    
    vtkSmartPointer<vtkUnstructuredGridWriter> ftmtreeEdgesWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    ftmtreeEdgesWriter->SetInputData(ftmTreeEdges);
    ftmtreeEdgesWriter->SetFileName((Directory+"/" + BaseFileName+"_FTM_arcs_ugrid.vtk").c_str());
    ftmtreeEdgesWriter->Write();
    
    
    
    vtkSmartPointer<vtkUnstructuredGridWriter> narcsWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    narcsWriter->SetInputConnection(ftmTree->GetOutputPort(1));
    narcsWriter->SetFileName((Directory+"/" + BaseFileName+"_FTM_arcs.vtk").c_str());
    narcsWriter->Write();
    
    
    
    double totalTime =  graphTimer.getElapsedTime();
    logger::mainlog << "Total time elapsed in the ftm tree module : " << totalTime << "(s)" << endl;
    
    return ftmTree;

}



void segmentor::accessibleVoidGraph(vtkSmartPointer <ttkFTMTree> ftmTree, double moleculeRadius, bool useAllCores){
    
    logger::mainlog << "\n\nSegmentor: Accessible Void graph module" << "\n" << flush;
    logger::mainlog << "Molecule Radius : " << moleculeRadius << endl << flush;
    
    ttk::Timer accessibleGraphTimer;
    
    double lowerThreshold = 0.0, upperThreshold = 0.0;
    std::string fileNameNodes = ""; std::string fileNameEdges = "";
    if (arrayName == "This is distance grid"){
        lowerThreshold = 1.0*moleculeRadius;
        upperThreshold = 9e9;
        fileNameNodes = Directory+"/" + BaseFileName+"_accessibleVoid_rad_" + to_string(moleculeRadius) + "_FTM_nodes.vtk";
        fileNameEdges = Directory+"/" + BaseFileName+"_accessibleVoid_rad_" + to_string(moleculeRadius) + "_FTM_edges.vtk";
    } else if (arrayName == "Potential Energy"){
        upperThreshold = 0.0;
        lowerThreshold = -9e9;
        fileNameNodes = Directory+"/" + BaseFileName+"_accessibleVoid_FTM_nodes.vtk";
        fileNameEdges = Directory+"/" + BaseFileName+"_accessibleVoid_FTM_edges.vtk";
    }
    
    vtkSmartPointer<vtkThreshold> accessibleGraphABC = vtkSmartPointer<vtkThreshold>::New();
    accessibleGraphABC->SetInputConnection(ftmTree->GetOutputPort(1));
    accessibleGraphABC->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"Scalar");
    accessibleGraphABC->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    accessibleGraphABC->SetLowerThreshold(lowerThreshold);
    accessibleGraphABC->SetUpperThreshold(upperThreshold);
    accessibleGraphABC->Update();

    vtkSmartPointer<vtkThreshold> accessibleGraphXYZ = vtkSmartPointer<vtkThreshold>::New();
    accessibleGraphXYZ->SetInputData(ftmTreeEdges);
    accessibleGraphXYZ->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"Scalar");
    accessibleGraphXYZ->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    accessibleGraphXYZ->SetLowerThreshold(lowerThreshold);
    accessibleGraphXYZ->SetUpperThreshold(upperThreshold);
    accessibleGraphXYZ->Update();
    
    vtkSmartPointer<vtkUnstructuredGridWriter> graphWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    graphWriter->SetInputConnection(accessibleGraphXYZ->GetOutputPort(0));
    graphWriter->SetFileName((fileNameEdges).c_str());
    graphWriter->Write();

    
    // save the graph in .nt2 format
    ofstream graphFile;
    std::string graphFileName = Directory + "/" + BaseFileName + "-voidGraph" + ".nt2";
    graphFile.open((graphFileName).c_str());
    assert(graphFile.is_open());
    
    // Dataset for the graph
    vtkSmartPointer<vtkUnstructuredGrid> ugrid = accessibleGraphABC->GetOutput();
    // We store all the vertices
    graphFile << "Nodes: " << "\n";
    for (size_t i = 0; i < ugrid->GetNumberOfPoints(); i++){
        
        double coordABC[3];
        ugrid->GetPoint(i,coordABC);
        double xt = coordABC[0]*ucVectors[0][0]+coordABC[1]*ucVectors[0][1]+coordABC[2]*ucVectors[0][2];
        double yt = coordABC[1]*ucVectors[1][1]+coordABC[2]*ucVectors[1][2];
        double zt = coordABC[2]*ucVectors[2][2];
        
        graphFile << i << " " << xt << " " << yt << " " << zt << "\n";
    }
    
    
    vtkIdType cellDims[3];
    double spacing[3];
    Grid->GetDimensions(cellDims);
    Grid->GetSpacing(spacing);
    double boxLength[3];
    for (size_t i = 0; i < 3; i++){
        boxLength[i] = spacing[i] * (cellDims[i]-1);
    }

    
    // Next we store all the edges
    graphFile << "Edges: " << "\n";
    vtkCellIterator *it = ugrid->NewCellIterator();
    for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextCell())
     {
         if (it->GetCellType() == VTK_LINE)
         {
             vtkIdList *pointIds = it->GetPointIds();
             
             double p1[3], p2[3];
             ugrid->GetPoint(pointIds->GetId(0),p1); // coordinates of birth point
             ugrid->GetPoint(pointIds->GetId(1),p2); // coordinates of death point
             
             int periodicity[3] = {0,0,0};
             double dp[3];for (size_t i = 0; i < 3; i++){
                 dp[i] = p2[i] - p1[i];
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
    double totalTime =  accessibleGraphTimer.getElapsedTime();
    logger::mainlog << "Total time elapsed in the accessible graph module : " << totalTime << "(s)" << endl;
    
    //Create a new graph just for visualization, this shows the nodes outside the periodic box.
    vtkSmartPointer<vtkUnstructuredGrid> vizGraph = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkNew<vtkPoints> allNodes;
    vizGraph->SetPoints(allNodes);
    it = ugrid->NewCellIterator();
    for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextCell())
     {
         if (it->GetCellType() == VTK_LINE)
         {
             vtkIdList *pointIds = it->GetPointIds();
             
             double p1[3], p2[3];
             ugrid->GetPoint(pointIds->GetId(0),p1); // coordinates of birth point
             ugrid->GetPoint(pointIds->GetId(1),p2); // coordinates of death point
             
             double pxyz1[3], pxyz2[3];
             abcToxyz(p1, pxyz1);
             abcToxyz(p2, pxyz2);
             
             vtkIdType id1 = allNodes->InsertNextPoint(pxyz1);
             vtkIdType id2 = allNodes->InsertNextPoint(pxyz2);
             
             
             int periodicity[3] = {0,0,0};
             double dp[3];for (size_t i = 0; i < 3; i++){
                 dp[i] = p2[i] - p1[i];
             }
             
             for (size_t i = 0; i < 3 ; i++){
                 if ( abs(dp[i]) > 0.5 * boxLength[i] )
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
                     
                     periodicNode2[i] = p2[i] + periodicity[i];
                     periodicNode1[i] = p1[i] - periodicity[i];
                 }
                 
                 double periodicNodeXYZ1[3], periodicNodeXYZ2[3];
                 abcToxyz(periodicNode1, periodicNodeXYZ1);
                 abcToxyz(periodicNode2, periodicNodeXYZ2);
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
    vizGraphWriter->SetFileName((Directory+"/" + BaseFileName+"_FTM_arcs_ugrid_viz_graph.vtk").c_str());
    vizGraphWriter->Write();
    
    
}



void segmentor::accessibleSolidGraph(vtkSmartPointer <ttkFTMTree> ftmTree, bool useAllCores){
    
    logger::mainlog << "\n\nSegmentor: Accessible Solid graph module" << "\n" << flush;
    ttk::Timer accessibleGraphTimer;
    
    double lowerThreshold = 0.0, upperThreshold = 0.0;
    std::string fileNameNodes = ""; std::string fileNameEdges = "";
    if (arrayName == "This is distance grid"){
        lowerThreshold = -9e9;
        upperThreshold = 0.0;
        fileNameNodes = Directory+"/" + BaseFileName+"_accessibleSolid" + "_FTM_nodes.vtk";
        fileNameEdges = Directory+"/" + BaseFileName+"_accessibleSolid" + "_FTM_edges.vtk";
        
    } else {
        logger::mainlog << "This module is to be used only with the distance function!" << endl;
        logger::errlog << "This module is to be used only with the distance function!" << endl;
        std::cout << "This module is to be used only with the distance function!" << endl;
        exit(1);
    }
    
    // updating tree type to get leaves of minima-saddle tree as this corresponds
    // to the solid in the israeli structures.
    ftmTree->SetTreeType(0);
    ftmTree->Update();
    
    vtkSmartPointer<vtkThreshold> accessibleGraph = vtkSmartPointer<vtkThreshold>::New();
    accessibleGraph->SetInputConnection(ftmTree->GetOutputPort(1));
    accessibleGraph->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"Scalar");
    accessibleGraph->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    accessibleGraph->SetLowerThreshold(lowerThreshold);
    accessibleGraph->SetUpperThreshold(upperThreshold);
    accessibleGraph->Update();

    
    vtkSmartPointer<vtkUnstructuredGridWriter> graphWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    graphWriter->SetInputConnection(accessibleGraph->GetOutputPort(0));
    graphWriter->SetFileName((fileNameEdges).c_str());
    graphWriter->Write();
    
    // save the graph in .nt2 format
    ofstream graphFile;
    std::string graphFileName = Directory + "/" + BaseFileName + "-voidGraph" + ".nt2";
    graphFile.open((graphFileName).c_str());
    assert(graphFile.is_open());
    
    // Dataset for the graph
    vtkSmartPointer<vtkUnstructuredGrid> ugrid = accessibleGraph->GetOutput();
    // We store all the vertices
    graphFile << "Nodes: " << "\n";
    for (size_t i = 0; i < ugrid->GetNumberOfPoints(); i++){
        
        double coords[3];
        ugrid->GetPoint(i,coords);
        graphFile << i << " " << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
    }
    
    
    vtkIdType cellDims[3];
    double spacing[3];
    Grid->GetDimensions(cellDims);
    Grid->GetSpacing(spacing);
    double boxLength[3];
    for (size_t i = 0; i < 3; i++){
        boxLength[i] = spacing[i] * (cellDims[i]-1);
    }

    
    // Next we store all the edges
    graphFile << "Edges: " << "\n";
    vtkCellIterator *it = ugrid->NewCellIterator();
    for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextCell())
     {
         if (it->GetCellType() == VTK_LINE)
         {
             vtkIdList *pointIds = it->GetPointIds();
             
             double p1[3], p2[3];
             ugrid->GetPoint(pointIds->GetId(0),p1); // coordinates of birth point
             ugrid->GetPoint(pointIds->GetId(1),p2); // coordinates of death point
             
             int periodicity[3] = {0,0,0};
             double dp[3];
             dp[0] = p2[0] - p1[0];
             dp[1] = p2[1] - p1[0];
             dp[2] = p2[2] - p1[2];
             
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
    double totalTime =  accessibleGraphTimer.getElapsedTime();
    logger::mainlog << "Total time elapsed in the accessible graph module : " << totalTime << "(s)" << endl;
    
}



vtkIdType segmentor::getNumberOfDescendingManifolds(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex){
        
    // Output the number of descending manifolds
    vtkSmartPointer<ttkExtract> descendingManifolds = vtkSmartPointer<ttkExtract>::New();
    //accessibleSpace->SetDebugLevel(1);
    descendingManifolds->SetUseAllCores(true);
    descendingManifolds->SetInputConnection(morseSmaleComplex->GetOutputPort(3));
    descendingManifolds->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"DescendingManifold");
    descendingManifolds->SetExtractionMode(3); //Array Values
    descendingManifolds->SetExtractUniqueValues(true);
    descendingManifolds->Update();

    auto descendingManifoldsDataset = vtkDataSet::SafeDownCast(descendingManifolds->GetOutputDataObject(0));
    auto descendingManifoldsID = descendingManifoldsDataset->GetFieldData()->GetAbstractArray("UniqueDescendingManifold");

    vtkIdType numberOfDescendingManifolds = descendingManifoldsID->GetNumberOfValues();
    
    return numberOfDescendingManifolds;
    
}




vtkIdType segmentor::getNumberOfAscendingManifolds(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex){
    
    // Output the number of descending manifolds
    vtkSmartPointer<ttkExtract> ascendingManifolds = vtkSmartPointer<ttkExtract>::New();
    //accessibleSpace->SetDebugLevel(1);
    ascendingManifolds->SetUseAllCores(true);
    ascendingManifolds->SetInputConnection(morseSmaleComplex->GetOutputPort(3));
    ascendingManifolds->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"AscendingManifold");
    ascendingManifolds->SetExtractionMode(3); //Array Values
    ascendingManifolds->SetExtractUniqueValues(true);
    ascendingManifolds->Update();

    auto ascendingManifoldsDataset = vtkDataSet::SafeDownCast(ascendingManifolds->GetOutputDataObject(0));
    auto ascendingManifoldsID = ascendingManifoldsDataset->GetFieldData()->GetAbstractArray("UniqueAscendingManifold");

    vtkIdType numberOfAscendingManifolds = ascendingManifoldsID->GetNumberOfValues();
    
    return numberOfAscendingManifolds;
    
}

/**
 * @brief Get the void space of the Morse Smale Complex Segmentation computed in the MSC function
 *
 * @param morseSmaleComplex Morse Smale Complex results from the MSC function
 * @param useAllCores Use all available cores in the computer
 
 */
void segmentor::voidSegmentation(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex, bool useAllCores)
{
    logger::mainlog << "\n\nSegmentor: Void Segmentation Module" << "\n" << flush;
    
    ttk::Timer VoidSegmentationTimer;
    //Writer of the .csv results file
    ofstream misDatos;
    misDatos.open((Directory + "/" + BaseFileName + "_Void_Segments.csv").c_str());
    assert(misDatos.is_open());
    misDatos << "regionID,x,y,z,Scalar,RegionMaxValue,isMaximum,isSaddle,numberOfPoints,numberOfConnections,xScaled,yScaled,zScaled" << "\n";
    
    ofstream materialInfo;
    
    if (arrayName == "Potential Energy"){
        materialInfo.open( (Directory+ "/"+ BaseFileName + "-materialInfo.csv").c_str());
        assert(materialInfo.is_open());
        materialInfo << "regionID,AverageMinimumEnergy,NumberOfPoints,NumberOfConexions,MaxMinEnergy,MinMinEnergy" << "\n";
    }
    
    //--------------------------------------------------------------------------------------------------------
    double maxMinima=-9e10;
    double minMinima=0;
    double acumulatedEnergy = 0.0;
    double averageMinimaEnergy = 0.0;
    if (arrayName == "Potential Energy")
    {
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
            if (currentEnergy < minMinima)
            {
                minMinima = currentEnergy;
            }
            if (currentEnergy > maxMinima)
            {
                maxMinima = currentEnergy;
            }
            acumulatedEnergy += currentEnergy;
            
        }

        averageMinimaEnergy = acumulatedEnergy/minimaDataSet->GetNumberOfValues();
        
    }
    //Compute cell dimensions of the input file
    //---------------------------------------------------------------------------------------------
    double cellDimensions[6];
    segmentation->GetCellBounds(0,cellDimensions);
    //Cell size of the current dataset
    double cellSize = cellDimensions[1] - cellDimensions[0];
    CellSize = cellSize;
    //---------------------------------------------------------------------------------------------


    //Segmentation corresponding to the void structure
    vtkSmartPointer<vtkThresholdPoints> voidSegmentation = vtkSmartPointer<vtkThresholdPoints>::New();
    voidSegmentation->SetInputData(segmentation);
    voidSegmentation->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    if (arrayName == "This is distance grid"){
        voidSegmentation->ThresholdBetween(0.0,9e9);
    } else if (arrayName == "Potential Energy"){
        voidSegmentation->ThresholdBetween(-9.9e10,0.0);
    }
    voidSegmentation->Update();
    
    auto currentVoidDataSet = vtkDataSet::SafeDownCast(voidSegmentation->GetOutputDataObject(0));
    // Set of ID list of the manifolds
    std::string ManifoldForAnalysis;
    if (arrayName == "This is distance grid"){
        ManifoldForAnalysis = "AscendingManifold";
    } else if (arrayName == "Potential Energy"){
        ManifoldForAnalysis = "DescendingManifold";
    }
    vtkSmartPointer<vtkAbstractArray> ascendingManifoldArray = currentVoidDataSet->GetPointData()->GetAbstractArray(ManifoldForAnalysis.c_str());
    
    std::set<int> ascendingManifoldIDList;
    for (size_t i = 0; i < ascendingManifoldArray->GetNumberOfValues(); i++ ){
        
        ascendingManifoldIDList.insert(ascendingManifoldArray->GetVariantValue(i).ToInt());
        
    }
    
    //Find the 2-saddle critical points
    vtkSmartPointer<vtkThresholdPoints> saddles = vtkSmartPointer<vtkThresholdPoints>::New();
    saddles->SetInputData(criticalPoints);
    saddles->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"CellDimension");
    
    if (arrayName == "This is distance grid"){
        saddles->ThresholdBetween(2,2);
    } else if (arrayName == "Potential Energy"){
        saddles->ThresholdBetween(1,1);
    }
    saddles->Update();
    //Find the 2-saddles on the void structure
    vtkSmartPointer<vtkThresholdPoints> positiveSaddles = vtkSmartPointer<vtkThresholdPoints>::New();
    positiveSaddles->SetInputConnection(saddles->GetOutputPort(0));
    positiveSaddles->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    
    if (arrayName == "This is distance grid"){
        positiveSaddles->ThresholdBetween(0.0,9e9);
    } else if (arrayName == "Potential Energy"){
        positiveSaddles->ThresholdBetween(-9e10,0.0);
    }
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
        pointLocator->FindPointsWithinRadius(sqrt(2.0) * CellSize,currentSaddleCoords,closestPoints);

        vector<int> closestRegionsToSaddle; //Closest Regions ID to the saddle
        //logger::mainlog << "Current Saddle ID: " << k << endl;
        for (size_t kk = 0; kk < closestPoints->GetNumberOfIds(); kk++)
        {
            auto currentClosestRegion = currentVoidDataSet->GetPointData()->GetAbstractArray(ManifoldForAnalysis.c_str())->GetVariantValue(closestPoints->GetId(kk)).ToInt();
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
        sectionID->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,ManifoldForAnalysis.c_str());
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
        
        if (arrayName == "Potential Energy"){
            materialInfo << currentRegion << "," << averageMinimaEnergy << "," << sectionIDDataset->GetNumberOfPoints()<< "," << numberOfConnections << "," << maxMinima << "," << minMinima << "\n";
        }

        
    }
    
    misDatos.close();
    if (arrayName == "Potential energy") materialInfo.close();
    
    double elapsedTime = VoidSegmentationTimer.getElapsedTime(); 
    logger::mainlog << "Time elapsed in the void segmentation module: " << elapsedTime << "(s)\n" << flush;
    
}

/**
 * @brief Get the accessible void space for a predefined molecule's radius in a material MSC segmentation. Then, writes a file with all the segments
 * information required for future postprocessing analysis.
 *
 * @param morseSmaleComplex Input MSC results from the MSC function
 * @param moleculeRadius Invited molecule radius
 * @param useAllCores Use All available cores to speed up computations(CAUTION: If this options is set to true when the program is being computed in parallel this could lead to
 * some computational crashes)
 
 */
void segmentor::accessibleVoidSpace(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex,double moleculeRadius, bool useAllCores)
{
    logger::mainlog << "\n\nSegmentor: Accessible Void Space Module" << "\n" << flush;
    
    if(arrayName == "Potential Energy") {
        logger::mainlog << "This module is to be used only with the distance function!" << endl;
        logger::errlog << "This module is to be used only with the distance function!" << endl;
        std::cout << "This module is to be used only with the distance function!" << endl;
        exit(1);
    }
    
    ttk::Timer VoidSpaceTimer;
    //Writer of the .csv results file
    ofstream segmentResults;
    segmentResults.open((Directory + "/" + BaseFileName + "-accessible-void-space-mRad-" + to_string(moleculeRadius) + ".csv").c_str());
    assert(segmentResults.is_open());
    segmentResults << "regionID,Scalar,Volume,NumberOfConexions" << "\n";

    //Compute cell dimensions of the input file
    //---------------------------------------------------------------------------------------------
    double cellDimensions[6];
    segmentation->GetCellBounds(0,cellDimensions);
    //Cell size of the current dataset
    double cellSize = cellDimensions[1] - cellDimensions[0];
    CellSize = cellSize;
    logger::mainlog << "CellSize = " << CellSize << endl;
    //---------------------------------------------------------------------------------------------

    //Volume of each tetrahedron. As we know the volume of an unit cubic cell and each
    //cubic cell is made of 6 tetrahedrons. We set their volume to be a sixth part of the total
    
    double unitCellVolume = determinant(GridResolution)/6.0;

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
        pointLocator->FindPointsWithinRadius(sqrt(2.0)*CellSize,currentSaddleCoords,closestPoints);

       
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

/**
 * @brief Get the solid segmentation of the MSC and creates a .csv file with the segments
 *
 * @param morseSmaleComplex
 * @return auto
 */
auto segmentor::solidSegmentation(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex)
{
    logger::mainlog << "\n\nSegmentor: Solid Segmentation Module" << "\n" << flush;
    
    ttk::Timer SolidSegmentationTimer;
    //Writer of the .csv results file
    ofstream misDatos;
    misDatos.open((Directory+"/"+ BaseFileName +"_Solid.csv").c_str());
    assert(misDatos.is_open());
    misDatos << "regionID,x,y,z,Scalar,RegionMinValue,isMinima,isSaddle,numberOfPoints,numberOfConnections,xScaled,yScaled,zScaled" << "\n";

    
    //Compute cell dimensions
    double dimensionesCelda[6];
    segmentation->GetCellBounds(0,dimensionesCelda);
    //Cell size of the current dataset
    double cellSize = dimensionesCelda[1] - dimensionesCelda[0];
    CellSize = cellSize;
    logger::mainlog << "Cell Size: " << cellSize << "\n" << flush;

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
        pointLocator->FindPointsWithinRadius(sqrt(2.0) * CellSize,currentSaddleCoords,closestPoints);

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
    
}



/**
 * @brief Reconstruct the segments , computes the Laplace Eigen Field of each segment and computes
 * the Persistence Diagram of this field
 *
 * @param morseSmaleComplex Output of the Morse Smale Complex Method
 * @param numberOfEigenFunctions Number of eigen vector to compute for the Spectral Decomposition
 * @param writeSegments Write an output file with the segment reconstructed and the file attached
 * @return auto
 */
void segmentor::eigenField(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex,int numberOfEigenFunctions, bool writeSegments,string scalar,bool useAllCores)
{
    logger::mainlog << "segmentor: Eigen Field Module " << "\n";
    logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";

   


    //Compute the outer box bounds of the material
    auto provisionalData = vtkDataSet::SafeDownCast(morseSmaleComplex->GetOutputDataObject(3));
    double materialBounds[6]; //Material outer box bounds
    provisionalData->GetBounds(materialBounds);

    //Segmentation corresponding to the void structure
    vtkSmartPointer<vtkThreshold> voidSegmentation = vtkSmartPointer<vtkThreshold>::New();
    voidSegmentation->SetInputConnection(morseSmaleComplex->GetOutputPort(3));
    voidSegmentation->SetAllScalars(1);
    voidSegmentation->SetInputArrayToProcess(0,0,0,0,"This is distance grid");
    voidSegmentation->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    voidSegmentation->SetLowerThreshold(-9e9);
    voidSegmentation->SetUpperThreshold(0.0);
    voidSegmentation->Update();

    //Same structure segmentation but with a Field Data added
    vtkSmartPointer<ttkExtract> descendingManifoldIDList = vtkSmartPointer<ttkExtract>::New();
    descendingManifoldIDList->SetDebugLevel(1);
    descendingManifoldIDList->SetUseAllCores(useAllCores);
    descendingManifoldIDList->SetInputConnection(voidSegmentation->GetOutputPort());
    descendingManifoldIDList->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"DescendingManifold");
    descendingManifoldIDList->SetExtractionMode(3); //Array Values
    descendingManifoldIDList->SetExtractUniqueValues(true);
    descendingManifoldIDList->Update();

    auto currentVoidDataSet = vtkDataSet::SafeDownCast(descendingManifoldIDList->GetOutputDataObject(0));
    auto uniqueDesSegIdDataSet = currentVoidDataSet->GetFieldData()->GetAbstractArray("UniqueDescendingManifold");

    //Stores the regions ID
    vector<int> regions;
    
    vector<int> acceptedRegions;
    //Stores all the regions to compute the Eigen Field using parallel computation
    vector<vtkSmartPointer<ttkEigenField>> eigenFields;

    for (size_t i = 0; i < uniqueDesSegIdDataSet->GetNumberOfValues(); i++) //For each of the void segments
    {
        //Current Region of the Descending Segmentation
        vtkSmartPointer<vtkThreshold> sectionID = vtkSmartPointer<vtkThreshold>::New();
        sectionID->SetInputConnection(voidSegmentation->GetOutputPort());
        sectionID->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"DescendingManifold");
        sectionID->SetAllScalars(1);
        int segmentID = uniqueDesSegIdDataSet->GetVariantValue(i).ToInt();
        sectionID->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
        sectionID->SetLowerThreshold(segmentID);
        sectionID->SetUpperThreshold(segmentID);
        sectionID->Update();
        regions.push_back(uniqueDesSegIdDataSet->GetVariantValue(i).ToInt());

        //DataSet of the specific region of the Descending Segmentation
        auto sectionIDDataset = vtkDataSet::SafeDownCast(sectionID->GetOutputDataObject(0));

        //Check for isolated regions in the segment
        vtkSmartPointer<vtkConnectivityFilter> segmentConnectivity = vtkSmartPointer<vtkConnectivityFilter>::New();
        segmentConnectivity->SetInputConnection(sectionID->GetOutputPort());
        segmentConnectivity->SetExtractionModeToAllRegions();
        segmentConnectivity->SetRegionIdAssignmentMode(0);
        segmentConnectivity->ColorRegionsOn();
        segmentConnectivity->Update();

        //Number of isolated fragments of the segment
        int numberOfRegions = segmentConnectivity->GetNumberOfExtractedRegions();
        if (numberOfRegions > 1) //If there is more than one isolated region
        {
            vtkSmartPointer<vtkAppendFilter> append = vtkSmartPointer<vtkAppendFilter>::New(); //Filter used to append several fragments of the region
            
            for (auto j = 0; j < numberOfRegions; j++) //For each of the isolated regions
            {
                vtkSmartPointer<vtkThreshold> currentIsolatedRegion = vtkSmartPointer<vtkThreshold>::New();
                currentIsolatedRegion->SetInputConnection(segmentConnectivity->GetOutputPort());
                currentIsolatedRegion->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"RegionId");
                currentIsolatedRegion->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
                currentIsolatedRegion->SetLowerThreshold(j);
                currentIsolatedRegion->SetUpperThreshold(j);
                currentIsolatedRegion->Update();
                auto currentIsolatedRegionDataSet = vtkDataSet::SafeDownCast(currentIsolatedRegion->GetOutputDataObject(0));
                //logger::mainlog << "Current Isolated Number Of Points : " << currentIsolatedRegionDataSet->GetNumberOfPoints() << endl;
                double currentIsolatedBounds[6];
                currentIsolatedRegionDataSet->GetBounds(currentIsolatedBounds);
                //logger::mainlog << "Current Isolated Region Bounds:" << endl;
                //logger::mainlog << currentIsolatedBounds[0] << " " << currentIsolatedBounds[1] << " " << currentIsolatedBounds[2] << " " << currentIsolatedBounds[3] << " " << currentIsolatedBounds[4] << " " << currentIsolatedBounds[5] << endl;
                
                
                //Move the isolated regions to a common location
                int xTransform = 0;
                int yTransform = 0;
                int zTransform = 0;

                if (currentIsolatedBounds[1] == materialBounds[1])
                {
                    xTransform -= materialBounds[1];
                }
                else if((currentIsolatedBounds[1] != materialBounds[1]) && (currentIsolatedBounds[0] != materialBounds[0]))
                {
                    xTransform -= materialBounds[1];
                }

                if (currentIsolatedBounds[3] == materialBounds[3])
                {
                    yTransform -= materialBounds[3];
                }
                else if((currentIsolatedBounds[3] != materialBounds[3]) && (currentIsolatedBounds[2] != materialBounds[2]))
                {
                    yTransform -= materialBounds[3];
                }

                if (currentIsolatedBounds[5] == materialBounds[5])
                {
                    zTransform -= materialBounds[5];
                }
                else if((currentIsolatedBounds[5] != materialBounds[5]) && (currentIsolatedBounds[4] != materialBounds[4]))
                {
                    zTransform -= materialBounds[5];
                }
                
                
                //Translate isolated fragments to a common place
                vtkSmartPointer<vtkTransform> aTransform = vtkSmartPointer<vtkTransform>::New();
                aTransform->Translate(xTransform,yTransform,zTransform);

                vtkSmartPointer<vtkTransformFilter> transform = vtkSmartPointer<vtkTransformFilter>::New();
                transform->SetInputConnection(currentIsolatedRegion->GetOutputPort());
                transform->SetTransform(aTransform);
                transform->Update();

                append->AddInputConnection(transform->GetOutputPort());
                append->MergePointsOn();
                append->SetTolerance(0.0);
                append->Update();
                

            }

            //Now we have all the isolated regions merged on a single one

            //logger::mainlog << "Get append data" << endl;
            auto appendDataSet = vtkDataSet::SafeDownCast(append->GetOutputDataObject(0));
            //logger::mainlog << "Remoce POINT REGION" << endl;
            appendDataSet->GetPointData()->RemoveArray("RegionId");
            //logger::mainlog << "Remove CELL REGION" << endl;
            appendDataSet->GetCellData()->RemoveArray("RegionId");
            //logger::mainlog << "Update" << endl;
            append->Update();

            //Check that all the fragments are merged
            vtkSmartPointer<vtkConnectivityFilter> segmentConnectivity2 = vtkSmartPointer<vtkConnectivityFilter>::New();
            segmentConnectivity2->SetInputConnection(append->GetOutputPort());
            //segmentConnectivity2->SetExtractionModeToAllRegions();
            segmentConnectivity2->SetExtractionModeToLargestRegion();
            segmentConnectivity2->SetRegionIdAssignmentMode(0);
            //segmentConnectivity2->ColorRegionsOn();
            segmentConnectivity2->Update();

            auto segmentConnectivity2DataSet = vtkDataSet::SafeDownCast(segmentConnectivity2->GetOutputDataObject(0));

            //auto appendDataSet = vtkDataSet::SafeDownCast(append->GetOutputDataObject(0));

            if (segmentConnectivity2DataSet->GetNumberOfPoints() > 300)
            {
                //Find the outer surface of the void space
                // vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
                // surfaceFilter->SetInputConnection(segmentConnectivity2->GetOutputPort());
                // surfaceFilter->SetPieceInvariant(true);
                // surfaceFilter->Update();

                //Check that all the fragments are merged
                vtkSmartPointer<vtkConnectivityFilter> segmentConnectivity3 = vtkSmartPointer<vtkConnectivityFilter>::New();
                segmentConnectivity3->SetInputConnection(segmentConnectivity2->GetOutputPort());
                segmentConnectivity3->SetExtractionModeToAllRegions();
                //segmentConnectivity3->SetExtractionModeToLargestRegion();
                segmentConnectivity3->SetRegionIdAssignmentMode(0);
                segmentConnectivity3->ColorRegionsOn();
                segmentConnectivity3->Update();
        
                vtkSmartPointer<vtkDataSetTriangleFilter> triangulation = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
                //triangulation->SetInputConnection(surfaceFilter->GetOutputPort());
                triangulation->SetInputConnection(segmentConnectivity3->GetOutputPort());
                triangulation->Update();


                eigenFields.push_back(vtkSmartPointer<ttkEigenField>::New());
                eigenFields.back()->SetInputConnection(triangulation->GetOutputPort());
                acceptedRegions.push_back(uniqueDesSegIdDataSet->GetVariantValue(i).ToInt());

                
            }
            
        }
        else
        {
            if (sectionIDDataset->GetNumberOfPoints() > 300)
            {
                //Find the outer surface of the void space
                // vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
                // surfaceFilter->SetInputConnection(sectionID->GetOutputPort());
                // surfaceFilter->SetPieceInvariant(true);
                // surfaceFilter->Update();

                vtkSmartPointer<vtkDataSetTriangleFilter> triangulation = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
                //triangulation->SetInputConnection(surfaceFilter->GetOutputPort());
                triangulation->SetInputConnection(sectionID->GetOutputPort());
                triangulation->Update();

                eigenFields.push_back(vtkSmartPointer<ttkEigenField>::New());
                eigenFields.back()->SetInputConnection(triangulation->GetOutputPort());

                acceptedRegions.push_back(uniqueDesSegIdDataSet->GetVariantValue(i).ToInt());
            }
            

            // vtkSmartPointer<vtkUnstructuredGridWriter> regionWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
            // regionWriter->SetInputConnection(triangulation->GetOutputPort());
            // //regionWriter->SetFileName((Directory + "/Region_" + i + ".vtk").c)
            // regionWriter->SetFileName((Directory+"/Region_" + to_string(uniqueDesSegIdDataSet->GetVariantValue(i).ToInt()) + ".vtk").c_str());
            // regionWriter->Write();
        }
                
    }

    //#pragma omp parallel for
    for (size_t i = 0; i < acceptedRegions.size(); i++)
    {
        logger::mainlog << acceptedRegions[i] << endl;
        
        eigenFields[i]->SetUseAllCores(useAllCores);
        eigenFields[i]->SetEigenNumber(numberOfEigenFunctions);
        eigenFields[i]->SetOutputFieldName("EigenFunctions");
        eigenFields[i]->Update();
        
        
        //Extract the last EigenFunction
        vtkSmartPointer<vtkArrayCalculator> extraction = vtkSmartPointer<vtkArrayCalculator>::New();
        extraction->SetInputConnection(eigenFields[i]->GetOutputPort());
        extraction->SetAttributeTypeToPointData();
        extraction->AddScalarArrayName("EigenFunctions",numberOfEigenFunctions-1);
        extraction->SetFunction("EigenFunctions");
        extraction->SetResultArrayName("maxEigenFunction");
        extraction->Update();

        if (writeSegments)
        {
            //logger::mainlog << "Writing segment" << acceptedRegions[i] << endl;
            
            vtkSmartPointer<vtkUnstructuredGridWriter> regionWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
            regionWriter->SetInputConnection(extraction->GetOutputPort());
            regionWriter->SetFileName((Directory+"/" + BaseFileName + "_EigenField_" + to_string(acceptedRegions[i]) + ".vtk").c_str());
            regionWriter->Write();
        }

        vtkSmartPointer<ttkPersistenceDiagram> persistenceDiagram = vtkSmartPointer<ttkPersistenceDiagram>::New();
        persistenceDiagram->SetUseAllCores(useAllCores);
        persistenceDiagram->SetInputConnection(extraction->GetOutputPort());
        persistenceDiagram->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS, scalar.c_str());
        persistenceDiagram->Update();
        
        
        

        //logger::mainlog << "Writing persistence diagram" << acceptedRegions[i] << endl;

        vtkSmartPointer<vtkUnstructuredGridWriter> regionWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
        regionWriter->SetInputConnection(persistenceDiagram->GetOutputPort());
        regionWriter->SetFileName(("../Results/PersistenceDiagrams/"+ BaseFileName + "_" + to_string(acceptedRegions[i]) + ".vtk").c_str());
        regionWriter->Write();
    }
    

   
    logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
}

/**
 * @brief   Get the eigen field of the whole stucture of the material. Then ,it writes a .grid file with the scalar stored and a Persistence Diagram as fingerprint
 *
 *
 * @param grid Input grid from a readers function
 * @param numberOfEigenFunctions Number of eigenfunctions used to capture the shape(the more eigenfunctions, the greater the accuracy but with worse computation efficiency)
 * @param useAllCores
 * @return auto
 */
void segmentor::eigenStructure(vtkSmartPointer<vtkImageData> grid, int numberOfEigenFunctions, bool useAllCores)
{
    
    //Get the name of the material and create a folder to store the results
    //---------------------------------------------------------------------------------------------
    std::string base_filename = BaseFileName.substr(BaseFileName.find_last_of("-") + 1);
    std::string::size_type const p(base_filename.find_last_of('.'));
    std::string file_without_extension = base_filename.substr(0, p); //Get the material filename without extension nor directory

    int state = system("mkdir -p ../Results/GridFiles"); //Create a directory to save the results
    if(!state)
    {
        logger::mainlog << "Grid Files Folder created" << "\n";

    }

    int state2 = system("mkdir -p ../Results/Done"); //CDirectory to write files when a material is completed
    //---------------------------------------------------------------------------------------------
    
    
    //Segmentation corresponding to the void structure
    vtkSmartPointer<vtkThreshold> voidSegmentation = vtkSmartPointer<vtkThreshold>::New();
    //voidSegmentation->SetInputConnection(grid->GetOutputPort(0));
    voidSegmentation->SetInputData(grid);
    voidSegmentation->SetAllScalars(1);
    voidSegmentation->SetInputArrayToProcess(0,0,0,0,"This is distance grid");
    voidSegmentation->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    voidSegmentation->SetLowerThreshold(0.0);
    voidSegmentation->SetUpperThreshold(9e9);
    voidSegmentation->Update();

    //Find the outer surface of the void space
    vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    //surfaceFilter->SetInputConnection(append->GetOutputPort());
    surfaceFilter->SetInputConnection(voidSegmentation->GetOutputPort());
    surfaceFilter->SetPieceInvariant(true);
    surfaceFilter->Update();

    vtkSmartPointer<vtkDataSetTriangleFilter> triangulation = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
    triangulation->SetInputConnection(surfaceFilter->GetOutputPort());
    triangulation->Update();

    vtkSmartPointer<ttkEigenField> eigen = vtkSmartPointer<ttkEigenField>::New();
    eigen->SetInputConnection(triangulation->GetOutputPort());
    eigen->SetUseAllCores(useAllCores);
    eigen->SetEigenNumber(numberOfEigenFunctions);
    eigen->SetOutputFieldName("EigenFunctions");
    eigen->Update();

    //Extract the last EigenFunction
    vtkSmartPointer<vtkArrayCalculator> extraction = vtkSmartPointer<vtkArrayCalculator>::New();
    extraction->SetInputConnection(eigen->GetOutputPort());
    extraction->SetAttributeTypeToPointData();
    extraction->AddScalarArrayName("EigenFunctions",numberOfEigenFunctions-1);
    extraction->SetFunction("EigenFunctions");
    extraction->SetResultArrayName("maxEigenFunction");
    extraction->Update();

    auto eigenDataSet = vtkDataSet::SafeDownCast(extraction->GetOutputDataObject(0))->GetPointData()->GetAbstractArray("maxEigenFunction");
    string gridName = "../Results/GridFiles/" + file_without_extension + ".grid";
    ofstream outputFile(gridName);

    for (size_t i = 0; i < eigenDataSet->GetNumberOfValues(); i++)
    {
        outputFile << eigenDataSet->GetVariantValue(i).ToDouble() << endl;
    }
    outputFile.close();
    

    vtkSmartPointer<ttkPersistenceDiagram> persistenceDiagram = vtkSmartPointer<ttkPersistenceDiagram>::New();
    persistenceDiagram->SetUseAllCores(useAllCores);
    persistenceDiagram->SetInputConnection(extraction->GetOutputPort());
    persistenceDiagram->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"maxEigenFunction");
    persistenceDiagram->Update();

 
    vtkSmartPointer<vtkThreshold> criticalPairs = vtkSmartPointer<vtkThreshold>::New();
    criticalPairs->SetInputConnection(persistenceDiagram->GetOutputPort());
    criticalPairs->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"PairIdentifier");
    criticalPairs->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    criticalPairs->SetLowerThreshold(0.0);
    criticalPairs->SetUpperThreshold(9e9);
    criticalPairs->Update();

    auto persistenceDiagramDataSet = vtkDataSet::SafeDownCast(criticalPairs->GetOutputDataObject(0));
    //auto persistenceData = persistenceDiagramDataSet->GetPointData()->GetAbstractArray("Persistence");
    

    ofstream persistenceDiagrams;
    persistenceDiagrams.open("../Results/PersistenceDiagrams/" + file_without_extension + ".csv");
    assert(persistenceDiagrams.is_open());
    persistenceDiagrams << "Birth,Death" << "\n";

    for (size_t i = 1; i < persistenceDiagramDataSet->GetNumberOfPoints(); i=i+2)
    {
        double currentPointCoords[3];
        persistenceDiagramDataSet->GetPoint(i,currentPointCoords);
        persistenceDiagrams << currentPointCoords[0] << "," << currentPointCoords[1]  << "\n";
    }
    persistenceDiagrams.close();

    ofstream doneFile;
    doneFile.open("../Results/Done/" + file_without_extension);
    doneFile << "Checking\n";
    doneFile.close();
    

}




void segmentor::persistencecurve(vtkSmartPointer<ttkTriangulationManager> grid, bool useAllCores)
{
    ttk::Timer percurveTimer;
    logger::mainlog << "\nSegmentor: Persistence Curve module" << "\n" << flush;
 
    // To create the persistence curve, first the persistence diagram is calculated 
    // which serves as input to the persistence curve
    vtkSmartPointer<ttkPersistenceDiagram> persistenceDiagram = vtkSmartPointer<ttkPersistenceDiagram>::New();
    //persistenceDiagram->SetDebugLevel(3);
    persistenceDiagram->SetUseAllCores(useAllCores);
    persistenceDiagram->SetInputConnection(grid->GetOutputPort());
    persistenceDiagram->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());

    vtkNew<ttkPersistenceCurve> curve{};
    curve->SetInputConnection(persistenceDiagram->GetOutputPort());
    curve->Update();
    
    /* Write the persistence curve in VTK format */
    vtkNew<vtkTableWriter> curveWriter0{};
    curveWriter0->SetInputConnection(curve->GetOutputPort(0));
    curveWriter0->SetFileName((Directory+"/" + BaseFileName+"_minSaddlePairs.vtk").c_str());
    curveWriter0->Write();
    
    vtkNew<vtkTableWriter> curveWriter1{};
    curveWriter1->SetInputConnection(curve->GetOutputPort(1));
    curveWriter1->SetFileName((Directory+"/" + BaseFileName+"_saddleSaddlePairs.vtk").c_str());
    curveWriter1->Write();
    
    vtkNew<vtkTableWriter> curveWriter2{};
    curveWriter2->SetInputConnection(curve->GetOutputPort(2));
    curveWriter2->SetFileName((Directory+"/" + BaseFileName+"_SaddleMaxPairs.vtk").c_str());
    curveWriter2->Write();
    
    vtkNew<vtkTableWriter> curveWriter3{};
    curveWriter3->SetInputConnection(curve->GetOutputPort(3));
    curveWriter3->SetFileName((Directory+"/" + BaseFileName+"_allPairs.vtk").c_str());
    curveWriter3->Write();
    
    
    double timeTakenForPerCurve = percurveTimer.getElapsedTime();
    logger::mainlog << "Time taken in persistence curve module: " << timeTakenForPerCurve << "(s)" << endl;

}



/**
 * @brief Get the separatrices corresponding to the void structure and create a .csv file with the different separatrices
 *
 * @param morseSmaleComplex
 * @return auto
 */
void segmentor::voidSeparatrices(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex)
{
    logger::mainlog << "Void Separatrices Module" << "\n";
    logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
    ofstream sepFile;
    sepFile.open((Directory+"/Separatrices.csv").c_str());
    sepFile << "x,y,z,pointCellId,separatrixID,SepOriginID,SepDestinationID,SepMinValue,isMinima,isSaddle,xScaled,yScaled,zScaled" << "\n";
    
    //Descending separatrices of the MSC
    vtkSmartPointer<vtkThreshold> descendingSeparatrices = vtkSmartPointer<vtkThreshold>::New();
    descendingSeparatrices->SetInputConnection(morseSmaleComplex->GetOutputPort(1));
    descendingSeparatrices->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"SeparatrixType");
    descendingSeparatrices->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    descendingSeparatrices->SetLowerThreshold(0.0);
    descendingSeparatrices->SetUpperThreshold(0.0);

    //Descending separatrices of the MSC corresponding to the void space
    vtkSmartPointer<vtkThreshold> voidDescendingSeparatrices = vtkSmartPointer<vtkThreshold>::New();
    voidDescendingSeparatrices->SetInputConnection(descendingSeparatrices->GetOutputPort());
    voidDescendingSeparatrices->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"SeparatrixFunctionMaximum");
    voidDescendingSeparatrices->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    voidDescendingSeparatrices->SetLowerThreshold(-9e9);
    voidDescendingSeparatrices->SetUpperThreshold(0.0);


    auto descendingSepDataSet = vtkDataSet::SafeDownCast(voidDescendingSeparatrices->GetOutputDataObject(0));
    auto descendingSepPointDataSet = descendingSepDataSet->GetPointData();
    auto descendingSepCellDataSet = descendingSepDataSet->GetCellData();

    //We get the unique separatrix IDs
    vtkNew<ttkExtract> separatrixUniqueIdList{};
    separatrixUniqueIdList->SetUseAllCores(true);
    separatrixUniqueIdList->SetInputConnection(voidDescendingSeparatrices->GetOutputPort());
    separatrixUniqueIdList->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"SeparatrixId");
    separatrixUniqueIdList->SetExtractionMode(3); //Array Values
    separatrixUniqueIdList->SetExtractUniqueValues(true);
    separatrixUniqueIdList->Update();
    auto separatrixIdDataSet = vtkDataSet::SafeDownCast(separatrixUniqueIdList->GetOutputDataObject(0))->GetFieldData()->GetArray("UniqueSeparatrixId");

    for (size_t i = 0; i < separatrixIdDataSet->GetNumberOfValues(); i++) //For each of the separatrix IDs
    {
        //We get the current separatrix by its ID
        vtkNew<vtkThreshold> currentSeparatrix{};
        currentSeparatrix->SetInputConnection(voidDescendingSeparatrices->GetOutputPort());
        currentSeparatrix->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"SeparatrixId");
        int sepID = separatrixIdDataSet->GetVariantValue(i).ToInt();
        currentSeparatrix->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
        currentSeparatrix->SetLowerThreshold(sepID);
        currentSeparatrix->SetUpperThreshold(sepID);
        currentSeparatrix->Update();

        auto currentSeparatrixDataSet = vtkDataSet::SafeDownCast(currentSeparatrix->GetOutputDataObject(0));
        auto currentSeparatrixCellSet = currentSeparatrixDataSet->GetCellData();
        auto currentSeparatrixPointSet = currentSeparatrixDataSet->GetPointData();

        //Separatrix ID
        auto separatrixID = currentSeparatrixCellSet->GetArray("SeparatrixId")->GetVariantValue(0).ToInt();
        //Destination critical point CellID
        auto destinationID = currentSeparatrixCellSet->GetArray("DestinationId")->GetVariantValue(0).ToInt();
        //Source critical point CellID
        auto sourceID = currentSeparatrixCellSet->GetArray("SourceId")->GetVariantValue(0).ToInt();
        //Minimimum scalar value in the separatrix
        auto separatrixFunctionMinimum = currentSeparatrixCellSet->GetArray("SeparatrixFunctionMinimum")->GetVariantValue(0).ToDouble();

        for (size_t j = 0; j < currentSeparatrixDataSet->GetNumberOfPoints(); j++) //For each of the points of the current separatrix
        {
            double pointCoordinates[3]; //Current point coordinates
            currentSeparatrixDataSet->GetPoint(j,pointCoordinates);
            auto pointCellId = currentSeparatrixDataSet->GetPointData()->GetArray("CellId")->GetVariantValue(j).ToInt();
            bool isMinima = false;
            bool isSaddle = false;
            if (pointCellId == destinationID)
            {
                isMinima = true;
            }
            if (pointCellId == sourceID)
            {
                isSaddle = true;
            }
            sepFile << pointCoordinates[0] << "," << pointCoordinates[1] << "," << pointCoordinates[2] << "," << pointCellId << "," << separatrixID << "," << sourceID << "," << destinationID << "," << separatrixFunctionMinimum << "," << isMinima << "," << isSaddle << "," << pointCoordinates[0] << "," << pointCoordinates[1] << "," <<  pointCoordinates[2]  << "\n";
            
        }

    }
    logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
    
}



/**
 * @brief Computes ascending segmentation of the previous computed Morse Smales Segmentation. It also computes the closest regions ID Besides, it creates three files:
 * 1) fileName_Completo.csv -> CSV file with the information of all the segments of the input file.
 * 2) fileName_Region_i.csv -> CSV file with the information of the i region of the input file.
 * 3) propertiesEvolutionPython.csv -> CSV file with the information of the regions for each of the stages of compression(if we have them).
 * @param morseSmaleComplex Container with all the Morse Smale Complex computation results.
 */
auto segmentor::evolutionFile2( vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex, vtkSmartPointer<vtkThreshold> previousSolid)
{
    logger::mainlog << "segmentor: Evolution Module" << "\n";
    logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
    
    //Output stream to create evolutionPython.csv
    ofstream evolutionFilePython;

    //Check if the propertiesEvolutionPython file exist in order to create a new one or append the data to the existing one
    ifstream inFile;
    inFile.open("../propertiesEvolutionPython.csv");
    if (inFile)
    {
        inFile.close();
        evolutionFilePython.open("../propertiesEvolutionPython.csv",std::ios_base::app);
    }
    else
    {
        evolutionFilePython.open("../propertiesEvolutionPython.csv");
        evolutionFilePython << "Stage,RegionID,maxXCoord,maxYCoord,maxxZCoord,NumberOfPoints,AveragePotentialEnergy,AverageSigmaXX,AverageSigmaYY,AverageSigmaZZ,AverageSigmaXY,AverageSigmaXZ,AverageSigmaYZ,ClosestRegionID,numberOfConnections" << "\n";
    }
    
    vtkDataSet * previousDataSet;
    if (previousSolid) //Check that we are not in the first file
    {
        logger::mainlog << "We have previous stages information" << "\n";
        previousDataSet = vtkDataSet::SafeDownCast(previousSolid->GetOutputDataObject(0));
        //logger::mainlog << "Previous data set number of points" << previousDataSet->GetNumberOfPoints() << "\n";
    }
    else
    {
        logger::mainlog << "We don't have previous stages information" << "\n";
    }
    //---------------------------------------------------------------------------------------------
    //We get the 2-saddle criticalPoints
    vtkSmartPointer<vtkThresholdPoints> saddles = vtkSmartPointer<vtkThresholdPoints>::New();
    saddles->SetInputConnection(morseSmaleComplex->GetOutputPort(0));
    saddles->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"CellDimension");
    saddles->ThresholdBetween(2,2);
    saddles->Update();

    auto prueba = vtkDataSet::SafeDownCast(saddles->GetOutputDataObject(0))->GetNumberOfPoints();
    //We get the 2-saddles with positve distance(solid)
    vtkSmartPointer<vtkThresholdPoints> positiveSaddles = vtkSmartPointer<vtkThresholdPoints>::New();
    positiveSaddles->SetInputConnection(saddles->GetOutputPort(0));
    positiveSaddles->SetInputArrayToProcess(0,0,0,0,"This is distance grid");
    positiveSaddles->ThresholdBetween(0,1e9);
    positiveSaddles->Update();

    auto positiveSaddlesDataSet = vtkDataSet::SafeDownCast(positiveSaddles->GetOutputDataObject(0));

    
    //------------------------------------------------------------------------------------------------
    
    //Segmentation corresponding to the solid structure
    vtkSmartPointer<vtkThreshold> solidSegmentation = vtkSmartPointer<vtkThreshold>::New();
    solidSegmentation->SetInputConnection(morseSmaleComplex->GetOutputPort(3));
    // solidSegmentation->SetInputArrayToProcess(0,0,0,0,"potentialEnergyAtom");
    
    solidSegmentation->SetInputArrayToProcess(0,0,0,0,"This is distance grid");
    solidSegmentation->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    solidSegmentation->SetLowerThreshold(-9.9e9);
    solidSegmentation->SetUpperThreshold(-1e-9);
    solidSegmentation->Update();


    //Solid Segmentation DataSet
    auto solidSegmentationDataSet = vtkDataSet::SafeDownCast(solidSegmentation->GetOutputDataObject(0));

    //------------------------------------------------------------------------------------------------
    // vtkSmartPointer<vtkIntArray> connectionNumber = vtkSmartPointer<vtkIntArray>::New();
    // connectionNumber->SetName("NumberOfConnections");
    // connectionNumber->SetNumberOfValues(solidSegmentationDataSet->GetNumberOfPoints());
    // solidSegmentationDataSet->GetPointData()->AddArray(connectionNumber);
    //------------------------------------------------------------------------------------------------

    //Array corresponding to the Ascending Segmentation of the Solid
    auto ascendingManifoldIDArray = solidSegmentationDataSet->GetPointData()->GetArray("AscendingManifold");
    //Same structure segmentation but with a Field Data added
    vtkSmartPointer<ttkExtract> ascendingManifoldIDList = vtkSmartPointer<ttkExtract>::New();
    ascendingManifoldIDList->SetDebugLevel(1);
    ascendingManifoldIDList->SetUseAllCores(true);
    ascendingManifoldIDList->SetInputConnection(solidSegmentation->GetOutputPort());
    ascendingManifoldIDList->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"AscendingManifold");
    ascendingManifoldIDList->SetExtractionMode(3); //Array Values
    ascendingManifoldIDList->SetExtractUniqueValues(true);
    ascendingManifoldIDList->Update();
    auto uniqueAscSegIdDataSet = vtkDataSet::SafeDownCast(ascendingManifoldIDList->GetOutputDataObject(0))->GetFieldData()->GetAbstractArray("UniqueAscendingManifold");
    //Writer of the .csv results file
    ofstream misDatos;
    misDatos.open((Directory+"/"+BaseFileName+"_Completo.csv").c_str());
    assert(misDatos.is_open());
    misDatos << "regionID,x,y,z,Scalar,RegionMaxValue,minValueID,isMaximum,numberOfPoints,numberOfConnections" << "\n";

    double averagePotentialEnergy; //Average Potential Energy of the Region
    double averageSigmaXXStress; //Average SigmaXX Stress of the Region
    double averageSigmaYYStress; //Average SigmaYY Stress of the Region
    double averageSigmaZZStress; //Average SigmaZZ Stress of the Region
    double averageSigmaXYStress; //Average SigmaXY Stress of the Region
    double averageSigmaXZStress; //Average SigmaXZ Stress of the Region
    double averageSigmaYZStress; //Average SigmaYZ Stress of the Region

    for (size_t i = 0; i < uniqueAscSegIdDataSet->GetNumberOfValues(); i++) //For each regoin of the Ascending Segmentation
    {
        
        
        
        ofstream regionData; //Stream to write segmentation data for each of the regions
        regionData.open((Directory+"/"+BaseFileName+"_Region"+to_string(i)+".csv").c_str()); //Create de file
        assert(regionData.is_open()); //Check if the file is open
        regionData << "regionID,x,y,z,Scalar,RegionMaxValue,maxValueID,isMaximum,numberOfPoints,numberOfConnections" << "\n"; //Write the column names

        //Current Region of the Ascending Segmentation
        vtkSmartPointer<vtkThreshold> sectionID = vtkSmartPointer<vtkThreshold>::New();
        sectionID->SetInputConnection(solidSegmentation->GetOutputPort());
        sectionID->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"AscendingManifold");
        int segID = uniqueAscSegIdDataSet->GetVariantValue(i).ToInt();
        sectionID->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
        sectionID->SetLowerThreshold(segID);
        sectionID->SetUpperThreshold(segID);
        sectionID->Update();
                
        //DataSet of the specific region of the Ascending Segmentation
        auto sectionIDDataset = vtkDataSet::SafeDownCast(sectionID->GetOutputDataObject(0));
        //-----------------------------------------------------------------------------------------
        int numberOfConnections = 0;
        if (sectionIDDataset->GetNumberOfPoints() > 0) //If the current region has more than one point
        {
            for (size_t j = 0; j < positiveSaddlesDataSet->GetNumberOfPoints(); j++) //For each of the saddles of the solid
            {
                double currentSaddleCoords[3];
                positiveSaddlesDataSet->GetPoint(j,currentSaddleCoords); //Save the saddle currentCoords
                int closestRegionPoint = sectionIDDataset->FindPoint(currentSaddleCoords); //Find the closest point in the current region
                double currentClosestCoords[3]; //Closest point coords
                sectionIDDataset->GetPoint(closestRegionPoint,currentClosestCoords);
                double distance = sqrt(pow(currentSaddleCoords[0]-currentClosestCoords[0],2)+pow(currentSaddleCoords[1]-currentClosestCoords[1],2)+pow(currentSaddleCoords[2]-currentClosestCoords[2],2));

                if (distance < (CellSize))
                {
                    ++numberOfConnections;
                }
                
            }
        }
        
        
        
        //-----------------------------------------------------------------------------------------
        //Array corresponding to the scalar values of the Region
        auto scalarValues = sectionIDDataset->GetPointData()->GetArray("This is distance grid");

        //Critical value of the maximum scalar value critical point
        double maximumValue = 0;
        //Coordinates of the maximum scalar value critical point
        double maximumCoords[3];
        //Index of the maximum scalar value critical point
        int maxID{};
        //Critical Pointds Region Maximum ID
        int criticalPointsRegionMaximumID = 0;

                
        double acumulatedPotentialEnergy = 0; //Acumulated Potential Energy of the Region
        double acumulatedSigmaXXStress = 0; //Acumulated SigmaXX Stress of the Region
        double acumulatedSigmaYYStress = 0; //Acumulated SigmaYY Stress of the Region
        double acumulatedSigmaZZStress = 0; //Acumulated SigmaZZ Stress of the Region
        double acumulatedSigmaXYStress = 0; //Acumulated SigmaXY Stress of the Region
        double acumulatedSigmaXZStress = 0; //Acumulated SigmaXZ Stress of the Region
        double acumulatedSigmaYZStress = 0; //Acumulated SigmaYZ Stress of the Region

        vector<int> closestPreviousStageRegionID; //Store all the closest regions to the points here

        for (size_t j = 0; j < scalarValues->GetNumberOfValues(); j++) //For each of the points of the region
        {
            
            if (previousSolid) //If we have information from previous stages
            {
                double currentCoords[3]; //Current Point coords
                sectionIDDataset->GetPoint(j,currentCoords); //Find current point coords
                auto previousPointID = previousDataSet->FindPoint(currentCoords); //Find the closest point in the previous data set
                closestPreviousStageRegionID.push_back(previousDataSet->GetPointData()->GetArray("AscendingManifold")->GetVariantValue(previousPointID).ToInt()); //Get the Region ID of the previous stage segmentation
            }
            
            
            double currentValue = scalarValues->GetVariantValue(j).ToDouble(); //Current scalar value of the point
            if (currentValue >= maximumValue) //Check if this point distance is greater than maximum
            {
                maximumValue = currentValue; //If greater, we set a new maximum value
                sectionIDDataset->GetPoint(j,maximumCoords); //Find maximum coords
                criticalPointsRegionMaximumID = vtkDataSet::SafeDownCast(morseSmaleComplex->GetOutputDataObject(0))->FindPoint(maximumCoords); //Find  maximum (critical point) ID
                maxID = j;
            }
            acumulatedPotentialEnergy += sectionIDDataset->GetPointData()->GetArray("potentialEnergyAtom")->GetVariantValue(j).ToDouble();
            acumulatedSigmaXXStress += sectionIDDataset->GetPointData()->GetArray("sigmaXX_Atom")->GetVariantValue(j).ToDouble();
            acumulatedSigmaYYStress += sectionIDDataset->GetPointData()->GetArray("sigmaYY_Atom")->GetVariantValue(j).ToDouble();
            acumulatedSigmaZZStress += sectionIDDataset->GetPointData()->GetArray("sigmaZZ_Atom")->GetVariantValue(j).ToDouble();
            acumulatedSigmaXYStress += sectionIDDataset->GetPointData()->GetArray("sigmaXY_Atom")->GetVariantValue(j).ToDouble();
            acumulatedSigmaXZStress += sectionIDDataset->GetPointData()->GetArray("sigmaXZ_Atom")->GetVariantValue(j).ToDouble();
            acumulatedSigmaYZStress += sectionIDDataset->GetPointData()->GetArray("sigmaYZ_Atom")->GetVariantValue(j).ToDouble();

        }

        //auto closestRegionToTheMaximum = previousDataSet->FindPoint(maximumCoords)
        
        //Is this point maximum?
        bool isMaxima = 0;
            
        for (size_t k = 0; k < sectionIDDataset->GetNumberOfPoints(); k++) //For each of the points of the region check if this point is maxima and write values
        {
            double pointCoords[3];
            sectionIDDataset->GetPoint(k,pointCoords);
            
            if(k == maxID)
            {
                isMaxima = 1;
            }
            else
            {
                isMaxima = 0;
            }
                
            misDatos << uniqueAscSegIdDataSet->GetVariantValue(i).ToInt() <<","<< pointCoords[0]<<","<<pointCoords[1]<<","<<pointCoords[2]<<","<<scalarValues->GetVariantValue(k).ToDouble()<< "," << maximumValue<<","<<criticalPointsRegionMaximumID<<","<< isMaxima <<","<< sectionIDDataset->GetNumberOfPoints()<<"," << numberOfConnections<<"\n";
            regionData << uniqueAscSegIdDataSet->GetVariantValue(i).ToInt() <<","<< pointCoords[0]<<","<<pointCoords[1]<<","<<pointCoords[2]<<","<<scalarValues->GetVariantValue(k).ToDouble()<< "," << maximumValue<<","<<criticalPointsRegionMaximumID<<","<< isMaxima <<","<< sectionIDDataset->GetNumberOfPoints()<<"," << numberOfConnections<<"\n";
                    
        }

        regionData.close();

        int closestRegion; //Closest region ID of the previous stage
        if (previousSolid)
        {
            closestRegion = findMostCommonValue(closestPreviousStageRegionID);
            //logger::mainlog << "Closest region: " << closestRegion << "\n";
            //logger::mainlog << "Closest region to the maxima: " << previousDataSet->GetPointData()->GetArray("AscendingManifold")->GetVariantValue(previousDataSet->FindPoint(maximumCoords)).ToInt() << "\n";

        }

        averagePotentialEnergy = acumulatedPotentialEnergy/scalarValues->GetNumberOfValues();
        averageSigmaXXStress = acumulatedSigmaXXStress/scalarValues->GetNumberOfValues();
        averageSigmaYYStress = acumulatedSigmaYYStress/scalarValues->GetNumberOfValues();
        averageSigmaZZStress = acumulatedSigmaZZStress/scalarValues->GetNumberOfValues();
        averageSigmaXYStress = acumulatedSigmaXYStress/scalarValues->GetNumberOfValues();
        averageSigmaXZStress = acumulatedSigmaXZStress/scalarValues->GetNumberOfValues();
        averageSigmaYZStress = acumulatedSigmaYZStress/scalarValues->GetNumberOfValues();

        if (previousSolid)
        {
            evolutionFilePython << BaseFileName.substr(BaseFileName.find_last_of(".")+1,BaseFileName.size()) << "," << uniqueAscSegIdDataSet->GetVariantValue(i).ToInt() << "," << maximumCoords[0] << "," << maximumCoords[1] << "," << maximumCoords[2] << "," << sectionIDDataset->GetNumberOfPoints() << "," << averagePotentialEnergy << "," << averageSigmaXXStress << "," << averageSigmaYYStress << "," << averageSigmaZZStress<< "," << averageSigmaXYStress << "," << averageSigmaXZStress << "," << averageSigmaYZStress << "," << closestRegion <<"," << numberOfConnections<< "\n";
        }
        else
        {
            evolutionFilePython << BaseFileName.substr(BaseFileName.find_last_of(".")+1,BaseFileName.size()) << "," << uniqueAscSegIdDataSet->GetVariantValue(i).ToInt() << "," << maximumCoords[0] << "," << maximumCoords[1] << "," << maximumCoords[2] << "," << sectionIDDataset->GetNumberOfPoints() << "," << averagePotentialEnergy << "," << averageSigmaXXStress << "," << averageSigmaYYStress << "," << averageSigmaZZStress<< "," << averageSigmaXYStress << "," << averageSigmaXZStress << "," << averageSigmaYZStress << "," << -1 <<"," << numberOfConnections<< "\n";
        }
        
        
           
    }
    misDatos.close();
    evolutionFilePython.close();
    
    logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
    return solidSegmentation;

}



/**
 * @brief Get the solid structure from a segmentation file obtained from a Morse Smale Complex
 * computation
 *
 * @param currentFile Name of the current segmentation file
 * @return auto Solid structure of the segmentation file
 */
auto segmentor::solidGetter(string currentFile)
{
    
    
    string inputFile = "../results_dump.stress." + currentFile + "/segmentation.vtk"; //Current stage input file
    vtkSmartPointer<vtkDataSetReader> currentStageReader = vtkSmartPointer<vtkDataSetReader>::New();
    currentStageReader->SetFileName((inputFile).c_str());
    currentStageReader->Update();

    vtkSmartPointer<vtkThresholdPoints> currentSolid = vtkSmartPointer<vtkThresholdPoints>::New();
    currentSolid->SetInputConnection(currentStageReader->GetOutputPort());
    currentSolid->SetInputArrayToProcess(0,0,0,0,"This is distance grid");
    currentSolid->ThresholdBetween(1e-9,1e9);
    currentSolid->Update();

    //Same structure segmentation but with a Field Data added
    vtkSmartPointer<ttkExtract> ascendingManifoldIDList = vtkSmartPointer<ttkExtract>::New();
    ascendingManifoldIDList->SetDebugLevel(1);
    ascendingManifoldIDList->SetUseAllCores(true);
    ascendingManifoldIDList->SetInputConnection(currentSolid->GetOutputPort());
    ascendingManifoldIDList->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"AscendingManifold");
    ascendingManifoldIDList->SetExtractionMode(3); //Array Values
    ascendingManifoldIDList->SetExtractUniqueValues(true);
    ascendingManifoldIDList->Update();

    
    auto currentSolidDataSet = vtkDataSet::SafeDownCast(ascendingManifoldIDList->GetOutputDataObject(0));
    //logger::mainlog << currentSolidDataSet->GetNumberOfPoints() << "\n";
    
    return ascendingManifoldIDList;


}

/**
 * @brief Find the 2-saddles located in the solid structure
 *
 * @param currentFile Name of the file currenly readed
 * @return auto 2-saddles located in the solid structure
 */
auto segmentor::saddlesGetter(string currentFile)
{
    string inputFile = "../results_dump.stress." + currentFile + "/criticalPoints.vtk"; //Current stage input file
    //Read critical points input file
    vtkSmartPointer<vtkPolyDataReader> currentStageReader = vtkSmartPointer<vtkPolyDataReader>::New();
    currentStageReader->SetFileName((inputFile).c_str());
    currentStageReader->Update();
    //Find the 2-saddle critical points
    vtkSmartPointer<vtkThresholdPoints> saddles = vtkSmartPointer<vtkThresholdPoints>::New();
    saddles->SetInputConnection(currentStageReader->GetOutputPort());
    saddles->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"CellDimension");
    saddles->ThresholdBetween(2,2);
    saddles->Update();
    //Find the 2-saddles on the solid structure
    vtkSmartPointer<vtkThresholdPoints> positiveSaddles = vtkSmartPointer<vtkThresholdPoints>::New();
    positiveSaddles->SetInputConnection(saddles->GetOutputPort(0));
    positiveSaddles->SetInputArrayToProcess(0,0,0,0,"This is distance grid");
    positiveSaddles->ThresholdBetween(1e-9,1e9);
    positiveSaddles->Update();

    return positiveSaddles;
}


/**
 * @brief Find, for each of the stages of compression, and for each regions(of the Ascending Manifolds) the closest region on the following and
 * the previous stage
 *
 * @param fileNames Name of the different regions that the function will compute
 * @return auto
 */
auto segmentor::stagesEvolution(vector<string> fileNames)
{
    logger::mainlog << "Stages Evolution Module" << "\n";
    logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
    ofstream evolutionData; //Stream to write segmentation data for each of the regions
    evolutionData.open("../stagesEvolution.csv"); //Opening the writer
    assert(evolutionData.is_open()); //Check if the file is open
    evolutionData << "Stage,RegionID,NextID,PrevID,NumberOfPoints,NumberOfConnections,AveragePotentialEnergy,AverageSigmaXX,AverageSigmaYY,AverageSigmaZZ,AverageSigmaXY,AverageSigmaXZ,AverageSigmaYZ" << "\n"; //Write the column names
    
    ofstream connectivityData; //Stream to write segmentation data for each of the regions
    connectivityData.open("../connectivityData.csv"); //Opening the writer
    assert(connectivityData.is_open()); //Check if the file is open
    connectivityData << "Stage,SaddleNumber,ConnectedRegion,ConnectedRegion,ConnectedRegion,ConnectedRegion,ConnectedRegion,ConnectedRegion,ConnectedRegion" << "\n"; //Write the column names
    

    for (size_t i = 0; i < fileNames.size(); i++) //For each of the stages of compression
    {
        string inputFile = "../results_dump.stress." + fileNames[i] + "/segmentation.vtk"; //Current stage input file
        vtkSmartPointer<vtkDataSetReader> currentStageReader = vtkSmartPointer<vtkDataSetReader>::New();
        currentStageReader->SetFileName((inputFile).c_str());
        currentStageReader->Update();

        auto provisionalDataSet = vtkDataSet::SafeDownCast(currentStageReader->GetOutputDataObject(0));

        double dimensionesCelda[6];
        provisionalDataSet->GetCellBounds(0,dimensionesCelda);
        //Cell size of the current dataset
        double cellSize = dimensionesCelda[1] - dimensionesCelda[0];
        logger::mainlog << "Tamaño celda " << cellSize << "\n";

        vtkSmartPointer<vtkThresholdPoints> currentSolid = vtkSmartPointer<vtkThresholdPoints>::New();
        currentSolid->SetInputConnection(currentStageReader->GetOutputPort());
        currentSolid->SetInputArrayToProcess(0,0,0,0,"This is distance grid");
        currentSolid->ThresholdBetween(1e-9,1e9);
        currentSolid->Update();

        //Same structure segmentation but with a Field Data added
        vtkSmartPointer<ttkExtract> currentStage = vtkSmartPointer<ttkExtract>::New();
        currentStage->SetDebugLevel(1);
        currentStage->SetUseAllCores(true);
        currentStage->SetInputConnection(currentSolid->GetOutputPort());
        currentStage->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"AscendingManifold");
        currentStage->SetExtractionMode(3); //Array Values
        currentStage->SetExtractUniqueValues(true);
        currentStage->Update();


        //auto currentStage = solidGetter(fileNames[i]);
        
        //DataSet corresponding to the current stage solid structure
        auto currentStageDataSet = vtkDataSet::SafeDownCast(currentStage->GetOutputDataObject(0));
       
        
        
        logger::mainlog << "Current Stage: " << fileNames[i] << "\n";
        logger::mainlog << "Current Stage:Number of Points: " << currentStageDataSet->GetNumberOfPoints() << "\n";
        auto segmentationIDS = currentStageDataSet->GetFieldData()->GetAbstractArray("UniqueAscendingManifold");
        //---------------------------------------------------------------------------

        double averagePotentialEnergy; //Average Potential Energy of the Region
        double averageSigmaXXStress; //Average SigmaXX Stress of the Region
        double averageSigmaYYStress; //Average SigmaYY Stress of the Region
        double averageSigmaZZStress; //Average SigmaZZ Stress of the Region
        double averageSigmaXYStress; //Average SigmaXY Stress of the Region
        double averageSigmaXZStress; //Average SigmaXZ Stress of the Region
        double averageSigmaYZStress; //Average SigmaYZ Stress of the Region

        
        vtkSmartPointer<ttkExtract> nextStageData;
        if (fileNames[i] != "300000") //If not on the last stage
        {
            nextStageData = solidGetter(fileNames[i+1]);
        }
        
        vtkSmartPointer<ttkExtract> prevStageData;
        //vtkDataSet * prevStageDataSet;
        if (fileNames[i] != "87690") //If not on the first stage
        {
            prevStageData = solidGetter(fileNames[i-1]);
        }
        
        //---------------------------------------------------------------------------
        //2-saddles corresponding to the solid structure of the current stage
        auto saddles = saddlesGetter(fileNames[i]);
        auto saddlesDataSet = vtkDataSet::SafeDownCast(saddles->GetOutputDataObject(0));
        
        vector<vector<int>> connectivityMatrix; //Matrix containing all the 2-saddles and the regions they are connected to
        connectivityMatrix.resize(saddlesDataSet->GetNumberOfPoints(),vector<int>(8,-1.0)); //Resize to make it a 2D vector
        for (size_t k = 0; k < connectivityMatrix.size(); k++)
        {
            connectivityMatrix[k][1] = k; //Critical Point Number
            connectivityMatrix[k][0] = i; //Stage number
        }
        
        for (size_t j = 0; j < segmentationIDS->GetNumberOfValues(); j++) //For each of the segmentation regions
        {
            
            
            //---------------------------------------------------------------------------
            //Current Region of the Ascending Segmentation
            vtkSmartPointer<vtkThresholdPoints> currentRegion = vtkSmartPointer<vtkThresholdPoints>::New();
            currentRegion->SetInputConnection(currentStage->GetOutputPort());
            currentRegion->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"AscendingManifold");
            currentRegion->ThresholdBetween(segmentationIDS->GetVariantValue(j).ToInt(),segmentationIDS->GetVariantValue(j).ToInt());
            currentRegion->Update();
            auto currentRegionDataSet = vtkDataSet::SafeDownCast(currentRegion->GetOutputDataObject(0));
            //logger::mainlog << "CurrentRegion:Number of Points: " << currentRegionDataSet->GetNumberOfPoints() << "\n";
            //---------------------------------------------------------------------------
            int numberOfConnections = 0; //Number of connections of the current region
            vector<int> regionsConnected; //ID of the regions connected to the current region
            for (size_t kk = 0; kk < 20; kk++)
            {
                regionsConnected.push_back(-1);
            }
            
            if (currentRegionDataSet->GetNumberOfPoints() > 0) //If not a void region
            {
                int contador = 0;
                for (size_t k = 0; k < saddlesDataSet->GetNumberOfPoints(); k++)
                {
                    double currentSaddleCoords[3]; //Coordinates of the current saddle
                    saddlesDataSet->GetPoint(k,currentSaddleCoords);
                    int closestRegionPoint = currentRegionDataSet->FindPoint(currentSaddleCoords); //ID of the closest point to the saddle from the region
                    double currentClosestCoords[3];
                    currentRegionDataSet->GetPoint(closestRegionPoint,currentClosestCoords); //Save the coordinates of the closest region points
                    double distance = sqrt(pow(currentSaddleCoords[0]-currentClosestCoords[0],2)+pow(currentSaddleCoords[1]-currentClosestCoords[1],2)+pow(currentSaddleCoords[2]-currentClosestCoords[2],2));

                    if (distance < (2*cellSize)) //If the distance is less than the cell size
                    {
                        ++numberOfConnections; //Increase the number of connections in one
                        connectivityMatrix[k][getIndex(connectivityMatrix[k],-1)] = segmentationIDS->GetVariantValue(j).ToInt();
                        
                    }

                }
                
                
            }
            



            //---------------------------------------------------------------------------

            
            double acumulatedPotentialEnergy = 0; //Acumulated Potential Energy of the Region
            double acumulatedSigmaXXStress = 0; //Acumulated SigmaXX Stress of the Region
            double acumulatedSigmaYYStress = 0; //Acumulated SigmaYY Stress of the Region
            double acumulatedSigmaZZStress = 0; //Acumulated SigmaZZ Stress of the Region
            double acumulatedSigmaXYStress = 0; //Acumulated SigmaXY Stress of the Region
            double acumulatedSigmaXZStress = 0; //Acumulated SigmaXZ Stress of the Region
            double acumulatedSigmaYZStress = 0; //Acumulated SigmaYZ Stress of the Region

            vector<int> closestNextStageRegionID; //Store the closest regions for each of the points
            vector<int> closestPrevStageRegionID; //Store the closest regions for each of the points
            
            for (size_t k = 0; k < currentRegionDataSet->GetNumberOfPoints(); k++) //For each of the points of the region
            {
                
                double currentPointCoords[3];
                currentRegionDataSet->GetPoint(k,currentPointCoords); //Save current point coordinates
                
                //logger::mainlog << currentPointCoords[0] <<","<< currentPointCoords[1] << "," << currentPointCoords[2] << "\n";

                acumulatedPotentialEnergy += currentRegionDataSet->GetPointData()->GetArray("potentialEnergyAtom")->GetVariantValue(k).ToDouble();
                acumulatedSigmaXXStress += currentRegionDataSet->GetPointData()->GetArray("sigmaXX_Atom")->GetVariantValue(k).ToDouble();
                acumulatedSigmaYYStress += currentRegionDataSet->GetPointData()->GetArray("sigmaYY_Atom")->GetVariantValue(k).ToDouble();
                acumulatedSigmaZZStress += currentRegionDataSet->GetPointData()->GetArray("sigmaZZ_Atom")->GetVariantValue(k).ToDouble();
                acumulatedSigmaXYStress += currentRegionDataSet->GetPointData()->GetArray("sigmaXY_Atom")->GetVariantValue(k).ToDouble();
                acumulatedSigmaXZStress += currentRegionDataSet->GetPointData()->GetArray("sigmaXZ_Atom")->GetVariantValue(k).ToDouble();
                acumulatedSigmaYZStress += currentRegionDataSet->GetPointData()->GetArray("sigmaYZ_Atom")->GetVariantValue(k).ToDouble();
                
                
                if (fileNames[i] != "300000") //If not on the last stage
                {
                    
                    auto nextStageDataSet = vtkDataSet::SafeDownCast(nextStageData->GetOutputDataObject(0));
                    //---------------------------------------------------------------------------
                    auto nextClosestPointID = nextStageDataSet->FindPoint(currentPointCoords); //Find the closest point ID in the following stage
                    double pruebaCoords[3];
                    nextStageDataSet->GetPoint(nextClosestPointID,pruebaCoords);
                    auto closestNextPoint =nextStageDataSet->GetPointData()->GetAbstractArray("AscendingManifold")->GetVariantValue(nextClosestPointID).ToInt();
                    closestNextStageRegionID.push_back(closestNextPoint);

                }

                if (fileNames[i] != "87690") //If not on the first stage
                {
                    auto prevStageDataSet = vtkDataSet::SafeDownCast(prevStageData->GetOutputDataObject(0));
                    //---------------------------------------------------------------------------
                    auto prevClosestPointID = prevStageDataSet->FindPoint(currentPointCoords); //Find the closest point ID
                    auto closestPrevPoint =prevStageDataSet->GetPointData()->GetAbstractArray("AscendingManifold")->GetVariantValue(prevClosestPointID).ToInt();
                    closestPrevStageRegionID.push_back(closestPrevPoint);

                }
                
            }
        
            
            

            averagePotentialEnergy = acumulatedPotentialEnergy/currentRegionDataSet->GetNumberOfPoints();
            averageSigmaXXStress = acumulatedSigmaXXStress/currentRegionDataSet->GetNumberOfPoints();
            averageSigmaYYStress = acumulatedSigmaYYStress/currentRegionDataSet->GetNumberOfPoints();
            averageSigmaZZStress = acumulatedSigmaZZStress/currentRegionDataSet->GetNumberOfPoints();
            averageSigmaXYStress = acumulatedSigmaXYStress/currentRegionDataSet->GetNumberOfPoints();
            averageSigmaXZStress = acumulatedSigmaXZStress/currentRegionDataSet->GetNumberOfPoints();
            averageSigmaYZStress = acumulatedSigmaYZStress/currentRegionDataSet->GetNumberOfPoints();

            int closestNext;
            if (fileNames[i] != "300000")
            {
                closestNext = findMostCommonValue(closestNextStageRegionID);
                //logger::mainlog << "Closest Next Region:" << closestNext << "\n";
            }
            else if(fileNames[i] == "300000")
            {
                closestNext = -1;
                //logger::mainlog << "Closest Next Region:" << closestNext << "\n";
            }

            int closestPrev;
            if (fileNames[i] != "87690")
            {
                closestPrev = findMostCommonValue(closestPrevStageRegionID);
                //logger::mainlog << "Closest Prev Region:" << closestPrev << "\n";
            }
            if (fileNames[i] == "87690")
            {
                closestPrev = -1;
                //logger::mainlog << "Closest Prev Region:" << closestPrev << "\n";
            }
            
            int currentReg = segmentationIDS->GetVariantValue(j).ToInt();
            //logger::mainlog << "Current Region:" << currentReg << "\n";
            evolutionData << fileNames[i] << "," << currentReg << "," << closestNext << "," << closestPrev << "," << currentRegionDataSet->GetNumberOfPoints() <<","<< numberOfConnections << "," << averagePotentialEnergy << "," << averageSigmaXXStress << "," << averageSigmaYYStress << "," << averageSigmaZZStress<< "," << averageSigmaXYStress << "," << averageSigmaXZStress << "," << averageSigmaYZStress << "\n";
            
            // if (regionsConnected.size() > 0)
            // {
            //     connectivityData << fileNames[i] << "," << currentReg << "," << numberOfConnections << ",";
            //     for (size_t m = 0; m < regionsConnected.size(); m++)
            //     {
            //         connectivityData << regionsConnected[m] << ",";
            //     }
            //     connectivityData << "\n";
                
            // }
           
        }
        for (size_t m = 0; m < connectivityMatrix.size(); m++)
        {
            for (size_t n = 0; n < connectivityMatrix[0].size(); n++)
            {
                if (n != connectivityMatrix[0].size())
                {
                    connectivityData << connectivityMatrix[m][n] << ",";
                }
                else if(n == connectivityMatrix[0].size())
                {
                    connectivityData << connectivityMatrix[m][n];
                }
                
                
            }
            connectivityData << "\n";
            
                
        }
        
        
        
        

    }
    evolutionData.close();
    connectivityData.close();
    
}

/**
 * @brief Find the most common element in a vector
 *
 * @param inputVector Input vector
 * @return int Most common element in the vector
 */
int segmentor::findMostCommonValue(vector<int> &inputVector)
{
    if (inputVector.empty())
    {
        return -1;
    }
    

    sort(inputVector.begin(),inputVector.end());
    
    auto lastInt = inputVector.front(); //Initialize in the beggining of the vector
    auto mostFreqInt = inputVector.front(); //Initialize in the beggining of the vector
    int maxFreq = 0, currentFreq = 0;

    for(const auto &i : inputVector)
    {
        if (i == lastInt) //Check that the current element is the same to the previous one
        {
            ++currentFreq; //Adds 1 to current frequency
        }
        else
        {
            if (currentFreq > maxFreq) //Check that the current frequency is the maximum frequency
            {
                maxFreq = currentFreq;
                mostFreqInt = lastInt;
            }
            lastInt = i; //Start again
            currentFreq = 1;
        }
    }
    if (currentFreq > maxFreq)
    {
        maxFreq = currentFreq;
        mostFreqInt = lastInt;
    }

    return mostFreqInt;
    
    
}

/**
 * @brief Check if the evolution file exists and delete it if it is true
 * @param fileName File to be checked
 */
void isFileExist(const char *fileName)
{
    ifstream checkFile;
    checkFile.open(fileName);
    if (checkFile)
    {
        checkFile.close();
        int status = remove(fileName);
        if (status == 0)
        {
            logger::mainlog << "Removed propertiesEvolutionPython.csv file from previous computations" << "\n";
        }
        else
        {
            logger::mainlog << "Error while trying to remove propertiesEvolutionPython.csv file from previous computations" << "\n";
        }
     
    }
}


/**
 * @brief Compute the persistence diagram of a energy field of an input material and writes it to a .grid file
 *
 * @param grid Input grid from a reader function
 * @param useAllCores Use all cores available
 */
void segmentor::energyDiagrams(vtkSmartPointer<vtkImageData> grid, bool useAllCores)
{
    logger::mainlog << "segmentor: Energy Diagrams Module" << "\n";
    logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
    
    //VTK Function to apply Periodic Boundary Conditions
    vtkSmartPointer<ttkTriangulationManager> periodGrid = vtkSmartPointer<ttkTriangulationManager>::New();;
    periodGrid->SetUseAllCores(useAllCores);
    periodGrid->SetInputData(grid);
    periodGrid->SetPeriodicity(true);
    periodGrid->Update();

    //Input Grid DataSet
    auto inputGridDataSet = vtkDataSet::SafeDownCast(periodGrid->GetOutputDataObject(0));
    //Potential Energy DataSet of the input grid
    auto energyDataSet = inputGridDataSet->GetPointData()->GetAbstractArray("Potential Energy");

    double minimumEnergy = 0.0; //Minimum energy value of the material
    double maxEnergy = 0.0;
    for (size_t i = 0; i < energyDataSet->GetNumberOfValues(); i++) //For each of the points in the energy DataSet
    {
        double currentEnergy = energyDataSet->GetVariantValue(i).ToDouble(); //Current energy value tested
        if (currentEnergy < minimumEnergy)
        {
            minimumEnergy = currentEnergy;
        }
        if (currentEnergy > maxEnergy)
        {
            maxEnergy = currentEnergy;
        }
        
    }

    for (size_t i = 0; i < energyDataSet->GetNumberOfValues(); i++)
    {
        double currentEnergy = energyDataSet->GetVariantValue(i).ToDouble(); //Current energy value tested
        if (currentEnergy > 0.0)
        {
            double scaledValue = (currentEnergy * abs(minimumEnergy)) / maxEnergy; //Scale values
            energyDataSet->SetVariantValue(i,scaledValue);
        }
        
    }
    //periodGrid->Update();

    //Persistence Diagram of the data
    vtkSmartPointer<ttkPersistenceDiagram> persistenceDiagram = vtkSmartPointer<ttkPersistenceDiagram>::New();
    persistenceDiagram->SetUseAllCores(useAllCores);
    persistenceDiagram->SetInputConnection(periodGrid->GetOutputPort());
    persistenceDiagram->SetInputArrayToProcess(0,0,0,0,"Potential Energy");
    
    //We delete the persistence pairs corresponding to the graph diagonal
    vtkSmartPointer<vtkThreshold> criticalPairs = vtkSmartPointer<vtkThreshold>::New();
    criticalPairs->SetInputConnection(persistenceDiagram->GetOutputPort());
    criticalPairs->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"PairIdentifier");
    criticalPairs->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    criticalPairs->SetLowerThreshold(-0.1);
    criticalPairs->SetUpperThreshold(9e9);
    criticalPairs->Update();

    // //We delete the persistence pairs corresponding to the graph diagonal
    // vtkSmartPointer<vtkThreshold> most = vtkSmartPointer<vtkThreshold>::New();
    // most->SetInputConnection(persistenceDiagram->GetOutputPort());
    // most->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"Persistence");
    // most->ThresholdBetween(0.01*abs(minimumEnergy),9e19);
    // most->Update();

    auto persistenceDiagramDataSet = vtkDataSet::SafeDownCast(criticalPairs->GetOutputDataObject(0));

    std::string base_filename = BaseFileName.substr(BaseFileName.find_last_of("-") + 1);

    ofstream persistenceDiagrams;
    persistenceDiagrams.open("../Results/PersistenceDiagrams/" + base_filename + ".csv");
    assert(persistenceDiagrams.is_open());
    persistenceDiagrams << "Birth,Death" << "\n";

    for (size_t i = 1; i < persistenceDiagramDataSet->GetNumberOfPoints(); i=i+2)
    {
        double currentPointCoords[3];
        persistenceDiagramDataSet->GetPoint(i,currentPointCoords);
        persistenceDiagrams << currentPointCoords[0] << "," << currentPointCoords[1]  << "\n";
    }
    persistenceDiagrams.close();


}

/**
 * @brief Comput the persistence diagram of a energy field of an input material. Then, it computes a persistence simplification to get
 * the most important features. Finally, it writes its outputs to a .grid file
 *
 * @param grid Input grid from a reader function
 * @param useAllCores Use all available cores
 */
void segmentor::energyDiagrams2(vtkSmartPointer<vtkImageData> grid, bool useAllCores)
{
    logger::mainlog << "Energy Diagrams Module" << "\n";
    logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
    
    //Apply Periodic Boundary Conditions
    vtkSmartPointer<ttkTriangulationManager> periodGrid = vtkSmartPointer<ttkTriangulationManager>::New();;
    periodGrid->SetUseAllCores(useAllCores);
    periodGrid->SetInputData(grid);
    periodGrid->SetPeriodicity(true);
    periodGrid->Update();

    //Given grid dataSet
    auto inputGridDataSet = vtkDataSet::SafeDownCast(periodGrid->GetOutputDataObject(0));
    //Potential Energy DataSet of the input grid
    auto energyDataSet = inputGridDataSet->GetPointData()->GetAbstractArray("Potential Energy");

    double minimumEnergy = 0.0; //Minimum energy value of the material
    double maxEnergy = 0.0;
    for (size_t i = 0; i < energyDataSet->GetNumberOfValues(); i++) //For each of the points in the energy DataSet
    {
        double currentEnergy = energyDataSet->GetVariantValue(i).ToDouble(); //Current energy value tested
        if (currentEnergy < minimumEnergy)
        {
            minimumEnergy = currentEnergy;
        }
        if (currentEnergy > maxEnergy)
        {
            maxEnergy = currentEnergy;
        }
        
    }

    for (size_t i = 0; i < energyDataSet->GetNumberOfValues(); i++)
    {
        double currentEnergy = energyDataSet->GetVariantValue(i).ToDouble(); //Current energy value tested
        if (currentEnergy > 0.0)
        {
            double scaledValue = (currentEnergy * abs(minimumEnergy)) / maxEnergy; //Scale values
            energyDataSet->SetVariantValue(i,scaledValue);
        }
        
    }
    

    //Persistence Diagram of the data
    vtkSmartPointer<ttkPersistenceDiagram> persistenceDiagram = vtkSmartPointer<ttkPersistenceDiagram>::New();
    persistenceDiagram->SetUseAllCores(useAllCores);
    persistenceDiagram->SetInputConnection(periodGrid->GetOutputPort());
    persistenceDiagram->SetInputArrayToProcess(0,0,0,0,"Potential Energy");
    
    //We delete the persistence pairs corresponding to the graph diagonal
    vtkSmartPointer<vtkThreshold> criticalPairs = vtkSmartPointer<vtkThreshold>::New();
    criticalPairs->SetInputConnection(persistenceDiagram->GetOutputPort());
    criticalPairs->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"PairIdentifier");
    criticalPairs->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    criticalPairs->SetLowerThreshold(-0.1);
    criticalPairs->SetUpperThreshold(9e9);
    criticalPairs->Update();

    //We capture only the features corresponding to saddles-minimum
    vtkSmartPointer<vtkThreshold> saddlesMin = vtkSmartPointer<vtkThreshold>::New();
    saddlesMin->SetInputConnection(criticalPairs->GetOutputPort());
    saddlesMin->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"CriticalType");
    saddlesMin->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    saddlesMin->SetLowerThreshold(0);
    saddlesMin->SetUpperThreshold(2);
    saddlesMin->Update();

    //We delete the persistence pairs corresponding to the graph diagonal
    vtkSmartPointer<vtkThreshold> most = vtkSmartPointer<vtkThreshold>::New();
    most->SetInputConnection(saddlesMin->GetOutputPort());
    most->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"Persistence");
    most->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    most->SetLowerThreshold(0.01*abs(minimumEnergy));
    most->SetUpperThreshold(9e21);
    most->Update();

    //auto persistenceDiagramDataSet = vtkDataSet::SafeDownCast(criticalPairs->GetOutputDataObject(0));
    auto persistenceDiagramDataSet = vtkDataSet::SafeDownCast(most->GetOutputDataObject(0));


    std::string base_filename = BaseFileName.substr(BaseFileName.find_last_of("-") + 1);

    ofstream persistenceDiagrams;
    persistenceDiagrams.open("../Results/PersistenceDiagrams/" + base_filename + ".csv");
    assert(persistenceDiagrams.is_open());
    persistenceDiagrams << "Birth,Death" << "\n";

    for (size_t i = 1; i < persistenceDiagramDataSet->GetNumberOfPoints(); i=i+2)
    {
        double currentPointCoords[3];
        persistenceDiagramDataSet->GetPoint(i,currentPointCoords);


        if ((currentPointCoords[0] < 0) && (currentPointCoords[1] < 0))
        {
            persistenceDiagrams << currentPointCoords[0] << "," << currentPointCoords[1]  << "\n";

        }
        
    }
    persistenceDiagrams.close();


}


/**
 * @brief Comput the persistence diagram of a distance field of an input material. Then, it computes a persistence simplification to get
 * the most important features. Finally, it writes its outputs to a .grid file
 *
 * @param grid Input grid from a reader function
 * @param useAllCores Use all available cores
 */
void segmentor::distanceDiagrams2(vtkSmartPointer<vtkImageData> grid, bool useAllCores)
{
    logger::mainlog << "Distance Diagrams Module" << "\n";
    logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
    
    //Apply Periodic Boundary Conditions
    vtkSmartPointer<ttkTriangulationManager> periodGrid = vtkSmartPointer<ttkTriangulationManager>::New();;
    periodGrid->SetUseAllCores(useAllCores);
    periodGrid->SetInputData(grid);
    periodGrid->SetPeriodicity(true);
    periodGrid->Update();


    //Persistence Diagram of the data
    vtkSmartPointer<ttkPersistenceDiagram> persistenceDiagram = vtkSmartPointer<ttkPersistenceDiagram>::New();
    persistenceDiagram->SetUseAllCores(useAllCores);
    //persistenceDiagram->SetDebugLevel(3);
    persistenceDiagram->SetInputConnection(periodGrid->GetOutputPort());
    persistenceDiagram->SetInputArrayToProcess(0,0,0,0,"This is distance grid");
    persistenceDiagram->Update();

    vtkSmartPointer<vtkThreshold> criticalPairs = vtkSmartPointer<vtkThreshold>::New();
    criticalPairs->SetInputConnection(persistenceDiagram->GetOutputPort());
    criticalPairs->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"PairIdentifier");
    criticalPairs->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    criticalPairs->SetLowerThreshold(-0.1);
    criticalPairs->SetUpperThreshold(9e9);
    criticalPairs->Update();
    
    

    auto persistenceDiagramDataSet = vtkDataSet::SafeDownCast(criticalPairs->GetOutputDataObject(0));
    //auto persistenceDiagramDataSet = vtkDataSet::SafeDownCast(persistenceDiagram->GetOutputDataObject(0));


    std::string base_filename = BaseFileName.substr(BaseFileName.find_last_of("-") + 1);

    ofstream persistenceDiagrams;
    persistenceDiagrams.open("../Results/PersistenceDiagrams/" + base_filename + ".csv");
    assert(persistenceDiagrams.is_open());
    persistenceDiagrams << "Birth,Death" << "\n";

    for (size_t i = 1; i < persistenceDiagramDataSet->GetNumberOfPoints(); i=i+2)
    {
        double currentPointCoords[3];
        persistenceDiagramDataSet->GetPoint(i,currentPointCoords);


        
        
        persistenceDiagrams << currentPointCoords[0] << "," << currentPointCoords[1]  << "\n";

        
        
    }
    persistenceDiagrams.close();


}

/**
 * @brief Get information from the distance and energy grid files for the same material in order to gather them in a single grid
 *
 * @param distanceGrid Input distance grid from a reader function
 * @param energyFile Input energy file to merge information
 * @param persistenceThreshold Persistence Threshold to compute the persistence simplification
 * @param useAllCores Use all cores available
 * @return auto
 */
auto segmentor::energyIncluder(vtkSmartPointer<vtkImageData> distanceGrid,string energyFile,double persistenceThreshold, bool useAllCores)
{
    logger::mainlog << "Energy Includer Module" << "\n";
    logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
    
    vtkSmartPointer<vtkGaussianCubeReader2> energyReader = vtkSmartPointer<vtkGaussianCubeReader2>::New();
    energyReader->SetFileName(energyFile.data()); //Set the input file
    energyReader->Update();
    //Image data output from the Gaussian Cube file
    vtkSmartPointer<vtkImageData> energyGrid = vtkSmartPointer<vtkImageData>::New();
    energyGrid = energyReader->GetGridOutput();
    
    
    //Apply Periodic Boundary Conditions
    vtkSmartPointer<ttkTriangulationManager> periodGrid = vtkSmartPointer<ttkTriangulationManager>::New();;
    periodGrid->SetUseAllCores(useAllCores);
    periodGrid->SetInputData(energyGrid);
    periodGrid->SetPeriodicity(true);
    periodGrid->Update();

    //Given grid dataSet
    auto inputGridDataSet = vtkDataSet::SafeDownCast(periodGrid->GetOutputDataObject(0));
    //Potential Energy DataSet of the input grid
    auto energyDataSet = inputGridDataSet->GetPointData()->GetAbstractArray("Potential Energy");

    double minimumEnergy = 0.0; //Minimum energy value of the material
    double maxEnergy = 0.0;
    for (size_t i = 0; i < energyDataSet->GetNumberOfValues(); i++) //For each of the points in the energy DataSet
    {
        double currentEnergy = energyDataSet->GetVariantValue(i).ToDouble(); //Current energy value tested
        if (currentEnergy < minimumEnergy)
        {
            minimumEnergy = currentEnergy;
        }
        if (currentEnergy > maxEnergy)
        {
            maxEnergy = currentEnergy;
        }
        
    }

    for (size_t i = 0; i < energyDataSet->GetNumberOfValues(); i++)
    {
        double currentEnergy = energyDataSet->GetVariantValue(i).ToDouble(); //Current energy value tested
        if (currentEnergy > 0.0)
        {
            double scaledValue = (currentEnergy * abs(minimumEnergy)) / maxEnergy; //Scale values
            energyDataSet->SetVariantValue(i,scaledValue);
        }
        
    }
    
    double energyPersistence = persistenceThreshold*abs(minimumEnergy);

    //Persistence Diagram of the data
    vtkSmartPointer<ttkPersistenceDiagram> persistenceDiagram = vtkSmartPointer<ttkPersistenceDiagram>::New();
    persistenceDiagram->SetUseAllCores(useAllCores);
    persistenceDiagram->SetInputConnection(periodGrid->GetOutputPort());
    persistenceDiagram->SetInputArrayToProcess(0,0,0,0,"Potential Energy");
    
    //We delete the persistence pairs corresponding to the graph diagonal
    vtkSmartPointer<vtkThreshold> criticalPairs = vtkSmartPointer<vtkThreshold>::New();
    criticalPairs->SetInputConnection(persistenceDiagram->GetOutputPort());
    criticalPairs->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"PairIdentifier");
    criticalPairs->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    criticalPairs->SetLowerThreshold(-0.1);
    criticalPairs->SetUpperThreshold(9e9);
    criticalPairs->Update();

    //Persistence threshold for future simplifications
    vtkSmartPointer<vtkThreshold> persistentPairs = vtkSmartPointer<vtkThreshold>::New();
    persistentPairs->SetInputConnection(criticalPairs->GetOutputPort());
    persistentPairs->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "Persistence");
    persistentPairs->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    persistentPairs->SetLowerThreshold(energyPersistence);
    persistentPairs->SetUpperThreshold(9.0e21);
    

    //Topological simplification from the persistence results
    vtkSmartPointer<ttkTopologicalSimplification> topologicalSimplification = vtkSmartPointer<ttkTopologicalSimplification>::New();
    //topologicalSimplification->SetDebugLevel(4);
    topologicalSimplification->SetUseAllCores(useAllCores);
    //topologicalSimplification->SetInputData(grid);
    topologicalSimplification->SetInputConnection(0,periodGrid->GetOutputPort());
    topologicalSimplification->SetInputArrayToProcess(0,0,0, 0,"Potential Energy");
    topologicalSimplification->SetInputConnection(1, persistentPairs->GetOutputPort());
    topologicalSimplification->Update();

    auto energyArray = vtkDataSet::SafeDownCast(topologicalSimplification->GetOutputDataObject(0))->GetPointData()->GetAbstractArray("Potential Energy");

    distanceGrid->GetPointData()->AddArray(energyArray);

    return distanceGrid;


}

/**
 * @brief From a materials name, it finds its grids in the given input directories in order to find the energy and distance values of each point in the grid.
 * Then, it writes an output file with these values
 *
 * @param inputFile
 * @param distanceDirectory
 * @param energyDirectory
 * @param distanceSimp
 * @param energySimp
 */
void segmentor::energyVsDistance(string inputFile, string distanceDirectory, string energyDirectory)
{
    logger::mainlog << "Energy vs Distance Module" << "\n";
    logger::mainlog << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
    
    int state = system("mkdir -p ../Results/GridFiles"); //Create a directory to save the results
    
    
    //Input materials' distance grid
    string distanceFile = distanceDirectory + inputFile;
    //Input materials' energy grid
    string energyFile = energyDirectory + inputFile;

    std::string base_filename = inputFile.substr(inputFile.find_last_of("-") + 1);
    std::string::size_type const p(base_filename.find_last_of('.'));
    std::string file_without_extension = base_filename.substr(0, p); //Get the material filename without extension nor directory
    
    BaseFileName = file_without_extension;

    //Read the distance grid
    vtkSmartPointer<vtkGaussianCubeReader2> distanceReader = vtkSmartPointer<vtkGaussianCubeReader2>::New();
    distanceReader->SetFileName(distanceFile.data()); //Set the input file
    distanceReader->Update();
    //Image data output from the Gaussian Cube file
    vtkSmartPointer<vtkImageData> distanceGrid = vtkSmartPointer<vtkImageData>::New();
    distanceGrid = distanceReader->GetGridOutput();

    //Read the energy grid
    vtkSmartPointer<vtkGaussianCubeReader2> energyReader = vtkSmartPointer<vtkGaussianCubeReader2>::New();
    energyReader->SetFileName(energyFile.data()); //Set the input file
    energyReader->Update();
    //Image data output from the Gaussian Cube file
    vtkSmartPointer<vtkImageData> energyGrid = vtkSmartPointer<vtkImageData>::New();
    energyGrid = energyReader->GetGridOutput();




    //---------------------------------------------------------------------------------------------

    

    auto distanceDataSet = distanceGrid->GetPointData()->GetAbstractArray("This is distance grid");
    auto energyDataSet = energyGrid->GetPointData()->GetAbstractArray("Potential Energy");

    ofstream gridFiles;
    gridFiles.open("../Results/GridFiles/" + BaseFileName + ".csv");
    gridFiles << "Energy,Distance" << endl;

    for (size_t i = 0; i < energyGrid->GetNumberOfPoints(); i++)
    {
        //Save the coords of the current point in the energy grid
        double energyPointCoords[3];
        energyGrid->GetPoint(i,energyPointCoords);

        //ID of the closest point in the distance grid of the same material
        int closestDistancePoint = distanceGrid->FindPoint(energyPointCoords);

        gridFiles << energyDataSet->GetVariantValue(i).ToDouble() << "," << distanceDataSet->GetVariantValue(closestDistancePoint).ToDouble() << endl;


    }

    gridFiles.close();
    








}

/**
 * @brief Given a list of persistence Diagrams. Compute the Wasserstein Distance for each of them to fin
 * the closest ones
 *
 * @param fileList  Vector storing the names of the Persistence Diagram files
 * @return auto
 */
auto segmentor::persistenceMatchings(bool useAllCores)
{
    //We create a directory to save the results
    std::string func = "mkdir -p "; //Commands needed to create the folder
    std::string directory = "../Results/VoidsComparative"; //Directory
    //directory.append(file_without_extension);
    func.append(directory);
    int status;
    status = system(func.c_str()); // Creating a directory
    if (status == -1)
    {
        cerr << "Error : " << strerror(errno) << endl;
    }
        
    else
    {
        logger::mainlog << "Voids Comparative Directories created" << endl;
    }

    system("rm -f diagramsList.txt"); //Delete previous one if existed

    auto function = system("ls ../Results/PersistenceDiagrams >> diagramsList.txt");

    ifstream myFile("diagramsList.txt");
    string line;
    vector<string> files;
    if (myFile.is_open())
    {
        while ( getline (myFile,line) )
        {
            files.push_back(line);
        }
        myFile.close();
    }
    logger::mainlog << files.size() << endl;


    #pragma omp parallel for
    for (size_t i = 0; i < files.size(); i++)
    {
        
        
        //We get the input file name
        std::string base_filename = files[i];
        std::string::size_type const p(base_filename.find_last_of('.'));
        //Input File name
        std::string file_without_extension = base_filename.substr(0, p);

        logger::mainlog << file_without_extension << endl;

        ofstream comparativeData; //Stream to write the comparatives of each region
        comparativeData.open("../Results/VoidsComparative/" + file_without_extension + ".csv"); //Opening the writer
        assert(comparativeData.is_open());
        comparativeData << "RegionToCompare,WassersteinDistance" << "\n";

        //Persistence Diagram Currently Analysing
        vtkSmartPointer<vtkUnstructuredGridReader> persistenceDiagram1 = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        persistenceDiagram1->SetFileName(("../Results/PersistenceDiagrams/"+files[i]).c_str());
        persistenceDiagram1->Update();


        for (size_t j = 0; j < files.size(); j++)
        {
            //Persistence Diagram to compare
            vtkSmartPointer<vtkUnstructuredGridReader> persistenceDiagram2 = vtkSmartPointer<vtkUnstructuredGridReader>::New();
            persistenceDiagram2->SetFileName(("../Results/PersistenceDiagrams/"+files[j]).c_str());
            persistenceDiagram2->Update();

            vtkSmartPointer<vtkMultiBlockDataGroupFilter> group = vtkSmartPointer<vtkMultiBlockDataGroupFilter>::New();
            group->AddInputConnection(persistenceDiagram1->GetOutputPort());
            group->AddInputConnection(persistenceDiagram2->GetOutputPort());
            group->Update();

            vtkSmartPointer<ttkPersistenceDiagramClustering> persistenceCluster = vtkSmartPointer<ttkPersistenceDiagramClustering>::New();
            persistenceCluster->SetUseAllCores(useAllCores);
            //persistenceCluster->SetDebugLevel(3);
            persistenceCluster->SetInputConnection(group->GetOutputPort());
            persistenceCluster->SetNumberOfClusters(1);
            persistenceCluster->Update();

            auto persistenceClusterDataSet = vtkMultiBlockDataSet::SafeDownCast(persistenceCluster->GetOutputDataObject(2));
            auto wassersteinDistance =persistenceClusterDataSet->GetBlock(0)->GetFieldData()->GetAbstractArray("WassersteinDistance")->GetVariantValue(0).ToDouble();

            comparativeData << files[j] << "," << wassersteinDistance << "\n";
        }

        comparativeData.close();

        
    }

    
    


    // #pragma omp parallel for
    // for(auto file1: fileList)
    // {
        
    //     //We get the input file name
    //     std::string base_filename = file1;
    //     std::string::size_type const p(base_filename.find_last_of('.'));
    //     //Input File name
    //     std::string file_without_extension = base_filename.substr(0, p);

    //     ofstream comparativeData; //Stream to write the comparatives of each region
    //     comparativeData.open("../Results/VoidsComparative/" + file_without_extension + ".csv"); //Opening the writer
    //     assert(comparativeData.is_open());
    //     comparativeData << "RegionToCompare,WassersteinDistance" << "\n";
        
    //     //Persistence Diagram Currently Analysing
    //     vtkSmartPointer<vtkUnstructuredGridReader> persistenceDiagram1 = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    //     persistenceDiagram1->SetFileName(("../Results/PersistenceDiagrams/"+file1).c_str());
    //     persistenceDiagram1->Update();

        

    //     for(auto file2:fileList)
    //     {
    //         //Persistence Diagram to compare
    //         vtkSmartPointer<vtkUnstructuredGridReader> persistenceDiagram2 = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    //         persistenceDiagram2->SetFileName(("../Results/PersistenceDiagrams/"+file2).c_str());
    //         //persistenceDiagram1->SetFileName("../PersistenceDiagrams/FAU_11.vtk");
    //         persistenceDiagram2->Update();

    //         vtkSmartPointer<vtkMultiBlockDataGroupFilter> group = vtkSmartPointer<vtkMultiBlockDataGroupFilter>::New();
    //         group->AddInputConnection(persistenceDiagram1->GetOutputPort());
    //         group->AddInputConnection(persistenceDiagram2->GetOutputPort());
    //         group->Update();

    //         vtkSmartPointer<ttkPersistenceDiagramClustering> persistenceCluster = vtkSmartPointer<ttkPersistenceDiagramClustering>::New();
    //         persistenceCluster->SetUseAllCores(1);
    //         persistenceCluster->SetDebugLevel(3);
    //         persistenceCluster->SetInputConnection(group->GetOutputPort());
    //         persistenceCluster->SetNumberOfClusters(1);
    //         persistenceCluster->Update();

    //         auto persistenceClusterDataSet = vtkMultiBlockDataSet::SafeDownCast(persistenceCluster->GetOutputDataObject(2));
    //         auto wassersteinDistance =persistenceClusterDataSet->GetBlock(0)->GetFieldData()->GetAbstractArray("WassersteinDistance")->GetVariantValue(0).ToDouble();

    //         comparativeData << file2 << "," << wassersteinDistance << "\n";
    //     }
    //     comparativeData.close();
    // }
}

auto segmentor::persistenceDiagramsWriter(bool useAllCores)
{
    //We create a directory to save the results
    std::string func = "mkdir -p "; //Commands needed to create the folder
    std::string directory = "../Results/DiagramsGrids"; //Directory
    //directory.append(file_without_extension);
    func.append(directory);
    int status;
    status = system(func.c_str()); // Creating a directory
    if (status == -1)
    {
        cerr << "Error : " << strerror(errno) << endl;
    }
        
    else
    {
        logger::mainlog << "Directories created" << endl;
    }

    system("rm -f diagramsList.txt"); //Delete previous one if existed

    auto function = system("ls ../Results/PersistenceDiagrams >> diagramsList.txt");

    ifstream myFile("diagramsList.txt");
    string line;
    vector<string> files;
    if (myFile.is_open())
    {
        while ( getline (myFile,line) )
        {
            files.push_back(line);
        }
        myFile.close();
    }
    logger::mainlog << files.size() << endl;

    #pragma omp parallel for
    for (size_t i = 0; i < files.size(); i++)
    {
        //We get the input file name
        std::string base_filename = files[i];
        std::string::size_type const p(base_filename.find_last_of('.'));
        //Input File name
        std::string file_without_extension = base_filename.substr(0, p);

        logger::mainlog << file_without_extension << endl;

        ofstream comparativeData; //Stream to write the comparatives of each region
        comparativeData.open("../Results/DiagramsGrids/" + file_without_extension + ".csv"); //Opening the writer
        assert(comparativeData.is_open());
        comparativeData << "Birth,Death" << "\n";

        //Persistence Diagram Currently Analysing
        vtkSmartPointer<vtkUnstructuredGridReader> persistenceDiagram1 = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        persistenceDiagram1->SetFileName(("../Results/PersistenceDiagrams/"+files[i]).c_str());
        persistenceDiagram1->Update();

        vtkSmartPointer<vtkThreshold> criticalPairs = vtkSmartPointer<vtkThreshold>::New();
        criticalPairs->SetInputConnection(persistenceDiagram1->GetOutputPort());
        criticalPairs->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"PairIdentifier");
        criticalPairs->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
        criticalPairs->SetLowerThreshold(-0.1);
        criticalPairs->SetUpperThreshold(9e9);
        criticalPairs->Update();

        auto persistenceDiagramDataSet = vtkDataSet::SafeDownCast(criticalPairs->GetOutputDataObject(0));

        for (size_t j = 1; j < persistenceDiagramDataSet->GetNumberOfPoints(); j=j+2)
        {
            double currentPointCoords[3];
            persistenceDiagramDataSet->GetPoint(j,currentPointCoords);


            
            
            comparativeData << currentPointCoords[0] << "," << currentPointCoords[1]  << "\n";

            
            
        }
        comparativeData.close();
    }
      
}


