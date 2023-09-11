/*********************************************************************

TDA-Segmentor    A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:         Aditya Vasudevan (adityavv.iitkgp@gmail.com)
                 Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
                 Maciek Haranczyk (maciej.haranczyk@imdea.org)

                 IMDEA Materiales Institute
 
**********************************************************************/

#include "segmentor.h"


segmentor::segmentor(const parameters &p):

isPersistenceDiagramComputed(false),
theMSC(vtkSmartPointer<ttkMorseSmaleComplex>::New()),
thePersistenceDiagram(vtkSmartPointer<ttkPersistenceDiagram>::New()),
thePersistenceCurve(vtkSmartPointer<ttkPersistenceCurve>::New())

{

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
    
    
}

/**
 * @brief segmentor Destructor class
 *
 */
segmentor::~segmentor()
{
    delete Grid;
    logger::mainlog << "segmentor: Closing class" << "\n";

}


/**
 * @brief Get the  grid(vtkImageData) from a Gaussian Cube file(.cube)
 *
 * @param inputFilePath Input File Path of the Nanoporous Material
 * @param gridResolution Grid resolution of the Gaussian Cube file
 * @param writeGridFile Write the results to an output file(OPTIONAL)
 * @return auto Grid from the Gaussian Cube file
 */
grid*  segmentor::readInputFile(const parameters &p, bool writeGridFile)
{
    logger::mainlog << "\n\nSegmentor: Reader Module" << "\n";
    logger::mainlog << "Reading " << fileName << endl;
    ttk::Timer readerTime;
    
    // We first find the array name based on the input file
    // and create either a distance grid or a potential energy grid
    getArrayName(arrayName);
    
    if (arrayName == "This is distance grid"){
        Grid = new distanceGrid(p);
    } else if (arrayName == "Potential Energy"){
        Grid = new PEgrid(p);
    }
    
    // Based on the extension, the properties of the grid are read and stored.
    
    if (extensionName == ".cube"){
        
        vtkSmartPointer<vtkGaussianCubeReader2> cubeReader = vtkSmartPointer<vtkGaussianCubeReader2>::New();
        cubeReader->SetFileName(fileName.c_str());
        cubeReader->Update();
        Grid->cubicGrid = cubeReader->GetGridOutput();
        
        getGridResolutionFromCubeFile(Grid->gridResolution);

        vtkIdType cellDims[3];
        Grid->cubicGrid->GetDimensions(cellDims);
        Grid->cubicGrid->SetSpacing(1.0/(cellDims[0]-1), 1.0/(cellDims[1]-1), 1.0/(cellDims[2]-1));
        
    } else if (extensionName == ".vti"){ // This .vti extension is for the input file that has the stresses from LAMMPS
                                        // (see postProcessing/generate-stress-grids.py)
        
        vtkSmartPointer<vtkXMLImageDataReader> dataReader = vtkSmartPointer<vtkXMLImageDataReader>::New();
        dataReader->SetFileName((fileName).c_str());
        dataReader->Update();
        Grid->cubicGrid = dataReader->GetOutput();
        double imageGridRes[3];
        Grid->cubicGrid->GetSpacing(imageGridRes);
        Grid->gridResolution[0][0] = imageGridRes[0];
        Grid->gridResolution[1][1] = imageGridRes[1];
        Grid->gridResolution[2][2] = imageGridRes[2];
        
    }
    
    
    vtkIdType cellDims[3];
    Grid->cubicGrid->GetDimensions(cellDims);
    Grid->nx = cellDims[0]; Grid->ny = cellDims[1]; Grid->nz = cellDims[2];
    
    Grid->defineUnitCellVectors();
    
    vtkIdType numberOfCells = Grid->cubicGrid->GetNumberOfCells();
    double unitCellVolume = determinant(Grid->gridResolution);
    double volume = unitCellVolume * numberOfCells;
    
    logger::mainlog << "Grid Vector (X)             : (" << Grid->gridResolution[0][0] << ", " << Grid->gridResolution[1][0] << ", " << Grid->gridResolution[2][0] << ")\n";
    logger::mainlog << "Grid Vector (Y)             : (" << Grid->gridResolution[0][1] << ", " << Grid->gridResolution[1][1] << ", " << Grid->gridResolution[2][1] << ")\n";
    logger::mainlog << "Grid Vector (Z)             : (" << Grid->gridResolution[0][2] << ", " << Grid->gridResolution[1][2] << ", " << Grid->gridResolution[2][2] << ")\n";
    logger::mainlog << "Number of points in the grid        : (" << cellDims[0] << " X " << cellDims[1] << " X "<< cellDims[2] << ")" << endl;
    logger::mainlog << "Volume of the unit cell             : "  << volume << " (A^o)^3" << endl;
    logger::mainlog << "Unit Cell Vector a :   " << Grid->unitCellVectors[0][0] << "    " << Grid->unitCellVectors[1][0] << "    " << Grid->unitCellVectors[2][0] << "\n";
    logger::mainlog << "Unit Cell Vector b :   " << Grid->unitCellVectors[0][1] << "    " << Grid->unitCellVectors[1][1] << "    " << Grid->unitCellVectors[2][1] << "\n";
    logger::mainlog << "Unit Cell Vector c :   " << Grid->unitCellVectors[0][2] << "    " << Grid->unitCellVectors[1][2] << "    " << Grid->unitCellVectors[2][2] << "\n";

    
    // Store all the locations of the grid points for a general triclinic lattice.
    double x = 0.0, y = 0.0, z = 0.0;
    for (unsigned int k = 0; k < cellDims[2]; k++){
        for (unsigned int j = 0; j < cellDims[1]; j++){
            for (unsigned int i = 0; i < cellDims[0]; i++){
                
              x = Grid->gridResolution[0][0]*i + Grid->gridResolution[0][1]*j + Grid->gridResolution[0][2] * k;
              y = Grid->gridResolution[1][0]*i + Grid->gridResolution[1][1]*j + Grid->gridResolution[1][2] * k;
              z = Grid->gridResolution[2][0]*i + Grid->gridResolution[2][1]*j + Grid->gridResolution[2][2] * k;
              Grid->gridPointsXYZ->InsertNextPoint(x, y, z);
          }
        }
      }

    
    if (writeGridFile == true)
    {
        vtkNew<vtkDoubleArray> pointValues;
        pointValues->SetNumberOfComponents(1);
        pointValues->SetNumberOfTuples(Grid->nx*Grid->ny*Grid->nz);
        
        for (size_t i = 0; i < (Grid->nx*Grid->ny*Grid->nz); ++i)
        {
          pointValues->SetValue(i, Grid->cubicGrid->GetPointData()->GetArray(arrayName.c_str())->GetVariantValue(i).ToDouble());
        }

        Grid->originalGrid->SetDimensions(static_cast<int>(Grid->nx), static_cast<int>(Grid->ny),
                                      static_cast<int>(Grid->nz));
        Grid->originalGrid->SetPoints(Grid->gridPointsXYZ);
        Grid->originalGrid->GetPointData()->SetScalars(pointValues);

       vtkNew<vtkStructuredGridWriter> strucGridWriter;
       strucGridWriter->SetInputData(Grid->originalGrid);
       strucGridWriter->SetFileName((Directory+"/"+BaseFileName+"_grid.vtk").c_str());
       strucGridWriter->Write();
        
    }
 
    double elapsedTime = readerTime.getElapsedTime();
    logger::mainlog << "Time elapsed in the reader module: " << elapsedTime << "(s)" << endl;
    
    return Grid;
}




void segmentor::getGridResolutionFromCubeFile(double  GridResolution[3][3]) {
    
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




void segmentor::getArrayName(std::string &nameOfArray){
    
    if (extensionName == ".cube"){
        
        if (nameOfArray.empty()) getArrayNameFromCubeFile(nameOfArray);
        logger::mainlog << "Array that is going to be used for TDA analysis: " << nameOfArray << endl;
        
    } else if (extensionName == ".vti"){
        
        if (nameOfArray.empty()) getArrayNameFromVTIFile(nameOfArray);
        logger::mainlog << "Array that is going to be used for TDA analysis: " << nameOfArray << endl;
        
    } else {
        
        logger::mainlog << "Extension type of the input file is neither .cube nor .vti" <<endl;
        logger::errlog << "Extension type of the input file neither .cube nor .vti" <<endl;
        exit(0);
    }
    
    
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
            if (lineNumber > 1) break;
        }
                    
    }
        
}




void segmentor::getArrayNameFromVTIFile(std::string &nameOfArray){
    
    vtkSmartPointer<vtkXMLImageDataReader> dataReader = vtkSmartPointer<vtkXMLImageDataReader>::New();
    dataReader->SetFileName((fileName).c_str());
    dataReader->Update();
    vtkSmartPointer<vtkImageData> data = vtkSmartPointer<vtkImageData>::New();
    data = dataReader->GetOutput();
    
    // The name of the first array in the grid is used as the array for segmentation.
    if (nameOfArray.empty()) {
        char * name = data->GetPointData()->GetAbstractArray(0)->GetName();
        nameOfArray = name;
    }
        
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
vtkSmartPointer<ttkTriangulationManager> segmentor::generatePeriodicGrid(vtkSmartPointer<vtkImageData> grid,bool periodicConditions, bool useAllCores)
{
    
    
    logger::mainlog << "\n\nSegmentor: generate periodic grid Module" << "\n";
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




void segmentor::computePersistenceDiagram(vtkSmartPointer<ttkTriangulationManager> grid, bool useAllCores){

    thePersistenceDiagram->SetUseAllCores(useAllCores);
    thePersistenceDiagram->SetInputConnection(grid->GetOutputPort());
    thePersistenceDiagram->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    
    isPersistenceDiagramComputed = true;
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
void segmentor::MSC(vtkSmartPointer<ttkTriangulationManager> grid,double persistenceThreshold, double saddlesaddleIncrement, bool writeOutputs, bool useAllCores)
{
    logger::mainlog << "\n\nSegmentor: Morse Smale Complex Module" << "\n" << flush;
    
    ttk::Timer MSCTimer;
    
    if (!isPersistenceDiagramComputed) computePersistenceDiagram(grid, useAllCores); // Persistence Diagram is stored in the variable thePersistenceDiagram of the class.
    
    //We delete the persistence pairs corresponding to the graph diagonal
    vtkSmartPointer<vtkThreshold> criticalPairs = vtkSmartPointer<vtkThreshold>::New();
    criticalPairs->SetInputConnection(thePersistenceDiagram->GetOutputPort());
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
    
    // If persistenceThreshold is not provided as input, then 1% of max is automatically taken.
    if (persistenceThreshold == 0.0) {
        persistenceThreshold = 0.01 * maximumPersistence;
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
    topologicalSimplification->SetUseAllCores(useAllCores);
    topologicalSimplification->SetInputConnection(0,grid->GetOutputPort());
    topologicalSimplification->SetInputArrayToProcess(0,0,0, 0,arrayName.c_str());
    topologicalSimplification->SetInputConnection(1, persistentPairs->GetOutputPort());

    //=============================================================================================
    //=============================================================================================
    //3.3 Morse Smale Complex Computation
    //Morse Smale Complex Computation
    theMSC->SetUseAllCores(useAllCores);
    theMSC->SetReturnSaddleConnectors(1);
    theMSC->SetSaddleConnectorsPersistenceThreshold(saddlesaddleIncrement*persistenceThreshold);
    theMSC->SetInputConnection(topologicalSimplification->GetOutputPort());
    theMSC->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    theMSC->SetComputeSaddleConnectors(false);
    theMSC->SetComputeAscendingSeparatrices1(false);
    theMSC->SetComputeAscendingSeparatrices2(false);
    theMSC->SetComputeDescendingSeparatrices1(false);
    theMSC->SetComputeDescendingSeparatrices2(false);
    theMSC->Update();
    
    double timeTakenForMSC = MSCTimer.getElapsedTime();
    MSCTimer.reStart();
    logger::mainlog << "Time taken for MSC creation: " << timeTakenForMSC << "(s)" << endl;
    
    vtkIdType numberOfDescendingManifolds = getNumberOfDescendingManifolds(theMSC);
    vtkIdType numberOfAscendingManifolds = getNumberOfAscendingManifolds(theMSC);
    
    logger::mainlog << "Total number of descending manifolds (typically solid segments) : " << numberOfDescendingManifolds << endl;
    logger::mainlog << "Total number of ascending manifolds (typically void segments) : " << numberOfAscendingManifolds << endl;

    if (writeOutputs)
    {
        // We write the critical points for a triclinic lattice:
        auto segmentationDataSet = vtkDataSet::SafeDownCast(theMSC->GetOutputDataObject(3));
        size_t numberOfArrays = segmentationDataSet->GetPointData()->GetNumberOfArrays();
        Grid->segmentation->SetDimensions((int) Grid->nx, (int) Grid->ny,(int) Grid->nz);
        Grid->segmentation->SetPoints(Grid->gridPointsXYZ);
        
        for (size_t i = 0; i < numberOfArrays; i++){
            vtkNew<vtkDoubleArray> pointValues;
            char * name = segmentationDataSet->GetPointData()->GetAbstractArray((int)i)->GetName();
            pointValues->SetName(name);
            pointValues->SetNumberOfComponents(1);
            pointValues->SetNumberOfTuples(Grid->nx*Grid->ny*Grid->nz);
            vtkIdType numberOfPoints = segmentationDataSet->GetNumberOfPoints();
            for (size_t j = 0; j < segmentationDataSet->GetNumberOfPoints(); j++){
                
                pointValues->SetValue(j, segmentationDataSet->GetPointData()->GetArray(name)->GetVariantValue(j).ToDouble());
                
            }
            
            Grid->segmentation->GetPointData()->AddArray(pointValues);
        }
        
        vtkNew<vtkStructuredGridWriter> segmentationWriter;
        segmentationWriter->SetInputData(Grid->segmentation);
        segmentationWriter->SetFileName((Directory+"/"+BaseFileName+"_Segmentation.vtk").c_str());
        segmentationWriter->Write();
        
        auto criticalPointsDataSet = vtkDataSet::SafeDownCast(theMSC->GetOutputDataObject(0));
        numberOfArrays = criticalPointsDataSet->GetPointData()->GetNumberOfArrays();
        
        vtkNew<vtkPoints> cPoints;
        for (size_t i = 0; i < criticalPointsDataSet->GetNumberOfPoints(); i++){
            double coordABC[3];
            criticalPointsDataSet->GetPoint(i,coordABC);
            double coordXYZ[3];
            abcToxyz(coordABC,coordXYZ,Grid->unitCellVectors);
            cPoints->InsertNextPoint(coordXYZ[0], coordXYZ[1], coordXYZ[2]);
            
        }
        
        Grid->criticalPoints->SetPoints(cPoints);
        
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
            
            Grid->criticalPoints->GetPointData()->AddArray(pointValues);
        }
        
        vtkSmartPointer<vtkPolyDataWriter> critPointWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
        critPointWriter->SetInputData(Grid->criticalPoints);
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
    
    
}




/*auto segmentor::ftmtree(vtkSmartPointer<ttkTriangulationManager> grid, double persistenceThreshold, bool useAllCores)
{
    
    logger::mainlog << "\n\nSegmentor: FTM tree module" << "\n" << flush;
    
    ttk::Timer graphTimer;
    //Persistence Diagram of the data
    if (!isPersistenceDiagramComputed) computePersistenceDiagram(grid);
    
    //We delete the persistence pairs corresponding to the graph diagonal
    vtkSmartPointer<vtkThreshold> criticalPairs = vtkSmartPointer<vtkThreshold>::New();
    criticalPairs->SetInputConnection(thePersistenceDiagram->GetOutputPort());
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

} */



/*void segmentor::accessibleVoidGraph(vtkSmartPointer <ttkFTMTree> ftmTree, double moleculeRadius, bool useAllCores){
    
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
    
}*/



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




void segmentor::persistencecurve(vtkSmartPointer<ttkTriangulationManager> grid, bool useAllCores)
{
    ttk::Timer percurveTimer;
    logger::mainlog << "\nSegmentor: Persistence Curve module" << "\n" << flush;
 
    // To create the persistence curve, first the persistence diagram is calculated 
    // which serves as input to the persistence curve
    if (!isPersistenceDiagramComputed) computePersistenceDiagram(grid, useAllCores);

    thePersistenceCurve->SetInputConnection(thePersistenceDiagram->GetOutputPort());
    thePersistenceCurve->Update();
    
    /* Write the persistence curve in VTK format */
    vtkNew<vtkTableWriter> curveWriter0{};
    curveWriter0->SetInputConnection(thePersistenceCurve->GetOutputPort(0));
    curveWriter0->SetFileName((Directory+"/" + BaseFileName+"_minSaddlePairs.vtk").c_str());
    curveWriter0->Write();
    
    vtkNew<vtkTableWriter> curveWriter1{};
    curveWriter1->SetInputConnection(thePersistenceCurve->GetOutputPort(1));
    curveWriter1->SetFileName((Directory+"/" + BaseFileName+"_saddleSaddlePairs.vtk").c_str());
    curveWriter1->Write();
    
    vtkNew<vtkTableWriter> curveWriter2{};
    curveWriter2->SetInputConnection(thePersistenceCurve->GetOutputPort(2));
    curveWriter2->SetFileName((Directory+"/" + BaseFileName+"_SaddleMaxPairs.vtk").c_str());
    curveWriter2->Write();
    
    vtkNew<vtkTableWriter> curveWriter3{};
    curveWriter3->SetInputConnection(thePersistenceCurve->GetOutputPort(3));
    curveWriter3->SetFileName((Directory+"/" + BaseFileName+"_allPairs.vtk").c_str());
    curveWriter3->Write();
    
    
    double timeTakenForPerCurve = percurveTimer.getElapsedTime();
    logger::mainlog << "Time taken in persistence curve module: " << timeTakenForPerCurve << "(s)" << endl;

}

