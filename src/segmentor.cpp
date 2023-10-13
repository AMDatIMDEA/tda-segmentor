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
thePersistenceCurve(vtkSmartPointer<ttkPersistenceCurve>::New()),
Grid(nullptr)

{

    fileName = p.inputfilename;
    BaseFileName = p.basefilename;
    extensionName = p.extensionname;
    arrayName = p.arrayName; 
    writeFractionalGrid = p.writeFractionalGrid; 
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




segmentor::~segmentor()
{
    delete Grid;
    logger::mainlog << "segmentor: Closing class" << "\n";

}


/**
 * @brief This reads the input file and generates the grid either of type distanceGrid or PEgrid.
 *  This function also sets the variables of the grid class, such as gridresolution, number of points,
 *  vtk grids, stored both in fractional (cubicGrid) and actual (originalGrid) coordinates.
 *
 * @param p  the parameter class variable that contains all the input parameters.
 * @param writeGridFile writes the results of the input grid.
 *
 * @return an instance of the grid.
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
        
        Grid->cubicGrid = readFromCubeFile();
        // the default spacing is a voxel, and is normalized with the size in each direction to
        // fit the grid on a unit cube.
        Grid->cubicGrid->SetSpacing(1.0/(Grid->nx-1), 1.0/(Grid->ny-1), 1.0/(Grid->nz-1));
        
    } else if (extensionName == ".vti"){ // This .vti extension is for the input file that has the stresses from LAMMPS
                                        // (see postProcessing/generate-stress-grids.py on how to generate them)
        
        vtkSmartPointer<vtkXMLImageDataReader> dataReader = vtkSmartPointer<vtkXMLImageDataReader>::New();
        dataReader->SetFileName((fileName).c_str());
        dataReader->Update();
        Grid->cubicGrid = dataReader->GetOutput();

        double imageGridRes[3];
        Grid->cubicGrid->GetSpacing(imageGridRes);
        Grid->gridResolution[0][0] = imageGridRes[0];
        Grid->gridResolution[1][1] = imageGridRes[1];
        Grid->gridResolution[2][2] = imageGridRes[2];
        vtkIdType cellDims[3];
        Grid->cubicGrid->GetDimensions(cellDims);
        Grid->nx = cellDims[0]; Grid->ny = cellDims[1]; Grid->nz = cellDims[2];
        Grid->cubicGrid->SetSpacing(1.0/(Grid->nx-1), 1.0/(Grid->ny-1), 1.0/(Grid->nz-1));
        
    }
    
    Grid->defineUnitCellVectors();
    
    vtkIdType numberOfCells = Grid->cubicGrid->GetNumberOfCells();
    double dV = determinant(Grid->gridResolution);//the volume of the voxel grid 
    double volume = dV * numberOfCells;
    
    logger::mainlog << "Grid Vector (X)                                              : (" << Grid->gridResolution[0][0] << ", "
                                                                                          << Grid->gridResolution[1][0] << ", "
                                                                                          << Grid->gridResolution[2][0] << ")\n";
    
    logger::mainlog << "Grid Vector (Y)                                              : (" << Grid->gridResolution[0][1] << ", "
                                                                                          << Grid->gridResolution[1][1] << ", "
                                                                                          << Grid->gridResolution[2][1] << ")\n";
    
    logger::mainlog << "Grid Vector (Z)                                              : (" << Grid->gridResolution[0][2] << ", "
                                                                                          << Grid->gridResolution[1][2] << ", "
                                                                                          << Grid->gridResolution[2][2] << ")\n";
    
    logger::mainlog << "Number of points in the grid                                 : (" << Grid->nx << " X " << Grid->ny << " X "<< Grid->nz << ")" << endl;
    logger::mainlog << "Volume of the unit cell                                      : "  << volume << " (A^o)^3" << endl;
    
    logger::mainlog << "Unit Cell Vector a                                           :   " << Grid->unitCellVectors[0][0] << "    "
                                                                                           << Grid->unitCellVectors[1][0] << "    "
                                                                                           << Grid->unitCellVectors[2][0] << "\n";
    logger::mainlog << "Unit Cell Vector b                                           :   " << Grid->unitCellVectors[0][1] << "    "
                                                                                           << Grid->unitCellVectors[1][1] << "    "
                                                                                           << Grid->unitCellVectors[2][1] << "\n";
    logger::mainlog << "Unit Cell Vector c                                           :   " << Grid->unitCellVectors[0][2] << "    "
                                                                                           << Grid->unitCellVectors[1][2] << "    "
                                                                                           << Grid->unitCellVectors[2][2] << "\n";

    
    // Store all the locations of the grid points for a general triclinic lattice.
    double x = 0.0, y = 0.0, z = 0.0;
    for (unsigned int k = 0; k < Grid->nz; k++){
        for (unsigned int j = 0; j < Grid->ny; j++){
            for (unsigned int i = 0; i < Grid->nx; i++){
                
              x = Grid->gridResolution[0][0]*i + Grid->gridResolution[0][1]*j + Grid->gridResolution[0][2] * k;
              y = Grid->gridResolution[1][0]*i + Grid->gridResolution[1][1]*j + Grid->gridResolution[1][2] * k;
              z = Grid->gridResolution[2][0]*i + Grid->gridResolution[2][1]*j + Grid->gridResolution[2][2] * k;
              Grid->gridPointsXYZ->InsertNextPoint(x, y, z);
          }
        }
      }

    /* As the input file is only voxel data, the coordinates of the points are computed based
       on the grid resolution vector, and the grid in actual coordinates of the nanoporous material
       is stored in Grid->orginalGrid */
    
    if (writeGridFile == true)
    {
        vtkNew<vtkDoubleArray> pointValues;
        pointValues->SetNumberOfComponents(1);
        pointValues->SetNumberOfTuples(Grid->nx*Grid->ny*Grid->nz);
        pointValues->SetName(arrayName.c_str());
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

        // It could be useful to dump the fractional grid too for debugging purposes,
        // for this use -writefractionalgrid at invocation

        if (writeFractionalGrid) {
            vtkNew<vtkXMLImageDataWriter> imageWriter;
            std::ostringstream oss;
            oss << Directory << "/" << BaseFileName << "_fractional_grid." << imageWriter->GetDefaultFileExtension();
            imageWriter->SetInputData(Grid->cubicGrid);
            imageWriter->SetFileName(oss.str().c_str());
            imageWriter->Write();
        }
        
    }
 
    double elapsedTime = readerTime.getElapsedTime();
    logger::mainlog << "Time elapsed in the reader module: " << elapsedTime << "(s)" << endl;
    
    return Grid;
}




/**
 @brief Reads the cube file line by line. The code branches into two for distance grids
 and PEgrids as zeo++ stores it in slightly different format with an additional line 1 1 before the start of data.
 @return the grid as an vtkImageData grid.
 */
vtkSmartPointer<vtkImageData> segmentor::readFromCubeFile(){
    
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    
    if (arrayName == "This is distance grid"){
        readDistanceGrid(imageData);
        
    } else if (arrayName == "Potential Energy"){
        readPEgrid(imageData);
    }
    
    return imageData;
    
}



/**
 @brief Reads the distance grid generated from zeo++
 @param imageData is a vtkImageData class, and is updated in this function
*/
void segmentor::readDistanceGrid(vtkSmartPointer<vtkImageData> imageData){
    
    ifstream inputFile;
    
    inputFile.open(fileName);
    
    if(!inputFile) {
        throw std::runtime_error("Failed to open file!!\n");
    }
    int nx, ny, nz;
    int natoms;
    if (inputFile.is_open())
    {
            
        string line;
        int lineNumber = 0;
        while (getline(inputFile,line))
        {
            if (lineNumber == 2){
                stringstream ss(line);
                ss >> natoms;
            }
            if (lineNumber == 3)
            {
                stringstream ss(line);
                ss >> Grid->nx;
                ss >> Grid->gridResolution[0][0]; ss >> Grid->gridResolution[1][0]; ss >> Grid->gridResolution[2][0];
            }
                
            if (lineNumber == 4)
            {
                stringstream ss(line);
                ss >> Grid->ny;
                ss >> Grid->gridResolution[0][1]; ss >> Grid->gridResolution[1][1]; ss >> Grid->gridResolution[2][1];
            }
                
            if (lineNumber == 5)
            {
                stringstream ss(line);
                ss >> Grid->nz;
                ss >> Grid->gridResolution[0][2]; ss >> Grid->gridResolution[1][2]; ss >> Grid->gridResolution[2][2];
            }
                
            lineNumber++;
        }
                    
    }
    inputFile.close();

    // We set the dimensions of the imageData grid
    imageData->SetDimensions((int) Grid->nx, (int) Grid->ny, (int) Grid->nz);
    vtkNew<vtkDoubleArray> pointValues;
    pointValues->SetName(arrayName.c_str());
    pointValues->SetNumberOfComponents(1);
    size_t numberOfPoints = Grid->nx * Grid->ny * Grid->nz;
    pointValues->SetNumberOfTuples(numberOfPoints);
    
    // We open the input file again, to read the grid
    // Note that in zeo++, the z loop is the innermost, then, y, then x,
    // However in VTK, the x loop is the innermost, then y, then z.
    size_t ix = 0, iy = 0, iz = 0;
    inputFile.open(fileName);
    if (inputFile.is_open())
    {
        string line;
        //Number of the line that is being read
        int lineNumber = 0;
        while (getline(inputFile,line))
        {
            if (lineNumber > 6 + natoms){
                
                stringstream ss(line);
                double value;
                while (ss >> value) {
                    size_t index = Grid->nx*(iy + Grid->ny*iz) + ix;
                    pointValues->SetValue(index, value);
                    iz++;
                }
                if (iz == Grid->nz){
                    iz = 0;
                    iy++;
                    if (iy == Grid->ny){
                        iy = 0;
                        ix++;
                    }
                }
                    
            }
            
            lineNumber++;
        }
    }
    
    inputFile.close();
    imageData->GetPointData()->AddArray(pointValues);
    
}



/**
 @brief Reads the PE  grid generated from PorousMaterials.jl repo. 
  Some data points can have Inf, and this is given a large value of 5e+20.
 @param imageData is a vtkImageData class, and is updated in this function
*/
void segmentor::readPEgrid(vtkSmartPointer<vtkImageData> imageData){
    
    ifstream inputFile;
    inputFile.open(fileName);
    if(!inputFile) {
        throw std::runtime_error("Failed to open file!!\n");
    }
    
    int nx, ny, nz;
    int natoms;
    size_t InfCount = 0;
    
    if (inputFile.is_open())
    {
            
        string line;
        //Number of the line that is being readed
        int lineNumber = 0;
        while (getline(inputFile,line))
        {
            if (lineNumber == 2){
                stringstream ss(line);
                ss >> natoms;
            }
            if (lineNumber == 3)
            {
                stringstream ss(line);
                ss >> Grid->nx;
                ss >> Grid->gridResolution[0][0]; ss >> Grid->gridResolution[1][0]; ss >> Grid->gridResolution[2][0];
            }
                
            if (lineNumber == 4)
            {
                stringstream ss(line);
                ss >> Grid->ny;
                ss >> Grid->gridResolution[0][1]; ss >> Grid->gridResolution[1][1]; ss >> Grid->gridResolution[2][1];
            }
                
            if (lineNumber == 5)
            {
                stringstream ss(line);
                ss >> Grid->nz;
                ss >> Grid->gridResolution[0][2]; ss >> Grid->gridResolution[1][2]; ss >> Grid->gridResolution[2][2];
            }
                
            lineNumber++;
        }
                    
    }
    inputFile.close();

    // We set the dimensions of the imageData grid
    imageData->SetDimensions((int) Grid->nx, (int) Grid->ny, (int) Grid->nz);
    vtkNew<vtkDoubleArray> pointValues;
    pointValues->SetName(arrayName.c_str());
    pointValues->SetNumberOfComponents(1);
    size_t numberOfPoints = Grid->nx * Grid->ny * Grid->nz;
    pointValues->SetNumberOfTuples(numberOfPoints);
    
    // We open the input file again, to read the grid
    // Note that in zeo++, the z loop is the innermost, then, y, then x,
    // However in VTK, the x loop is the innermost, then y, then z.
    size_t ix = 0, iy = 0, iz = 0;
    inputFile.open(fileName);
    if (inputFile.is_open())
    {
        string line;
        //Number of the line that is being readed
        int lineNumber = 0;
        while (getline(inputFile,line))
        {
            if (lineNumber > 5 + natoms){
                
                stringstream ss(line);
                double value;
                std::string temp;
                while (ss >> temp) {
                    if (temp == "Inf"){
                        value = 5.0e+20;
                        InfCount++;
                    } else {
                        value = ::atof(temp.c_str());
                    }
                    size_t index = Grid->nx*(iy + Grid->ny*iz) + ix;
                    pointValues->SetValue(index, value);
                    iz++;
                }
                if (iz == Grid->nz){
                    iz = 0;
                    iy++;
                    if (iy == Grid->ny){
                        iy = 0;
                        ix++;
                    }
                }
                    
            }
            
            lineNumber++;
        }
    }
    
    inputFile.close();
    imageData->GetPointData()->AddArray(pointValues);
    logger::mainlog << "Number of Inf encountered while reading PEgrid               : " << InfCount << endl<<flush;
    
}



/**
 @brief The arrayname that is used for TDA is automatically obtained from the input file
  For file of type .cube, this is the string in the second line of the input file
  For file of type .vti, the first array of the grid (which is the distance grid) is the array for TDA.
 @param nameofArray which is updated in this function.
 */
void segmentor::getArrayName(std::string &nameOfArray){
    
    if (extensionName == ".cube"){
        
        if (nameOfArray.empty()) {
            getArrayNameFromCubeFile(nameOfArray);
        } else {
            std::string temp;
            getArrayNameFromCubeFile(temp);
            if (temp.compare(nameOfArray) != 0){
                std::cout << "(ERROR) Array name in the .cube file is " << temp << ", while that provided as input is " << nameOfArray << endl;
                exit(0);
            }
        } 
        logger::mainlog << "Array that is going to be used for TDA analysis              : " << nameOfArray << endl;
        
    } else if (extensionName == ".vti"){
        
        if (nameOfArray.empty()) getArrayNameFromVTIFile(nameOfArray);
        logger::mainlog << "Array that is going to be used for TDA analysis              : " << nameOfArray << endl;
        
    } else {
        
        logger::mainlog << "Extension type of the input file is neither .cube nor .vti" <<endl;
        logger::errlog << "Extension type of the input file neither .cube nor .vti" <<endl;
        exit(0);
    }
    
    
}




void segmentor::getArrayNameFromCubeFile(std::string &nameOfArray){
    
    ifstream inputFile;
    inputFile.open(fileName);
    if(!inputFile) {
        throw std::runtime_error("Failed to open file!!\n");
    }
    
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
    // If PorousMaterials.jl is used to generate the energy grids,
    // then the second line is Loop order - this is replaced with Potential Energy
    // for the code to read better .
    if (nameOfArray == "Loop order: x, y, z") nameOfArray = "Potential Energy";
    
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
 * @brief  Set periodic conditions for the smoothed grid.
 * @param smoothgrid Input grid
 * @param periodicConditions  Periodic Boundary Conditions
 * @param useAllCores
 * @return a pointer to a ttk variable of type ttkTriangulationManager.
 */
vtkSmartPointer<ttkTriangulationManager> segmentor::generatePeriodicGrid(vtkSmartPointer<ttkScalarFieldSmoother> smoothgrid, bool periodicConditions, bool useAllCores)
{

    logger::mainlog << "\n\nSegmentor: generate periodic grid Module"
                    << "\n";
    ttk::Timer periodicTimer;

    // VTK function used to set Periodic Boundary Conditions
    vtkSmartPointer<ttkTriangulationManager> periodGrid = vtkSmartPointer<ttkTriangulationManager>::New();
    periodGrid->SetUseAllCores(useAllCores);
    periodGrid->SetInputConnection(smoothgrid->GetOutputPort());
    periodGrid->SetPeriodicity(periodicConditions);
    periodGrid->Update();

    double elapsedTime = periodicTimer.getElapsedTime();
    logger::mainlog << "Time elapsed in periodic condition setter module: " << elapsedTime << "\n"
                    << flush;

    return periodGrid;
}



/**
 * @brief  Smooths the scalar grid by taking average of the neighbors for niterations.
 * @param grid Input grid
 * @param niterations number of iterations for smoothing. 
 * @param useAllCores
 * @return a pointer to a ttk variable of type ttkScalarFieldSmoother.
 */
vtkSmartPointer<ttkScalarFieldSmoother>  segmentor::smoothInputGrid(vtkSmartPointer<vtkImageData> grid, int niterations, bool useAllCores)
{
    logger::mainlog << "\n\nSegmentor: Smooth Input Grid module"
                    << "\n";
    ttk::Timer smootherTimer;

    vtkSmartPointer<ttkScalarFieldSmoother> smoother = vtkSmartPointer<ttkScalarFieldSmoother>::New();
    smoother->AddInputData(grid);
    smoother->SetNumberOfIterations(niterations);
    smoother->SetInputArrayToProcess(0, 0, 0, 0, arrayName.c_str());
    smoother->SetUseAllCores(useAllCores);
    smoother->Update();

    double elapsedTime = smootherTimer.getElapsedTime();
    logger::mainlog << "Time elapsed in in smoothing the input grid: " << elapsedTime << "\n"
                    << flush;

    return smoother;
}


/**
 @brief computes the persistence diagram for the periodic grid. 
 @param grid - Input grid
 @param useAllCores

 The variable thePersistenceDiagram is updated storing the persistencediagram and sets 
 the flag isPersistenceDiagramComputed to TRUE. 
*/
void segmentor::computePersistenceDiagram(vtkSmartPointer<ttkTriangulationManager> grid, bool useAllCores){

    logger::mainlog << "Computing persistence diagram ...";
    thePersistenceDiagram->SetUseAllCores(useAllCores);
    thePersistenceDiagram->SetInputConnection(grid->GetOutputPort());
    thePersistenceDiagram->SetInputArrayToProcess(0,0,0,0,arrayName.c_str());
    thePersistenceDiagram->Update();
    
    isPersistenceDiagramComputed = true;
    logger::mainlog << ". done" << endl;
}



/**
 * @brief Computes a topological simplification based on the persistence of the scalar field
 * attached to the input grid. Besides, it computes the Morse Smale Complex Segmentation. It writes 2 files: 1) Critical points file.
 * 2) Segmentation file. 
 * @param grid  Input grid to analyse
 * @param persistenceThreshold Persistence threshold to discard noisy maxima's and minima's.
 * @param saddlesaddleIncrement Persistence threshold increment(if needed) for the
 *  simplification of the sadde-saddle connectors
 * @param writeOutputs Write MSC results to external files
 * @param useAllCores Use all cores available to speed up computations
 *
 * @return nothing is returned, but the results of the MSC in the cubic grid is stored in theMSC,
 * the persistence diagram is stored in the thePersitenceDiagram, and in the actual coordinates,
 * the segmentation data is stored in the Grid->segmentation and the critical points in Grid->criticalPoints.
 * Moreover, critical points are also stored in native C++ data structures without the VTK wrappers. 
 */
void segmentor::MSC(vtkSmartPointer<ttkTriangulationManager> grid,double persistenceThreshold, double saddlesaddleIncrement, bool writeOutputs, bool useAllCores)
{
    logger::mainlog << "\n\nSegmentor: Morse Smale Complex Module" << "\n" << flush;
    
    ttk::Timer MSCTimer;
    
    if (!isPersistenceDiagramComputed) computePersistenceDiagram(grid, useAllCores); 
    
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

    //Morse Smale Complex Computation
    logger::mainlog << "Computing Morse Smale Complex ...";
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
    logger::mainlog << ". done" << endl;
    
    double timeTakenForMSC = MSCTimer.getElapsedTime();
    MSCTimer.reStart();
    logger::mainlog << "Time taken for MSC creation: " << timeTakenForMSC << "(s)" << endl;
    
    vtkIdType numberOfDescendingManifolds = getNumberOfDescendingManifolds();
    vtkIdType numberOfAscendingManifolds = getNumberOfAscendingManifolds();
    
    logger::mainlog << "Total number of descending manifolds (typically solid segments) : " << numberOfDescendingManifolds << endl;
    logger::mainlog << "Total number of ascending manifolds (typically void segments) : " << numberOfAscendingManifolds << endl;

    if (writeOutputs)
    {
        // We store the critical points and the segmentation data in actual coordinates.
        // These are stored in Grid->criticalPoints and Grid->segmentation
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

        /* Convert the segmentation data that is stored in a VTK grid to a native C++ structure
           that can be easily used for postprocessing */
        
        bool nativeData = true;
        if (nativeData){
            
            ttk::Timer nativeDataTimer;
            logger::mainlog << "Generating native C++ data without VTK dependency ...\n";
            Grid->segmentationData->setDimensions(Grid->nx, Grid->ny, Grid->nz);
            Grid->criticalPointsData.resize(criticalPointsDataSet->GetNumberOfPoints());
            
            // Storing the critical points as a vector of critical points class without VTK //
            for (size_t i = 0; i < cPoints->GetNumberOfPoints(); i++){
                double p[3];
                cPoints->GetPoint(i, p);
                Grid->criticalPointsData[i].setCoordinates(p);
                int critType = criticalPointsDataSet->GetPointData()->GetArray("CellDimension")->GetVariantValue(i).ToInt();
                double scalarVal = criticalPointsDataSet->GetPointData()->GetArray((arrayName).c_str())->GetVariantValue(i).ToDouble();
                Grid->criticalPointsData[i].setCriticalType((size_t) critType);
                Grid->criticalPointsData[i].setScalarValue(scalarVal);
                
            }

            // Storing the segmentation data without VTK //
            for (size_t i = 0; i < segmentationDataSet->GetNumberOfPoints(); i++){
                double p[3];
                Grid->gridPointsXYZ->GetPoint(i, p);
                int ascendingManifoldID = segmentationDataSet->GetPointData()->GetArray("AscendingManifold")->GetVariantValue(i).ToInt();
                int descendingManifoldID = segmentationDataSet->GetPointData()->GetArray("DescendingManifold")->GetVariantValue(i).ToInt();
                double scalarVal = segmentationDataSet->GetPointData()->GetArray(arrayName.c_str())->GetVariantValue(i).ToDouble();

                Grid->segmentationData->setCoordinates(i, p);
                Grid->segmentationData->setAscendingManifoldID(i, (size_t) ascendingManifoldID);
                Grid->segmentationData->setDescendingManifoldID(i, (size_t) descendingManifoldID);
                Grid->segmentationData->setScalarValue(i, scalarVal);
                                
            }
            
            double timeTakenForNativeData = nativeDataTimer.getElapsedTime();
            logger::mainlog << "Time taken to create native segmentation data: " << timeTakenForNativeData << endl;

        }
        
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



/**
 @brief Output the total number of Descending Manifolds. 
*/
vtkIdType segmentor::getNumberOfDescendingManifolds(){
        
    
    vtkSmartPointer<ttkExtract> descendingManifolds = vtkSmartPointer<ttkExtract>::New();
    descendingManifolds->SetUseAllCores(true);
    descendingManifolds->SetInputConnection(theMSC->GetOutputPort(3));
    descendingManifolds->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,"DescendingManifold");
    descendingManifolds->SetExtractionMode(3); //Array Values
    descendingManifolds->SetExtractUniqueValues(true);
    descendingManifolds->Update();

    auto descendingManifoldsDataset = vtkDataSet::SafeDownCast(descendingManifolds->GetOutputDataObject(0));
    auto descendingManifoldsID = descendingManifoldsDataset->GetFieldData()->GetAbstractArray("UniqueDescendingManifold");

    vtkIdType numberOfDescendingManifolds = descendingManifoldsID->GetNumberOfValues();
    
    return numberOfDescendingManifolds;
    
}



/**
 @brief Output the total number of Ascending Manifolds. 
*/
vtkIdType segmentor::getNumberOfAscendingManifolds(){
    
    vtkSmartPointer<ttkExtract> ascendingManifolds = vtkSmartPointer<ttkExtract>::New();
    ascendingManifolds->SetUseAllCores(true);
    ascendingManifolds->SetInputConnection(theMSC->GetOutputPort(3));
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
 @brief for the input periodic grid, this plots the persistence curve and write all the four pairs of persistence curves:
 @write 1) Minumum-saddle Pairs
      2) Saddle-saddle Pairs
      3) Saddle-Maximum Pairs
      4) All-pairs.
 These are output as tables in text files, and to plot one can use postProcessing/plotPersistenceCurve.py BaseFileName
 */
void segmentor::persistencecurve(vtkSmartPointer<ttkTriangulationManager> grid, bool useAllCores)
{
    ttk::Timer percurveTimer;
    logger::mainlog << "\nSegmentor: Persistence Curve module" << "\n" << flush;
 
    // To create the persistence curve, first the persistence diagram is calculated 
    // which serves as input to the persistence curve
    if (!isPersistenceDiagramComputed) computePersistenceDiagram(grid, useAllCores);

    thePersistenceCurve->SetInputConnection(thePersistenceDiagram->GetOutputPort());
    thePersistenceCurve->Update();

    /* Write Minimum Saddle Pairs */
    ofstream pcurveResults1;
    pcurveResults1.open((Directory + "/" + BaseFileName + "-minSaddlePairs.txt").c_str());
    assert(pcurveResults1.is_open());
    pcurveResults1 << "   #persistence     minSaddlePairs  " << "\n";
    
    vtkTable* outputTable = thePersistenceCurve->GetOutput(0);
    vtkIdType nrows = outputTable->GetNumberOfRows();
    vtkAbstractArray * output1 = outputTable->GetColumn(0);
    vtkAbstractArray * output2 = outputTable->GetColumn(1);

    for (size_t i = 0; i < nrows; i++){
        pcurveResults1 << scientific << setw(18) << output1->GetVariantValue(i).ToDouble() << setw(18) << output2->GetVariantValue(i).ToDouble() << "\n";
    }


    /* Write Saddle Saddle Pairs */
    ofstream pcurveResults2;
    pcurveResults2.open((Directory + "/" + BaseFileName + "-SaddleSaddlePairs.txt").c_str());
    assert(pcurveResults2.is_open());
    pcurveResults2 << "   #persistence   SaddleSaddlePairs " << "\n";
    
    outputTable = thePersistenceCurve->GetOutput(1);
    nrows = outputTable->GetNumberOfRows();
    output1 = outputTable->GetColumn(0);
    output2 = outputTable->GetColumn(1);

    for (size_t i = 0; i < nrows; i++){
        pcurveResults2 << scientific << setw(18) << output1->GetVariantValue(i).ToDouble() << setw(18) << output2->GetVariantValue(i).ToDouble() << "\n";
    }

    /* Write maximum Saddle Pairs */
    ofstream pcurveResults3;
    pcurveResults3.open((Directory + "/" + BaseFileName + "-maxSaddlePairs.txt").c_str());
    assert(pcurveResults3.is_open());
    pcurveResults3 << "   #persistence     MaxSaddlePairs  " << "\n";
    
    outputTable = thePersistenceCurve->GetOutput(2);
    nrows = outputTable->GetNumberOfRows();
    output1 = outputTable->GetColumn(0);
    output2 = outputTable->GetColumn(1);

    for (size_t i = 0; i < nrows; i++){
        pcurveResults3 << scientific << setw(18) << output1->GetVariantValue(i).ToDouble() << setw(18) << output2->GetVariantValue(i).ToDouble() << "\n";
    }

    /* Write all Pairs */
    ofstream pcurveResults4;
    pcurveResults4.open((Directory + "/" + BaseFileName + "-allPairs.txt").c_str());
    assert(pcurveResults4.is_open());
    pcurveResults4 << "   #persistence        AllPairs     " << "\n";
    
    outputTable = thePersistenceCurve->GetOutput(3);
    nrows = outputTable->GetNumberOfRows();
    output1 = outputTable->GetColumn(0);
    output2 = outputTable->GetColumn(1);

    for (size_t i = 0; i < nrows; i++){
        pcurveResults4 << scientific << setw(18) << output1->GetVariantValue(i).ToDouble() << setw(18) << output2->GetVariantValue(i).ToDouble() << "\n";
    }

    pcurveResults1.close();
    pcurveResults2.close();
    pcurveResults3.close();
    pcurveResults4.close();    
    
    double timeTakenForPerCurve = percurveTimer.getElapsedTime();
    logger::mainlog << "Time taken in persistence curve module: " << timeTakenForPerCurve << "(s)" << endl;

}

