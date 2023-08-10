/*********************************************************************

TDA-Segmentor    A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:         Aditya Vasudevan (adityavv.iitkgp@gmail.com)J
                 Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
                 Maciek Haranczyk (maciej.haranczyk@imdea.org)
                 IMDEA Materiales Institute
 
**********************************************************************/

#ifndef segmentor_h
#define segmentor_h

#include <stdio.h>
#include "headers.h"
#include "logger.h"

#define DEBUG 1

class segmentor
{
private:
    

    
    
public:
    //Nanoporous Material name
    string BaseFileName, fileName, extensionName;
    string arrayName; 
    //Directory where the results will be store
    string Directory;
    //Input grid
    vtkImageData * Grid = vtkImageData::New();
    //Input grid resolution
    double GridResolution;
    //Cell Size(default to zero)
    double CellSize = 0;

    //Cleaned functions
    segmentor(const parameters &p);
    ~segmentor();
    auto superCell(vtkSmartPointer<vtkImageData> grid);

    void getVolume(string inputFolder, double resolution);
    vtkIdType getNumberOfDescendingManifolds(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex);
    vtkIdType getNumberOfAscendingManifolds(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex);

    void getGridResolutionFromCubeFile(double gridres[3]);
    void getArrayNameFromCubeFile(std::string &arrayname);
    void getCubeVolume(string inputFolder, double resolution);
    auto reader(string inputFilePath,int gridResolution, bool writeGridFile);
    auto readInputFile(bool writeGridFile = false);
    auto readerCombined(string inputFilePath1,string inputFilePath2,double gridResolution,bool writeGridFile);
    auto inputPrecondition2(vtkSmartPointer<vtkImageData> grid,bool periodicConditions, bool computeDistanceField, bool writeFile);
    auto inputPrecondition(vtkSmartPointer<vtkImageData> grid, bool changeValues,bool periodicConditions, bool useAllCores);
    auto inputPrecondition_E(vtkSmartPointer<vtkImageData> grid,bool periodicConditions);
    auto MSC(vtkSmartPointer<ttkTriangulationManager> grid,double persistenceThreshold, double saddlesaddleIncrement, bool writeOutputs, bool useAllCores);
    auto MSC_E(vtkSmartPointer<ttkTriangulationManager> grid,double persistencePercentage, double saddlesaddleIncrement, bool writeOutputs, bool useAllCores);
    void voidSegmentation(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex, bool useAllCores);
    void accessibleVoidSpace(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex,double moleculeRadius, bool useAllCores);
    void voidSegmentation_E(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex);
    auto solidSegmentation(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex);
    void eigenField(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex,int numberOfEigenFunctions, bool writeSegments,string scalar, bool useAllCores);
    void eigenStructure(vtkSmartPointer<vtkImageData> grid, int numberOfEigenFunctions, bool useAllCores);
    void persistencecurve(vtkSmartPointer<ttkTriangulationManager> grid, bool useAllCores);



    //-------------------------------------------------------------------------------

   
    void gridFileCreator(string scalarName, string inputFilePath, double persistencePercentage, bool useAllCores);
    

    void oneGridFileCreator(string scalarName, string inputFilePath, double persistencePercentage, bool useAllCores);
   
    auto evolutionFile2( vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex,vtkSmartPointer<vtkThreshold> previousSolid);


     
    
     
    void voidSeparatrices(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex);
    auto stagesEvolution(vector<string> fileNames);
    auto getIndex(vector<int> v, int K);
    auto solidGetter(string currentFile);
    auto saddlesGetter(string currentFile);
  
    int findMostCommonValue(vector<int> &inputVector);
  

    void energyDiagrams(vtkSmartPointer<vtkImageData> grid, bool useAllCores);
    void energyDiagrams2(vtkSmartPointer<vtkImageData> grid, bool useAllCores);
    void distanceDiagrams2(vtkSmartPointer<vtkImageData> grid, bool useAllCores);


    void energyVsDistance(string inputFile, string distanceDirectory, string energyDirectory);

    auto persistenceMatchings(bool useAllCores);
    auto persistenceDiagramsWriter(bool useAllCores);
    auto energyIncluder(vtkSmartPointer<vtkImageData> distanceGrid,string energyFile,double persistenceThreshold, bool useAllCores);
    
    auto segmentsShapes(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex,int numberOfEigenFunctions, bool writeSegments,string scalar, bool useAllCores);
    auto segmentsShapes2(vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex,int numberOfEigenFunctions, bool writeSegments,string scalar, bool useAllCores);
    auto segmentSelection(string inputFile, int numberOfEigenFunctions, bool writeOutputs,bool useAllCores);


};

#include "segmentor.cpp"


#endif 

