/*********************************************************************

TDA-Segmentor    A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:         Aditya Vasudevan (adityavv.iitkgp@gmail.com)
                 Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
                 Maciek Haranczyk (maciej.haranczyk@imdea.org)
                 IMDEA Materiales Institute
 
**********************************************************************/

#include "segmentor.h" 
#include "instructions.h"
#include "parameters.h"
#include "logger.h"

int main(int argc, char ** argv){
    
    ttk::Timer programTimer;
    
    parameters param;
    param.parser(argc,argv);
    logger::openLogfiles(param);
    param.printinvocation(argc, argv);
    param.writetoLogFile();
        
    
    segmentor * analysis = new segmentor(param);
    
    vtkSmartPointer<vtkImageData> grid = analysis->readInputFile();
    
    if (param.useSuperCell) grid = analysis->superCell(grid);
    
    
    auto periodicGrid = analysis->inputPrecondition(grid,true,true,param.useAllCores);

    vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex;
    vtkSmartPointer<ttkFTMTree> ftmTree;
    if (param.segmentationFlag){
        morseSmaleComplex = analysis->MSC(periodicGrid,param.persistenceThreshold,1.0,true,param.useAllCores);
    }
    if (param.ftmTreeFlag){
        ftmTree = analysis->ftmtree(periodicGrid, param.persistenceThreshold, param.useAllCores); 
    }

    for (size_t i = 0; i < param.moduleNames.size(); i++ )
    {
        if (param.moduleNames[i] == "accessiblevoidspace")
            analysis->accessibleVoidSpace(morseSmaleComplex,param.probeRadius,param.useAllCores);
        else if (param.moduleNames[i] == "voidsegmentation")
            analysis->voidSegmentation(morseSmaleComplex,0);
        else if (param.moduleNames[i] == "persistencecurve")
        {
            analysis->persistencecurve(periodicGrid, param.useAllCores);
        } else if (param.moduleNames[i] == "solidsegmentation"){
            analysis->solidSegmentation(morseSmaleComplex);
        } else if (param.moduleNames[i] == "accessiblevoidgraph"){
            analysis->accessibleVoidGraph(ftmTree, param.probeRadius, param.useAllCores);
        } else if (param.moduleNames[i] == "accessiblesolidgraph") {
            analysis->accessibleSolidGraph(ftmTree, param.useAllCores); 
        }
        
    }
     

    //Computations' finish time(TDA)
    double elapsedTime = programTimer.getElapsedTime();
    logger::mainlog << "\nTotal Time taken for all the modules: " << elapsedTime << " seconds" <<endl;

    delete analysis;
    logger::closeLogfiles();

    return 0;
}
