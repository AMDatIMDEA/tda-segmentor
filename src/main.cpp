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

    grid* inputGrid = analysis->readInputFile(param);
    
    if (param.useSuperCell) inputGrid->convertToSuperCell(); // inputGrid is now updated to a supercell 

    // Data is smoothed by NITERATIONS (set in headers.h) and periodicity is applied
    auto smoothGrid = analysis->smoothInputGrid(inputGrid->cubicGrid, NITERATIONS, param.useAllCores);
    auto periodicGrid = analysis->generatePeriodicGrid(smoothGrid, true, param.useAllCores);

    vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex;
    if (param.segmentationFlag){
         analysis->MSC(periodicGrid,param.persistenceThreshold,1.0,true,param.useAllCores);
    }

    for (size_t i = 0; i < param.moduleNames.size(); i++ )
    {
        if (param.moduleNames[i] == module::ACCESSIBLE_VOID_SPACE)
            inputGrid->accessibleVoidSpace(param.probeRadius,param.useAllCores);
        else if (param.moduleNames[i] == module::VOID_SEGMENTATION)
            inputGrid->voidSegmentation();
        else if (param.moduleNames[i] == module::PERSISTENCE_CURVE)
        {
            analysis->persistencecurve(periodicGrid, param.useAllCores);
        } else if (param.moduleNames[i] == module::SOLID_SEGMENTATION){
            inputGrid->solidSegmentation();
        } else if (param.moduleNames[i] == module::ACCESSIBLE_VOID_GRAPH){
            inputGrid->accessibleVoidGraph(param.probeRadius, param.useAllCores);
        } else if (param.moduleNames[i] == module::ACCESSIBLE_SOLID_GRAPH) {
            inputGrid->accessibleSolidGraph(param.probeRadius, param.useAllCores);
        }
        
    }
     

    double elapsedTime = programTimer.getElapsedTime();
    logger::mainlog << "\nTotal Time taken for all the modules: " << elapsedTime << " seconds" <<endl;

    delete analysis;
    logger::closeLogfiles();

    return 0;
}
