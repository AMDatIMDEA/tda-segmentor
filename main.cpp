/*********************************************************************

TDA-Segmentor     -     A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:                       Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
                 Aditya Vasudevan (adityavv.iitkgp@gmail.com)
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
    logger::openLogfiles(param.inputfilename, param.saveLogFile);
    param.printinvocation(argc, argv);
    param.writetoLogFile();
        
    
    segmentor * analysis = new segmentor(param.inputfilename);
    

     auto grid = analysis->reader();
     auto periodicGrid = analysis->inputPrecondition(grid,true,true,false);
     
    vtkSmartPointer<ttkMorseSmaleComplex> morseSmaleComplex;
    
    if (param.segmentationFlag){
        morseSmaleComplex = analysis->MSC(periodicGrid,param.persistenceThreshold,1.0,true,false);
    }

    for (size_t i = 0; i < param.moduleNames.size(); i++ )
    {
        if (param.moduleNames[i] == "accessiblevoidspace")
            analysis->accessibleVoidSpace(morseSmaleComplex,param.probeRadius,false);
        else if (param.moduleNames[i] == "voidsegmentation")
            analysis->voidSegmentation(morseSmaleComplex,0);
        else if (param.moduleNames[i] == "persistencecurve")
        {
            // To be completed
        } else if (param.moduleNames[i] == "solidsegmentation"){
            // to be completed
        }
        
    }
     

    //Computations' finish time(TDA)
    double elapsedTime = programTimer.getElapsedTime();
    logger::mainlog << "\nTotal Time taken for all the files: " << elapsedTime << " seconds" <<endl;

    delete analysis;
    logger::closeLogfiles();

    return 0;
}
