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


int main(int argc, char * argv[]){
       

     ttk::Timer programTimer;
    
    if (argc != 4) {
        printInstructions();
        exit(0);
    }
    
    double persistenceValue = stod(argv[1]);
    double proberadius = stod(argv[2]);
    string inputfilename = argv[3];


     segmentor * analysis = new segmentor(inputfilename);
     auto grid = analysis->reader();
     auto periodicGrid = analysis->inputPrecondition(grid,true,true,false);
     auto morseSmale = analysis->MSC(periodicGrid,persistenceValue,1.0,true,false);
     analysis->accessibleVoidSpace(morseSmale,proberadius,false);
     analysis->voidSegmentation(morseSmale,0);


    //Computations' finish time(TDA)
    double elapsedTime = programTimer.getElapsedTime();
    segmentor::mainlog << "\nTotal Time taken for all the files: " << elapsedTime << " seconds" <<endl;

    delete analysis;

    return 0;
}
