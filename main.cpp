/*********************************************************************

TDA-Segmentor     -     A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:                       Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
                 Aditya Vasudevan (adityavv.iitkgp@gmail.com)
                 Maciek Haranczyk (maciej.haranczyk@imdea.org)
                 IMDEA Materiales Institute
 
**********************************************************************/

#include "segmentor.h" 

int main(int argc, char * argv[]){
       

     ttk::Timer programTimer;

     string inputfilename = argv[1];


     segmentor * analysis = new segmentor(inputfilename);
     auto grid = analysis->reader(0.15,0);
     auto period = analysis->inputPrecondition(grid,true,true,false);
     auto morseSmale = analysis->MSC(period,0.12,1.0,true,false);
     analysis->accessibleVoidSpace(morseSmale,1.6,false);
     analysis->voidSegmentation(morseSmale,0);


    //Computations' finish time(TDA)
    double elapsedTime = programTimer.getElapsedTime();
    segmentor::mainlog << "\nTotal Time taken for all the files: " << elapsedTime << " seconds" <<endl;

    delete analysis;

    return 0;
}
