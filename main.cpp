/*********************************************************************

TDA-Segmentor     -     A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:                       Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
                 Aditya Vasudevan (adityavv.iitkgp@gmail.com)
                 Maciek Haranczyk (maciej.haranczyk@imdea.org)
                 IMDEA Materiales Institute
 
**********************************************************************/

#include "poretda.h" 

int main(int argc, char * argv[]){
       

     auto startTotal = high_resolution_clock::now();

     string inputfilename = argv[1];


     poretda * analysis = new poretda(inputfilename);
     auto grid = analysis->reader(inputfilename,0.15,0);
     auto period = analysis->inputPrecondition(grid,true,true,false);
     auto morseSmale = analysis->MSC(period,0.12,1.0,true,false);
     analysis->accessibleVoidSpace(morseSmale,1.6,false);
     analysis->voidSegmentation(morseSmale,0);


    //Computations' finish time(TDA)
    auto stopTotal = high_resolution_clock::now();
    auto durationTotal = duration_cast<seconds>(stopTotal - startTotal);
    poretda::mainlog << "Total Time taken for all the files: " << durationTotal.count() << " seconds" <<endl;

    delete analysis;

    return 0;
}
