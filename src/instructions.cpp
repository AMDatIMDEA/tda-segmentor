/*********************************************************************

TDA-Segmentor     -     A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:                       Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
                 Aditya Vasudevan (adityavv.iitkgp@gmail.com)
                 Maciek Haranczyk (maciej.haranczyk@imdea.org)
                 IMDEA Materiales Institute
 
**********************************************************************/

#include "instructions.h"

using namespace std;

void printInstructions()
{
    
    cout << "\ntda-segmentor command-line syntax is as follows : " << "\n\n" ;
    
    cout << "tda-segmentor [-module [name-of-module]]" << "\n";
    cout << "                [-persistencethreshold [persistence-val]]" << "\n";
    cout << "                [-proberadius [probe-radius-val]]" << "\n";
    cout << "                [-usesupercell] " << "\n";
    cout << "                [-usetbb] " << "\n";
    cout << "                [-savelogfile] " << "\n";
    cout << "                INPUTFILE" << "\n";
    
    
    
    cout << "\nImplemented module names : " << "\n";
    cout << "                             -module segmentation" << "\n";
    cout << "                             -module accessiblevoidspace" << "\n";
    cout << "                             -module voidsegmentation" << "\n";
    cout << "                             -module solidsegmentation" << "\n";
    cout << "                             -module persistencecurve" << "\n";
    
    
    cout << "\n Example of typical invocation: " << "\n";
    cout << " poretda -module segmentation -persistencethreshold 0.12 -module accessiblevoidspace -proberadius 1.6 FAU.cube" << "\n";
    
}
