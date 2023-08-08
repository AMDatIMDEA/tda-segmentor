/*********************************************************************

TDA-Segmentor    A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:         Aditya Vasudevan (adityavv.iitkgp@gmail.com)
                 Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)              
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
    cout << "                [-scalar [name-of-array]]" << "\n";
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
    
    cout << "\nPossible options for scalar input are : " << "\n";
    cout << "                             -scalar distance" << "\n";
    cout << "                             -scalar energy" << "\n";
    cout << "                             -scalar sigmaXX" << "\n";
    cout << "                             -scalar sigmaYY" << "\n";
    cout << "                             -scalar sigmaZZ" << "\n";
    cout << "                             -scalar sigmaXY" << "\n";
    cout << "                             -scalar sigmaXZ" << "\n";
    cout << "                             -scalar sigmaYZ" << "\n";



    cout << "\n Example of typical invocation: " << "\n";
    cout << " tda-segmentor -module segmentation -persistencethreshold 0.12 -module accessiblevoidspace -proberadius 1.6 FAU.cube" << "\n";
    
}
