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
    
    cout << "\n\ntda-segmentor command-line syntax is as follows : " << "\n\n" ;
    
    cout << "tda-segmentor [-persistency [persistency_val]]" << "\n";
    cout << "                [-proberad [rad_val]]" << "\n";
    cout << "                INPUTFILE" << "\n";
    
    
}
