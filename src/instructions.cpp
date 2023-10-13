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
    
    cout << "tda-segmentor [-msc    [persistence-threshold]]" << "\n";
    cout << "              [-vs     [persistence-threshold]]" << "\n";
    cout << "              [-ss     [persistence-threshold]]" << "\n";
    cout << "              [-avs    [persistence-threshold] [probe-radius]]" << "\n";
    cout << "              [-avg    [persistence-threshold] [probe-radius]]" << "\n";
    cout << "              [-asg    [persistence-threshold]]" << "\n";
    cout << "              [-scalar [array-name]]" << "\n";
    cout << "              [-pc]" << "\n";
    cout << "              [-usesupercell] " << "\n";
    cout << "              [-useallcores] " << "\n";
    cout << "              [-writefractionalgrid] " << "\n";
    cout << "              [-savelogfile] " << "\n";
    cout << "                INPUTFILE" << "\n";
    
    
    cout << "Detailed explanation of commands: " << "\n";

    cout << "                        [-msc [persistence-threshold]]                 - Generates the morse-smale-complex (segmentation) of the grid." << "\n";
    cout << "                                                                         It takes a persistence threshold as an optional input." << "\n";
    cout << "                        [-vs [persistence-threshold]]                  - Saves the segmentation corresponding to the void structure." << "\n";
    cout << "                                                                         It takes a persistence threshold as an optional input." << "\n";
    cout << "                        [-ss [persistence-threshold]]                  - Saves the segmentation corresponding to the solid structure." << "\n";
    cout << "                                                                         It takes a persistence threshold as an optional input." << "\n";
    cout << "                        [-avs [persistence-threshold] [probe-radius]]  - Saves the segmentation of the void space accessible to a probe." << "\n";
    cout << "                                                                         It takes first, a persistence threshold as an optional input, and next, radius of the probe atom." << "\n";
    cout << "                        [-avg [persistence-threshold] [probe-radius]]  - Saves the graph corresponding to the void space accessible to a probe." << "\n";
    cout << "                                                                         It takes first, a persistence threshold as an optional input, and next, radius of the probe atom." << "\n";
    cout << "                        [-asg [persistence-threshold]]                 - Saves the graph corresponding to the solid space." << "\n";
    cout << "                                                                         It takes a persistence threshold as an optional input." << "\n";
    cout << "                        [-pc ]                                         - It saves the persistence curve for the input grid." << "\n";
    cout << "                                                                         This option does not have any input." << "\n";
    cout << "                        [-usesupercell]                                - This flag generates a 2 X 2 X 2 super cell for analysis." << "\n";
    cout << "                                                                         This flag does not have any input." << "\n";
    cout << "                        [-useallcores]                                 - This switches on the multi-core capabilities of TTK." << "\n";
    cout << "                                                                         For this to work, TTK must be compiled with MPI. See installation instructions. " << "\n";
    cout << "                        [-writefractionalgrid]                         - Saves the input grid stored in fractional coordinates. " << "\n";
    cout << "                        [-savelogfile]                                 - Saves the previous log file by renaming it. " << "\n";
    cout << "                        [-scalar  [array-name]]                        - Sets the name of the array that will be used for analysis. " << "\n";
    cout << "                                                                         Possible options of array-name are : " << "\n";
    cout << "                                                                                 -scalar distance" << "\n";
    cout << "                                                                                 -scalar energy" << "\n";
    cout << "                                                                                 -scalar sigmaXX" << "\n";
    cout << "                                                                                 -scalar sigmaYY" << "\n";
    cout << "                                                                                 -scalar sigmaZZ" << "\n";
    cout << "                                                                                 -scalar sigmaXY" << "\n";
    cout << "                                                                                 -scalar sigmaXZ" << "\n";
    cout << "                                                                                 -scalar sigmaYZ" << "\n";


    cout << "\n Example of typical invocation: " << "\n";
    cout << " tda-segmentor -msc 0.04 -vs -avg 1.6 -useallcores -writefractionalgrid FAU.cube" << "\n";
    
}
