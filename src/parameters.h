/*********************************************************************

TDA-Segmentor    A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:         Aditya Vasudevan (adityavv.iitkgp@gmail.com)
                 Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)              
                 Maciek Haranczyk (maciej.haranczyk@imdea.org)
                 IMDEA Materiales Institute
 
**********************************************************************/


#ifndef inputhandler_hpp
#define inputhandler_hpp

#include <stdio.h>
#include "headers.h"

class parameters
{
public:

    parameters();
    
    void                               parser(int nargs, char** args);
    void                               printinvocation(int nargs, char** args);

    void                               writetoLogFile();
    
    
    std::vector <std::string>          moduleNames;
    std::string                        inputfilename, basefilename, extensionname, Directory;
    std::string                        arrayName;
    bool                               useSuperCell;
    bool                               useAllCores;
    bool                               saveLogFile;
    bool                               segmentationFlag;
    bool                               ftmTreeFlag;
    double                             persistenceThreshold, probeRadius;
    
private:
    
    bool                               checkinputfile();
    

    
} ;

#endif /* inputhandler_hpp */
