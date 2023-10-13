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

/* This class parses all the input information such as 
    - name of the input file,
    - the type of analysis that is to be performed,
    - setting optional flags.
*/


enum class module {

    SEGMENTATION,
    VOID_SEGMENTATION,
    SOLID_SEGMENTATION,
    ACCESSIBLE_VOID_SPACE,
    ACCESSIBLE_VOID_GRAPH,
    ACCESSIBLE_SOLID_GRAPH,
    PERSISTENCE_CURVE
};


class parameters
{
public:

    parameters();
    
    void                               parser(int nargs, char** args);
    void                               printinvocation(int nargs, char** args);
    std::string                        moduleNameAsString(const module modulename);                               

    void                               writetoLogFile();
    
    
    std::vector <module>               moduleNames;
    std::string                        inputfilename, basefilename, extensionname, Directory;
    std::string                        arrayName;
    bool                               useSuperCell;
    bool                               useAllCores;
    bool                               writeFractionalGrid;
    bool                               saveLogFile;
    bool                               segmentationFlag;
    double                             persistenceThreshold, probeRadius;
    
private:
    
    bool                               checkinputfile();
    

    
} ;

#endif /* inputhandler_hpp */
