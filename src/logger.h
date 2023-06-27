/*********************************************************************

TDA-Segmentor     -     A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:                       Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
                 Aditya Vasudevan (adityavv.iitkgp@gmail.com)
                 Maciek Haranczyk (maciej.haranczyk@imdea.org)
                 IMDEA Materiales Institute
 
**********************************************************************/

#ifndef logger_hpp
#define logger_hpp

#include <stdio.h>
#include "headers.h"


class logger {
    
public:
    static  std::ofstream          mainlog;
    static  std::ofstream          errlog;

    static void                    openLogfiles(std::string inputfileName, bool saveLogFile);
    static void                    closeLogfiles(); 
};

#endif /* logger_hpp */
