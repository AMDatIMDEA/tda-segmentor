/*********************************************************************

TDA-Segmentor    A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:         Aditya Vasudevan (adityavv.iitkgp@gmail.com)
                 Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)              
                 Maciek Haranczyk (maciej.haranczyk@imdea.org)
                 IMDEA Materiales Institute
 
**********************************************************************/

#ifndef logger_hpp
#define logger_hpp

#include <stdio.h>
#include "headers.h"
#include "parameters.h"


/*
  The logger class has two streams mainlog and errlog on which information
  is continuously written as the code is executed. 
*/
class logger {
    
public:
    static  std::ofstream          mainlog;
    static  std::ofstream          errlog;

    static void                    openLogfiles(const parameters &p);
    static void                    closeLogfiles(); 
};

#endif /* logger_hpp */
