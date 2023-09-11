/*********************************************************************

TDA-Segmentor    A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:         Aditya Vasudevan (adityavv.iitkgp@gmail.com)
          Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
          Maciek Haranczyk (maciej.haranczyk@imdea.org)
          
          IMDEA Materiales Institute
 
**********************************************************************/
#ifndef PEgrid_hpp
#define PEgrid_hpp

#include <stdio.h>
#include "grid.h"

/*
 This is the Potential Energy grid instance
 */
class PEgrid : public grid {
    
public:
    
                               PEgrid(const parameters& p);
    virtual                    ~PEgrid();
    
    virtual void               voidSegmentation();
    virtual void               solidSegmentation();
    virtual void               accessibleVoidSpace(double moleculeRadius, bool useAllCores);

};


#endif /* PEgrid_hpp */
