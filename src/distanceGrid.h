/*********************************************************************

TDA-Segmentor    A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:         Aditya Vasudevan (adityavv.iitkgp@gmail.com)
          Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
          Maciek Haranczyk (maciej.haranczyk@imdea.org)
          
          IMDEA Materiales Institute
 
**********************************************************************/
#ifndef distanceGrid_hpp
#define distanceGrid_hpp

#include <stdio.h>
#include "grid.h"


/*
 This is the distance grid instance.
 */

class distanceGrid : public grid {
    
public:
    
                             distanceGrid(const parameters& p);
    virtual                  ~distanceGrid();
    
    
    virtual void             voidSegmentation();
    virtual void             solidSegmentation();
    virtual void             accessibleVoidSpace(double moleculeRadius, bool useAllCores);
    virtual void             accessibleVoidGraph(double moleculeRadius, bool useAllCores);
    virtual void             accessibleSolidGraph(bool useAllCores);
    
};



#endif /* distanceGrid_hpp */
