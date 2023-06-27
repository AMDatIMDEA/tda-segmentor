/*********************************************************************

TDA-Segmentor     -     A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:                       Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
                 Aditya Vasudevan (adityavv.iitkgp@gmail.com)
                 Maciek Haranczyk (maciej.haranczyk@imdea.org)
                 IMDEA Materiales Institute
 
**********************************************************************/

#include "logger.h"


ofstream      logger::mainlog;
ofstream      logger::errlog;


void logger::openLogfiles(const parameters &p) {
    
    ostringstream logstrm;
    ostringstream errstrm;
    logstrm << p.basefilename << ".log";
    errstrm << p.basefilename << ".err";
    string logFilename = logstrm.str();
    string errFilename = errstrm.str();
    mainlog.open(logFilename.c_str());
    errlog.open(errFilename.c_str());
    
    logger::mainlog << "\n";
    logger::mainlog << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << "\n";
    logger::mainlog << "                                                                                      " << "\n";
    logger::mainlog << "                T   D   A    -    S   E   G   M   E   N   T   O   R                   " << "\n";
    logger::mainlog << "                                                                                      " << "\n";
    logger::mainlog << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << "\n\n"<< flush;
    
    logger::mainlog << "Base file name  :                         " << p.basefilename << "\n" << flush;
    logger::mainlog << "Extension       :                         " << p.extensionname << "\n\n" << flush;
    
    if (!mainlog){
        std::cout << "Error opening mainlog file!" << endl;
        exit(0);
    } else if (!errlog) {
        std::cout << "Error opening error log file!" << endl;
    }
    
    
}




void logger::closeLogfiles(){
    
    
    if (mainlog.good())
    {
        mainlog << "\n\nFinished Analysis!" << endl;
        mainlog.close();
    }
    
    if (errlog.good())
    {
        errlog << "\n\nFinished Analysis!" << endl;
        errlog.close();
    }
}

