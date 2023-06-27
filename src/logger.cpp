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


void logger::openLogfiles(std::string inputfileName, bool saveLogFile) {
    
    std::string::size_type const p(inputfileName.find_last_of('.'));
    std::string BaseFileName = inputfileName.substr(0, p);
    std::string extensionName = inputfileName.substr(p,inputfileName.length());
    
    
    ostringstream logstrm;
    ostringstream errstrm;
    logstrm << BaseFileName << ".log";
    errstrm << BaseFileName << ".err";
    string logFilename = logstrm.str();
    string errFilename = errstrm.str();
    mainlog.open(logFilename.c_str());
    errlog.open(errFilename.c_str());
    
    
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

