/*********************************************************************

TDA-Segmentor    A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:         Aditya Vasudevan (adityavv.iitkgp@gmail.com)
                 Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
                 Maciek Haranczyk (maciej.haranczyk@imdea.org)
                 IMDEA Materiales Institute
 
**********************************************************************/

#include "logger.h"


ofstream      logger::mainlog;
ofstream      logger::errlog;


void logger::openLogfiles(const parameters &p) {
    
    DIR *dp;
    struct dirent *dirp;
    string curdir = ".";
    if ( (dp = opendir(curdir.c_str())) == NULL)
    {
        cout << "Error(" << errno << ") opening " << curdir << endl;
    }

    vector<string> files;
    while ((dirp = readdir(dp)) != NULL)
    {
        files.push_back(string(dirp->d_name));
    }
    closedir(dp);


    // see if there are any iris-save.* files
    if (p.saveLogFile)
    {
        int num = -1;
        for (unsigned a=0; a<files.size(); a++)
        {
            if ( files[a].substr(0,14) == "save-segmentor")
            {
                string t = files[a].substr(15, files[a].size());
                num = std::max<int>(num, stoi(t));
            }
        }

        for (unsigned a=0; a<files.size(); a++)
        {
            if ( files[a] == p.basefilename + ".log" )
            {
                ostringstream strm;
                strm << "save-segmentor." << num+1;
                if( std::rename((p.basefilename + ".log").c_str(), strm.str().c_str()) != 0)
                    cout << "\n Unable to save log file";
            }
        }
    }
    
    
    
    
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

