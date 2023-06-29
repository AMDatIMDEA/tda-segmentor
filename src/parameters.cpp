/*********************************************************************

TDA-Segmentor     -     A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:                       Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
                 Aditya Vasudevan (adityavv.iitkgp@gmail.com)
                 Maciek Haranczyk (maciej.haranczyk@imdea.org)
                 IMDEA Materiales Institute
 
**********************************************************************/

#include "parameters.h"
#include "logger.h"
#include "string.h"
#include "instructions.h"


parameters::parameters() :
moduleNames(),
inputfilename(),
useSuperCell(false),
useTBB(false),
saveLogFile(false),
segmentationFlag(false),
persistenceThreshold(0.0),
probeRadius(0.0)
{
    
    
};




void parameters::printinvocation(int nargs, char **args)
{
    int i;
    logger::mainlog << "Input received: " << endl;
    for ( i = 0; i < nargs; ++i )
    {
        logger::mainlog << args[i] << " " ;
    }
    logger::mainlog << "\n" ;
}




bool parameters::checkinputfile(){
    
    std::string::size_type const p(inputfilename.find_last_of('.'));
    if (p != std::string::npos){
        basefilename = inputfilename.substr(0, p);
        extensionname = inputfilename.substr(p,inputfilename.length());
    }

    // More extension names need to be added here for energy and density calculations
    if (extensionname != ".cube" || extensionname.empty()) {
        std::cout << "\n\nInput file check failed - check invocation syntax" << endl;
        return false;
    }
        
    return true;
    
}



void parameters::parser(int nargs, char **args)
{
    int i, k, j = 0;
    double dummy;
    
    inputfilename = args[nargs-1];
    if (!checkinputfile()){
        
        std::cout << "Input file needs to be the last argument!!" << endl;
        std::cout << "inputfilename received is : " << inputfilename << endl;
        std::cout << "\n\n";
        printInstructions();
        exit(0);
    }
    
    for ( i = 1; i < nargs ; i++ )
    {
        k = 0;
        while ( args[i][k] )
        {
            args[i][k] = tolower( args[i][k] );
            k++;
        }
        
        if ( strcmp(args[i] , "-module") == 0 && i < nargs - 1)
        {
            size_t nModules = moduleNames.size();
            moduleNames.resize(nModules+1);
            char ch[256];
            sscanf(args[i+1], "%s", ch);
            moduleNames[nModules] = ch;
            i++;j++;
        }
        
        if ( strcmp(args[i] , "-persistencethreshold") == 0 && i < nargs - 1)
        {
            sscanf(args[i+1], "%lf", &persistenceThreshold);
            i++; j++;
        }
        if ( strcmp(args[i] , "-proberadius") == 0 && i < nargs - 1)
        {
            sscanf(args[i+1], "%lf", &probeRadius);
            i++; j++;
        }
        
        if ( strcmp(args[i] , "-usesupercell") == 0 && i < nargs - 1)
        {
            useSuperCell = true;
            j++;
        }
        
        if ( strcmp(args[i] , "-usetbb") == 0 && i < nargs - 1)
        {
            useSuperCell = true;
            j++;
        }
        
        if ( strcmp(args[i] , "-savelogfile") == 0 && i < nargs - 1)
        {
            saveLogFile = true;
            j++;
        }
        
    }
    
    for (size_t i = 0; i < moduleNames.size(); i++){
        
        if (moduleNames[i] == "segmentation"
            || moduleNames[i] == "accessiblevoidspace"
            || moduleNames[i] == "voidsegmentation"
            || moduleNames[i] == "solidsegmentation")
        {
            segmentationFlag = true;
        }
        
    }
    
}




void parameters::writetoLogFile() {
    
    logger::mainlog << "Parsed input :" << endl;
    logger::mainlog << "Number of modules: " << moduleNames.size() << endl;
    for (size_t i = 0; i < moduleNames.size(); i++) {
        logger::mainlog << "Module " << (i+1) << ": "<< moduleNames[i] << endl;
        
        if (moduleNames[i] != "segmentation"
            and moduleNames[i] != "accessiblevoidspace"
            and moduleNames[i] != "voidsegmentation"
            and moduleNames[i] != "solidsegmentation"
            and moduleNames[i] != "persistencecurve")
        {
            logger::mainlog << "Module " << (i+1) << ": "<< moduleNames[i] << " is not implemented!!" << endl;
            logger::errlog << "Module " << (i+1) << ": "<< moduleNames[i] << " is not implemented!!" << endl;
            printInstructions();
            exit(0);
        }
    }
    
    for (size_t i = 0; i < moduleNames.size(); i++) {
        if (moduleNames[i] == "segmentation"){
            if (persistenceThreshold == 0.0)
            {
                logger::mainlog << "Persistence Threshold is not given and will be chosen automatically!" << endl;
                persistenceThreshold = 0.1;
            } else
            {
                logger::mainlog << "Persistence Threshold : " << persistenceThreshold << endl;
            }
            
        }
        if (moduleNames[i] == "accessiblevoidspace"){
            if (persistenceThreshold == 0.0)
            {
                logger::mainlog << "Persistence Threshold is not given and will be chosen automatically!" << endl;
                persistenceThreshold = 0.1;
            } else {
                logger::mainlog << "Persistence Threshold : " << persistenceThreshold << endl;
            }
            
            if (probeRadius == 0.0)
            {
                logger::mainlog << "Radius of the probe atom is not given, and is chosen to be zero!" << endl;
            } else
            {
                logger::mainlog << "Probe radius          :  " << probeRadius << endl;
            }
        }
        
    }
    
    std::string dummy("False");
    if (useSuperCell) dummy = "True";
    logger::mainlog << "Use Super Cell : " << dummy << endl;
    dummy = "False";
    if (useTBB) dummy = "True";
    logger::mainlog << "Use TBB : " << dummy << endl;
    dummy = "False";
    if (saveLogFile) dummy = "True";
    logger::mainlog << "Save log file : " << dummy << endl;
    
    
}