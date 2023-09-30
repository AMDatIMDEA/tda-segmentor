/*********************************************************************

TDA-Segmentor    A segmentation tool for porous structures using the topology
                 toolkit (https://topology-tool-kit.github.io/)

Authors:         Aditya Vasudevan (adityavv.iitkgp@gmail.com)
                 Jorge Zorrilla Prieto (jorge.zorrilla.prieto@gmail.com)
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
useAllCores(false),
saveLogFile(false),
segmentationFlag(false),
ftmTreeFlag(false),
persistenceThreshold(0.0),
probeRadius(0.0),
arrayName()
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
    
    extensionname = vtksys::SystemTools::GetFilenameLastExtension(inputfilename);
    basefilename = vtksys::SystemTools::GetFilenameWithoutExtension(inputfilename);
    
    // More extension names need to be added here for energy and density calculations
    if ((extensionname != ".cube" && extensionname != ".vti") || extensionname.empty()) {
        std::cout << "Input file check failed - check invocation syntax" << endl;
        return false;
    }
    
    std::stringstream str;
    str << "segmentor-" << basefilename << ".results";
    Directory = str.str();

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
        
        if ( strcmp(args[i] , "-scalar") == 0 && i < nargs - 1)
        {
            char ch[256];
            sscanf(args[i+1], "%s", ch);
            if (strcmp(ch,"distance") == 0) {
                arrayName = "This is distance grid"; 
            } else if (strcmp(ch,"energy") == 0) {
                arrayName = "Potential Energy"; 
            } else {
                arrayName = ch; 
            }
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
        
        if ( strcmp(args[i] , "-useallcores") == 0 && i < nargs - 1)
        {
            useAllCores = true;
            j++;
        }
        
        if ( strcmp(args[i] , "-savelogfile") == 0 && i < nargs - 1)
        {
            saveLogFile = true;
            j++;
        }
        
    }
    
    // by default run the segmentation module
    if (moduleNames.empty()) moduleNames.push_back("segmentation");
    
    for (size_t i = 0; i < moduleNames.size(); i++){
        
        if (moduleNames[i] == "segmentation"
            || moduleNames[i] == "accessiblevoidspace"
            || moduleNames[i] == "voidsegmentation"
            || moduleNames[i] == "solidsegmentation"
            || moduleNames[i] == "accessiblevoidgraph"
            || moduleNames[i] == "accessiblesolidgraph")
        {
            segmentationFlag = true;
        }

        if (moduleNames[i] == "graph" || 
            moduleNames[i] == "accessiblevoidgraph" ||
            moduleNames[i] == "accessiblesolidgraph")
        {
            ftmTreeFlag = true;
        }
        
    }
    
}




void parameters::writetoLogFile() {
    
    logger::mainlog << "\n\nParsed input :" << endl;
    logger::mainlog << "Number of modules                                            : " << moduleNames.size() << endl;
    for (size_t i = 0; i < moduleNames.size(); i++) {
        logger::mainlog << "Module " << (i+1) << "                                                     : "<< moduleNames[i] << endl;
        
        if (moduleNames[i] != "segmentation"
            and moduleNames[i] != "accessiblevoidspace"
            and moduleNames[i] != "voidsegmentation"
            and moduleNames[i] != "solidsegmentation"
            and moduleNames[i] != "persistencecurve"
            and moduleNames[i] != "graph"
            and moduleNames[i] != "accessiblevoidgraph"
            and moduleNames[i] != "accessiblesolidgraph")
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
            } else
            {
                logger::mainlog << "Persistence Threshold for segmentation module                : " << persistenceThreshold << endl;
            }
            
        }
        if (moduleNames[i] == "accessiblevoidspace" || moduleNames[i] == "accessiblevoidgraph"){
            if (persistenceThreshold == 0.0)
            {
                logger::mainlog << "Persistence Threshold is not given and will be chosen automatically!" << endl;
            } else {
                logger::mainlog << "Persistence Threshold for accessible void space module       : " << persistenceThreshold << endl;
            }
            
            if (probeRadius == 0.0)
            {
                logger::mainlog << "Radius of the probe atom is not given, and is chosen to be zero!" << endl;
            } else
            {
                logger::mainlog << "Probe radius                                                 : " << probeRadius << endl;
            }
        }
        
    }
    
    std::string givenArrayName;
    if (arrayName.empty()) {
        givenArrayName = "Chosen Automatically";
    } else {
        givenArrayName = arrayName;
    }
    logger::mainlog << "Scalar array given for segmentation                          : " << givenArrayName << endl;
    
    std::string dummy("False");
    if (useSuperCell) dummy = "True";
    logger::mainlog << "Use Super Cell                                               : " << dummy << endl;
    dummy = "False";
    if (useAllCores) dummy = "True";
    logger::mainlog << "Use all Cores                                                : " << dummy << endl;
    dummy = "False";
    if (saveLogFile) dummy = "True";
    logger::mainlog << "Save log file                                                : " << dummy << endl;
    
}
