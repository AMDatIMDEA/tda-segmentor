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



parameters::parameters() : moduleNames(),
                           inputfilename(),
                           useSuperCell(false),
                           useAllCores(false),
                           saveLogFile(false),
                           segmentationFlag(false),
                           writeFractionalGrid(false),
                           persistenceThreshold(0.0),
                           probeRadius(0.0),
                           arrayName(){

                           };




void parameters::printinvocation(int nargs, char **args)
{
    int i;
    logger::mainlog << "Input received: " << endl;
    for (i = 0; i < nargs; ++i)
    {
        logger::mainlog << args[i] << " ";
    }
    logger::mainlog << "\n";
}




bool parameters::checkinputfile()
{

    extensionname = vtksys::SystemTools::GetFilenameLastExtension(inputfilename);
    basefilename = vtksys::SystemTools::GetFilenameWithoutExtension(inputfilename);

    // More extension names need to be added here for energy and density calculations
    if ((extensionname != ".cube" && extensionname != ".vti") || extensionname.empty())
    {
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
    // The last argument is the input file
    inputfilename = args[nargs - 1];

    if (!checkinputfile())
    {
        std::cout << "Input file needs to be the last argument!!" << endl;
        std::cout << "Last argument received is : " << inputfilename << endl;
        std::cout << "\n\n";
        printInstructions();
        exit(0);
    }

    // store in a vector each of the commands starting with -
    vector<vector<string>> commands = vector<vector<string>>();
    vector<string> currentCommand = vector<string>();

    int argcount = 1;
    while (argcount < nargs - 1)
    {
        if (args[argcount][0] == '-')
        {
            if (currentCommand.size() != 0)
                commands.push_back(currentCommand);
            currentCommand = vector<string>();
        }

        currentCommand.push_back(args[argcount]);
        argcount++;
    }
    if (currentCommand.size() != 0)
        commands.push_back(currentCommand);

    
    int numCommands = commands.size();

    bool error = true;

    for (int i = 0; i < numCommands; i++)
    {
        vector<string> command = commands[i];

        if (command[0].compare("-msc") == 0)
        {
            moduleNames.push_back(module::SEGMENTATION);
            segmentationFlag = true;
            error = false;

            if (command.size() < 2)
            {
                std::cout << "Persistence threshold is not provided as an input and will be chosen automatically as 1 percent of maximum persistence." << endl;
                std::cout << "However it is strongly recommended to provide persistence threshold as an input, and can be provided as -msc followed by a value." << endl;
                std::cout << "(HINT) Try plotting the persistence curve first:  tda-segmentor -pc " << inputfilename << " to get an idea of persistence threshold." << endl;

                persistenceThreshold = 0.0;

            } else if (command.size() == 2) {

                persistenceThreshold = std::strtod(command[1].data(), NULL);

            } else if (command.size() > 2){

                std::cout <<"(ERROR) The option -msc only accepts one optional parameter - the persistence threshold and you have entered " << (command.size() - 1 ) << " parameters." << endl; 
                std::cout << "(HINT) Try plotting the persistence curve first:  tda-segmentor -pc " << inputfilename << " to get an idea of persistence threshold." << endl;
                exit(0);
            }

        }


        if (command[0].compare("-vs") == 0)
        {
            moduleNames.push_back(module::VOID_SEGMENTATION);
            segmentationFlag = true;
            error = false;

           
            if (command.size() < 2)
            {
                if (persistenceThreshold == 0.0) {
 
                    std::cout << "Persistence threshold is not provided as an input and will be chosen automatically as 1 percent of maximum persistence." << endl;
                    std::cout << "However it is strongly recommended to provide persistence threshold as an input, and can be provided as -vs followed by a value." << endl;
                    std::cout << "(HINT) Try plotting the persistence curve first:  tda-segmentor -pc " << inputfilename << " to get an idea of persistence threshold." << endl;

                } 


            } else if (command.size() == 2) {

                persistenceThreshold = std::strtod(command[1].data(), NULL);

            } else if (command.size() > 2){

                std::cout <<"(ERROR) The option -vs only accepts one optional parameter - the persistence threshold and you have entered " << (command.size() - 1 ) << " parameters." << endl; 
                std::cout << "(HINT) Try plotting the persistence curve first:  tda-segmentor -pc " << inputfilename << " to get an idea of persistence threshold." << endl;
                exit(0);
            }
        }


        if (command[0].compare("-ss") == 0)
        {
            moduleNames.push_back(module::SOLID_SEGMENTATION);
            segmentationFlag = true;
            error = false;
            
            if (command.size() < 2)
            {
                if (persistenceThreshold == 0.0) {
 
                    std::cout << "Persistence threshold is not provided as an input and will be chosen automatically as 1 percent of maximum persistence." << endl;
                    std::cout << "However it is strongly recommended to provide persistence threshold as an input, and can be provided as -ss followed by a value." << endl;
                    std::cout << "(HINT) Try plotting the persistence curve first:  tda-segmentor -pc " << inputfilename << " to get an idea of persistence threshold." << endl;

                } 

            } else if (command.size() == 2) {

                persistenceThreshold = std::strtod(command[1].data(), NULL);

            } else if (command.size() > 2){

                std::cout <<"(ERROR) The option -ss only accepts one optional parameter - the persistence threshold and you have entered " << (command.size() - 1 ) << " parameters." << endl; 
                std::cout << "(HINT) Try plotting the persistence curve first:  tda-segmentor -pc " << inputfilename << " to get an idea of persistence threshold." << endl;

                exit(0);
            }
        }


        if (command[0].compare("-avs") == 0)
        {
            moduleNames.push_back(module::ACCESSIBLE_VOID_SPACE);
            segmentationFlag = true;
            error = false;

            if (command.size() < 2 || command.size() > 3){

                std::cout << "(ERROR) The command -avs accepts two parameters, an optional persistence threshold first, and a mandatory probe radius next and " << (command.size()-1) << " was provided." << endl;
                std::cout << "(HINT) Try plotting the persistence curve first:  tda-segmentor -pc " << inputfilename << " to get an idea of persistence threshold." << endl;

                exit(0);

            } else if (command.size() == 2 ){
                if (persistenceThreshold == 0.0){

                    std::cout << "The command -avs accepts two parameters, an optional persistence threshold first, and a mandatory probe radius next and " << (command.size() - 1) << " was provided." << endl;
                    std::cout << "Persistence threshold is chosen automatically as 1 percent of maximum persistence." << endl;
                    std::cout << "However it is strongly recommended to provide persistence threshold as an input, and can be given as -avs followed by the threshold and then the probe radius." << endl;
                    std::cout << "(HINT) Try plotting the persistence curve first:  tda-segmentor -pc " << inputfilename << " to get an idea of persistence threshold." << endl;

                }

                probeRadius = std::strtod(command[1].data(), NULL);

            } else if (command.size() == 3) {
                persistenceThreshold = std::strtod(command[1].data(), NULL);
                probeRadius = std::strtod(command[2].data(), NULL);

            }

        }


        if (command[0].compare("-avg") == 0)
        {
            moduleNames.push_back(module::ACCESSIBLE_VOID_GRAPH);
            segmentationFlag = true;
            error = false;

            if (command.size() < 2 || command.size() > 3){

                std::cout << "(ERROR) The command -avg accepts two parameters, an optional persistence threshold first, and a mandatory probe radius next and " << (command.size()-1) << " was provided." << endl;
                std::cout << "(HINT) Try plotting the persistence curve first:  tda-segmentor -pc " << inputfilename << " to get an idea of persistence threshold." << endl;

                exit(0);

            } else if (command.size() == 2 ){

                if (persistenceThreshold == 0.0){

                    std::cout << "The command -avg accepts two parameters, an optional persistence threshold first, and a mandatory probe radius next and " << (command.size() - 1) << " was provided." << endl;
                    std::cout << "Persistence threshold is chosen automatically as 1 percent of maximum persistence." << endl;
                    std::cout << "However it is strongly recommended to provide persistence threshold as an input, and can be given as -avg followed by the threshold and then the probe radius." << endl;
                    std::cout << "(HINT) Try plotting the persistence curve first:  tda-segmentor -pc " << inputfilename << " to get an idea of persistence threshold." << endl;
                }

                probeRadius = std::strtod(command[1].data(), NULL);


            } else if (command.size() == 3) {
                persistenceThreshold = std::strtod(command[1].data(), NULL);
                probeRadius = std::strtod(command[2].data(), NULL);

            }

        }


        if (command[0].compare("-asg") == 0)
        {
            moduleNames.push_back(module::ACCESSIBLE_SOLID_GRAPH);
            segmentationFlag = true;
            error = false;

            if (command.size() < 2)
            {

                if (persistenceThreshold == 0.0)
                {
                    std::cout << "Persistence threshold is not provided as an input and will be chosen automatically as 1 percent of maximum persistence." << endl;
                    std::cout << "However it is strongly recommended to provide persistence threshold as an input, and can be provided as -asg followed by a value." << endl;
                    std::cout << "(HINT) Try plotting the persistence curve first:  tda-segmentor -pc " << inputfilename << " to get an idea of persistence threshold." << endl;
                }

            } else if (command.size() == 2) {

                persistenceThreshold = std::strtod(command[1].data(), NULL);

            } else if (command.size() > 2){

                std::cout <<"(ERROR) The option -asg only accepts one optional parameter - the persistence threshold and you have entered " << (command.size() - 1 ) << " parameters." << endl; 
                std::cout << "(HINT) Try plotting the persistence curve first:  tda-segmentor -pc " << inputfilename << " to get an idea of persistence threshold." << endl;

                exit(0);
            }
        }


        if (command[0].compare("-pc") == 0)
        {
            moduleNames.push_back(module::PERSISTENCE_CURVE);
            segmentationFlag = false;
            error = false;
        }


        if (command[0].compare("-scalar") == 0)
        {
            if (strcmp(command[1].data(), "distance") == 0)
            {
                arrayName = "This is distance grid";
            }
            else if (strcmp(command[1].data(), "energy") == 0)
            {
                arrayName = "Potential Energy";
            }
            else
            {
                arrayName = command[1].data();
            }
        }
        

        if (command[0].compare("-usesupercell") == 0)
        {
            useSuperCell = true;
        }


        if (command[0].compare("-writefractionalgrid") == 0)
        {
            writeFractionalGrid = true;
        }


        if (command[0].compare("-useallcores") == 0)
        {
            useAllCores = true;
        }


        if (command[0].compare("-savelogfile") == 0)
        {
            saveLogFile = true;
        }


        if (command[0].compare("-help") == 0){
            printInstructions();
            exit(0);
        }

    }

    // by default run the segmentation module
    if (moduleNames.empty()) {
        moduleNames.push_back(module::SEGMENTATION);
        error = false;
        persistenceThreshold = 0.0; // If set to 0, then 1% of the maximum persistence is taken. 
        segmentationFlag = true;

        std::cout << "No options are provided, and by default -msc is taken with the persistence threshold chosen 1 percent of the maximum." << endl;
        std::cout << "However it is strongly recommended to generate segmentation using -msc followed by the persistence threshold." << endl;
        std::cout << "(HINT) Try plotting the persistence curve first:  tda-segmentor -pc " << inputfilename << " to get an idea of persistence threshold." << endl;

    }
    
    
    if (error){
        std::cout << "There has been an error with the invocation command !!" << endl;
        printInstructions();
        exit(0);
    }



}

void parameters::writetoLogFile()
{

    logger::mainlog << "\n\nParsed input :" << endl;
    logger::mainlog << "Number of modules                                            : " << moduleNames.size() << endl;


    for (size_t i = 0; i < moduleNames.size(); i++)
    {
        std::string name = moduleNameAsString(moduleNames[i]);
        logger::mainlog << "Module " << (i+1) << "                                                     : " << name << endl;

    }

    for (size_t i = 0; i < moduleNames.size(); i++)
    {
        if (moduleNames[i] == module::SEGMENTATION)
        {
            if (persistenceThreshold == 0.0)
            {
                logger::mainlog << "Persistence Threshold is not given and will be chosen automatically!" << endl;
            }
            else
            {
                logger::mainlog << "Persistence Threshold for segmentation module                : " << persistenceThreshold << endl;
            }
        }
        if (moduleNames[i] == module::ACCESSIBLE_VOID_SPACE || moduleNames[i] == module::ACCESSIBLE_VOID_GRAPH)
        {
            if (persistenceThreshold == 0.0)
            {
                logger::mainlog << "Persistence Threshold is not given and will be chosen automatically!" << endl;
            }
            else
            {
                logger::mainlog << "Persistence Threshold for accessible void space module       : " << persistenceThreshold << endl;
            }

                logger::mainlog << "Probe radius                                                 : " << probeRadius << endl;
        }
    }

    std::string givenArrayName;
    if (arrayName.empty())
    {
        givenArrayName = "Chosen Automatically";
    }
    else
    {
        givenArrayName = arrayName;
    }
    logger::mainlog << "Scalar array given for segmentation                          : " << givenArrayName << endl;

    std::string dummy("False");
    if (useSuperCell)
        dummy = "True";
    logger::mainlog << "Use Super Cell                                               : " << dummy << endl;
    dummy = "False";
    if (useAllCores)
        dummy = "True";
    logger::mainlog << "Use all Cores                                                : " << dummy << endl;
    dummy = "False";
    if (saveLogFile)
        dummy = "True";
    logger::mainlog << "Save log file                                                : " << dummy << endl;
    if (writeFractionalGrid)
        dummy = "True";
    logger::mainlog << "Write fractional grid                                        : " << dummy << endl;
}




std::string  parameters::moduleNameAsString(const module modulename){


    std::string stringVal;
    std::vector < std::pair<module,string> > namesOfModules {{module::SEGMENTATION, "segmentation"}, 
                                                             {module::VOID_SEGMENTATION, "void segmentation"},
                                                             {module::SOLID_SEGMENTATION, "solid segmentation"},
                                                             {module::ACCESSIBLE_VOID_SPACE, "accessible void space"},
                                                             {module::ACCESSIBLE_VOID_GRAPH, "accessible void graph"},
                                                             {module::ACCESSIBLE_SOLID_GRAPH, "accessible solid graph"},
                                                             {module::PERSISTENCE_CURVE, "persistence curve"}  };

    
    for (size_t j = 0; j < namesOfModules.size(); j++){
            if (namesOfModules[j].first == modulename){
                stringVal = namesOfModules[j].second;
            }

        }

    return stringVal;

}