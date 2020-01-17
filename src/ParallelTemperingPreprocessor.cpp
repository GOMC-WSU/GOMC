/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include "ParallelTemperingPreprocessor.h"

#if GOMC_LIB_MPI
ParallelTemperingPreprocessor::ParallelTemperingPreprocessor( int argc, 
                                                              char *argv[]){
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);                                                              
  
  //CHECK IF ARGS/FILE PROVIDED IN CMD LINE
  if (argc < 2) {
    std::cout << "Error: Input parameter file (*.dat or *.conf) not specified on command line!\n";
    MPI_Finalize();
    exit(EXIT_FAILURE);
  } else {
    if(argc == 2) {
      //FIRST PARAMETER WILL BE FILE NAME
      inputFileStringMPI = argv[1];
    } else {
      //SECOND PARAMETER WILL BE FILE NAME
      inputFileStringMPI = argv[2];

      if(argv[1][0] == '+' && argv[1][1] == 'p') {
      // placeholder
      } else {
        std::cout << "Error: Undefined command to set number of threads!\n";
        std::cout << "Use +p# command to set number of threads.\n";
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
    }
    //OPEN FILE
    inputFileReaderMPI.open(inputFileStringMPI.c_str(), ios::in | ios::out);

    //CHECK IF FILE IS OPENED...IF NOT OPENED EXCEPTION REASON FIRED
    if (!inputFileReaderMPI.is_open()) {
      std::cout << "Error: Cannot open/find " << inputFileStringMPI <<
                " in the directory provided!\n";
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }

    //CLOSE FILE TO NOW PASS TO SIMULATION
    inputFileReaderMPI.close();

    if(checkIfParallelTempering(inputFileStringMPI.c_str())){
      checkIfValid(inputFileStringMPI.c_str());
      if(worldRank >= getNumberOfReplicas(inputFileStringMPI.c_str())){
        std::cout << "You may not request more processes (" << worldSize
          << ") than there are replicas(" << getNumberOfReplicas(inputFileStringMPI.c_str()) 
          << ")! Exiting this process[" << worldRank << "]!\n";
          MPI_Finalize();
          exit(EXIT_FAILURE);
      } else {
      #if ENSEMBLE == NPT
        pathToReplicaDirectory = setupReplicaDirectoriesAndRedirectSTDOUTToFile  (  getMultiSimFolderName(inputFileStringMPI.c_str()),
                                                                                    getTemperature(inputFileStringMPI.c_str(), worldRank),
                                                                                    "",
                                                                                    getPressure(inputFileStringMPI.c_str(), worldRank));                                                                            
      #elif ENSEMBLE == GCMC
        pathToReplicaDirectory = setupReplicaDirectoriesAndRedirectSTDOUTToFile  (  getMultiSimFolderName(inputFileStringMPI.c_str()),
                                                                                    getTemperature(inputFileStringMPI.c_str(), worldRank),
                                                                                    getChemicalPotential(inputFileStringMPI.c_str(), worldRank),
      #elif ENSEMBLE == GEMC
      if(isGEMCNPT){
        pathToReplicaDirectory = setupReplicaDirectoriesAndRedirectSTDOUTToFile  (  getMultiSimFolderName(inputFileStringMPI.c_str()),
                                                                                    getTemperature(inputFileStringMPI.c_str(), worldRank),
                                                                                    "",
                                                                                    getPressure(inputFileStringMPI.c_str(), worldRank)); 
      } else {
        pathToReplicaDirectory = setupReplicaDirectoriesAndRedirectSTDOUTToFile  (  getMultiSimFolderName(inputFileStringMPI.c_str()),
                                                                                    getTemperature(inputFileStringMPI.c_str(), worldRank),
                                                                                    "",
                                                                                    ""; 
      #else
        pathToReplicaDirectory = setupReplicaDirectoriesAndRedirectSTDOUTToFile  (  getMultiSimFolderName(inputFileStringMPI.c_str()),
                                                                                    getTemperature(inputFileStringMPI.c_str(), worldRank),
                                                                                    "",
                                                                                    ""); 
      #endif
      } 
    }
  }
}

bool ParallelTemperingPreprocessor::checkIfValidRank(){
  if (worldRank >= getNumberOfReplicas(inputFileStringMPI.c_str())){
    return false;
  } else {
    return true;
  }
}


bool ParallelTemperingPreprocessor::checkIfParallelTempering(const char *fileName){
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(fileName);
  isParallelTemperingInTemperature = false;
  isParallelTemperingInChemicalPotential = false;
  isParallelTemperingInFugacity = false;
  //isParallelTemperingInFreeEnergyCoulomb = false;
  //isParallelTemperingInFreeEnergyVDW = false;
  isParallelTemperingInPressure = false;
  isGEMCNPT = false;

  while(reader.readNextLine(line)) {
    if(line.size() == 0)
      continue;
    if(CheckString(line[0], "Temperature")) {
      if (line.size() > 2)
        isParallelTemperingInTemperature = true;
    } 
    #if ENSEMBLE == GCMC
    else if (CheckString(line[0], "ChemPot")) {
      if (line.size() > 3)
        isParallelTemperingInChemicalPotential = true;
    } else if (CheckString(line[0], "Fugacity")) {
      if (line.size() > 3)
        isParallelTemperingInFugacity = true;
    }
    #endif
    #if ENSEMBLE == GEMC || ENSEMBLE == NPT
    else if (CheckString(line[0], "Pressure")) {
      if (line.size() > 3)
        isParallelTemperingInPressure = true;
    } 
    #endif
    #if ENSEMBLE == GEMC
    else if(CheckString(line[0], "GEMC")) {
      if(CheckString(line[1], "NPT")) {      
        isGEMCNPT = true;
      }
    }
    #endif
    /*else if (CheckString(line[0], "LambdaCoulomb")) {
      if (line.size() > 2)
        isParallelTemperingInFreeEnergyCoulomb = true;
    } else if (CheckString(line[0], "LambdaVDW")) {
      if (line.size() > 2)
        isParallelTemperingInFreeEnergyVDW = true;
    }*/
    // Clear and get ready for the next line
    line.clear();
  }

  #if ENSEMBLE == GEMC
  isParallelTemperingInPressure = isParallelTemperingInPressure && isGEMCNPT;
  #endif
/*
  return isParallelTemperingInTemperature || isParallelTemperingInChemicalPotential || isParallelTemperingInPressure || 
          isParallelTemperingInFugacity || isParallelTemperingInFreeEnergyCoulomb || isParallelTemperingInFreeEnergyVDW;
*/  
  return isParallelTemperingInTemperature || isParallelTemperingInChemicalPotential || isParallelTemperingInPressure || 
          isParallelTemperingInFugacity;
}

void ParallelTemperingPreprocessor::checkIfValid(const char *fileName){
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(fileName);
  int numberOfTemperatures = 0;
  int numberOfPressures = 0;
  vector < int > numberOfChemPots;
  vector < int > numberOfFugacity;
  //int numberOfLambdaCoulombs = 0;
  //int numberOfLambdaVDWs = 0;

  while(reader.readNextLine(line)) {
    if(line.size() == 0)
      continue;
    if(CheckString(line[0], "Temperature")) {
      numberOfTemperatures = line.size() - 1;
    } else if (CheckString(line[0], "Pressure")){
      numberOfPressures = line.size() - 1;
    } else if (CheckString(line[0], "ChemPot")) {
      numberOfChemPots.push_back(line.size() - 2);
    } else if (CheckString(line[0], "Fugacity")) {
      numberOfFugacity.push_back(line.size() - 2);
    }/* else if (CheckString(line[0], "LambdaCoulomb")) {
      numberOfLambdaCoulombs = line.size() - 1;
    }  else if (CheckString(line[0], "LambdaVDW")) {
      numberOfLambdaVDWs = line.size() - 1;
    }*/
    // Clear and get ready for the next line
    line.clear();
  }

  #if ENSEMBLE == GCMC
  if(isParallelTemperingInChemicalPotential){
    int numRepl = getNumberOfReplicas(inputFileStringMPI.c_str());
    for( vector < int >::iterator it = numberOfChemPots.begin(); it != numberOfChemPots.end(); ++it ){
      if (*it > 1 && *it != numRepl){
        std::cout << "Error: Unequal number of Chemical Potentials between different Residues!\n";
        std::cout << "Please provide either the same number of Chemical Potentials for each residue\n";
        std::cout << "or provide only one value for a residue.\n";
        std::cout << "Number of Chemical Potentials provided: " << numRepl << "\n";
        std::cout << "Number of Chemical Potentials provided for another residue: " << *it << "\n";
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
    }
  }

  if(isParallelTemperingInFugacity){
    int numRepl = getNumberOfReplicas(inputFileStringMPI.c_str());
    for( vector < int >::iterator it = numberOfFugacity.begin(); it != numberOfFugacity.end(); ++it ){
      if (*it > 1 && *it != numRepl){
        std::cout << "Error: Unequal number of fugacities between different Residues!\n";
        std::cout << "Please provide either the same number of fugacities for each residue\n";
        std::cout << "or provide only one value for a residue.\n";
        std::cout << "Number of fugacities provided: " << numRepl << "\n";
        std::cout << "Number of fugacities provided for another residue: " << *it << "\n";
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
    }
  }

  if (isParallelTemperingInChemicalPotential  && isParallelTemperingInTemperature){
    for( vector < int >::iterator it = numberOfChemPots.begin(); it != numberOfChemPots.end(); ++it ){
      if (*it > 1 && numberOfTemperatures > 1 && *it != numberOfTemperatures){
        std::cout << "Error: Unequal number of temperatures and chemical potentials in Multicanonical!\n";
        std::cout << "If you only want to only sample mu-space or temperature-space\n";
        std::cout << "provide only one temperature or only one chemical potential.\n";
        std::cout << "Number of temperatures provided: " << numberOfTemperatures << "\n";
        std::cout << "Number of chemical potentials provided: " << *it << "\n";
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
    }
  }

  if (isParallelTemperingInFugacity  && isParallelTemperingInTemperature){
    for( vector < int >::iterator it = numberOfFugacity.begin(); it != numberOfFugacity.end(); ++it ){
      if (*it > 1 && numberOfTemperatures > 1 && *it != numberOfTemperatures){
        std::cout << "Error: Unequal number of temperatures and fugacities in Multicanonical!\n";
        std::cout << "If you only want to only sample mu-space or temperature-space\n";
        std::cout << "provide only one temperature or only one chemical potential.\n";
        std::cout << "Number of temperatures provided: " << numberOfTemperatures << "\n";
        std::cout << "Number of fugacities provided: " << *it << "\n";
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
    }
  }

  #endif

  #if ENSEMBLE == NPT || ENSEMBLE == GEMC

  if (isParallelTemperingInTemperature && isParallelTemperingInPressure){
    if (numberOfTemperatures != numberOfPressures){
        std::cout << "Error: Unequal number of temperatures and pressures in Multicanonical!\n";
        std::cout << "If you only want to only sample pressure-space or temperature-space\n";
        std::cout << "provide only one temperature or only one pressure.\n";
        std::cout << "Number of temperatures provided: " << numberOfTemperatures << "\n";
        std::cout << "Number of pressures provided: " << numberOfPressures << "\n";
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
  }

  #endif
/*
  if ( numberOfLambdaCoulombs != numberOfLambdaVDWs){
      std::cout << "Error: Unequal number of LambdaCoulombs and LambdaVDWs in Free Energy calculation!\n";
      std::cout << "Number of temperatures provided: " << numberOfLambdaCoulombs << "\n";
      std::cout << "Number of temperatures provided: " << numberOfLambdaVDWs << "\n";
      MPI_Finalize();
      exit(EXIT_FAILURE);
  }
*/
}

int ParallelTemperingPreprocessor::getNumberOfReplicas(const char *fileName){
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(fileName);
  int numberOfTemperatures = 0;
  vector < int > numberOfChemPots;
  //int numberOfLambdaCoulombs = 0;
  //int numberOfLambdaVDWs = 0;
  int numberOfPressures = 0;
  vector < int > numberOfFugacity;

  int numberOfReplicas = 0;

  while(reader.readNextLine(line)) {
    if(line.size() == 0)
      continue;
    if(CheckString(line[0], "Temperature")) {
      numberOfTemperatures = line.size() - 1;
    } if(CheckString(line[0], "Pressure")) {
      numberOfPressures = line.size() - 1;
    } else if (CheckString(line[0], "ChemPot")) {
      numberOfChemPots.push_back(line.size() - 2);
    } else if (CheckString(line[0], "Fugacity")) {
      numberOfFugacity.push_back(line.size() - 2);
    } /*else if (CheckString(line[0], "LambdaCoulomb")) {
      numberOfLambdaCoulombs = line.size() - 1;
    }  else if (CheckString(line[0], "LambdaVDW")) {
      numberOfLambdaVDWs = line.size() - 1;
    }*/
    #if ENSEMBLE == GEMC
    else if(CheckString(line[0], "GEMC")) {
      if(CheckString(line[1], "NPT")) {      
        isGEMCNPT = true;
      }
    }
    #endif
    // Clear and get ready for the next line
    line.clear();
  }

  #if ENSEMBLE == GCMC
  for( vector < int >::iterator it = numberOfChemPots.begin(); it != numberOfChemPots.end(); ++it ){
    numberOfReplicas = std::max(numberOfReplicas, *it);
  }

  for( vector < int >::iterator it = numberOfFugacity.begin(); it != numberOfFugacity.end(); ++it ){
    numberOfReplicas = std::max(numberOfReplicas, *it);
  }
  #endif

  numberOfReplicas = std::max(numberOfReplicas, numberOfTemperatures);
  
  #if ENSEMBLE == NPT 
  numberOfReplicas = std::max(numberOfReplicas, numberOfPressures);
  #endif 
  
  #if ENSEMBLE == GEMC
  if(isGEMCNPT){
    numberOfReplicas = std::max(numberOfReplicas, numberOfPressures);
  }
  #endif  
  //numberOfReplicas = std::max(numberOfReplicas, numberOfLambdaCoulombs);
  //numberOfReplicas = std::max(numberOfReplicas, numberOfLambdaVDWs);

  return numberOfReplicas;
}

std::string ParallelTemperingPreprocessor::getMultiSimFolderName(const char *fileName){
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(fileName);
  string folderName;

  while(reader.readNextLine(line)) {
    if(line.size() == 0){
      continue;
    } else if(line[0] == "MultiSimFolderName"){
        std::stringstream ss;
        for (int i = 1; i < line.size(); i++){
            if (line[i] == " ") {
                ss << '_';
            } else {    
                ss << line[i];
                if (i+1 != line.size())
                    ss << "_";
            }
        }
        folderName = ss.str();
    } 
    // Clear and get ready for the next line
    line.clear();
  }
  if (folderName.empty()){
    folderName = "MultiSimFolderName";
  }
  return folderName;
}

std::string ParallelTemperingPreprocessor::getTemperature(const char *fileName, int worldRank){
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(fileName);
  string temperature;


  while(reader.readNextLine(line)) {
    if(line.size() == 0){
      continue;
    } else if(line[0] == "Temperature"){
        if (line.size() > 2){
          temperature = "temp_" + line[worldRank+1];
        } else {
          //temperature = line[1];
        }
    } 
    // Clear and get ready for the next line
    line.clear();
  }
  return temperature;
}

std::string ParallelTemperingPreprocessor::getPressure(const char *fileName, int worldRank){
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(fileName);
  string pressure;


  while(reader.readNextLine(line)) {
    if(line.size() == 0){
      continue;
    } else if(line[0] == "Pressure"){
        if (line.size() > 2){
          pressure = "_pres_" + line[worldRank+1];
        } else {
          //pressure = "_pres_" + line[1];
        }
    } 
    // Clear and get ready for the next line
    line.clear();
  }
  return pressure;
}

std::string ParallelTemperingPreprocessor::getChemicalPotential(const char *fileName, int worldRank){
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(fileName);
  std::stringstream chemPotStream;
  std::string resName;
  std::string val;

  while(reader.readNextLine(line)) {
    if(line.size() == 0){
      continue;
    } else if(CheckString(line[0], "ChemPot") || CheckString(line[0], "Fugacity")) {
      if (line.size() > 3){
        resName = line[1];
        val = line[2 + worldRank];
        chemPotStream << "_" << resName << "_" << val;
      } else if(line.size() != 3) {
        std::cout << "Error: Chemical potential / fugacity parameters are not specified!\n";
        MPI_Finalize();
        exit(EXIT_FAILURE);
      } else {
        //resName = line[1];
        //val = line[2];
        //chemPotStream << "_" << resName << "_" << val;
      }
    } 
    // Clear and get ready for the next line
    line.clear();
  }
  return chemPotStream.str();
}

 std::string ParallelTemperingPreprocessor::setupReplicaDirectoriesAndRedirectSTDOUTToFile(std::string multiSimTitle, std::string temperature, std::string chemPot, std::string pressure){   
    std::stringstream replicaTemp;
    replicaTemp << temperature << chemPot << pressure;
    std::string replicaDirectory = replicaTemp.str();
    return mkdirWrapper(multiSimTitle, replicaDirectory);
 }
 
std::string ParallelTemperingPreprocessor::mkdirWrapper(std::string multisimDirectoryName, string replicaDirectoryName){
    std::stringstream replicaStream;

    replicaStream << "." << OS_SEP << multisimDirectoryName << OS_SEP 
                << replicaDirectoryName << OS_SEP;
    std::string replicaDirectoryPath = replicaStream.str();
     
    //printf("Creating directory : %s\n", multisimDirectoryName.c_str());
    mkdir(multisimDirectoryName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(replicaDirectoryPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
    std::string pathToReplicaDirectory = replicaStream.str();
    replicaStream << "ConsoleOut.dat";
    std::string pathToReplicaLogFile = replicaStream.str();
    if(worldRank == 0){
      std::cout << "Monitor progress of your simulation by navigating to a replica output directory and issuing:\n"
        << "\t$ tail -f \"YourUniqueFileName\".console" << std::endl; 
    }
    freopen(pathToReplicaLogFile.c_str(),"w",stdout); 
    return pathToReplicaDirectory;
    
 }

bool ParallelTemperingPreprocessor::CheckString(string str1, string str2)
{
  for(int k = 0; k < str1.length(); k++) {
    str1[k] = toupper(str1[k]);
  }

  for(int j = 0; j < str2.length(); j++) {
    str2[j] = toupper(str2[j]);
  }

  return (str1 == str2);
}

MultiSim::MultiSim(ParallelTemperingPreprocessor & pt) : 
  worldSize(pt.worldSize), worldRank(pt.worldRank), pathToReplicaDirectory(pt.pathToReplicaDirectory)
{}

#endif