/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/

#include "ParallelTemperingPreprocessor.h"

#if GOMC_LIB_MPI
ParallelTemperingPreprocessor::ParallelTemperingPreprocessor(int argc,
                                                             char *argv[]) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

  // CHECK IF ARGS/FILE PROVIDED IN CMD LINE
  if (argc < 2) {
    std::cout << "Error: Input parameter file (*.dat or *.conf) not specified "
                 "on command line!\n";
    MPI_Finalize();
    exit(EXIT_FAILURE);
  } else {
    if (argc == 2) {
      // FIRST PARAMETER WILL BE FILE NAME
      inputFileStringMPI = argv[1];
    } else {
      // SECOND PARAMETER WILL BE FILE NAME
      inputFileStringMPI = argv[2];

      if (argv[1][0] == '+' && argv[1][1] == 'p') {
        // placeholder
      } else {
        std::cout << "Error: Undefined command to set number of threads!\n";
        std::cout << "Use +p# command to set number of threads.\n";
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
    }
    // OPEN FILE
    inputFileReaderMPI.open(inputFileStringMPI.c_str(),
                            std::ios::in | std::ios::out);

    // CHECK IF FILE IS OPENED...IF NOT OPENED EXCEPTION REASON FIRED
    if (!inputFileReaderMPI.is_open()) {
      std::cout << "Error: Cannot open/find " << inputFileStringMPI
                << " in the directory provided!\n";
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }

    // CLOSE FILE TO NOW PASS TO SIMULATION
    inputFileReaderMPI.close();

    if (checkIfExpandedEnsemble(inputFileStringMPI.c_str())) {
      checkIfValid(inputFileStringMPI.c_str());
      if (worldRank >= getNumberOfReplicas(inputFileStringMPI.c_str())) {
        std::cout << "You may not request more processes (" << worldSize
                  << ") than there are replicas("
                  << getNumberOfReplicas(inputFileStringMPI.c_str())
                  << ")! Exiting this process[" << worldRank << "]!\n";
        MPI_Finalize();
        exit(EXIT_FAILURE);
      } else {
#if ENSEMBLE == GCMC
        setupReplicaDirectoriesAndRedirectSTDOUTToFile(
            inputFileStringMPI,
            getTemperature(inputFileStringMPI.c_str(), worldRank),
            getChemicalPotential(inputFileStringMPI.c_str(), worldRank));
#else
        setupReplicaDirectoriesAndRedirectSTDOUTToFile(
            inputFileStringMPI,
            getTemperature(inputFileStringMPI.c_str(), worldRank));
#endif
        restart = checkIfRestart(inputFileStringMPI.c_str());
        restartFromCheckpoint =
            checkIfRestartFromCheckpoint(inputFileStringMPI.c_str());
        parallelTemperingEnabled =
            checkIfParallelTemperingEnabled(inputFileStringMPI.c_str());
      }
    }
  }
}

bool ParallelTemperingPreprocessor::checkIfValidRank() {
  if (worldRank >= getNumberOfReplicas(inputFileStringMPI.c_str())) {
    return false;
  } else {
    return true;
  }
}

bool ParallelTemperingPreprocessor::checkIfExpandedEnsemble(
    const char *fileName) {
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(fileName);
  bool isParallelTemperingInTemperature = false;
  bool isParallelTemperingInChemicalPotential = false;
  bool isParallelTemperingInFreeEnergyCoulomb = false;
  bool isParallelTemperingInFreeEnergyVDW = false;

  while (reader.readNextLine(line)) {
    if (line.size() == 0)
      continue;
    if (checkString(line[0], "Temperature")) {
      if (line.size() > 2)
        isParallelTemperingInTemperature = true;
    } else if (checkString(line[0], "ChemPot")) {
      if (line.size() > 3)
        isParallelTemperingInChemicalPotential = true;
    } else if (checkString(line[0], "LambdaCoulomb")) {
      if (line.size() > 2)
        isParallelTemperingInFreeEnergyCoulomb = true;
    } else if (checkString(line[0], "LambdaVDW")) {
      if (line.size() > 2)
        isParallelTemperingInFreeEnergyVDW = true;
    }
    // Clear and get ready for the next line
    line.clear();
  }
  return isParallelTemperingInTemperature ||
         isParallelTemperingInChemicalPotential ||
         isParallelTemperingInFreeEnergyCoulomb ||
         isParallelTemperingInFreeEnergyVDW;
}

void ParallelTemperingPreprocessor::checkIfValid(const char *fileName) {
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(fileName);
  int numberOfTemperatures = 0;
  std::vector<int> numberOfChemPots;
  int numberOfLambdaCoulombs = 0;
  int numberOfLambdaVDWs = 0;

  while (reader.readNextLine(line)) {
    if (line.size() == 0)
      continue;
    if (checkString(line[0], "Temperature")) {
      numberOfTemperatures = line.size() - 1;
    } else if (checkString(line[0], "ChemPot")) {
      numberOfChemPots.push_back(line.size() - 2);
    } else if (checkString(line[0], "LambdaCoulomb")) {
      numberOfLambdaCoulombs = line.size() - 1;
    } else if (checkString(line[0], "LambdaVDW")) {
      numberOfLambdaVDWs = line.size() - 1;
    }
    // Clear and get ready for the next line
    line.clear();
  }

  for (std::vector<int>::iterator it = numberOfChemPots.begin();
       it != numberOfChemPots.end(); ++it) {
    if (*it > 1 && numberOfTemperatures > 1 && *it != numberOfTemperatures) {
      std::cout << "Error: Unequal number of temperatures and chemical "
                   "potentials in Multicanonical!\n";
      std::cout
          << "If you only want to only sample mu-space or temperature-space\n";
      std::cout
          << "provide only one temperature or only one chemical potential.\n";
      std::cout << "Number of temperatures provided: " << numberOfTemperatures
                << "\n";
      std::cout << "Number of chemical potentials provided: " << *it << "\n";
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
  }

  if (numberOfLambdaCoulombs != numberOfLambdaVDWs) {
    std::cout << "Error: Unequal number of LambdaCoulombs and LambdaVDWs in "
                 "Free Energy calculation!\n";
    std::cout << "Number of temperatures provided: " << numberOfLambdaCoulombs
              << "\n";
    std::cout << "Number of temperatures provided: " << numberOfLambdaVDWs
              << "\n";
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
}

int ParallelTemperingPreprocessor::getNumberOfReplicas(const char *fileName) {
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(fileName);
  int numberOfTemperatures = 0;
  std::vector<int> numberOfChemPots;
  int numberOfLambdaCoulombs = 0;
  int numberOfLambdaVDWs = 0;

  int numberOfReplicas = 0;

  while (reader.readNextLine(line)) {
    if (line.size() == 0)
      continue;
    if (checkString(line[0], "Temperature")) {
      numberOfTemperatures = line.size() - 1;
    } else if (checkString(line[0], "ChemPot")) {
      numberOfChemPots.push_back(line.size() - 2);
    } else if (checkString(line[0], "LambdaCoulomb")) {
      numberOfLambdaCoulombs = line.size() - 1;
    } else if (checkString(line[0], "LambdaVDW")) {
      numberOfLambdaVDWs = line.size() - 1;
    }
    // Clear and get ready for the next line
    line.clear();
  }

  for (std::vector<int>::iterator it = numberOfChemPots.begin();
       it != numberOfChemPots.end(); ++it) {
    numberOfReplicas = std::max(numberOfReplicas, *it);
  }
  numberOfReplicas = std::max(numberOfReplicas, numberOfTemperatures);
  numberOfReplicas = std::max(numberOfReplicas, numberOfLambdaCoulombs);
  numberOfReplicas = std::max(numberOfReplicas, numberOfLambdaVDWs);

  return numberOfReplicas;
}

std::string
ParallelTemperingPreprocessor::getInputFolderName(const char *fileName) {
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(fileName);
  std::string folderName;

  while (reader.readNextLine(line)) {
    if (line.size() == 0) {
      continue;
    } else if (line[0] == "InputFolderName") {
      std::stringstream ss;
      for (int i = 1; i < line.size(); i++) {
        if (line[i] == " ") {
          ss << '_';
        } else {
          ss << line[i];
          if (i + 1 != line.size())
            ss << "_";
        }
      }
      folderName = ss.str();
    }
    // Clear and get ready for the next line
    line.clear();
  }
  // This will default to the input files being in a flat directory
  if (folderName.empty()) {
    folderName = "";
  }
  return folderName;
}

std::string
ParallelTemperingPreprocessor::getOutputFolderName(const char *fileName) {
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(fileName);
  std::string folderName;

  while (reader.readNextLine(line)) {
    if (line.size() == 0) {
      continue;
    } else if (line[0] == "OutputFolderName") {
      std::stringstream ss;
      for (int i = 1; i < line.size(); i++) {
        if (line[i] == " ") {
          ss << '_';
        } else {
          ss << line[i];
          if (i + 1 != line.size())
            ss << "_";
        }
      }
      folderName = ss.str();
    }
    // Clear and get ready for the next line
    line.clear();
  }

  // This will default to the output files being in a directory tree structure
  // as follows
  // ./OutputFolder/temp_xxx ./OutputFolder/temp_yyy
  if (folderName.empty()) {
    folderName = "OutputFolder";
  }
  return folderName;
}

std::string ParallelTemperingPreprocessor::getTemperature(const char *fileName,
                                                          int worldRank) {
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(fileName);
  std::string temperature;

  while (reader.readNextLine(line)) {
    if (line.size() == 0) {
      continue;
    } else if (line[0] == "Temperature") {
      if (line.size() > 2) {
        temperature = line[worldRank + 1];
      } else {
        temperature = line[1];
      }
    }
    // Clear and get ready for the next line
    line.clear();
  }
  return temperature;
}

std::string
ParallelTemperingPreprocessor::getChemicalPotential(const char *fileName,
                                                    int worldRank) {
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(fileName);
  std::stringstream chemPotStream;
  std::string resName;
  std::string val;

  while (reader.readNextLine(line)) {
    if (line.size() == 0) {
      continue;
    } else if (checkString(line[0], "ChemPot")) {
      if (line.size() > 3) {
        resName = line[1];
        val = line[2 + worldRank];
        chemPotStream << "_" << resName << "_" << val;
      } else if (line.size() != 3) {
        std::cout
            << "Error: Chemical potential parameters are not specified!\n";
        MPI_Finalize();
        exit(EXIT_FAILURE);
      } else {
        resName = line[1];
        val = line[2];
        chemPotStream << "_" << resName << "_" << val;
      }
    }
    // Clear and get ready for the next line
    line.clear();
  }
  return chemPotStream.str();
}

bool ParallelTemperingPreprocessor::checkIfParallelTemperingEnabled(
    const char *fileName) {
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(fileName);

  bool parallelTemperingEnabled = false;

  while (reader.readNextLine(line)) {
    if (line.size() == 0)
      continue;

    if (checkString(line[0], "ParallelTemperingFreq")) {
      parallelTemperingEnabled = checkBool(line[1]);
    }
    // Clear and get ready for the next line
    line.clear();
  }
  return parallelTemperingEnabled;
}

bool ParallelTemperingPreprocessor::checkIfRestart(const char *fileName) {
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(fileName);

  bool restart = false;

  while (reader.readNextLine(line)) {
    if (line.size() == 0)
      continue;

    if (checkString(line[0], "Restart")) {
      restart = checkBool(line[1]);
      return restart;
    }

    // Clear and get ready for the next line
    line.clear();
  }
  return restart;
}

bool ParallelTemperingPreprocessor::checkIfRestartFromCheckpoint(
    const char *fileName) {
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(fileName);

  bool restartFromCheckpoint = false;

  while (reader.readNextLine(line)) {
    if (line.size() == 0)
      continue;

    if (checkString(line[0], "RestartCheckpoint")) {
      restartFromCheckpoint = checkBool(line[1]);
      return restartFromCheckpoint;
    }

    // Clear and get ready for the next line
    line.clear();
  }
  return restartFromCheckpoint;
}

void ParallelTemperingPreprocessor::
    setupReplicaDirectoriesAndRedirectSTDOUTToFile(std::string multiSimTitle,
                                                   std::string temperature) {
  std::stringstream replicaTemp;
  replicaTemp << "temp_" << temperature;
  std::string replicaDirectory = replicaTemp.str();
  mkdirWrapper(multiSimTitle, replicaDirectory);
}

void ParallelTemperingPreprocessor::
    setupReplicaDirectoriesAndRedirectSTDOUTToFile(std::string multiSimTitle,
                                                   std::string temperature,
                                                   std::string chemPot) {
  std::stringstream replicaTemp;
  replicaTemp << "temp_" << temperature << chemPot;
  std::string replicaDirectory = replicaTemp.str();
  mkdirWrapper(multiSimTitle, replicaDirectory);
}

void ParallelTemperingPreprocessor::mkdirWrapper(
    std::string multisimDirectoryName, std::string replicaDirectoryName) {
  std::string multiSimInputFolder;
  std::string multiSimOutputFolder;

  multiSimInputFolder = getInputFolderName(multisimDirectoryName.c_str());
  multiSimOutputFolder = getOutputFolderName(multisimDirectoryName.c_str());

  std::stringstream replicaInputStream;
  std::stringstream replicaOutputStream;
  std::stringstream replicaErrorStream;

  replicaInputStream << multiSimInputFolder << OS_SEP << replicaDirectoryName
                     << OS_SEP;

  replicaOutputStream << multiSimOutputFolder << OS_SEP << replicaDirectoryName
                      << OS_SEP;

  replicaErrorStream << multiSimOutputFolder << OS_SEP << replicaDirectoryName
                     << OS_SEP;

  // If InputFolderName not provided, use empty string.
  if (multiSimInputFolder.empty()) {
    replicaInputDirectoryPath = multiSimInputFolder;
  } else {
    // Else use the provided InputFolderName
    replicaInputDirectoryPath = replicaInputStream.str();
  }
  replicaOutputDirectoryPath = replicaOutputStream.str();

  // printf("Creating directory : %s\n", multisimDirectoryName.c_str());
  mkdir(multiSimOutputFolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir(replicaOutputDirectoryPath.c_str(),
        S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  replicaOutputStream << "ConsoleOut.dat";

  replicaErrorStream << "ErrorsMessages.dat";

  std::string pathToReplicaLogFile = replicaOutputStream.str();
  std::string pathToReplicaErrorLogFile = replicaErrorStream.str();
  if (worldRank == 0) {
    std::cout << "Monitor progress of your simulation by navigating to a "
                 "replica output directory and issuing:\n"
              << "\t$ tail -f \"YourUniqueFileName\".console" << std::endl;
  }

  stdOut = freopen(pathToReplicaLogFile.c_str(), "w", stdout);
  stdErr = freopen(pathToReplicaErrorLogFile.c_str(), "w", stderr);
}

bool ParallelTemperingPreprocessor::checkString(std::string str1,
                                                std::string str2) {
  for (int k = 0; k < str1.length(); k++) {
    str1[k] = toupper(str1[k]);
  }

  for (int j = 0; j < str2.length(); j++) {
    str2[j] = toupper(str2[j]);
  }

  return (str1 == str2);
}

bool ParallelTemperingPreprocessor::checkBool(std::string str) {
  int k;
  // capitalize string
  for (k = 0; k < str.length(); k++) {
    str[k] = toupper(str[k]);
  }

  if (str == "ON" || str == "TRUE" || str == "YES")
    return true;
  else if (str == "OFF" || str == "FALSE" || str == "NO")
    return false;
  std::cout << "Error: " << str << "couldn't be recognized!" << std::endl;
  exit(EXIT_FAILURE);
}

MultiSim::MultiSim(ParallelTemperingPreprocessor &pt)
    : worldSize(pt.worldSize), worldRank(pt.worldRank),
      replicaInputDirectoryPath(pt.replicaInputDirectoryPath),
      replicaOutputDirectoryPath(pt.replicaOutputDirectoryPath),
      restart(pt.restart), restartFromCheckpoint(pt.restartFromCheckpoint),
      parallelTemperingEnabled(pt.parallelTemperingEnabled) {
  std::string filename = replicaOutputDirectoryPath + "ParallelTempering.dat";
  fplog = fopen(filename.c_str(), "w");
}

#endif