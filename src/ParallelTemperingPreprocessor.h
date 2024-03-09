
/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef ParallelTemperingPreprocessor_H
#define ParallelTemperingPreprocessor_H
#include <sys/stat.h>
#include <sys/types.h>

#include <cstdlib>
#include <iostream>
#include <sstream>

#include "GOMC_Config.h" //For version number
#include "InputFileReader.h"
#if GOMC_LIB_MPI
#include <mpi.h>
#endif

#ifdef WIN32
#define OS_SEP '\\'
#else
#define OS_SEP '/'
#endif

class ParallelTemperingPreprocessor {
public:
#if GOMC_LIB_MPI

  explicit ParallelTemperingPreprocessor(int argc, char *argv[]);
  bool checkIfValidRank();
  bool checkIfExpandedEnsemble(const char *fileName);
  void checkIfValid(const char *fileName);
  bool checkIfParallelTemperingEnabled(const char *fileName);
  bool checkIfRestart(const char *fileName);
  bool checkIfRestartFromCheckpoint(const char *fileName);
  int getNumberOfReplicas(const char *fileName);
  std::string getInputFolderName(const char *fileName);
  std::string getOutputFolderName(const char *fileName);
  std::string getTemperature(const char *fileName, int worldRank);
  std::string getChemicalPotential(const char *fileName, int worldRank);
  void setupReplicaDirectoriesAndRedirectSTDOUTToFile(std::string multiSimTitle,
                                                      std::string temperature);
  void setupReplicaDirectoriesAndRedirectSTDOUTToFile(std::string multiSimTitle,
                                                      std::string temperature,
                                                      std::string chemPot);
  void mkdirWrapper(std::string multisimDirectoryName,
                    std::string replicaDirectoryName);
  bool checkString(std::string str1, std::string str2);
  bool checkBool(std::string str);

private:
  friend class MultiSim;
  std::string inputFileStringMPI;
  std::fstream inputFileReaderMPI;
  std::string getMultiSimFolderName;
  int worldSize, worldRank;
  bool restartFromCheckpoint = false;
  bool restart = false;
  bool parallelTemperingEnabled = false;
  std::string replicaInputDirectoryPath;
  std::string replicaOutputDirectoryPath;
  FILE *stdOut;
  FILE *stdErr;

#endif
};

class MultiSim {
public:
  explicit MultiSim(ParallelTemperingPreprocessor &pt);
  const int worldSize, worldRank;
  const std::string replicaInputDirectoryPath;
  const std::string replicaOutputDirectoryPath;
  bool restart, restartFromCheckpoint;
  bool parallelTemperingEnabled;
  FILE *fplog;

private:
};
#endif /*ParallelTemperingPreprocessor_H*/
