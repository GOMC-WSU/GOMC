
/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.51
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef ParallelTemperingPreprocessor_H
#define ParallelTemperingPreprocessor_H
#include "InputFileReader.h"
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <sstream>
#include <cstdlib>

#include "GOMC_Config.h"    //For version number
#if GOMC_LIB_MPI
#include <mpi.h>
#include <experimental/filesystem>
//#include "filesystem.hpp"
#endif

#ifdef _WIN32
#define OS_SEP '\\'
#else
#define OS_SEP '/'
#endif

class ParallelTemperingPreprocessor
{
public:

#if GOMC_LIB_MPI

  explicit ParallelTemperingPreprocessor( int argc,
                                          char *argv[]);
  bool checkIfValidRank();
  bool checkIfParallelTempering(const char *fileName);
  void checkIfValid(const char *fileName);
  bool checkIfRestart(const char *fileName);
  bool checkIfRestartFromCheckpoint(const char *fileName);
  int getNumberOfReplicas(const char *fileName);
  std::string getMultiSimFolderName(const char *fileName);
  std::string getTemperature(const char *fileName, int worldRank);
  std::string getChemicalPotential(const char *fileName, int worldRank);
  std::string setupReplicaDirectoriesAndRedirectSTDOUTToFile(std::string multiSimTitle, std::string temperature);
  std::string setupReplicaDirectoriesAndRedirectSTDOUTToFile(std::string multiSimTitle, std::string temperature, std::string chemPot);
  std::string mkdirWrapper(std::string multisimDirectoryName, std::string replicaDirectoryName);
  bool checkString(string str1, string str2);
  bool checkBool(string str);

private:
  friend class MultiSim;
  std::string inputFileStringMPI;
  fstream inputFileReaderMPI;
  std::string pathToReplicaDirectory;
  int worldSize, worldRank;
  bool restartFromCheckpoint;
  bool restart;

#endif

};

class MultiSim
{
public:
  explicit MultiSim(ParallelTemperingPreprocessor & pt);
  const int worldSize, worldRank;
  const std::string pathToReplicaDirectory;
  bool restart, restartFromCheckpoint;
private:
};
#endif /*ParallelTemperingPreprocessor_H*/
