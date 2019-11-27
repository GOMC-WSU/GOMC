
/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef ParallelTemperingPreprocessor_H
#define ParallelTemperingPreprocessor_H
#include "InputFileReader.h"
#include <dirent.h>
#include <iostream>
#include <sys/types.h> 
#include <sys/stat.h>
#include <sstream>
#include <cstdlib>
#include <mpi.h>


#ifdef WIN32
#define OS_SEP '\\'
#else
#define OS_SEP '/'
#endif

class ParallelTemperingPreprocessor {
public:

explicit ParallelTemperingPreprocessor( int argc, 
                                        char *argv[]);
bool checkIfValidRank();
bool checkIfParallelTempering(const char *fileName);
void checkIfValid(const char *fileName);
int getNumberOfReplicas(const char *fileName);
std::string getMultiSimFolderName(const char *fileName);
std::string getTemperature(const char *fileName, int worldRank);
std::string getChemicalPotential(const char *fileName, int worldRank);
std::string setupReplicaDirectoriesAndRedirectSTDOUTToFile(std::string multiSimTitle, std::string temperature);   
std::string setupReplicaDirectoriesAndRedirectSTDOUTToFile(std::string multiSimTitle, std::string temperature, std::string chemPot);   
std::string mkdirWrapper(std::string multisimDirectoryName, std::string replicaDirectoryName);
bool CheckString(string str1, string str2);
private:
  friend class MultiSim;
  std::string inputFileStringMPI;
  fstream inputFileReaderMPI;
  std::string pathToReplicaDirectory;
  int worldSize, worldRank;
};

class MultiSim {
public:
explicit MultiSim(ParallelTemperingPreprocessor & pt);
const int worldSize, worldRank; 
const std::string pathToReplicaDirectory;
private:
};
#endif /*ParallelTemperingPreprocessor_H*/