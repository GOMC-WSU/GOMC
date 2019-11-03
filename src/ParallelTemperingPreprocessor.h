
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

#ifdef WIN32
#define OS_SEP '\\'
#else
#define OS_SEP '/'
#endif

class ParallelTemperingPreprocessor {
public:
explicit ParallelTemperingPreprocessor();
bool checkIfParallelTempering(std::string inputFileString);
void checkIfValid(std::string inputFileString);
int getNumberOfReplicas(std::string inputFileString);
std::string getMultiSimFolderName(std::string inputFileString);
std::string getTemperature(std::string inputFileString, int world_rank);
std::string getChemicalPotential(std::string inputFileString, int world_rank);
std::string setupReplicaDirectoriesAndRedirectSTDOUTToFile(std::string multiSimTitle, std::string temperature);   
std::string setupReplicaDirectoriesAndRedirectSTDOUTToFile(std::string multiSimTitle, std::string temperature, std::string chemPot);   
std::string mkdirWrapper(std::string multisimDirectoryName, string replicaDirectoryName);
bool CheckString(string str1, string str2);
private:
};

class MultiSim {
public:
explicit MultiSim(int worldSize, int worldRank, std::string pathToReplicaDirectory);
int worldSize, worldRank; 
std::string pathToReplicaDirectory;


private:
};
#endif /*ParallelTemperingPreprocessor_H*/