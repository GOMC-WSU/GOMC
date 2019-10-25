
/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef ParallelTemperingPreprocessor_H
#define ParallelTemperingPreprocessor_H

#include "InputFileReader.h"
class ParallelTemperingPreprocessor {
public:
explicit ParallelTemperingPreprocessor();
bool checkIfParallelTempering(std::string inputFileString);
void checkIfValid(std::string inputFileString);
int getNumberOfReplicas(std::string inputFileString);
bool CheckString(string str1, string str2);
};
#endif /*ParallelTemperingPreprocessor_H*/