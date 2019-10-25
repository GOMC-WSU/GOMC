/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include "ParallelTemperingPreprocessor.h"
ParallelTemperingPreprocessor::ParallelTemperingPreprocessor(){}

bool ParallelTemperingPreprocessor::checkIfParallelTempering(std::string inputFileString){
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(inputFileString);
  bool isParallelTemperingInTemperature = false;
  bool isParallelTemperingInChemicalPotential = false;
  bool isParallelTemperingInFreeEnergyCoulomb = false;
  bool isParallelTemperingInFreeEnergyVDW = false;

  while(reader.readNextLine(line)) {
    if(line.size() == 0)
      continue;
    if(CheckString(line[0], "Temperature")) {
      if (line.size() > 2)
        isParallelTemperingInTemperature = true;
    } else if (CheckString(line[0], "ChemPot")) {
      if (line.size() > 3)
        isParallelTemperingInChemicalPotential = true;
    } else if (CheckString(line[0], "LambdaCoulomb")) {
      if (line.size() > 2)
        isParallelTemperingInFreeEnergyCoulomb = true;
    }  else if (CheckString(line[0], "LambdaVDW")) {
      if (line.size() > 2)
        isParallelTemperingInFreeEnergyVDW = true;
    }
    // Clear and get ready for the next line
    line.clear();
  }
  return isParallelTemperingInTemperature || isParallelTemperingInChemicalPotential || 
          isParallelTemperingInFreeEnergyCoulomb || isParallelTemperingInFreeEnergyVDW;
}

void ParallelTemperingPreprocessor::checkIfValid(std::string inputFileString){
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(inputFileString);
  int numberOfTemperatures = 0;
  vector < int > numberOfChemPots;
  int numberOfLambdaCoulombs = 0;
  int numberOfLambdaVDWs = 0;

  while(reader.readNextLine(line)) {
    if(line.size() == 0)
      continue;
    if(CheckString(line[0], "Temperature")) {
      numberOfTemperatures = line.size() - 1;
    } else if (CheckString(line[0], "ChemPot")) {
      numberOfChemPots.push_back(line.size() - 2);
    } else if (CheckString(line[0], "LambdaCoulomb")) {
      numberOfLambdaCoulombs = line.size() - 1;
    }  else if (CheckString(line[0], "LambdaVDW")) {
      numberOfLambdaVDWs = line.size() - 1;
    }
    // Clear and get ready for the next line
    line.clear();
  }

  for( vector < int >::iterator it = numberOfChemPots.begin(); it != numberOfChemPots.end(); ++it ){
    if (*it > 1 && numberOfTemperatures > 1 && *it != numberOfTemperatures){
      std::cout << "Error: Unequal number of temperatures and chemical potentials in Multicanonical!\n";
      std::cout << "If you only want to only sample mu-space or temperature-space\n";
      std::cout << "provide only one temperature or only one chemical potential.\n";
      std::cout << "Number of temperatures provided: " << numberOfTemperatures << "\n";
      std::cout << "Number of chemical potentials provided: " << *it << "\n";
      exit(EXIT_FAILURE);
    }
  }

  if ( numberOfLambdaCoulombs != numberOfLambdaVDWs){
      std::cout << "Error: Unequal number of LambdaCoulombs and LambdaVDWs in Free Energy calculation!\n";
      std::cout << "Number of temperatures provided: " << numberOfLambdaCoulombs << "\n";
      std::cout << "Number of temperatures provided: " << numberOfLambdaVDWs << "\n";
      exit(EXIT_FAILURE);
  }
}

int ParallelTemperingPreprocessor::getNumberOfReplicas(std::string inputFileString){
  InputFileReader reader;
  std::vector<std::string> line;
  reader.Open(inputFileString);
  int numberOfTemperatures = 0;
  vector < int > numberOfChemPots;
  int numberOfLambdaCoulombs = 0;
  int numberOfLambdaVDWs = 0;

  int numberOfReplicas = 0;

  while(reader.readNextLine(line)) {
    if(line.size() == 0)
      continue;
    if(CheckString(line[0], "Temperature")) {
      numberOfTemperatures = line.size() - 1;
    } else if (CheckString(line[0], "ChemPot")) {
      numberOfChemPots.push_back(line.size() - 2);
    } else if (CheckString(line[0], "LambdaCoulomb")) {
      numberOfLambdaCoulombs = line.size() - 1;
    }  else if (CheckString(line[0], "LambdaVDW")) {
      numberOfLambdaVDWs = line.size() - 1;
    }
    // Clear and get ready for the next line
    line.clear();
  }

  for( vector < int >::iterator it = numberOfChemPots.begin(); it != numberOfChemPots.end(); ++it ){
    numberOfReplicas = std::max(numberOfReplicas, *it);
  }
  numberOfReplicas = std::max(numberOfReplicas, numberOfTemperatures);
  numberOfReplicas = std::max(numberOfReplicas, numberOfLambdaCoulombs);
  numberOfReplicas = std::max(numberOfReplicas, numberOfLambdaVDWs);

  return numberOfReplicas;
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
