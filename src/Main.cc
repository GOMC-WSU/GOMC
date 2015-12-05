/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.0 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Simulation.h"
#include "GOMC_Config.h"    //For version number
#include <iostream>
#include <ctime>

//find and include appropriate files for getHostname
#ifdef _WIN32
#include <Winsock2.h>
#define HOSTNAME
#elif defined(__linux__) || defined(__apple__) || defined(__FreeBSD__)
#include <unistd.h>
#define HOSTNAME
#endif

namespace{
    std::ostream& PrintTime(std::ostream& stream);
    std::ostream& PrintHostname(std::ostream& stream);
    std::ostream& PrintVersion(std::ostream& stream);
    void PrintSimulationHeader();
    void PrintSimulationFooter();
    void PrintDebugMode();
    bool CheckAndPrintEnsemble();
	void _PAUSE(const char* _MSG);
}

int main(void)
{
	//const char * nm = "in.dat";//REMOVED TO ALLOW PASSAGE OF ANY NAME
	PrintSimulationHeader();
	//Only run if valid ensemble was detected.
	if (CheckAndPrintEnsemble())
	{   
#ifndef NDEBUG
		PrintDebugMode();
#endif
		//FOLLOWING LINES ADDED TO OBTAIN INPUT PARAMETER FILE
		string inputFileString;
		fstream inputFileReader;
		bool fileFound;

		do
		{
		  try {
			  cout << "Please provide location and name of input parameter file (*.dat): ";
			  cin >> inputFileString;
			  remove(inputFileString.begin(), inputFileString.end(), ' ');//REMOVE ANY TRAILING SPACE
			  if (inputFileString.length() == 0)
				  throw invalid_argument("No filename provided!");

			  inputFileReader.open(inputFileString.c_str(), ios::in | ios::out); //OPEN FILE
			  if (!inputFileReader.is_open()) //CHECK IF FILE IS OPENED...IF NOT OPENED EXCEPTION REASON FIRED
				  throw invalid_argument("Cannot open/find file in the directory provided!");

			  inputFileReader.close(); //CLOSE FILE TO NOW PASS TO SIMULATION
			  fileFound = true; //SET TRUE FILE FOUND
		  }
		  //EXCEPTION HANDLING
		  catch (invalid_argument &ex){
			  cout << "ERROR: " << ex.what() << endl;
			  _PAUSE("Please press ENTER to continue..."); //PAUSE FOR USER TO ACKNOWLEDGE ERROR
			  fileFound = false; //SET FALSE FILE FOUND
		  }
		} while (!fileFound);
		//ONCE FILE FOUND PASS STRING TO SIMULATION CLASS TO READ AND HANDLE PDB|PSF FILE
		Simulation sim(inputFileString.c_str());
		sim.RunSimulation();
		PrintSimulationFooter();
	}
	return 0;
}


namespace {

	void PrintSimulationHeader()
   {
      std::cout << PrintVersion << '\n'
		<< "Started at: " << PrintTime
#ifdef HOSTNAME
		<< "On hostname: " << PrintHostname
#endif
		<< "\n\n";
   }

	bool CheckAndPrintEnsemble()
   {
      bool healthy = true;
      std::cout << "------------------------------------------------------"
		<< std::endl
		<< "This code was compiled to support the ";
#if ENSEMBLE == NVT
      std::cout << "canonical (NVT)";
#elif ENSEMBLE == GEMC
      std::cout << "Gibbs";
#elif ENSEMBLE == GCMC
      std::cout << "grand canonical";
#else
      std::cerr << "CRITICAL ERROR! Preprocessor value ENSEMBLE is "
		<< "invalid or undefined." << std::endl
		<< "Code will exit.";
      healthy = false;
#endif
      std::cout << " ensemble." << std::endl
		<< "------------------------------------------------------"
		<< std::endl << std::endl;
      return healthy;
   }

	void PrintDebugMode()
	{
	  std::cout << "#########################################################\n";
	  std::cout << "################# RUNNING IN DEBUG MODE #################\n";
	  std::cout << "#########################################################\n\n";
	}

	void PrintSimulationFooter()
	{
		std::cout << PrintVersion << '\n'
			<< "Completed at: " << PrintTime
			<< "On hostname: " << PrintHostname
			<< '\n';
	}

	std::ostream& PrintVersion(std::ostream& stream)
	{
		stream << "GOMC Serial Version " << GOMC_VERSION_MAJOR 
			<< '.' << GOMC_VERSION_MINOR;
		return stream;
	}

	std::ostream& PrintTime(std::ostream& stream)
	{
		time_t timer;
		time(&timer);
		stream << asctime(localtime(&timer));
		return stream;
	}

	std::ostream& PrintHostname(std::ostream& stream)
	{
	#ifdef HOSTNAME
	#ifdef _WIN32
		//setup WINSOCK
		WSADATA wsaData;
		WSAStartup(MAKEWORD(2,0), &wsaData);
	#endif

		const int maxNameLength = 80;
		char hostname[maxNameLength];
		gethostname(hostname, maxNameLength);
		//gethostname does not guarantee null termination
		hostname[maxNameLength - 1] = '\0';
		stream << hostname;

	#ifdef _WIN32
		//teardown WINSOCK
		WSACleanup();
	#endif
	#else
		stream << "Hostname Unavailable";
	#endif
		return stream;
	}

	void _PAUSE(const char* _MSG)
	{
	#ifdef max	//	windows.h library defines a function-like macro with the name of max()
	#define _TEMP_MACRO_ max	//	store the predefined macro in a new one
	#undef max	//	undefine the problamatic macro.
	#endif
		cout << _MSG;
		cin.ignore(numeric_limits< streamsize >::max(), '\n');
		if (cin.get() != '\n')
			cin.ignore(numeric_limits< streamsize >::max(), '\n');
		cout << endl;
	#ifdef _Temp_MACRO_
	#define max _TEMP_MACRO_	// restore the max() macro.
	#undef _TEMP_MACRO_	// undefine the temporary macro.
	#endif 
	}
}


