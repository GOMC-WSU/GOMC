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
#include <sys/sysinfo.h>
#include <sys/utsname.h>
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
  uint ReadNum(char *argv);
}

void PrintHardwareInfo();

int main(int argc, char *argv[])
{
   PrintSimulationHeader();
   PrintHardwareInfo();
   //Only run if valid ensemble was detected.
   if (CheckAndPrintEnsemble())
   {   
#ifndef NDEBUG
      PrintDebugMode();
#endif

      //FOLLOWING LINES ADDED TO OBTAIN INPUT PARAMETER FILE
      string inputFileString;
      fstream inputFileReader;
      uint numThreads;

      //CHECK IF ARGS/FILE PROVIDED IN CMD LINE
      if (argc < 2)
      {
	 std::cout<<"Input parameter file (*.dat or *.conf) not specified on command line!\n";
	 exit(0);
      }
      else
      {
	if(argc ==2)
	{
	   //FIRST PARAMETER WILL BE FILE NAME
	  inputFileString = argv[1];
	  numThreads = 1;
	}
	else
	{
	  //SECOND PARAMETER WILL BE FILE NAME
	  inputFileString = argv[2];
	  
	  if(argv[1][0] == '+' && argv[1][1] == 'p')
	  {
	    numThreads = ReadNum(argv[1]);
	  }
	  else
	  {
	    std::cout<<"Error: Undefined command to set number of threads!\n";
	    std::cout<< "Use +p# command to set number of threads.\n";
	    exit(0);
	  }
	  
	} 
      }

      //SET NUMBER OF THREADS
#ifdef _OPENMP
      omp_set_num_threads(numThreads);
      std::cout << "Number of threads has been set to: " <<
	numThreads << std::endl;
#endif

      //OPEN FILE
      inputFileReader.open(inputFileString.c_str(), ios::in | ios::out);
 
      //CHECK IF FILE IS OPENED...IF NOT OPENED EXCEPTION REASON FIRED
      if (!inputFileReader.is_open()) 
      {
	std::cout<<"Cannot open/find " << inputFileString <<
	  " in the directory provided!\n";
	 exit(0);
      }

      //CLOSE FILE TO NOW PASS TO SIMULATION
      inputFileReader.close(); 

      //ONCE FILE FOUND PASS STRING TO SIMULATION CLASS TO READ AND
      //HANDLE PDB|PSF FILE
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
#elif ENSEMBLE == NPT
      std::cout << "isobaric-isothermal";
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
    std::cout << "#########################################################\n";
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

  uint ReadNum(char *argv)
  {
    uint thread = 0;
    
    for(uint i=2; argv[i] != 0; i++)
      {
	thread = thread * 10 + (argv[i]-'0');
      }

    return thread;
  }
}
  
void PrintHardwareInfo()
{
#ifndef _WIN32
  struct sysinfo mem;
  const double megabyte = 1024 * 1024;
  struct utsname name;
  uname(&name);
  std::cout << std::setprecision(1) << std::fixed;
  std::cout << "Total number of CPUs: " << get_nprocs() << std::endl;
  std::cout << "Total number of CPUs available: " << sysconf(_SC_NPROCESSORS_ONLN) << std::endl;
  std::cout << "Model name:" << std::flush;
  if(!system("awk -F: '/model name/ {print $2;exit}' /proc/cpuinfo"))
    std::cout << "Couldn't retrieve CPU information" << std::endl;
  std::cout << std::endl;
  std::cout << "System name: " << name.sysname << std::endl;
  std::cout << "Release: " << name.release << std::endl;
  std::cout << "Version: " << name.version << std::endl;
  std::cout << "Kernel Architecture: " << name.machine << std::endl;
  if(sysinfo(&mem)==0)
  {
    std::cout << "Total Ram: " << mem.totalram / megabyte << "MB" << std::endl;
    std::cout << "Used Ram: " << mem.totalram / megabyte - mem.freeram / megabyte << "MB" << std::endl;
  }
  std::cout << "Working in the current directory " << get_current_dir_name() << std::endl;
  std::cout << std::endl;
#endif
}

