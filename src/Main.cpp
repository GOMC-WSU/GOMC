/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Simulation.h"
#include "GOMC_Config.h"    //For version number
#ifdef GOMC_CUDA
#include "cuda.h"
#include <cuda_runtime_api.h>
#endif
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

namespace
{
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
#ifdef GOMC_CUDA
void PrintGPUHardwareInfo();
#endif

int main(int argc, char *argv[])
{
#ifndef NDEBUG
  PrintDebugMode();
#endif
  PrintSimulationHeader();
  PrintHardwareInfo();
#ifdef GOMC_CUDA
  PrintGPUHardwareInfo();
#endif
  //Only run if valid ensemble was detected.
  if (CheckAndPrintEnsemble()) {
    //FOLLOWING LINES ADDED TO OBTAIN INPUT PARAMETER FILE
    string inputFileString;
    fstream inputFileReader;
    uint numThreads;

    //CHECK IF ARGS/FILE PROVIDED IN CMD LINE
    if (argc < 2) {
      std::cout << "Error: Input parameter file (*.dat or *.conf) not specified on command line!\n";
      exit(EXIT_FAILURE);
    } else {
      if(argc == 2) {
        //FIRST PARAMETER WILL BE FILE NAME
        inputFileString = argv[1];
        numThreads = 1;
      } else {
        //SECOND PARAMETER WILL BE FILE NAME
        inputFileString = argv[2];

        if(argv[1][0] == '+' && argv[1][1] == 'p') {
          numThreads = ReadNum(argv[1]);
        } else {
          std::cout << "Error: Undefined command to set number of threads!\n";
          std::cout << "Use +p# command to set number of threads.\n";
          exit(EXIT_FAILURE);
        }

      }
    }

    //SET NUMBER OF THREADS
#ifdef _OPENMP
    omp_set_num_threads(numThreads);
    printf("%-40s %-d \n", "Info: Number of threads", numThreads);
#else
    printf("%-40s %-d \n", "Info: Number of threads", 1);
#endif

    //OPEN FILE
    inputFileReader.open(inputFileString.c_str(), ios::in | ios::out);

    //CHECK IF FILE IS OPENED...IF NOT OPENED EXCEPTION REASON FIRED
    if (!inputFileReader.is_open()) {
      std::cout << "Error: Cannot open/find " << inputFileString <<
                " in the directory provided!\n";
      exit(EXIT_FAILURE);
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


namespace
{

void PrintSimulationHeader()
{
  std::cout << PrintVersion << '\n'
            << "Info: Start Time: " << PrintTime
#ifdef HOSTNAME
            << "Info: Host Name: " << PrintHostname
#endif
            << "\n";
}

bool CheckAndPrintEnsemble()
{
  bool healthy = true;
  std::cout << "Info: GOMC COMPILED TO RUN ";
#if ENSEMBLE == NVT
  std::cout << "CANONICAL (NVT)";
#elif ENSEMBLE == GEMC
  std::cout << "GIBBS";
#elif ENSEMBLE == GCMC
  std::cout << "GRAND CANONICAL";
#elif ENSEMBLE == NPT
  std::cout << "ISOBARIC-ISOTHERMAL";
#else
  std::cerr << "CRITICAL ERROR! Preprocessor value ENSEMBLE is "
            << "invalid or undefined." << std::endl
            << "Code will exit.";
  healthy = false;
#endif
  std::cout << " ENSEMBLE." << std::endl;
  return healthy;
}

void PrintDebugMode()
{
  std::cout << "################################################################################\n";
  std::cout << "########################## RUNNING GOMC IN DEBUG MODE ##########################\n";
  std::cout << "################################################################################\n";
}

void PrintSimulationFooter()
{
  std::cout << PrintVersion << '\n'
            << "Info: Completed at: " << PrintTime
            << "Info: On hostname: " << PrintHostname
            << '\n';
}

std::ostream& PrintVersion(std::ostream& stream)
{
  stream << "Info: GOMC Version " << GOMC_VERSION_MAJOR
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
  WSAStartup(MAKEWORD(2, 0), &wsaData);
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
  stream << "Info: Hostname Unavailable";
#endif
  return stream;
}

uint ReadNum(char *argv)
{
  uint thread = 0;

  for(uint i = 2; argv[i] != 0; i++) {
    thread = thread * 10 + (argv[i] - '0');
  }

  return thread;
}
}

void PrintHardwareInfo()
{
#ifdef __linux__
  struct sysinfo mem;
  const real megabyte = 1024 * 1024;
  struct utsname name;
  uname(&name);
  printf("CPU information:\n");
  std::cout << std::setprecision(1) << std::fixed;
  std::cout << "Info: Total number of CPUs: " << get_nprocs() << std::endl;
  std::cout << "Info: Total number of CPUs available: " << sysconf(_SC_NPROCESSORS_ONLN) << std::endl;
  std::cout << "Info: Model name:" << std::flush;
  if(system("awk -F: '/model name/ {print $2;exit}' /proc/cpuinfo") == -1)
    std::cout << "Error: Couldn't retrieve CPU information" << std::endl;
  std::cout << "Info: System name: " << name.sysname << std::endl;
  std::cout << "Info: Release: " << name.release << std::endl;
  std::cout << "Info: Version: " << name.version << std::endl;
  std::cout << "Info: Kernel Architecture: " << name.machine << std::endl;
  if(sysinfo(&mem) == 0) {
    std::cout << "Info: Total Ram: " << mem.totalram / megabyte << "MB" << std::endl;
    std::cout << "Info: Used Ram: " << mem.totalram / megabyte - mem.freeram / megabyte << "MB" << std::endl;
  }
  std::cout << "Info: Working in the current directory: " << get_current_dir_name();
  std::cout << std::endl;
#endif
}

#ifdef GOMC_CUDA
void PrintGPUHardwareInfo()
{
  int nDevices;
  int fast = 0;
  int fastIndex = 0;

  cudaGetDeviceCount(&nDevices);
  if(nDevices <= 4) {
    printf("GPU information:\n");
    for (int i = 0; i < nDevices; i++) {
      cudaDeviceProp prop;
      cudaGetDeviceProperties(&prop, i);
      printf("Info: Device Number: %d\n", i);
      printf("Info: Device name: %s\n", prop.name);
      printf("Info: Memory Clock Rate (KHz): %d\n", prop.memoryClockRate);
      printf("Info: Memory Bus Width (bits): %d\n", prop.memoryBusWidth);
      printf("Info: Peak Memory Bandwidth (GB/s): %f\n",
             2.0 * prop.memoryClockRate * (prop.memoryBusWidth / 8) / 1.0e6);
      if( prop.memoryClockRate > fast) {
        fast =  prop.memoryClockRate;
        fastIndex = i;
      }
    }
  }
  cudaSetDevice(fastIndex);
  printf("Info: Using Device Number: %d\n", fastIndex);
}
#endif
