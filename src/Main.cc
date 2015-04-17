/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Simulation.h"
//#include "GOMC_Config.h"    //For version number
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

//namespace{
//    std::ostream& PrintTime(std::ostream& stream);
//    std::ostream& PrintHostname(std::ostream& stream);
//    std::ostream& PrintVersion(std::ostream& stream);
//    void PrintSimulationHeader();
//    void PrintSimulationFooter();
//}

int main(void)
{
   const char * nm = "in.dat";
 //  PrintSimulationHeader();
   Simulation sim(nm);
   sim.RunSimulation();
   //PrintSimulationFooter();
   return 0;
}


//namespace {
//
//void PrintSimulationHeader()
//{
//    std::cout << PrintVersion << '\n'
//        << "Started at: " << PrintTime
//#ifdef HOSTNAME
//        << "On hostname: " << PrintHostname
//#endif
//        << "\n\n";
//}
//
//void PrintSimulationFooter()
//{
//    std::cout << PrintVersion << '\n'
//        << "Completed at: " << PrintTime
//        << "On hostname: " << PrintHostname
//        << '\n';
//}
//
//std::ostream& PrintVersion(std::ostream& stream)
//{
//    stream << "GOMC Serial Version " << GOMC_VERSION_MAJOR 
//        << '.' << GOMC_VERSION_MINOR;
//    return stream;
//}
//
//std::ostream& PrintTime(std::ostream& stream)
//{
//    time_t timer;
//    time(&timer);
//    stream << asctime(localtime(&timer));
//    return stream;
//}
//
//std::ostream& PrintHostname(std::ostream& stream)
//{
//#ifdef HOSTNAME
//#ifdef _WIN32
//    //setup WINSOCK
//    WSADATA wsaData;
//    WSAStartup(MAKEWORD(2,0), &wsaData);
//#endif
//
//    const int maxNameLength = 80;
//    char hostname[maxNameLength];
//    gethostname(hostname, maxNameLength);
//    //gethostname does not guarantee null termination
//    hostname[maxNameLength - 1] = '\0';
//    stream << hostname;
//
//#ifdef _WIN32
//    //teardown WINSOCK
//    WSACleanup();
//#endif
//#else
//    stream << "Hostname Unavailable";
//#endif
//    return stream;
//}
//
//}


