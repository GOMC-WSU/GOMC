/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#ifndef REPLDIRSETUP_H
#define REPLDIRSETUP_H

#include <dirent.h>
#include <iostream>
#include <sys/types.h> 
#include <sys/stat.h>
#include "StaticVals.h"
#include "System.h"


class ReplDirSetup{
public:
    
StaticVals * statV;
System * sys;
string path_to_replica_log_file;
string path_to_replica_directory;
std::stringstream replica_temp;
std::stringstream replica_stream;

ReplDirSetup(StaticVals * staticVals, System * system, int temperature, ReplicaExchangeParameters replExParams){
   statV = staticVals;
   sys = system;
   ReplDirSetup::setupReplicaDirectories(temperature, replExParams);
}

  void setupReplicaDirectories(int temperature, ReplicaExchangeParameters replExParams){   
     
   #if ENSEMBLE == GCMC
      replica_temp << "temp_" << temperature;
      for (uint kind = 0; kind < sys->molLookup.GetNumKind(); kind++){
         replica_temp  << "_" << statV->mol.kinds[kind].name << "_" << statV->mol.kinds[kind].chemPot;
      }   
   #else
    replica_temp << "temp_" << temperature;
   #endif
    
    std::string replica_directory = replica_temp.str();
    mkdirWrapper(replExParams.multiSimTitle, replica_directory);
 }
 
 void mkdirWrapper(std::string multisim_directory_name, string replica_directory_name){

    replica_stream << multisim_directory_name << "/" 
                << replica_directory_name << "/";
    std::string replica_directory_path = replica_stream.str();
     
    printf("Creating directory : %s\n", multisim_directory_name.c_str());
    mkdir(multisim_directory_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(replica_directory_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
    replica_stream << "/";
    path_to_replica_directory = replica_stream.str();
    replica_stream << "repl_log.txt";
    path_to_replica_log_file = replica_stream.str();    
 }

 

};
#endif /* REPLDIRSETUP_H */

