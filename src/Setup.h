/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef SETUP_H
#define SETUP_H

#include <cereal/archives/binary.hpp>
#include <cstdlib>
#include <string> //for filename
#include <chrono>

#include "ConfigSetup.h"
#include "FFSetup.h"
#include "GOMC_Config.h" //For PT
#include "MolSetup.h"
#include "PDBSetup.h"
#include "PRNGSetup.h"
#include "ParallelTemperingPreprocessor.h"
class Setup {
public:
  // Read order follows each item
  ConfigSetup config; // 1
  PDBSetup pdb;       // 2
  FFSetup ff;         // 4
  PRNGSetup prng;     // 5
#if GOMC_LIB_MPI
  PRNGSetup prngParallelTemp; // 4
#endif
  MolSetup mol; // 5

  void Init(char const *const configFileName, MultiSim const *const &multisim) {
    // Read in all config data
    auto start = std::chrono::steady_clock::now();
    config.Init(configFileName, multisim);
    auto end = std::chrono::steady_clock::now();
    std::cout << "Config reading time: "
              << std::chrono::duration<double, std::milli>(end - start).count()
              << "ms" << std::endl;

    start = std::chrono::steady_clock::now();
    ff.Init(config.in.files.param, config.in.ffKind.isCHARMM); //read in ff data
    end = std::chrono::steady_clock::now();
    std::cout << "FF reading time: "
              << std::chrono::duration<double, std::milli>(end - start).count()
              << "ms" << std::endl;

    start = std::chrono::steady_clock::now();
    pdb.Init(config.in.restart, config.in.files.pdb.name); //read in PDB data
    end = std::chrono::steady_clock::now();
    std::cout << "PDB reading time: "
              << std::chrono::duration<double, std::milli>(end - start).count()
              << "ms" << std::endl;

    prng.Init(config.in.restart, config.in.prng, config.in.files.seed.name); //read in PRNG data
#if GOMC_LIB_MPI
    if (multisim->parallelTemperingEnabled)
      prngParallelTemp.Init(config.in.restart, config.in.prngParallelTempering,
                            config.in.files.seed.name);
#endif
    // Load the MolSetup from checkpoint
    if (!config.in.restart.restartFromCheckpoint) {
      // Read molecule data from psf
      auto molStart = std::chrono::steady_clock::now();
      if (mol.Init(config.in.files.psf.name, config.in.files.psf.defined,
                  pdb.atoms) != 0) {
        exit(EXIT_FAILURE);
      }
      auto molEnd = std::chrono::steady_clock::now();
      std::cout << "MolSetup reading time: "
                << std::chrono::duration<double, std::milli>(molEnd - molStart).count()
                << "ms" << std::endl;
      mol.AssignKinds(mol.molVars, ff);
    }
  }
};

#endif /*SETUP_H*/
