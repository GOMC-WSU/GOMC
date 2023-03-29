/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef SETUP_H
#define SETUP_H

#include <cereal/archives/binary.hpp>
#include <cstdlib>
#include <string> //for filename

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
    config.Init(configFileName, multisim);
    // Read in FF data.
    ff.Init(config.in.files.param, config.in.ffKind.isCHARMM);
    // Read in PDB data
    pdb.Init(config.in.restart, config.in.files.pdb.name);
    // Initialize PRNG
    prng.Init(config.in.restart, config.in.prng, config.in.files.seed.name);
#if GOMC_LIB_MPI
    if (multisim->parallelTemperingEnabled)
      prngParallelTemp.Init(config.in.restart, config.in.prngParallelTempering,
                            config.in.files.seed.name);
#endif
    // Load the MolSetup from checkpoint
    if (!config.in.restart.restartFromCheckpoint) {
      // Read molecule data from psf
      if (mol.Init(config.in.files.psf.name, config.in.files.psf.defined,
                   pdb.atoms) != 0) {
        exit(EXIT_FAILURE);
      }
      mol.AssignKinds(mol.molVars, ff);
    }
  }
};

#endif
