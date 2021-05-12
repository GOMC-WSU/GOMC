#include <gtest/gtest.h>
#include "PDBSetup.h"
#include "MolSetup.h"
#include "BondAdjacencyList.h"
#include "ConfigSetup.h"
#include "FFSetup.h"        //For geometry kinds
#include "FFConst.h"
#include "Reader.h"
#include "InputFileReader.h"
#include "AlphaNum.h"
#include "CheckpointSetup.h"


TEST(ConsistentTrajectoryTest, CheckCheckpointVersusPSF) {

    config_setup::RestartSettings rs;
    std::string pdbnames[2];

    pdbnames[0] = "./test/input/ConsistentTrajectory/BOX_0.pdb";
    pdbnames[1] = "./test/input/ConsistentTrajectory/BOX_1.pdb";

    PDBSetup pdb;
    config_setup::InFiles if2;

    pdb.Init(rs, if2, pdbnames);

    std::string psfnames[2];
    bool psfdefined[2];

    MolSetup ms;

    /*
        Let GOMC parse dialanine
    */

    psfnames[0] = "./test/input/ConsistentTrajectory/BOX_0.psf";
    psfnames[1] = "./test/input/ConsistentTrajectory/BOX_1.psf";

    psfdefined[0] = true;
    psfdefined[1] = true;

    ms.Init(psfnames, psfdefined, pdb.atoms);

    CheckpointSetup csp("./test/input/ConsistentTrajectory/Checkpoint.chk");
    csp.ReadAll();

    
    typedef std::vector<std::string>::const_iterator segmentNameIterator;
    typedef std::vector<uint>::const_iterator checkpointIterator;
    std::pair<segmentNameIterator, checkpointIterator> itPair(ms.molVars.moleculeSegmentNames.cbegin(), csp.originalMoleculeIndicesVec.cbegin());
    AlphaNum seg2Uint;
    for (; itPair.second != csp.originalMoleculeIndicesVec.cend(); ++itPair.first, ++itPair.second){
        EXPECT_EQ(*itPair.first, *itPair.second);
    }
}