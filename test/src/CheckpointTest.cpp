#include "BondAdjacencyList.h"
#include "ConfigSetup.h"
#include "FFConst.h"
#include "FFSetup.h" //For geometry kinds
#include "InputFileReader.h"
#include "MolSetup.h"
#include "MoveSettings.h"
#include "PDBSetup.h"
#include "Reader.h"
#include "Simulation.h"
#include <gtest/gtest.h>
#include <unistd.h>

TEST(CheckpointTest, Check_PEN_HEX) {

  ulong base_runsteps, Continued_runsteps;
  ulong Continued_true_step;
#if !GOMC_CUDA
  int result = chdir("./test/input/Systems/PEN_HEX/Base/");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  Simulation base("in.conf");
  base.RunSimulation();
  result = chdir("../Continued");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  Simulation Continued("in.conf");
  result = chdir("../SingleRun");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  Simulation SingleRun("in100.conf");
  SingleRun.RunSimulation();

  MoleculeLookup &Continued_ml = Continued.GetMolLookup();
  MoleculeLookup &SingleRun_ml = SingleRun.GetMolLookup();
  MoveSettings &Continued_ms = Continued.GetMoveSettings();
  MoveSettings &SingleRun_ms = SingleRun.GetMoveSettings();
  Coordinates &Continued_coords = Continued.GetCoordinates();
  Coordinates &SingleRun_coords = SingleRun.GetCoordinates();

  EXPECT_EQ(Continued_ml.operator==(SingleRun_ml), true);
  EXPECT_EQ(Continued_ms.operator==(SingleRun_ms), true);
  EXPECT_EQ(Continued_coords.operator==(SingleRun_coords), true);

  result = chdir("../../../../..");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = chdir("./test/input/Systems/PEN_HEX");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = system("exec rm -r ./Base/Base_*");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = system("exec rm -r ./Continued/Continued_*");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = system("exec rm -r ./SingleRun/SingleRun_*");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = chdir("../../../..");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
#endif
}

TEST(CheckpointTest, Check_BPTI_TIP3) {

  ulong base_runsteps, Continued_runsteps;
  ulong Continued_true_step;
#if !GOMC_CUDA

  int result = chdir("./test/input/Systems/BPTI_TIP3/Base/");

  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  Simulation base("in.conf");
  base.RunSimulation();
  result = chdir("../Continued");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  Simulation Continued("in.conf");
  result = chdir("../SingleRun");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }

  Simulation SingleRun("in100.conf");
  SingleRun.RunSimulation();

  MoleculeLookup &Continued_ml = Continued.GetMolLookup();
  MoleculeLookup &SingleRun_ml = SingleRun.GetMolLookup();
  MoveSettings &Continued_ms = Continued.GetMoveSettings();
  MoveSettings &SingleRun_ms = SingleRun.GetMoveSettings();
  Coordinates &Continued_coords = Continued.GetCoordinates();
  Coordinates &SingleRun_coords = SingleRun.GetCoordinates();
  Velocity &Continued_velocs = Continued.GetVelocities();
  Velocity &SingleRun_velocs = SingleRun.GetVelocities();
  BoxDimensions &Continued_boxDim = Continued.GetBoxDim();
  BoxDimensions &SingleRun_boxDim = SingleRun.GetBoxDim();
  PRNG &Continued_PRNG = Continued.GetPRNG();
  PRNG &SingleRun_PRNG = SingleRun.GetPRNG();
  Molecules &Continued_mols = Continued.GetMolecules();
  Molecules &SingleRun_mols = SingleRun.GetMolecules();

  EXPECT_EQ(Continued_ml.operator==(SingleRun_ml), true);
  EXPECT_EQ(Continued_ms.operator==(SingleRun_ms), true);
  EXPECT_EQ(Continued_coords.operator==(SingleRun_coords), true);
  EXPECT_EQ(Continued_velocs.operator==(SingleRun_velocs), true);
  EXPECT_EQ(Continued_boxDim.operator==(SingleRun_boxDim), true);
  EXPECT_EQ(Continued_PRNG.operator==(SingleRun_PRNG), true);
  EXPECT_EQ(Continued_mols.operator==(SingleRun_mols), true);

  result = chdir("../../../../..");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = chdir("./test/input/Systems/BPTI_TIP3");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = system("exec rm -r ./Base/Base_*");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = system("exec rm -r ./Continued/Continued_*");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }

  result = system("exec rm -r ./SingleRun/SingleRun_*");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = chdir("../../../..");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
#endif
}

TEST(CheckpointTest, Check_K_CHANNEL_TIP3) {

  ulong base_runsteps, Continued_runsteps;
  ulong Continued_true_step;
#if !GOMC_CUDA

  int result = chdir("./test/input/Systems/K_CHANNEL_TIP3/Base/");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  Simulation base("in.conf");
  base.RunSimulation();
  result = chdir("../Continued");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  Simulation Continued("in.conf");
  result = chdir("../SingleRun");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }

  Simulation SingleRun("in100.conf");
  SingleRun.RunSimulation();

  MoleculeLookup &Continued_ml = Continued.GetMolLookup();
  MoleculeLookup &SingleRun_ml = SingleRun.GetMolLookup();
  MoveSettings &Continued_ms = Continued.GetMoveSettings();
  MoveSettings &SingleRun_ms = SingleRun.GetMoveSettings();
  Coordinates &Continued_coords = Continued.GetCoordinates();
  Coordinates &SingleRun_coords = SingleRun.GetCoordinates();
  Velocity &Continued_velocs = Continued.GetVelocities();
  Velocity &SingleRun_velocs = SingleRun.GetVelocities();
  BoxDimensions &Continued_boxDim = Continued.GetBoxDim();
  BoxDimensions &SingleRun_boxDim = SingleRun.GetBoxDim();
  PRNG &Continued_PRNG = Continued.GetPRNG();
  PRNG &SingleRun_PRNG = SingleRun.GetPRNG();
  Molecules &Continued_mols = Continued.GetMolecules();
  Molecules &SingleRun_mols = SingleRun.GetMolecules();

  EXPECT_EQ(Continued_ml.operator==(SingleRun_ml), true);
  EXPECT_EQ(Continued_ms.operator==(SingleRun_ms), true);
  EXPECT_EQ(Continued_coords.operator==(SingleRun_coords), true);
  EXPECT_EQ(Continued_velocs.operator==(SingleRun_velocs), true);
  EXPECT_EQ(Continued_boxDim.operator==(SingleRun_boxDim), true);
  EXPECT_EQ(Continued_PRNG.operator==(SingleRun_PRNG), true);
  EXPECT_EQ(Continued_mols.operator==(SingleRun_mols), true);

  result = chdir("../../../../..");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = chdir("./test/input/Systems/K_CHANNEL_TIP3");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = system("exec rm -r ./Base/Base_*");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = system("exec rm -r ./Continued/Continued_*");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }

  result = system("exec rm -r ./SingleRun/SingleRun_*");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = chdir("../../../..");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
#endif
}

TEST(CheckpointTest, Check_BPTI_TIP3_FORCE_SWAP) {

  ulong base_runsteps, Continued_runsteps;
  ulong Continued_true_step;
#if !GOMC_CUDA

  int result = chdir("./test/input/Systems/BPTI_TIP3/Base/");

  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  Simulation base("forceSwap100.conf");
  base.RunSimulation();
  result = chdir("../Continued");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  Simulation Continued("forceSwap100.conf");
  result = chdir("../SingleRun");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }

  Simulation SingleRun("forceSwap100.conf");
  SingleRun.RunSimulation();

  MoleculeLookup &Continued_ml = Continued.GetMolLookup();
  MoleculeLookup &SingleRun_ml = SingleRun.GetMolLookup();
  MoveSettings &Continued_ms = Continued.GetMoveSettings();
  MoveSettings &SingleRun_ms = SingleRun.GetMoveSettings();
  Coordinates &Continued_coords = Continued.GetCoordinates();
  Coordinates &SingleRun_coords = SingleRun.GetCoordinates();
  Velocity &Continued_velocs = Continued.GetVelocities();
  Velocity &SingleRun_velocs = SingleRun.GetVelocities();
  BoxDimensions &Continued_boxDim = Continued.GetBoxDim();
  BoxDimensions &SingleRun_boxDim = SingleRun.GetBoxDim();
  PRNG &Continued_PRNG = Continued.GetPRNG();
  PRNG &SingleRun_PRNG = SingleRun.GetPRNG();
  Molecules &Continued_mols = Continued.GetMolecules();
  Molecules &SingleRun_mols = SingleRun.GetMolecules();

  EXPECT_EQ(Continued_ml.operator==(SingleRun_ml), true);
  EXPECT_EQ(Continued_ms.operator==(SingleRun_ms), true);
  EXPECT_EQ(Continued_coords.operator==(SingleRun_coords), true);
  EXPECT_EQ(Continued_velocs.operator==(SingleRun_velocs), true);
  EXPECT_EQ(Continued_boxDim.operator==(SingleRun_boxDim), true);
  EXPECT_EQ(Continued_PRNG.operator==(SingleRun_PRNG), true);
  EXPECT_EQ(Continued_mols.operator==(SingleRun_mols), true);

  result = chdir("../../../../..");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = chdir("./test/input/Systems/BPTI_TIP3");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = system("exec rm -r ./Base/Base_*");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = system("exec rm -r ./Continued/Continued_*");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }

  result = system("exec rm -r ./SingleRun/SingleRun_*");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = chdir("../../../..");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
#endif
}

TEST(CheckpointTest, Check_K_CHANNEL_TIP3_FORCE_SWAP) {

  ulong base_runsteps, Continued_runsteps;
  ulong Continued_true_step;
#if !GOMC_CUDA

  int result = chdir("./test/input/Systems/K_CHANNEL_TIP3/Base/");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  Simulation base("forceSwap100.conf");
  base.RunSimulation();
  result = chdir("../Continued");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  Simulation Continued("forceSwap100.conf");
  result = chdir("../SingleRun");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }

  Simulation SingleRun("forceSwap100.conf");
  SingleRun.RunSimulation();

  MoleculeLookup &Continued_ml = Continued.GetMolLookup();
  MoleculeLookup &SingleRun_ml = SingleRun.GetMolLookup();
  MoveSettings &Continued_ms = Continued.GetMoveSettings();
  MoveSettings &SingleRun_ms = SingleRun.GetMoveSettings();
  Coordinates &Continued_coords = Continued.GetCoordinates();
  Coordinates &SingleRun_coords = SingleRun.GetCoordinates();
  Velocity &Continued_velocs = Continued.GetVelocities();
  Velocity &SingleRun_velocs = SingleRun.GetVelocities();
  BoxDimensions &Continued_boxDim = Continued.GetBoxDim();
  BoxDimensions &SingleRun_boxDim = SingleRun.GetBoxDim();
  PRNG &Continued_PRNG = Continued.GetPRNG();
  PRNG &SingleRun_PRNG = SingleRun.GetPRNG();
  Molecules &Continued_mols = Continued.GetMolecules();
  Molecules &SingleRun_mols = SingleRun.GetMolecules();

  EXPECT_EQ(Continued_ml.operator==(SingleRun_ml), true);
  EXPECT_EQ(Continued_ms.operator==(SingleRun_ms), true);
  EXPECT_EQ(Continued_coords.operator==(SingleRun_coords), true);
  EXPECT_EQ(Continued_velocs.operator==(SingleRun_velocs), true);
  EXPECT_EQ(Continued_boxDim.operator==(SingleRun_boxDim), true);
  EXPECT_EQ(Continued_PRNG.operator==(SingleRun_PRNG), true);
  EXPECT_EQ(Continued_mols.operator==(SingleRun_mols), true);

  result = chdir("../../../../..");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = chdir("./test/input/Systems/K_CHANNEL_TIP3");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = system("exec rm -r ./Base/Base_*");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = system("exec rm -r ./Continued/Continued_*");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }

  result = system("exec rm -r ./SingleRun/SingleRun_*");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
  result = chdir("../../../..");
  if (result) {
    std::cout << "System call failed!" << std::endl;
    exit(1);
  }
#endif
}
