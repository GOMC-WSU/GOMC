#include "BondAdjacencyList.h"
#include "ConfigSetup.h"
#include "FFConst.h"
#include "FFSetup.h" //For geometry kinds
#include "InputFileReader.h"
#include "MolSetup.h"
#include "PDBSetup.h"
#include "Reader.h"
#include <gtest/gtest.h>

TEST(PSFParserTest, CheckProtAndWaterTest) {

  config_setup::RestartSettings rs;
  rs.recalcTrajectory = false;
  std::string pdbnames[2];

  pdbnames[0] = "./test/input/PSFParser/BOX_0.pdb";
  pdbnames[1] = "./test/input/PSFParser/BOX_1.pdb";

  PDBSetup pdb;

  pdb.Init(rs, pdbnames);

  std::string psfnames[2];
  bool psfdefined[2];

  MolSetup ms;

  /*
      Let GOMC parse dialanine
  */

  psfnames[0] = "./test/input/PSFParser/BOX_0.psf";
  psfnames[1] = "./test/input/PSFParser/BOX_1.psf";

  psfdefined[0] = true;
  psfdefined[1] = true;

  ms.Init(psfnames, psfdefined, pdb.atoms);

  /*
      Manually build dialanine
  */

  mol_setup::MolMap kindMap;

  kindMap["PROTA"] = mol_setup::MolKind();
  kindMap["PROTA"].isMultiResidue = true;

  /*
      Build the 3 atoms of SPCE
  */

  std::vector<mol_setup::Atom> SPCE;

  mol_setup::Atom SPCE_1_atom_1 =
      mol_setup::Atom("O1", "SPCE", 1, "0000", "OT", -0.847600, 15.9990);
  mol_setup::Atom SPCE_1_atom_2 =
      mol_setup::Atom("H1", "SPCE", 1, "0000", "HT", 0.423800, 1.0080);
  mol_setup::Atom SPCE_1_atom_3 =
      mol_setup::Atom("H2", "SPCE", 1, "0000", "HT", 0.423800, 1.0080);

  SPCE.push_back(SPCE_1_atom_1);
  SPCE.push_back(SPCE_1_atom_2);
  SPCE.push_back(SPCE_1_atom_3);

  /*
                  HB2
                   |
            HB1---CB---HB3      O
                   |           | |
                   |           | |        HT1
      OT2```       |           | |   HA    |
          ``     `.CA-.`     `. C.`  |     N --HT2
           .- C :  |   `` N.` .`   `.CA.-  |
             | |   HA      |         |    HT3
             | |          HN         |
             | |                     |
             OT1               HB1---CB---HB3
                                     |
                                     HB2
  */

  /*
      Manually Build the 23 atoms of dialanine
  */
  /*                              PSF Col    4       3    2    1       5 6
   * 7  */

  std::vector<mol_setup::Atom> DIALA;

  mol_setup::Atom DIALA_atom_1 =
      mol_setup::Atom("N", "ALA", 2, "0003", "NH3", -0.300000, 14.0070);
  mol_setup::Atom DIALA_atom_2 =
      mol_setup::Atom("HT1", "ALA", 2, "0003", "HC", 0.330000, 1.0080);
  mol_setup::Atom DIALA_atom_3 =
      mol_setup::Atom("HT2", "ALA", 2, "0003", "HC", 0.330000, 1.0080);
  mol_setup::Atom DIALA_atom_4 =
      mol_setup::Atom("HT3", "ALA", 2, "0003", "HC", 0.330000, 1.0080);
  mol_setup::Atom DIALA_atom_5 =
      mol_setup::Atom("CA", "ALA", 2, "0003", "CT1", 0.210000, 12.0110);
  mol_setup::Atom DIALA_atom_6 =
      mol_setup::Atom("HA", "ALA", 2, "0003", "HB1", 0.100000, 1.0080);
  mol_setup::Atom DIALA_atom_7 =
      mol_setup::Atom("CB", "ALA", 2, "0003", "CT3", -0.270000, 12.0110);
  mol_setup::Atom DIALA_atom_8 =
      mol_setup::Atom("HB1", "ALA", 2, "0003", "HA3", 0.090000, 1.0080);
  mol_setup::Atom DIALA_atom_9 =
      mol_setup::Atom("HB2", "ALA", 2, "0003", "HA3", 0.090000, 1.0080);
  mol_setup::Atom DIALA_atom_10 =
      mol_setup::Atom("HB3", "ALA", 2, "0003", "HA3", 0.090000, 1.0080);
  mol_setup::Atom DIALA_atom_11 =
      mol_setup::Atom("C", "ALA", 2, "0003", "C", 0.510000, 12.0110);
  mol_setup::Atom DIALA_atom_12 =
      mol_setup::Atom("O", "ALA", 2, "0003", "O", -0.510000, 15.9990);
  mol_setup::Atom DIALA_atom_13 =
      mol_setup::Atom("N", "ALA", 3, "0003", "NH1", -0.470000, 14.0070);
  mol_setup::Atom DIALA_atom_14 =
      mol_setup::Atom("HN", "ALA", 3, "0003", "H", 0.310000, 1.0080);
  mol_setup::Atom DIALA_atom_15 =
      mol_setup::Atom("CA", "ALA", 3, "0003", "CT1", 0.070000, 12.0110);
  mol_setup::Atom DIALA_atom_16 =
      mol_setup::Atom("HA", "ALA", 3, "0003", "HB1", 0.090000, 1.0080);
  mol_setup::Atom DIALA_atom_17 =
      mol_setup::Atom("CB", "ALA", 3, "0003", "CT3", -0.270000, 12.0110);
  mol_setup::Atom DIALA_atom_18 =
      mol_setup::Atom("HB1", "ALA", 3, "0003", "HA3", 0.090000, 1.0080);
  mol_setup::Atom DIALA_atom_19 =
      mol_setup::Atom("HB2", "ALA", 3, "0003", "HA3", 0.090000, 1.0080);
  mol_setup::Atom DIALA_atom_20 =
      mol_setup::Atom("HB3", "ALA", 3, "0003", "HA3", 0.090000, 1.0080);
  mol_setup::Atom DIALA_atom_21 =
      mol_setup::Atom("C", "ALA", 3, "0003", "CC", 0.340000, 12.0110);
  mol_setup::Atom DIALA_atom_22 =
      mol_setup::Atom("OT1", "ALA", 3, "0003", "OC", -0.670000, 15.9990);
  mol_setup::Atom DIALA_atom_23 =
      mol_setup::Atom("OT2", "ALA", 3, "0003", "OC", -0.670000, 15.9990);

  DIALA.push_back(DIALA_atom_1);
  DIALA.push_back(DIALA_atom_2);
  DIALA.push_back(DIALA_atom_3);
  DIALA.push_back(DIALA_atom_4);
  DIALA.push_back(DIALA_atom_5);
  DIALA.push_back(DIALA_atom_6);
  DIALA.push_back(DIALA_atom_7);
  DIALA.push_back(DIALA_atom_8);
  DIALA.push_back(DIALA_atom_9);
  DIALA.push_back(DIALA_atom_10);
  DIALA.push_back(DIALA_atom_11);
  DIALA.push_back(DIALA_atom_12);
  DIALA.push_back(DIALA_atom_13);
  DIALA.push_back(DIALA_atom_14);
  DIALA.push_back(DIALA_atom_15);
  DIALA.push_back(DIALA_atom_16);
  DIALA.push_back(DIALA_atom_17);
  DIALA.push_back(DIALA_atom_18);
  DIALA.push_back(DIALA_atom_19);
  DIALA.push_back(DIALA_atom_20);
  DIALA.push_back(DIALA_atom_21);
  DIALA.push_back(DIALA_atom_22);
  DIALA.push_back(DIALA_atom_23);

  /* Compare GOMC Parser's PROTA vs our Manual Dialanine */

  typedef std::vector<mol_setup::Atom>::const_iterator atomIterator;
  std::pair<atomIterator, atomIterator> itPair1(
      ms.kindMap["PROTA"].atoms.cbegin(), DIALA.cbegin());
  for (; itPair1.second != DIALA.cend(); ++itPair1.first, ++itPair1.second) {
    EXPECT_EQ(*itPair1.first == *itPair1.second, true);
    if (*itPair1.first == *itPair1.second) {

    } else {
      std::cout << (*itPair1.first).name << " " << (*itPair1.first).mass
                << " not equal to " << (*itPair1.second).name << " "
                << (*itPair1.second).mass << std::endl;
    }
  }

  std::pair<atomIterator, atomIterator> itPair2(
      ms.kindMap["SPCE"].atoms.cbegin(), SPCE.cbegin());
  for (; itPair2.second != SPCE.cend(); ++itPair2.first, ++itPair2.second) {
    EXPECT_EQ(*itPair2.first == *itPair2.second, true);
    if (*itPair2.first == *itPair2.second) {

    } else {
      std::cout << (*itPair2.first).name << " " << (*itPair2.first).mass
                << " not equal to " << (*itPair2.second).name << " "
                << (*itPair2.second).mass << std::endl;
    }
  }

  EXPECT_EQ(ms.kindMap["PROTA"].isMultiResidue, true);
  EXPECT_EQ(ms.kindMap["SPCE"].isMultiResidue, false);
}
