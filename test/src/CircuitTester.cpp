#include "CircuitFinder.h"
#include "FloydWarshallCycle.h"
#include <gtest/gtest.h>

#include "BondAdjacencyList.h"
#include "ConfigSetup.h"
#include "FFConst.h"
#include "FFSetup.h" //For geometry kinds
#include "MolSetup.h"
#include "PDBSetup.h"
#include "Reader.h"

TEST(CircuitTester, DialaTest) {

  config_setup::RestartSettings rs2;

  std::string pdbnames2[2];

  pdbnames2[0] = "./test/input/CircuitFinder/START_BOX_0.pdb";
  pdbnames2[1] = "./test/input/CircuitFinder/START_BOX_1.pdb";

  PDBSetup pdb2;

  pdb2.Init(rs2, pdbnames2);

  std::string psfnames2[2];
  bool psfdefined2[2];

  MolSetup ms2;

  /*
      Let GOMC parse dialanine
  */

  psfnames2[0] = "./test/input/CircuitFinder/START_BOX_0.psf";
  psfnames2[1] = "./test/input/CircuitFinder/START_BOX_1.psf";

  psfdefined2[0] = true;
  psfdefined2[1] = true;

  ms2.Init(psfnames2, psfdefined2, pdb2.atoms);

  CircuitFinder cf(ms2.kindMap["PROTA"].atoms.size());
  FloydWarshallCycle fw(ms2.kindMap["PROTA"].atoms.size());

  typedef std::vector<mol_setup::Bond>::const_iterator bondIterator;
  for (bondIterator it = ms2.kindMap["PROTA"].bonds.cbegin();
       it != ms2.kindMap["PROTA"].bonds.cend(); ++it) {
    cf.addEdge(it->a0, it->a1);
    cf.addEdge(it->a1, it->a0);
    fw.AddEdge(it->a0, it->a1);
    // FW is undirected
    // fw.AddEdge(it->a1, it->a0);
  }

  std::vector<std::vector<int>> cyclicCFAtoms;
  std::vector<std::vector<int>> cyclicFWAtoms;

  cyclicCFAtoms = cf.GetAllCommonCycles();
  cyclicFWAtoms = fw.GetAllCommonCycles();

  /* FW leaves empty vectors when unioning */
  for (int i = 0; i < (int)cyclicFWAtoms.size();) {
    if (cyclicFWAtoms[i].size() == 0) {
      cyclicFWAtoms.erase(cyclicFWAtoms.begin() + i);
    } else
      ++i;
  }

  /* Same number of cycles */
  EXPECT_EQ(cyclicCFAtoms.size() == cyclicFWAtoms.size(), true);

  /* The two algorithms don't guaruntee ordering of cyclic atoms */
  /* First we sort within cycles */
  for (int i = 0; i < (int)cyclicFWAtoms.size(); i++) {
    std::sort(cyclicCFAtoms[i].begin(), cyclicCFAtoms[i].end());
    std::sort(cyclicFWAtoms[i].begin(), cyclicFWAtoms[i].end());
  }
  /* Second we sort the cycles */
  std::sort(cyclicCFAtoms.begin(), cyclicCFAtoms.end(),
            std::greater<std::vector<int>>());
  std::sort(cyclicFWAtoms.begin(), cyclicFWAtoms.end(),
            std::greater<std::vector<int>>());

  /* Cycles are the same size */
  for (int i = 0; i < (int)cyclicCFAtoms.size(); i++) {
    EXPECT_EQ(cyclicCFAtoms[i].size() == cyclicFWAtoms[i].size(), true);
  }

  typedef std::vector<std::vector<int>>::const_iterator cycleIterator;
  std::pair<cycleIterator, cycleIterator> itPair1(cyclicCFAtoms.cbegin(),
                                                  cyclicFWAtoms.cbegin());
  for (; itPair1.first != cyclicCFAtoms.cend();
       ++itPair1.first, ++itPair1.second) {
    typedef std::vector<int>::const_iterator atomIterator;
    std::pair<atomIterator, atomIterator> itPair2((*itPair1.first).cbegin(),
                                                  (*itPair1.second).cbegin());
    for (; itPair2.second != (*itPair1.second).cend();
         ++itPair2.first, ++itPair2.second) {

      /* Cycles contain same atoms */
      EXPECT_EQ(*itPair2.first == *itPair2.second, true);
      if (*itPair2.first == *itPair2.second) {

      } else {
        std::cout << *itPair2.first << " not equal to " << *itPair2.second
                  << std::endl;
      }
    }
  }
}
