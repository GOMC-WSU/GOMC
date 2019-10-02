#pragma once

#include "CellList.h"
#include "Coordinates.h"
#include "CalculateEnergy.h"
#include <vector>

class VerletList
{
public:
  VerletList(CellList & cellList, Coordinates & coordinates,
             CalculateEnergy & calcEnergy);
  ~VerletList() {}

  void Init(const pdb_setup::Atoms& atomData);
  void UpdateNeighborList();
  std::vector<uint>& GetNeighborList(uint box);

private:
  const CellList& cellList;
  const Coordinates& coordinates;
  const CalculateEnergy& calcEnergy;
  std::vector< std::vector<uint> > neighbors;
  std::vector<uint> box;
};