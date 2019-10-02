#include "VerletList.h"

VerletList::VerletList(CellList & cellList, Coordinates & coordinates,
                       CalculateEnergy & calcEnergy) : 
cellList(cellList), coordinates(coordinates), calcEnergy(calcEnergy)
{ }

// Let's copy the box info so we would know each atom belongs to which box
void VerletList::Init(const pdb_setup::Atoms& atomData)
{
  this->box = atomData.box;
  neighbors.resize(box.size());
}

void VerletList::UpdateNeighborList()
{
  int lengthOfCoordinates = coordinates.Count();
  assert(box.size() == lengthOfCoordinates);
  assert(box.size() != 0);
  int i, b;

  // Clear neighbors list first
  for(i=0; i<lengthOfCoordinates; i++)
    neighbors[i].clear();

  // For each atom find the neigbors and then add them to the neighbors list
  // If atom 0 has 4 neighbors then the following will be pushed to list:
  // 4 1 2 9 14
  // First number '4' is the number of neighbors and the following 4 numbers
  // are the indeces of those neighbors
  // Atom 1 follows right after atom 0
  // e.g. 4 1 2 9 14 4 2 3 4 0  ...
  //      ========== =========
  //        Atom 0     Atom 1   ...
  for(i=0; i<lengthOfCoordinates; i++) {
    neighbors[i].clear();
    b = box[i];

    CellList::Neighbors n = cellList.EnumerateLocal(coordinates[i], b);

    //store atom index in neighboring cell
    while (!n.Done()) {
      if(*n != i && !calcEnergy.SameMolecule(i, *n) && i < *n)
        neighbors[i].push_back(*n);
      n.Next();
    }
  }
}

std::vector<uint> VerletList::GetNeighborList(uint atom) const
{
  return neighbors[atom];
}