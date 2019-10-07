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

  // For each atom find the neigbors and then add them to the neighbors list
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

std::vector<uint> & VerletList::GetNeighborList(uint atom)
{
  return neighbors[atom];
}