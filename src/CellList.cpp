/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.51
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "CellList.h"
#include "BoxDimensions.h"
#include "BoxDimensionsNonOrth.h"
#include "Molecules.h"
#include "XYZArray.h"
#include "MoleculeLookup.h"

#include <algorithm>

const int CellList::END_CELL;

CellList::CellList(const Molecules& mols,  BoxDimensions& dims)
  : mols(&mols)
{
  dimensions = &dims;
  isBuilt = false;
  for(uint b = 0; b < BOX_TOTAL; b++) {
    edgeCells[b][0] = edgeCells[b][1] = edgeCells[b][2] = 0;
  }
}

void CellList::SetCutoff()
{
  for(uint b = 0; b < BOX_TOTAL; b++) {
    cutoff[b] = dimensions->rCut[b];
  }
}

bool CellList::IsExhaustive() const
{
  std::vector<int> particles(list);
  for(int b = 0; b < BOX_TOTAL; ++b) {
    particles.insert(particles.end(), head[b].begin(), head[b].end());
  }
  particles.erase(std::remove(particles.begin(), particles.end(), -1), particles.end());
  std::sort(particles.begin(), particles.end());
  for(int i = 0; i < particles.size(); ++i) {
    if (i != particles[i]) return false;
  }
  return true;
}

void CellList::RemoveMol(const int molIndex, const int box, const XYZArray& pos)
{
  // For each atom in molecule
  int p = mols->MolStart(molIndex);
  int end = mols->MolEnd(molIndex);
  while(p != end) {
    int cell = PositionToCell(pos[p], box);
    int at = head[box][cell];

    //If particle we're looking for is at head of list assign its
    //pointed at index (should be -1) to head.... this is the case
    //for removing the head molecule in the cell, which is often when
    //we're removing the last molecule/particle from a particular cell.
    //
    //If particle isn't at the head of the list, traverse links to find it,
    //relinking once found.
    if (at == p) {
      head[box][cell] = list[p];
    } else {
      while(at != END_CELL) {
        if (list[at] == p) {
          list[at] = list[p];
          break;
        }
        at = list[at];
      }
    }
    ++p;
  }
}

void CellList::AddMol(const int molIndex, const int box, const XYZArray& pos)
{
  // For each atom in molecule
  int p = mols->MolStart(molIndex);
  int end = mols->MolEnd(molIndex);

  //Note: GridAll assigns everthing to END_CELL
  // so list should point to that
  // if this is the first particle in a particular cell.
  while(p != end) {
    int cell = PositionToCell(pos[p], box);
    //Make the current head index the index the new head points at.
    list[p] = head[box][cell];
    //Assign the new head as our particle index
    head[box][cell] = p;
    ++p;
  }
}

// Resize all boxes to match current axes
void CellList::ResizeGrid(const BoxDimensions& dims)
{
  for(uint b = 0; b < BOX_TOTAL; ++b) {
    XYZ sides = dims.axis[b];
    bool rebuild = false;
    int* eCells = edgeCells[b];
    int oldCells = eCells[0];
    eCells[0] = std::max((int)floor(sides.x / cutoff[b]), 3);
    cellSize[b].x = sides.x / eCells[0];
    rebuild |= (!isBuilt || (oldCells != eCells[0]));

    oldCells = eCells[1];
    eCells[1] = std::max((int)floor(sides.y / cutoff[b]), 3);
    cellSize[b].y = sides.y / eCells[1];
    rebuild |= (!isBuilt || (oldCells != eCells[1]));

    oldCells = eCells[2];
    eCells[2] = std::max((int)floor(sides.z / cutoff[b]), 3);
    cellSize[b].z = sides.z / eCells[2];
    rebuild |= (!isBuilt || (oldCells != eCells[2]));

    if (rebuild) {
      RebuildNeighbors(b);
    }
  }
  isBuilt = true;
}

// Resize one boxes to match current axes
void CellList::ResizeGridBox(const BoxDimensions& dims, const uint b)
{
  XYZ sides = dims.axis[b];
  bool rebuild = false;
  int* eCells = edgeCells[b];
  int oldCells = eCells[0];
  eCells[0] = std::max((int)floor(sides.x / cutoff[b]), 3);
  cellSize[b].x = sides.x / eCells[0];
  rebuild |= (!isBuilt || (oldCells != eCells[0]));

  oldCells = eCells[1];
  eCells[1] = std::max((int)floor(sides.y / cutoff[b]), 3);
  cellSize[b].y = sides.y / eCells[1];
  rebuild |= (!isBuilt || (oldCells != eCells[1]));

  oldCells = eCells[2];
  eCells[2] = std::max((int)floor(sides.z / cutoff[b]), 3);
  cellSize[b].z = sides.z / eCells[2];
  rebuild |= (!isBuilt || (oldCells != eCells[2]));

  if (rebuild) {
    RebuildNeighbors(b);
  }
  isBuilt = true;
}

void CellList::RebuildNeighbors(int b)
{
  int* eCells = edgeCells[b];
  int nCells = eCells[0] * eCells[1] * eCells[2];
  head[b].resize(nCells);
  neighbors[b].resize(nCells);
  for (int i = 0; i < nCells; ++i) {
    neighbors[b][i].clear();
  }

  for (int x = 0; x < eCells[0]; ++x) {
    for (int y = 0; y < eCells[1]; ++y) {
      for (int z = 0; z < eCells[2]; ++z) {
        int cell = x * eCells[2] * eCells[1] + y * eCells[2] + z;
        for (int dx = -1; dx <= 1; ++dx) {
          for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
              // Cache adjacent cells, wrapping if needed
              neighbors[b][cell].push_back(
                ((x + dx + eCells[0]) % eCells[0]) *
                eCells[2] * eCells[1] +
                ((y + dy + eCells[1]) % eCells[1]) *
                eCells[2] +
                ((z + dz + eCells[2]) % eCells[2]));
            }
          }
        }
      }
    }
  }
}

void CellList::GridAll(BoxDimensions& dims, const XYZArray& pos,
                       const MoleculeLookup& lookup)
{
  dimensions = &dims;
  list.resize(pos.Count());
  ResizeGrid(dims);
  for (int b = 0; b < BOX_TOTAL; ++b) {
    head[b].assign(edgeCells[b][0] * edgeCells[b][1] *
                   edgeCells[b][2],
                   END_CELL);
    MoleculeLookup::box_iterator it = lookup.BoxBegin(b),
                                 end = lookup.BoxEnd(b);

    // For each molecule per box
    while (it != end) {
      AddMol(*it, b, pos);
      ++it;
    }
  }
}

void CellList::GridBox(BoxDimensions& dims, const XYZArray& pos,
                       const MoleculeLookup& lookup, const uint b)
{
  dimensions = &dims;
  list.resize(pos.Count());
  ResizeGridBox(dims, b);
  head[b].assign(edgeCells[b][0] * edgeCells[b][1] *
                 edgeCells[b][2], END_CELL);
  MoleculeLookup::box_iterator it = lookup.BoxBegin(b),
                               end = lookup.BoxEnd(b);

  // For each molecule per box
  while (it != end) {
    AddMol(*it, b, pos);
    ++it;
  }
}


CellList::Pairs CellList::EnumeratePairs(int box) const
{
  return CellList::Pairs(*this, box);
}

void CellList::GetCellListNeighbor(uint box, int coordinateSize,
    vector<int> &cellVector, vector<int> &cellStartIndex,
    vector<int> &mapParticleToCell) const
{
  cellVector.resize(coordinateSize);
  cellStartIndex.resize(head[box].size());
  mapParticleToCell.resize(coordinateSize);
  int vector_index = 0;
  for(int cell = 0; cell < head[box].size(); cell++) {
    cellStartIndex[cell] = vector_index;
    int particleIndex = head[box][cell];
    while(particleIndex != END_CELL) {
      cellVector[vector_index] = particleIndex;
      mapParticleToCell[particleIndex] = cell;
      vector_index++;
      particleIndex = list[particleIndex];
    }
    // we are going to sort particles in each cell for better memor access
    sort(cellVector.begin()+cellStartIndex[cell], cellVector.begin()+vector_index);
  }
}

std::vector<std::vector<int> > CellList::GetNeighborList(uint box) const
{
  return neighbors[box];
}

