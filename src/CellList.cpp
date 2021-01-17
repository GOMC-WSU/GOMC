/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
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
#include <set>

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

CellList::CellList(const CellList & other) : mols(other.mols)
{
  dimensions = other.dimensions;
  isBuilt = true;
  for(uint b = 0; b < BOX_TOTAL; b++) {
    edgeCells[b][0] = other.edgeCells[b][0];
    edgeCells[b][1] = other.edgeCells[b][1];
    edgeCells[b][2] = other.edgeCells[b][2];
  }

  for(uint b = 0; b < BOX_TOTAL; b++) {
    RebuildNeighbors(b);
  }

  list.resize(other.list.size());

  for (size_t i = 0; i < other.list.size(); i++) {
    list[i] = other.list[i];
  }

  for(uint b = 0; b < BOX_TOTAL; b++) {
    for (size_t i = 0; i < other.neighbors[b].size(); i++) {
      neighbors[b][i] = other.neighbors[b][i];
    }
  }

  for(uint b = 0; b < BOX_TOTAL; b++) {
    for (size_t i = 0; i < other.head[b].size(); i++) {
      head[b][i] = other.head[b][i];
    }
  }
  //neighbors(other.neighbors);
  //head(other.head);
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
  for(int i = 0; i < (int) particles.size(); ++i) {
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

  //Note: GridAll assigns everything to END_CELL
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
                                   std::vector<int> &cellVector, std::vector<int> &cellStartIndex,
                                   std::vector<int> &mapParticleToCell) const
{
  cellVector.resize(coordinateSize);
  cellStartIndex.resize(head[box].size());
  mapParticleToCell.resize(coordinateSize);
  int vector_index = 0;
  for(size_t cell = 0; cell < head[box].size(); cell++) {
    cellStartIndex[cell] = vector_index;
    int particleIndex = head[box][cell];
    while(particleIndex != END_CELL) {
      cellVector[vector_index] = particleIndex;
      mapParticleToCell[particleIndex] = cell;
      vector_index++;
      particleIndex = list[particleIndex];
    }
    // we are going to sort particles in each cell for better memory access
    std::sort(cellVector.begin() + cellStartIndex[cell], cellVector.begin() + vector_index);
  }
  // push one last cellStartIndex for the last cell
  cellStartIndex.push_back(vector_index);

  // in case there are two boxes we need to remove the extra space allocated here
  cellVector.resize(vector_index);
}

std::vector<std::vector<int> > CellList::GetNeighborList(uint box) const
{
  return neighbors[box];
}

bool CellList::CompareCellList(CellList & other, int coordinateSize)
{
  std::vector<int> cellVector, cellStartIndex, mapParticleToCell;
  std::vector<int> otherCellVector, otherCellStartIndex, otherMapParticleToCell;

  for(uint box = 0; box < BOX_TOTAL; box++) {
    cellVector.resize(coordinateSize);
    cellStartIndex.resize(head[box].size());
    mapParticleToCell.resize(coordinateSize);

    otherCellVector.resize(coordinateSize);
    otherCellStartIndex.resize(head[box].size());
    otherMapParticleToCell.resize(coordinateSize);
  }

  for(uint box = 0; box < BOX_TOTAL; box++) {
    int vector_index = 0;
    for(size_t cell = 0; cell < head[box].size(); cell++) {
      cellStartIndex[cell] = vector_index;
      int particleIndex = head[box][cell];
      while(particleIndex != END_CELL) {
        cellVector[vector_index] = particleIndex;
        mapParticleToCell[particleIndex] = cell;
        vector_index++;
        particleIndex = list[particleIndex];
      }
    }
  }

  for(uint box = 0; box < BOX_TOTAL; box++) {
    int vector_index = 0;
    for(size_t cell = 0; cell < other.head[box].size(); cell++) {
      otherCellStartIndex[cell] = vector_index;
      int particleIndex = other.head[box][cell];
      while(particleIndex != END_CELL) {
        otherCellVector[vector_index] = particleIndex;
        otherMapParticleToCell[particleIndex] = cell;
        vector_index++;
        particleIndex = other.list[particleIndex];
      }
    }
  }

  if(list.size() == other.list.size()) {
    for(size_t i = 0; i < list.size(); i++) {
      if (list[i] != other.list[i])
        std::cout << "List objects are different" << std::endl;
    }
  }

  for (size_t i = 0; i < mapParticleToCell.size(); i++) {
    if (mapParticleToCell[i] != otherMapParticleToCell[i])
      return false;
  }

  std::cout << "CellList objects have equal states" << std::endl;

  return true;
}

void CellList::PrintList()
{
  for(size_t i = 0; i < list.size(); i++)
    std::cout << list[i] << std::endl;

  std::cout << "head vector" << std::endl;
  for(int i = 0; i < BOX_TOTAL; i++) {
    for( size_t j = 0; j < head[i].size(); j++) {
      std::cout << head[i][j] << std::endl;
    }
  }
}

void CellList::GetRangeDimension(int start_cell, int end_cell, int maximum, std::vector<std::pair<int, int>> & range) const
{
  if(start_cell < 0 && end_cell >= maximum) { // wrap around from both end, not sure if this even happens
    range.push_back(std::pair<int, int>(0, maximum-1));
  } else if(start_cell < 0) { // wrap around from start
    range.push_back(std::pair<int, int>(0, end_cell));
    range.push_back(std::pair<int, int>(start_cell + maximum, maximum-1));
  } else if(end_cell >= maximum) {
    range.push_back(std::pair<int, int>(start_cell, maximum-1));
    range.push_back(std::pair<int, int>(0, end_cell - maximum));
  } else { // no wrap around
    range.push_back(std::pair<int, int>(start_cell, end_cell));
  }
}

std::vector<int> CellList::CellsInCavity(XYZ const& center, XYZ const& cavDim, uint box) const
{
  // using std::set to avoid duplicates
  std::set<int> cells;
  // final return vector
  std::vector<int> cellsInsideCavity;
  // range of cells for each dimension
  std::vector<std::pair<int, int>> all_x_range, all_y_range, all_z_range;

  XYZ start_corner = center - cavDim / 2.0;
  XYZ end_corner   = center + cavDim / 2.0;

  int start_cell_x = PositionToCellDimension(start_corner, box, 0);
  int end_cell_x   = PositionToCellDimension(end_corner, box, 0);
  GetRangeDimension(start_cell_x, end_cell_x, edgeCells[box][0], all_x_range);

  int start_cell_y = PositionToCellDimension(start_corner, box, 1);
  int end_cell_y   = PositionToCellDimension(end_corner, box, 1);
  GetRangeDimension(start_cell_y, end_cell_y, edgeCells[box][1], all_y_range);

  int start_cell_z = PositionToCellDimension(start_corner, box, 2);
  int end_cell_z   = PositionToCellDimension(end_corner, box, 2);
  GetRangeDimension(start_cell_z, end_cell_z, edgeCells[box][2], all_z_range);

  // first loop over the possible ranges in each dimension
  // in most cases there are only 1 or 2 possible ranges.
  for(auto & x_range : all_x_range) {
    for(auto & y_range : all_y_range) {
      for(auto & z_range : all_z_range) {
        // loop over each individual range
        for(int cell_x_idx = x_range.first; cell_x_idx <= x_range.second; cell_x_idx++) {
          for(int cell_y_idx = y_range.first; cell_y_idx <= y_range.second; cell_y_idx++) {
            for(int cell_z_idx = z_range.first; cell_z_idx <= z_range.second; cell_z_idx++) {
              // once we have each cell's index in each dimension,
              // we can calculate the cell index and add it to ``cells'' set
              int cell_idx = cell_x_idx * edgeCells[box][1] * edgeCells[box][2] +
                             cell_y_idx * edgeCells[box][2] +
                             cell_z_idx;
              cells.insert(cell_idx);
            }
          }
        }
      }
    }
  }

  std::copy(cells.begin(), cells.end(), cellsInsideCavity.begin());

  return cellsInsideCavity;
}

// returns all the particles that exist in the list of cells
std::vector<int> CellList::GetParticlesInCells(std::vector<int> cells, uint box) const
{
  std::vector<int> particlesInsideCells;

  for(auto & cell : cells) {
    int particleIndex = head[box][cell];
    while(particleIndex != END_CELL) {
      particlesInsideCells.push_back(particleIndex);
      particleIndex = list[particleIndex];
    }
  }

  return particlesInsideCells;
}