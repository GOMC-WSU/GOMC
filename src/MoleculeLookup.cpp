/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "MoleculeLookup.h" //Header spec.
#include "EnsemblePreprocessor.h" //For box total
#include "PRNG.h" //For index selection
#include "PDBSetup.h" //For init.
#include "Molecules.h" //For init.
#include <algorithm>
#include <utility>
#include <iostream>
#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include "CUDAMemoryManager.cuh"
#include "VariablesCUDA.cuh"
#endif

void MoleculeLookup::Init(const Molecules& mols,
                          const pdb_setup::Atoms& atomData,
                          Forcefield &ff,
                          bool restartFromCheckpoint)
{
  this->restartFromCheckpoint = restartFromCheckpoint;

  numKinds = mols.GetKindsCount();

  molLookup = new uint[mols.count];
  molLookupCount = mols.count;
  // beta has same size as total number of atoms
  molIndex = new int[atomData.beta.size()];
  atomIndex = new int[atomData.beta.size()];
  molKind = new int[atomData.beta.size()];
  atomKind = new int[atomData.beta.size()];
  atomCharge = new double[atomData.beta.size()];

  //+1 to store end value
  boxAndKindStart = new uint[numKinds * BOX_TOTAL + 1];
  boxAndKindSwappableCounts = new uint[numKinds * BOX_TOTAL];

  boxAndKindStartCount = numKinds * BOX_TOTAL + 1;

  // vector[box][kind] = list of mol indices for kind in box
  std::vector<std::vector<std::vector<uint> > > indexVector;
  indexVector.resize(BOX_TOTAL);
  fixedMolecule.resize(mols.count);

  for (int i = 0; i < (int) numKinds * BOX_TOTAL; i++){
    boxAndKindSwappableCounts[i] = 0;
  }


  for (uint b = 0; b < BOX_TOTAL; ++b) {
    indexVector[b].resize(numKinds);
  }

  int counter = 0;
  for(uint m = 0; m < mols.count; ++m) {
    uint box = atomData.box[mols.start[m]];
    uint kind = mols.kIndex[m];
    indexVector[box][kind].push_back(m);
    const MoleculeKind& mk = mols.GetKind(m);
    for(uint a = 0; a < mk.NumAtoms(); ++a) {
      molIndex[counter] = int(m);
      atomIndex[counter] = int(a);
      molKind[counter] = int(kind);
      atomKind[counter] = int(mk.AtomKind(a));
      atomCharge[counter] = mk.AtomCharge(a);
      ++counter;
    }



    /* We don't currently support hybrid molecules - part fixed part flexible
      so we get a consensus based on the precendent of betas defined in this method */
    uint pStart = 0, pEnd = 0;
    mols.GetRangeStartStop(pStart, pEnd, m);
    fixedMolecule[m] = GetConsensusMolBeta(pStart, pEnd, atomData.beta, m, box, mols.kinds[mols.kIndex[m]].name);

    //Find the kind that can be swap(beta == 0) or move(beta == 0 or 2)
    if(fixedMolecule[m] == 0) {
      if(std::find(canSwapKind.begin(), canSwapKind.end(), kind) ==
          canSwapKind.end())
        canSwapKind.push_back(kind);

      if(std::find(canMoveKind.begin(), canMoveKind.end(), kind) ==
          canMoveKind.end())
        canMoveKind.push_back(kind);

      boxAndKindSwappableCounts[box * numKinds + kind]++;

    } else if(fixedMolecule[m] == 2) {
      if(std::find(canMoveKind.begin(), canMoveKind.end(), kind) ==
          canMoveKind.end())
        canMoveKind.push_back(kind);

    }
  }


  uint* progress = molLookup;
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    for (uint k = 0; k < numKinds; ++k) {
      boxAndKindStart[b * numKinds + k] = progress - molLookup;
      progress = std::copy(indexVector[b][k].begin(),
                           indexVector[b][k].end(), progress);
    }
  }

  boxAndKindStart[numKinds * BOX_TOTAL] = mols.count;

  /* originalMoleculeIndices have 2 sources
    if a new run, they are depedent on the originalMolInds set below
    if a checkpointed run, they are the originalInds permuted through mol transfers */
  if (!restartFromCheckpoint){
    originalMoleculeIndices = new uint[mols.count];
    permutedMoleculeIndices = new uint[mols.count];
    uint molCounter = 0, b, k, kI, countByKind, molI;
    for (b = 0; b < BOX_TOTAL; ++b) {
      for (k = 0; k < mols.kindsCount; ++k) {
        countByKind = NumKindInBox(k, b);
        for (kI = 0; kI < countByKind; ++kI) {
          molI = GetMolNum(kI, k, b);
          originalMoleculeIndices[molCounter] = molI;
          /* This allows us to use input files that aren't of the form
          box 0 kind 0 
          box 0 kind 1
          etc.
          without an n^3 search operation.
          Since checkpointed outputs are sorted in the above form,
          we use origMolInds only when rest from checkpoint.
          */
          permutedMoleculeIndices[molCounter] = molCounter;
          ++molCounter;
        }
      }
    }
  } else {
        permutedMoleculeIndices = new uint[mols.count];

    for (uint i = 0; i < molLookupCount; ++i)
      permutedMoleculeIndices[i] = i;
  }

// allocate and set gpu variables
#ifdef GOMC_CUDA
  VariablesCUDA *cudaVars = ff.particles->getCUDAVars();
  int numMol = mols.count + 1;
  // allocate memory to store molecule start atom index
  CUMALLOC((void**) &cudaVars->gpu_startAtomIdx, numMol * sizeof(int));
  // copy start atom index
  cudaMemcpy(cudaVars->gpu_startAtomIdx, mols.start, numMol * sizeof(int), cudaMemcpyHostToDevice);
#endif

}

uint MoleculeLookup::NumInBox(const uint box) const
{
  return boxAndKindStart[(box + 1) * numKinds]
         - boxAndKindStart[box * numKinds];
}

void MoleculeLookup::TotalAndDensity
(uint * numByBox, uint * numByKindBox, double * molFractionByKindBox,
 double * densityByKindBox, double const*const volInv) const
{
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    numByBox[b] = 0.0;
    for (uint k = 0; k < numKinds; ++k) {
      uint numMK = NumKindInBox(k, b);
      uint mkIdx = k + numKinds * b;
      numByKindBox[mkIdx] = numMK;
      densityByKindBox[mkIdx] = numMK * volInv[b];
      numByBox[b] += numMK;
    }

    //Calculate mol fractions
    if (numKinds > 1) {
      for (uint k = 0; k < numKinds; ++k) {
        uint mkIdx = k + numKinds * b;
        if (numByBox[b] > 0) {
          molFractionByKindBox[mkIdx] = (double)numByKindBox[mkIdx] /
                                        (double)numByBox[b];
        } else {
          molFractionByKindBox[mkIdx] = 0.0;
        }
      }
    }

  }

}

#ifdef VARIABLE_PARTICLE_NUMBER


bool MoleculeLookup::ShiftMolBox(const uint mol, const uint currentBox,
                                 const uint intoBox, const uint kind)
{
  uint index = std::find(
                 molLookup + boxAndKindStart[currentBox * numKinds + kind],
                 molLookup + boxAndKindStart[currentBox * numKinds + kind + 1], mol)
               - molLookup;
  assert(index != boxAndKindStart[currentBox * numKinds + kind + 1]);
  assert(molLookup[index] == mol);
  Shift(index, currentBox, intoBox, kind);
  return true;
}

void MoleculeLookup::Shift(const uint index, const uint currentBox,
                           const uint intoBox, const uint kind)
{
  uint oldIndex = index;
  uint newIndex;
  uint section = currentBox * numKinds + kind;

  boxAndKindSwappableCounts[section]--;
  boxAndKindSwappableCounts[intoBox * numKinds + kind]++;

  if(currentBox >= intoBox) {
    while (section != intoBox * numKinds + kind) {
      newIndex = boxAndKindStart[section]++;
      std::swap(molLookup[oldIndex], molLookup[newIndex]);
      std::swap(oldIndex, newIndex);
      if (!restartFromCheckpoint)
      std::swap(permutedMoleculeIndices[oldIndex], permutedMoleculeIndices[newIndex]);
      --section;
    }
  } else {
    while (section != intoBox * numKinds + kind) {
      newIndex = --boxAndKindStart[++section];
      std::swap(molLookup[oldIndex], molLookup[newIndex]);
      std::swap(oldIndex, newIndex);
      if (!restartFromCheckpoint)
      std::swap(permutedMoleculeIndices[oldIndex], permutedMoleculeIndices[newIndex]);
    }
  }
}

#endif /*ifdef VARIABLE_PARTICLE_NUMBER*/

uint MoleculeLookup::GetConsensusMolBeta( const uint pStart,
                                          const uint pEnd,
                                          const std::vector<double> & betas,
                                          const uint m,
                                          const uint box,
                                          const std::string & name){
  double firstBeta = 0.0;
  double consensusBeta = 0.0;
  for (uint p = pStart; p < pEnd; ++p) {
    if (p == pStart){
      firstBeta = betas[p];
      consensusBeta = firstBeta;
    } 
    if (firstBeta != betas[p]){
      std::cout << 
      "\nWarning: Conflict between the beta values of atoms in the same molecule" <<
      "\nName : " << name << 
      "\nMol Index : " << m <<
      "\nConflicting Indices" <<
      "\nStarting Index : " << pStart << 
      "\nConflicting Index : " << p  << std::endl;
      /* A beta == 1 functions like multiplying any number  by 0,
        even if there is only 1 atom with beta == 1, 
        the entire molecule will be fixed */
      if (firstBeta == 1.0 || betas[p] == 1.0){
        std::cout << "Beta == 1.0 takes precedent, so " << 
        "\nName : " << name << 
        "\nMol Index : " << m <<
        "\nis fixed in place."  << std::endl;
        consensusBeta = 1.0;
        break;
      } else {
        std::cout << "Beta == 2.0 takes precedent, so " << 
        "\nName: " << name << 
        "\nMol Index: " << m <<
        "\nis fixed in box " << box << std::endl;
        consensusBeta = 2.0;
      }
    }
  }
  return (uint)consensusBeta;
}

MoleculeLookup::box_iterator MoleculeLookup::box_iterator::operator++(int)
{
  box_iterator tmp = *this;
  ++(*this);
  return tmp;
}


MoleculeLookup::box_iterator::box_iterator(uint* _pLook, uint* _pSec)
  : pIt(_pLook + * _pSec) {}


MoleculeLookup::box_iterator MoleculeLookup::BoxBegin(const uint box) const
{
  return box_iterator(molLookup, boxAndKindStart + box * numKinds);
}

MoleculeLookup::box_iterator MoleculeLookup::BoxEnd(const uint box) const
{
  return box_iterator(molLookup, boxAndKindStart + (box + 1) * numKinds);
}
