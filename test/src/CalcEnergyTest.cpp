#include <gtest/gtest.h>
#include "Simulation.h"
#include<unistd.h> 

TEST(CalcEnergyTest, CheckISOPEN) {
  double tempREn = 0.0, tempLJEn = 0.0;
  std::vector<int> cellVector, cellStartIndex, mapParticleToCell;
  std::vector< std::vector<int> > neighborList;
  cellList.GetCellListNeighbor(box, currentCoords.Count(),
                               cellVector, cellStartIndex, mapParticleToCell);
  neighborList = cellList.GetNeighborList(box);

#ifdef GOMC_CUDA
  //update unitcell in GPU
  UpdateCellBasisCUDA(forcefield.particles->getCUDAVars(), box,
                      boxAxes.cellBasis[box].x, boxAxes.cellBasis[box].y,
                      boxAxes.cellBasis[box].z);

  if(!boxAxes.orthogonal[box]) {
    //In this case, boxAxes is really an object of type BoxDimensionsNonOrth,
    // so cast and copy the additional data to the GPU
    const BoxDimensionsNonOrth *NonOrthAxes = static_cast<const BoxDimensionsNonOrth*>(&boxAxes);
    UpdateInvCellBasisCUDA(forcefield.particles->getCUDAVars(), box,
                           NonOrthAxes->cellBasis_Inv[box].x,
                           NonOrthAxes->cellBasis_Inv[box].y,
                           NonOrthAxes->cellBasis_Inv[box].z);
  }

  CallBoxInterGPU(forcefield.particles->getCUDAVars(), cellVector, cellStartIndex,
                  neighborList, coords, boxAxes, electrostatic, particleCharge,
                  particleKind, particleMol, tempREn, tempLJEn, forcefield.sc_coul,
                  forcefield.sc_sigma_6, forcefield.sc_alpha,
                  forcefield.sc_power, box);
  double sumREn = 0.0, sumLJEn = 0.0;
  double sumREnMol = 0.0, sumLJEnMol = 0.0;
  double sumREnMolAbb = 0.0, sumLJEnMolAbb = 0.0;
  MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box);
  MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);
  while(thisMol!=end){
      uint mol = *thisMol;
      std::cout << "Box " << box <<  " Mol " << mol << std::endl;
      double tempREnMol = 0.0, tempLJEnMol = 0.0;
      double tempREnMolAbb = 0.0, tempLJEnMolAbb = 0.0;
      uint length = mols.GetKind(mol).NumAtoms();
      uint start = mols.MolStart(mol);
                        
      CallMolInterGPU(forcefield.particles->getCUDAVars(), 
                  start, length, cellVector, cellStartIndex,
                  neighborList, currentCoords, currentAxes, electrostatic, particleCharge,
                  particleKind, particleMol, tempREnMol, tempLJEnMol, forcefield.sc_coul,
                  forcefield.sc_sigma_6, forcefield.sc_alpha,
                  forcefield.sc_power, box);

    XYZArray molCoords(1);
    std::vector<int> molCoordsToCell(1);
    CallMolInterGPU(forcefield.particles->getCUDAVars(), 
                  start, length, cellVector, cellStartIndex,
                  neighborList, currentCoords, molCoords, mapParticleToCell, molCoordsToCell, currentAxes, electrostatic, particleCharge,
                  particleKind, particleMol, tempREnMolAbb, tempLJEnMolAbb, forcefield.sc_coul,
                  forcefield.sc_sigma_6, forcefield.sc_alpha,
                  forcefield.sc_power, box);
                  
      sumREn += tempREnMol;
      sumLJEn += tempLJEnMol;
      sumREnMolAbb += tempREnMolAbb;
      sumLJEnMolAbb += tempLJEnMolAbb;
      thisMol++;
  }
  /*
  if (fabs(sumREn - tempREn) >= std::numeric_limits<double>::epsilon()){
    std::cout << "sumREn not eq tempRen" << sumREn << " " << tempREn << std::endl;
    exit(1);
  } else {
    std::cout << "sumREn eq tempRen" << sumREn << " " << tempREn << std::endl;
  }
  if (fabs(sumLJEn - tempLJEn) >= std::numeric_limits<double>::epsilon()){
    std::cout << "sumLJEn not eq tempLJEn " << sumLJEn << " " << tempLJEn << std::endl;
    exit(1);
  } else {
    std::cout << "sumLJEn eq tempLJEn" << sumLJEn << " " << tempLJEn << std::endl;
  }
  if (fabs(sumREnMolAbb - tempREn) >= std::numeric_limits<double>::epsilon()){
    std::cout << "sumREn not eq tempRen" << sumREn << " " << tempREn << std::endl;
    exit(1);
  } else {
    std::cout << "sumREn eq tempRen" << sumREn << " " << tempREn << std::endl;
  }
  if (fabs(sumLJEnMolAbb - tempLJEn) >= std::numeric_limits<double>::epsilon()){
    std::cout << "sumLJEn not eq tempLJEn " << sumLJEn << " " << tempLJEn << std::endl;
    exit(1);
  } else {
    std::cout << "sumLJEn eq tempLJEn" << sumLJEn << " " << tempLJEn << std::endl;
  }
  */
  if (fabs(sumREnMolAbb - sumREn) >= 0.001){
    std::cout << "sumREn not eq sumREnMolAbb" << sumREn << " " << sumREnMolAbb << std::endl;
    exit(1);
  } else {
    std::cout << "sumREn eq sumREnMolAbb" << sumREn << " " << sumREnMolAbb << std::endl;
  }
  if (fabs(sumLJEnMolAbb - sumLJEn) >= 0.001){
    std::cout << "sumLJEn not eq sumLJEnMolAbb " << sumLJEn << " " << sumLJEnMolAbb << std::endl;
    exit(1);
  } else {
    std::cout << "sumLJEn eq sumLJEnMolAbb" << sumLJEn << " " << sumLJEnMolAbb << std::endl;
  }
}