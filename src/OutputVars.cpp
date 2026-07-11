/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "OutputVars.h"

#include "MoleculeLookup.h"
#include "MoveSettings.h"
#include "System.h"
#include "UnitConst.h" //For unit conversion factors

#if ENSEMBLE == GEMC
#include "MoveConst.h" //For box constants, if we're calculating Hv
#endif

OutputVars::OutputVars(System &sys, StaticVals const &statV,
                       const std::vector<std::string> &molKindNames)
    : T_in_K(statV.forcefield.T_in_K), calc(sys.calcEnergy),
      molKindNames(molKindNames) {
  InitRef(sys, statV);
  for (int b = 0; b < BOXES_WITH_U_NB; ++b) {
    compressibility[b] = 0.0;
    enthalpy[b] = 0.0;
  }
#if ENSEMBLE == GEMC
  liqBox = 0;
  vapBox = 0;
  heatOfVap = 0.0;
  for (int b = 0; b < BOXES_WITH_U_NB; ++b) {
    heatOfVap_energy_term_box[b] = 0.0;
    heatOfVap_density_term_box[b] = 0.0;
  }
#endif
}

void OutputVars::InitRef(System &sys, StaticVals const &statV) {
  volumeRef = sys.boxDimRef.volume;
  axisRef = &sys.boxDimRef.axis;
  volInvRef = sys.boxDimRef.volInv;
  energyTotRef = &sys.potential.totalEnergy;
  virialTotRef = &sys.potential.totalVirial;
  energyRef = sys.potential.boxEnergy;
  virialRef = sys.potential.boxVirial; //
  kindsRef = statV.mol.kinds;
  molLookupRef = &sys.molLookupRef;
  moveSetRef = &sys.moveSettings;
  movePercRef = statV.movePerc;
  pCalcFreq = statV.simEventFreq.pCalcFreq;
  pressureCalc = statV.simEventFreq.pressureCalc;
  virial = new Virial[BOXES_WITH_U_NB];
}

bool OutputVars::Performed(uint moveKind) {
  return (movePercRef[moveKind] > 0.0);
}

uint OutputVars::GetTries(uint box, uint sub) {
  return moveSetRef->GetTrialTot(box, sub);
}

uint OutputVars::GetAccepted(uint box, uint sub) {
  return moveSetRef->GetAcceptTot(box, sub);
}

double OutputVars::GetScale(uint box, uint sub) {
  return moveSetRef->GetScaleTot(box, sub);
}

double OutputVars::GetAcceptPercent(uint box, uint sub) {
  if (GetTries(box, sub) == 0)
    return 0.0;
  return (double)(GetAccepted(box, sub)) / (double)(GetTries(box, sub)) * 100.0;
}

void OutputVars::Init() {
  // Init vals.
  numKinds = molLookupRef->GetNumKind();

  // Allocate arrays,
  uint kTot = BOX_TOTAL * numKinds;
  numByBox = new uint[BOX_TOTAL];
  numByKindBox = new uint[kTot];
  densityByKindBox = new double[kTot];
  if (numKinds > 1)
    molFractionByKindBox = new double[kTot];
}

OutputVars::~OutputVars(void) {
  if (numByBox != NULL)
    delete[] numByBox;
  if (numByKindBox != NULL)
    delete[] numByKindBox;
  if (numKinds > 1 && molFractionByKindBox != NULL)
    delete[] molFractionByKindBox;
  if (densityByKindBox != NULL)
    delete[] densityByKindBox;
  if (virial != NULL)
    delete[] virial;
}

void OutputVars::CalcAndConvert(ulong step) {
  double rawPressure[BOXES_WITH_U_NB];

  molLookupRef->TotalAndDensity(numByBox, numByKindBox, molFractionByKindBox,
                                densityByKindBox, volInvRef);

  for (uint b = 0; b < BOX_TOTAL; b++) {
    densityTot[b] = 0.0;
    for (uint k = 0; k < numKinds; k++) {
      double density = densityByKindBox[k + numKinds * b];

      // Convert density to g/ml (which is equivalent to g/cm3)
      // To get kg/m3, multiply output densities by 1000.
      density *= unit::MOLECULES_PER_A3_TO_MOL_PER_CM3 * kindsRef[k].molMass;
      densityTot[b] += density;
    }
    densityTot[b] *= 1000;
  }

  for (uint b = 0; b < BOXES_WITH_U_NB; b++) {
    // Account for dimensionality of virial (raw "virial" is actually a
    // multiple of the true virial, based on the dimensions stress is exerted
    // in)

    if (pressureCalc) {
      if ((step + 1) % pCalcFreq == 0 || step == 0) {
        if (step != 0) {
          virialRef[b] = calc.VirialCalc(b);
          *virialTotRef += virialRef[b];
        }
        // calculate surface tension in mN/M
        surfaceTens[b] = (virialRef[b].totalTens[2][2] -
                          0.5 * (virialRef[b].totalTens[0][0] +
                                 virialRef[b].totalTens[1][1]));
        surfaceTens[b] *=
            unit::K_TO_MN_PER_M / (2.0 * axisRef->Get(b).x * axisRef->Get(b).y);

        // save the pressure tensor for print
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            // convert K to bar
            pressureTens[b][i][j] = virialRef[b].totalTens[i][j] *
                                    unit::K_MOLECULE_PER_A3_TO_BAR /
                                    volumeRef[b];
          }
        }

        virial[b] = virialRef[b];
        virial[b] /= unit::DIMENSIONALITY;
        virial[b] /= volumeRef[b];
        virial[b].Total();
        // Calculate the pressure
        rawPressure[b] = 0.0;
        for (uint k = 0; k < numKinds; k++) {
          // Instead of converting to mass first
          // (an alternate route to calculate the ideal gas pressure)
          // the form of the Boltzmann constant that kcal/mol/K is used
          // such that a single conversion factor can be applied to both
          // the ideal and virial components of the pressure.
          rawPressure[b] += densityByKindBox[k + numKinds * b];
        }
        // Finish ideal component
        rawPressure[b] *= T_in_K;
        // Add the virial component
        rawPressure[b] += virial[b].total;

        // Convert to desired units
        // ( starting: K * molecule / Angstrom^3 )
        pressure[b] = rawPressure[b];
        pressure[b] *= unit::K_MOLECULE_PER_A3_TO_BAR;
        if (numByBox[b] != 0) {
          compressibility[b] = (pressure[b]) * (volumeRef[b]) / numByBox[b] /
                               (T_in_K) /
                               (UNIT_CONST_H::unit::K_MOLECULE_PER_A3_TO_BAR);
          enthalpy[b] = (energyRef[b].total / numByBox[b] +
                         rawPressure[b] * volumeRef[b] / numByBox[b]) *
                        UNIT_CONST_H::unit::K_TO_KJ_PER_MOL;
        } else {
          compressibility[b] = 0.0;
          enthalpy[b] = 0.0;
        }
      }
    }
  }
#if ENSEMBLE == GEMC
  if (pressureCalc) {
    if ((step + 1) % pCalcFreq == 0 || step == 0) {
      // Determine which box is liquid for purposes of heat of vap.
      if (densityTot[mv::BOX1] >= densityTot[mv::BOX0]) {
        vapBox = mv::BOX0;
        liqBox = mv::BOX1;
      } else {
        vapBox = mv::BOX1;
        liqBox = mv::BOX0;
      }
      // delta Hv = (Uv-Ul) + P(Vv-Vl)
      if (numByBox[vapBox] != 0) {
        heatOfVap_energy_term_box[vapBox] =
            energyRef[vapBox].total / numByBox[vapBox];
        heatOfVap_density_term_box[vapBox] =
            volumeRef[vapBox] / numByBox[vapBox];
      } else {
        heatOfVap_energy_term_box[vapBox] = 0.0;
        heatOfVap_density_term_box[vapBox] = 0.0;
      }

      if (numByBox[liqBox] != 0) {
        heatOfVap_energy_term_box[liqBox] =
            energyRef[liqBox].total / numByBox[liqBox];
        heatOfVap_density_term_box[liqBox] =
            volumeRef[liqBox] / numByBox[liqBox];
      } else {
        heatOfVap_energy_term_box[liqBox] = 0.0;
        heatOfVap_density_term_box[liqBox] = 0.0;
      }

      heatOfVap =
          heatOfVap_energy_term_box[vapBox] - heatOfVap_energy_term_box[liqBox];
      heatOfVap += rawPressure[vapBox] * (heatOfVap_density_term_box[vapBox] -
                                          heatOfVap_density_term_box[liqBox]);
      heatOfVap *= unit::K_TO_KJ_PER_MOL;
    }
  }
#endif
}
