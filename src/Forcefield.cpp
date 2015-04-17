/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Forcefield.h" //Header spec.
//Setup partner classes
#include "Setup.h"
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <math.h>

const double BOLTZMANN = 0.0019872041;

void Forcefield::Init(const Setup& set)
{
   InitBasicVals(set.config.sys.ff);
   particles.Init(set.ff.mie, rCut);
   bonds.Init(set.ff.bond);
   angles.Init(set.ff.angle);
   dihedrals.Init(set.ff.dih);
}

void Forcefield::InitBasicVals(config_setup::FFValues const& val)
{
   useLRC = val.doTailCorr;
   T_in_K = val.temperature; 
   rCut = val.cutoff; 
   rCutSq = rCut * rCut;
   rCutOver2 = rCut / 2.0;
   scl_14 = val.oneFourScale;
   beta = 1/T_in_K;
}

