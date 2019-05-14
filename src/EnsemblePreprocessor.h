/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef ENSEMBLE_PREPROCESSOR_H
#define ENSEMBLE_PREPROCESSOR_H

//Value table of ensembles, for preprocessor conditionals
#define NVT 1
#define GEMC 2
#define GCMC 3
#define NPT 4

#ifndef ENSEMBLE
//The choice of ensemble.
#define ENSEMBLE GCMC
#endif

//Ensemble specific defines, such as what data members are variable
//... will be used later to decide what moves to include.
#if ENSEMBLE == GEMC
#define VARIABLE_PARTICLE_NUMBER
#define VARIABLE_VOLUME
#define VARIABLE_DENSITY
#define BOX_TOTAL 2
//Do not do vol. trans. if molecule is over the boundary
//Towhee does this for volume transfers
//#define NO_VOL_TRANS_IF_OVER_PERIODIC_BOUND
#elif ENSEMBLE == GCMC
#define VARIABLE_PARTICLE_NUMBER
#define VARIABLE_DENSITY
#define BOX_TOTAL 2
#elif ENSEMBLE == NVT
#define VARIABLE_PARTICLE_NUMBER
#define BOX_TOTAL 1
#elif ENSEMBLE == NPT
#define VARIABLE_PARTICLE_NUMBER
#define VARIABLE_VOLUME
#define VARIABLE_DENSITY
#define BOX_TOTAL 1
#endif

#endif /*ENSEMBLE_PREPROCESSOR_H*/
