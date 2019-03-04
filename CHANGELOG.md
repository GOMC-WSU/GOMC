# Change Log
All notable changes to this project will be documented in this file.

## [2.40] - 3/15/2019
+ Added support for Cyclic molecules.
+ Added Inter Molecular Exchange Monte Carlo move in GCMC and GEMC simulation.
+ Added Intra Molecular Exchange Monte Carlo move.
+ Added Crankshaft move.
+ Added separate cutoff for short range electrostatic, for each simulation box.
+ Added reporting the molecule density value for each species.
+ Added support to process the configuration file's keyword with case-insensitive.
+ Changed printing format to scientific mode.
+ Fixed the overlap detection with hard cutoff value.
+ Fixed a bug where, cutoff value was greater than half of the box length and no error was generated.
+ Fixed to the bug in the updating adjustable value in MoveSetting, when we have more than one component.
+ Fixed to the bug in GPU pressure calculation for NPT simulation.

## [2.31] - 5/21/2018
+ Compiling problem fixed on CYGWIN
+ Fix to the error checking of volume exchange if simulation volume became negative
+ Fix the equation to impose fix angle if the angle is less than 90. Generate error if constrained angle is not possible.
+ Fixed a bug where move timings were 0 in windows.
+ Fixed the issue where input files imported from NVT to NPT had zero volume, default value of 500 will be assigned in those cases.

## [2.30] - 5/10/2018
+ Added Regrowth move
+ Cleaned up the code and removed any memory leaks and bugs.
+ Redesigned Ewald, EwaldCached, and NoEwald to be more efficient and allocate less memory.
+ Added an error when multiple atoms have zero coordinates.
+ Added an error when resorvoir is empty and should print an error message and exit.
+ Fixed some typo in output.
+ Fix a bug where seed number output and input was not the same size, so we changed int to uint while reading to avoid overflow. They both use uint now.
+ Increase memory allocation for NPT simulation in Ewald.
+ Fixed a bug in PickWeighted() function where it was return a value larger than the array.
+ Printing individual timing information for each move.
+ To avoid getting large energy we are recalculating the total system at equilbrium step.
+ Print CBMC information including first atom trial, secondary atom trials, angle trials, and dihedral trials.
+ Fix a bug where in windows the clock() function was returning a wrong value
+ Fix compiler warning for macOS in ConfigSetup file

## [2.20] - 1/2/2018
+ Non orthogonal implementation
+ Removed CUB library from our source code and we will automatically download CUB library before compiling (newest version all the time!).
+ Removed compute_20 and compute_21 since they are depricated and CUDA 9 will generate fatal error. 
+ Fixed a problem where the input reader will skip random_seed because there was a "info" print in an "else if" and would block anything after that.
+ Fixed a bug in MoveSettings where two for loops where not initialized.
+ Fixed a bug where in simulation of polar molecules, if any atom of the molecule has no VDW parameter (e.g. water), there is a chance that during insertion, atom of opposite charge will overlap and generate a very large energy value. Current mechanism of avoiding overlap is not working properly
+ Fixed a bug where STDOUT was printing out garbage because input reader was passing empty vector.
+ Fixed a problem where we were receiving seg fault at the end of the simulation (removed [] from destructor in Static.cpp). 

## [2.11] - 9/25/2017
+ Bug fixes

## [2.1] - 9/4/2017
+ Adsorption
+ Fixed the bug where GPU and CPU total energy calculations where slighty different.
+ Removed some unused variables
+ Removed some warnings
+ Fixed pow function ambiguity in some compilers like clang
+ Fixed compiling bug for Clock.h when using clang (mac users)
+ Set a maximum of 9999999999 for output energies

## [2.0] - 6/11/2017
+ NPT added
+ GPU implementation has been added to the project

## [1.9] - 12/12/2016
+ Revision on output. We now only generate one output file.
+ Bug fixes
+ Printing some hardware information
+ Changes to example files which can be found at GOMC_Examples repository

## [1.8] - 10/12/2016
+ Parallelizing CBMC branched algorithm using OpenMP
+ Bug Fixes

## [1.7] - 04/21/2016
+ Fix the bug of free Ewald class, which returns memory segmentation fault if ewald is not activated.
TO DO: IntraSwap does not work properly. The bug is going to be fixed.


##[1.6] - 04/01/2016
+ IntraSwap is added into the code.
+ Blockoutput, fluctoutput are modified to print out ewald energies.
+ Tests of this modification are not done yet.

## [1.7] - 03/24/2016
+ I/O fixes

## [1.6] - 03/24/2016
+ Code cleaning

## [1.5] - 03/24/2016
+ Support for Ewald Summation

##[1.00] - 3/15/2016
+ Fix the bug of GEMC simulation. Now all NVT, GCMC, and GEMC are working.
+ Enable the Cache version of GEMC to save time; however, the scalability is limited.

##[1.00] - 01/25/2016
+ Fix the coherence issue on two boxes simulations, including NPT, GEMC, and two boxes GCMC simulations.

+ To fix the issue, Calp() and SetupRecip() have to be recalculated everytime before and after the volume transfer; RecipSinSum and RecipCosSum arrays have to be synchronized in two boxes simulations; GEMC_NVT is different from GEMC_NPT, the difference between GEMC_NVT and GEMC_NPT requires "if" statement.

? GEMC is returning good result from computation; however, it does not return correct results. A potential bug inside of the computation logic.

## [1.1] - 01/20/2016
+ Supporting Martini forcefield

##[1.00] - 12/22/2015
+ Update CPUSide.cpp, move hist into the ifdef of GCMC. Target at the floating operation issue returned by NVT and GEMC simulations.

##[1.00] Ewald branch - 12/8/2015
+ Improve the initialization of Ewald's parameters, so that non-ewald simulation can also run on ewald code without crashing.
+ CMakeList.txt file is added.

##[1.00] - 12/4/2015
+ Added Examples of water and DME for NVT and GCMC
+ Added bin/linux, which contains all executable files for linux

## [1.0] - 10/18/2015
+ Changes for input file system.
+ Added Cell list optimization to the GPU and the Serial code.

## [0.972] - 5/11/2015
+ Fixes for PDB output for NVT in the GPU code.

## [0.971] - 5/6/2015
+ Added missing CMAKE files for GOMC serial.
+ Updates to the test systems to be compatible with the new input formats.

## [0.97] - 4/11/2015

+ Added support for grand canonical ensemble.
+ Fixed calculation of angular weights for branched molecules coupled-decoupled configuration bias algorithm.
+ Improved move adjustments for better targeting of desired acceptance rate.
+ Various minor bug fixes for fringe conditions.
+ Improvements to I/O handling, inclusion of new output types relating to grand canonical ensemble (fluctuations, energy/particle samples, and distribution of molecule counts).

