# GOMC - GPU Optimized Monte Carlo

Current Release: 2.40 (5/9/2019)

[![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/GOMC_WSU/Lobby?utm_source=share-link&utm_medium=link&utm_campaign=share-link)
[![Build Status](https://travis-ci.org/GOMC-WSU/GOMC.svg?branch=master)](https://travis-ci.org/GOMC-WSU/GOMC)

We recommend the [GOMC Project Website](http://gomc.eng.wayne.edu/ "GOMC Website") and the [user manual](https://gomc-wsu.github.io/Manual/ "User Manual") for further information and examples.

To cite GOMC project, please use [GOMC SoftwareX paper](https://www.sciencedirect.com/science/article/pii/S2352711018301171?via%3Dihub "SoftwareX").

## Building GOMC on GNU/Linux, macOS, or Cygwin:

  1. Clone or download our code from GitHub:
      ```bash
      git clone https://github.com/GOMC-WSU/GOMC.git
      ```
  2. Go into the GOMC directory: 
      ```bash
      cd GOMC
      ```
  3. Give execution permission: 
      ```bash
      chmod u+x metamake.sh
      ```
  4. Run metamake file:
      ```bash
      ./metamake.sh
      ```
  5. Step 4 should generate all the executables in ```bin``` directory

  You can set the number of the threads using the +pN argument, where N is the number of threads.
  For example:
  ```bash
  ./GOMC_<CPU|GPU>_XXXX +p4 in.conf
  ```

  Which will run 4 threads and reads input file "in.conf".

  NOTES:
  Building GOMC requires cmake, available at http://www.cmake.org and in most Linux package repositories (as cmake).
  If you wish to utilize NVIDIA graphic cards you will need to install NVIDIA toolkit before compiling. The metamake file will automatically detect the location of CUDA installation. (More info in Manual)

## BUILDING GOMC ON WINDOWS:
  1. Open the Windows-compatible CMake GUI.
  2. Set the Source Folder to the GOMC root folder.
  3. Set the build Folder to your Build Folder.
  4. Click configure, select your compiler/environment
  5. Wait for CMake to finish the configuration.
  6. Click configure again and click generate.
  7. Open the CMake-generated project/solution etc. to the desired IDE (e.g Visual Studio).
  8. Using the solution in the IDE of choice build GOMC per the IDE's standard release compilation/executable generation methods.

   NOTES:
      You can also use CMake from the Windows command line if its directory is
      added to the PATH environment variable.
