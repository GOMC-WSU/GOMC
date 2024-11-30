# GOMC - GPU Optimized Monte Carlo

Current Release: 2.75 (6/21/2022)

[![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/GOMC_WSU/Lobby?utm_source=share-link&utm_medium=link&utm_campaign=share-link)
[![Build Status](https://travis-ci.org/GOMC-WSU/GOMC.svg?branch=master)](https://travis-ci.org/GOMC-WSU/GOMC)

We recommend the [GOMC Project Website](https://gomc-wsu.org/ "GOMC Website") and the [user manual](https://gomc-wsu.github.io/Manual/ "User Manual") for further information and examples.

To cite GOMC project, please cite the following papers:
1.  [Y. Nejahi, M. Soroush Barhaghi,  G. Schwing, L. Schwiebert, J. Potoff. SoftwareX, 13, 100627 (2021). doi: 10.1016/j.softx.2020.100627.](https://www.sciencedirect.com/science/article/pii/S235271102030340X)
2.  [Y. Nejahi, M. Soroush Barhaghi, J. Mick, B. Jackman, K. Rushaidat, Y. Li, L. Schwiebert, J. Potoff. SoftwareX, 9, 20-27 (2019). doi: 10.1016/j.softx.2018.11.005.](https://www.sciencedirect.com/science/article/pii/S2352711018301171?via%3Dihub "SoftwareX")

## Building GOMC on GNU/Linux, macOS, or Cygwin:

1.  Clone or download our code from GitHub:
     ```bash
     git clone https://github.com/GOMC-WSU/GOMC.git
     ```
2.  Go into the GOMC directory: 
     ```bash
     cd GOMC
     ```
3.  Give execution permission: 
     ```bash
     chmod u+x metamake.sh
     ```
4.  Run metamake file:
     ```bash
     ./metamake.sh
     ```
5.  Step 4 will place all the executables in the ```bin``` directory.

  `./metamake.sh` accepts a list of which ensembles to compile. Default behavior, listing no ensembles, is to compile all CPU ensembles and, if CUDA is available, all GPU ensembles. Multiple ensemble names must be separated by spaces. Current accepted values are: `CPU` to compile all CPU ensembles, `GPU` to compile all GPU ensembles, or you can compile ensembles individually by using any of the following keywords:
  `NVT`, `NPT`, `GCMC`, `GEMC`, `GPU_NVT`, `GPU_NPT`, `GPU_GCMC`, `GPU_GEMC`.

> NOTES: Building GOMC requires [CMake](https://cmake.org/) version 3.18 or newer. CMake is available in most Linux package repositories (as cmake). If you wish to utilize NVIDIA graphics cards you will need to install the NVIDIA toolkit before compiling. The metamake file will automatically detect the location of your CUDA installation. More detailed info can be found in the [user manual](https://gomc-wsu.github.io/Manual/) "User Manual".

## Building GOMC on Windows:
1.  Open the Windows-compatible CMake GUI.
2.  Set the Source Folder to the GOMC root folder.
3.  Set the Build Folder to your build folder.
4.  Click Configure, select your compiler/environment.
5.  Wait for CMake to finish creating the configuration.
6.  Click Generate.
7.  If building GPU executables and the CUDA version is older than CUDA 11, download the [CUB library](https://nvlabs.github.io/cub/download_cub.html).
8.  If building GPU executables and the CUDA version is older than CUDA 11, extract the CUB library and copy the "cub" folder from the CUB library into the "lib" folder inside the GOMC directory.
9.  Open the CMake-generated project/solution file, located in your Build Folder, in the desired IDE (e.g., Visual Studio).
10. Using the solution in the IDE, build GOMC per the IDE's standard release compilation/executable generation methods.

> NOTES: You can also use CMake from the Windows command line if its directory is added to the PATH environment variable.

## Executing GOMC:
  You can set the number of CPU threads using the +pN argument, where N is the number of threads.
  For example:
  ```bash
  ./GOMC_GPU_GEMC +p4 in.conf
  ```

  will run a simulation with the Gibbs ensemble on the GPU using 4 CPU threads and loads configuration settings from the file "in.conf".
