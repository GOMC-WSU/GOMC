# GOMC - GPU Optimized Monte Carlo

GOMC is an acronym for GPU Optimized Monte Carlo. GOMC is a parallel molecular simulation code designed for high-performance simulation of large systems that runs on both CPUs, with or without OpenMP, and NVIDIA GPUs.

Current Release: 2.80 (12/03/2025)

The latest stable release is available on the main branch. The development branch is used for staging changes that should be stable but have not been thoroughly tested. If you encounter a problem running GOMC using the main branch, please try the development branch to see if your problem has already been fixed. If not, please post a GitHub issue, with the requested documents and details, so we can investigate the problem.

You can find many example configurations for running GOMC in our [GOMC Examples repository](https://github.com/GOMC-WSU/GOMC_Examples "GOMC Examples"). We recommend the [GOMC Project Website](https://gomc-wsu.org/ "GOMC Website") and the [user manual](https://gomc-wsu.github.io/Manual/ "User Manual") for additional information and examples.

If using GOMC in your work, we ask that you cite the following article:
[Y. Nejahi, M. Soroush Barhaghi, J. Mick, B. Jackman, K. Rushaidat, Y. Li, L. Schwiebert, J. Potoff. SoftwareX, 9, 20-27 (2019). doi: 10.1016/j.softx.2018.11.005.](https://www.sciencedirect.com/science/article/pii/S2352711018301171?via%3Dihub "SoftwareX 2019")

If you are using some more recent features of GOMC, such as the MultiParticle move, please consider also citing the following article:
[Y. Nejahi, M. Soroush Barhaghi,  G. Schwing, L. Schwiebert, J. Potoff. SoftwareX, 13, 100627 (2021). doi: 10.1016/j.softx.2020.100627.](https://www.sciencedirect.com/science/article/pii/S235271102030340X "SoftwareX 2021")

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

> NOTE: Building GOMC requires [CMake](https://cmake.org/) version 3.18 or newer. CMake is available in most Linux distributions (as cmake). If you wish to utilize NVIDIA graphics cards you will need to install the NVIDIA toolkit before compiling. The metamake file will automatically detect the location of your CUDA installation. More detailed info can be found in the [user manual](https://gomc-wsu.github.io/Manual/) "User Manual".

## Building GOMC on Windows:

1.  If building GPU executables and the CUDA version is older than CUDA 11, download the [CUB library](https://nvlabs.github.io/cub/download_cub.html).
2.  If building GPU executables and the CUDA version is older than CUDA 11, extract the CUB library and copy the "cub" folder from the CUB library into the "lib" folder inside the GOMC directory.
3.  Open the Windows-compatible CMake GUI.
4.  Set the Source Folder to the GOMC root folder.
5.  Set the Build Folder to your build folder.
6.  Click `Configure`, select your compiler/environment.
8.  Click `Generate` after CMake finishes configurating the project.
9.  Click `Open Project` after CMake finishes generating the project.
10. Using the solution in the IDE, build GOMC per the IDE's standard release compilation/executable generation methods.

> NOTE: You can also use CMake from the Windows command line if its directory is added to the PATH environment variable.

## Executing GOMC:

  You can set the number of OpenMP threads using the +pN argument, where N is the number of threads.
  For example:
  ```bash
  ./GOMC_GPU_GEMC +p4 in.conf
  ```

runs a simulation with the Gibbs ensemble on the GPU using 4 OpenMP threads and loads configuration settings from the file "in.conf".
