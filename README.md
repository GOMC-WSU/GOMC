# GOMC - GPU Optimized Monte Carlo

Current Release: 2.75 (6/21/2022)

[![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/GOMC_WSU/Lobby?utm_source=share-link&utm_medium=link&utm_campaign=share-link)
[![Build Status](https://travis-ci.org/GOMC-WSU/GOMC.svg?branch=master)](https://travis-ci.org/GOMC-WSU/GOMC)

We recommend the [GOMC Project Website](http://gomc.eng.wayne.edu/ "GOMC Website") and the [user manual](https://gomc-wsu.github.io/Manual/ "User Manual") for further information and examples.

To cite GOMC project, please use cite the following papers:
1.  [Y. Nejahi, M. Soroush Barhaghi,  G. Schwing, L. Schwiebert, J. Potoff. SoftwareX, 13, 100627 (2021)](https://www.sciencedirect.com/science/article/pii/S235271102030340X)
2.  [Y. Nejahi, M. Soroush Barhaghi, J. Mick, B. Jackman, K. Rushaidat, Y. Li, L. Schwiebert, J. Potoff. SoftwareX, 9, 20-27 (2019)](https://www.sciencedirect.com/science/article/pii/S2352711018301171?via%3Dihub "SoftwareX")

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
  5. Step 4 should generate all the executables in ```bin``` directory.

  `./metamake.sh` accepts flags which indicates which ensembles to compile. Default behavior with no flag will compile all CPU compilers and if CUDA available, all GPU ensembles. Multiple flags can be used by separating with a space. Current accepted flags are: `CPU` to compile all CPU ensembles, `GPU` to compile all GPU ensembles, or you can compile ensembles individually by using any of the following flags:
  `NVT`, `NPT`, `GCMC`, `GEMC`, `GPU_NVT`, `GPU_NPT`, `GPU_GCMC`, `GPU_GEMC`.

> NOTES: Building GOMC requires cmake, available at http://www.cmake.org and in most Linux package repositories (as cmake). If you wish to utilize NVIDIA graphic cards you will need to install NVIDIA toolkit before compiling. The metamake file will automatically detect the location of CUDA installation. (More info in Manual)

## Building GOMC in Docker (Linux)

  1. Clone or download our code from GitHub:
      ```bash
      git clone https://github.com/GOMC-WSU/GOMC.git
      ```
  2. Go into the GOMC directory: 
      ```bash
      cd GOMC
      ```
  3. Download/Install Docker:
      ```bash
      sudo apt-get update && sudo apt-get install -y docker.io
      ```
  4. Issue docker build command:
      ```bash
      docker build -t gomc/gomc:cpu -f dockerfiles/GOMC_CPU.dockerfile
      ```

  5. Run docker images:
      ```bash
      sudo docker run -i -t gomc/gomc:cpu
      ```

## Building GOMC on Windows:
  1. Open the Windows-compatible CMake GUI.
  2. Set the Source Folder to the GOMC root folder.
  3. Set the build Folder to your Build Folder.
  4. Click configure, select your compiler/environment
  5. Wait for CMake to finish the configuration.
  6. Click configure again and click generate.
  7. Download [CUB library](https://nvlabs.github.io/cub/download_cub.html)
  8. Extract CUB library and copy the "cub" folder from CUB library into "lib" folder inside GOMC directory.
  9. Open the CMake-generated project/solution etc. to the desired IDE (e.g Visual Studio).
  10. Using the solution in the IDE of choice build GOMC per the IDE's standard release compilation/executable generation methods.

> NOTES: You can also use CMake from the Windows command line if its directory is added to the PATH environment variable.

## Executing GOMC:
  You can set the number of the threads using the +pN argument, where N is the number of threads.
  For example:
  ```bash
  ./GOMC_<CPU|GPU>_XXXX +p4 in.conf
  ```

  Which will run 4 threads and reads input file "in.conf".

