#!/bin/bash
# Check if nvcc is available
use_cuda=0
use_profiler=0
use_gtest=0
use_gcc=0
use_mpi=0
use_debug=0
ENSEMBLES=""
CMAKEARGS=""


# Check if nvcc is available
if command -v nvcc &> /dev/null
then
	use_cuda=1
	nvcc_version=($(python scripts/get_cuda_version.py))
	if [ -z "$nvcc_version" ]
	then
		echo "python command not on path. trying python3"
		nvcc_version=($(python3 scripts/get_cuda_version.py))
		if [ -z "$nvcc_version" ]
		then
			echo "python3 command not on path. Please install/add to path. Exiting..."
			exit 1
		else
			echo "python3 found"
		fi
	fi
	# Check cuda version, if less than 11 then download CUB, otherwise skip
	if [[ "$nvcc_version" < "11" ]]
	then
		# Check if ./lib/cub exists
		if [ ! -d "./lib/cub" ]; then
			cd lib
			mkdir -p temp
			cd temp
			echo "==== GOMC needs CUB library to run..."
			echo "==== Finding latest CUB library..."

			# Find the link
			LINK=$(curl -s https://api.github.com/repos/NVIDIA/cub/releases/latest | grep "zip" | cut -d : -f 2,3 | tr -d \",\ )
			echo "==== Link found at ${LINK}"

			# download the zip file 
			echo "==== Downloading the CUB library... (Shouldn't take too long)"
			wget -O cub.zip "${LINK}" > /dev/null 2>&1

			#unzip
			echo "==== Extracting the CUB library..."
			unzip cub.zip > /dev/null 2>&1

			# move the cub directory to and remove the rest
			mv */cub ../cub
			rm -rf NVIDIA*
			rm cub.zip
			cd ..
			rmdir temp
			cd ..
		else
			echo "==== cub library already exists. Skipping..."
		fi
	else
		echo "CUDA version is 11.0 or higher, no need to download CUB library! Skipping..."
		rm -rf ./lib/cub/
	fi
fi

while getopts 'mptgd' opt; do
    case "$opt" in
        p)
            use_profiler=1;;
        m)
            use_mpi=1
            CMAKEARGS+="-DGOMC_MPI=on ";;
        g)
            use_gcc=1;;
        t)
            use_gtest=1;;
        d)
            use_debug=1;;
        *)  echo 'Error in command line options' >&2
            echo "Available options are: "
            echo "-p (NVTX tags),"
            echo "-t (disables Intel compiler to allow GTests to compile),"
            echo "-m, enables MPI support (Required for Parallel Tempering)"
            echo "-d, enables Debug Mode compilation"
            echo "For combined usage: -ptmg"
            exit 1
    esac
done

shift "$(( OPTIND - 1 ))"

while [ "$#" -ne 0 ]; do
	if [[ "$1" == 'CPU' ]]; then
		ENSEMBLES+="NVT NPT GCMC GEMC "
		shift
		continue
	fi
	if [[ "$1" == 'GPU' ]]; then
		ENSEMBLES+="GPU_NVT GPU_NPT GPU_GCMC GPU_GEMC "
		shift
		continue
	fi
    case "$1" in
        NVT|NPT|GCMC|GEMC|GPU_NVT|GPU_NPT|GPU_GCMC|GPU_GEMC)                   # or just:  -t|--t*)
            ENSEMBLES+="$1 ";;
         [!-]*) echo 'Error in Ensembles' >&2
                echo 'Valid Options: {NVT|NPT|GCMC|GEMC|GPU_NVT|GPU_NPT|GPU_GCMC|GPU_GEMC}' >&2
                exit 1;;# non-option detected, break out of loop
        --)    shift; break ;; # explicit end of options detected, break
        *)  echo 'Error in command line options' >&2
            exit 1
    esac
    shift
done

mkdir -p bin
cd bin

if (( !use_gtest )); then
    if (( !use_gcc )); 
    then
        ICC_PATH="$(which icc 2> /dev/null)"
        ICPC_PATH="$(which icpc 2> /dev/null)"
        if [ -z "$ICC_PATH" ]
        then
            export CC="$(which gcc 2> /dev/null)"
            export CXX="$(which g++ 2> /dev/null)"
        else
            export CC=${ICC_PATH}
            export CXX=${ICPC_PATH}
        fi
    else
        export CC="$(which gcc 2> /dev/null)"
        export CXX="$(which g++ 2> /dev/null)"
    fi
else
    if (( use_mpi )); 
    then
        ENSEMBLES+="GOMC_NVT_MPI_Test "
		ENSEMBLES+="GOMC_NPT_MPI_Test "
		ENSEMBLES+="GOMC_GCMC_MPI_Test "
		ENSEMBLES+="GOMC_GEMC_MPI_Test "
		if(( use_cuda ))
		then
        	ENSEMBLES+="GOMC_GPU_NVT_MPI_Test "
        	ENSEMBLES+="GOMC_GPU_NPT_MPI_Test "
        	ENSEMBLES+="GOMC_GPU_GCMC_MPI_Test "
        	ENSEMBLES+="GOMC_GPU_GEMC_MPI_Test "
		fi
        CMAKEARGS+="-DGOMC_GTEST_MPI=on "
    else
        ENSEMBLES+="GOMC_NVT_Test "
        ENSEMBLES+="GOMC_NPT_Test "
        ENSEMBLES+="GOMC_GCMC_Test "
        ENSEMBLES+="GOMC_GEMC_Test "
		if(( use_cuda ))
		then
        	ENSEMBLES+="GOMC_GPU_NVT_Test "
        	ENSEMBLES+="GOMC_GPU_NPT_Test "
        	ENSEMBLES+="GOMC_GPU_GCMC_Test "
        	ENSEMBLES+="GOMC_GPU_GEMC_Test "
		fi
        CMAKEARGS+="-DGOMC_GTEST=on "
    fi
    export CC="$(which gcc 2> /dev/null)"
    export CXX="$(which g++ 2> /dev/null)"
fi

# If user hasn't specified any ensemble, cmake automatically compiles all ensembles.
# This will ensure we don't print empty for ensembles.
if [ -z "$ENSEMBLES" ];
then
	ENSEMBLES="NVT NPT GCMC GEMC"
	if (( use_cuda ))
	then
		ENSEMBLES+=" GPU_NVT GPU_NPT GPU_GCMC GPU_GEMC"
	fi
fi
echo "Ensembles To Compile: $ENSEMBLES"

if (( use_profiler )); then
    if (( use_cuda )); then
      	echo "Enabling NVTX profiling for CUDA "
	    CMAKEARGS+="-DGOMC_NVTX_ENABLED=1 "
    else
      	echo "Warning: Cannot enable NVTX profiling without CUDA enabled."
    fi
fi

if (( use_debug )); then
	echo "Enabling Debug Compilation "
	CMAKEARGS+="-DCMAKE_BUILD_TYPE=Debug "
fi

cmake .. $CMAKEARGS
make -j8 $ENSEMBLES
