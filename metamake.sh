#!/bin/bash

# Initialize the commandline options and flags
use_cuda=0
use_profiler=0
use_gtest=0
use_gcc=0
use_clang=0
use_mpi=0
use_asan=0
use_opt=1
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

while getopts 'acdglmnpt' opt; do
    case "$opt" in
        a)
            use_asan=1;;
        c)
            CMAKEARGS+="-DGOMC_TIDY=on ";;
        d)
            use_debug=1;;
        g)
            use_gcc=1;;
        l)
            use_clang=1;;
        m)
            use_mpi=1
            CMAKEARGS+="-DGOMC_MPI=on ";;
        n)
            use_opt=0;;
        p)
            use_profiler=1;;
        t)
            use_gtest=1
            use_gcc=1;;
        *)  echo 'Error in command line options' >&2
            echo "Available options are: "
            echo "-a, enables address sanitizer runtime checking"
            echo "-c, enables clang-tidy source code checks"
            echo "-d, enables Debug Mode compilation"
            echo "-g, use the GNU compiler"
			echo "-l, use the Clang compiler"
            echo "-m, enables MPI support (Required for Parallel Tempering)"
            echo "-n, disables most compiler optimization flags"
            echo "-p enables GPU code profiling (NVTX tags)"
            echo "-t disables Intel compiler to allow GTests to compile"
            echo "For combined usage, concatenate flags, e.g.: -ptmg"
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

mkdir -p bin
cd bin

if (( !use_gtest )); then
    if (( !use_gcc && !use_clang ));
    then
# comment out this check until CUDA supports the newer Intel Compiler
#        ICC_PATH="$(which icx 2> /dev/null)"
#        ICPC_PATH="$(which icpx 2> /dev/null)"
#        if [ -z "$ICC_PATH" ]
#        then
            ICC_PATH="$(which icc 2> /dev/null)"
            ICPC_PATH="$(which icpc 2> /dev/null)"
#		fi
        if [ -z "$ICC_PATH" ]
        then
            export CC="$(which gcc 2> /dev/null)"
            export CXX="$(which g++ 2> /dev/null)"
        else
            if (( use_asan )); then
				echo "Warning: Address sanitizer unset. Not compatible with the Intel compiler."
				use_asan=0
			fi
            export CC=${ICC_PATH}
            export CXX=${ICPC_PATH}
        fi
	elif (( use_clang )); then
        CLANG_PATH="$(which clang 2> /dev/null)"
        CLANGXX_PATH="$(which clang++ 2> /dev/null)"
        if [ -z "$CLANG_PATH" ]
        then
            export CC="$(which gcc 2> /dev/null)"
            export CXX="$(which g++ 2> /dev/null)"
	    else
            export CC=${CLANG_PATH}
            export CXX=${CLANGXX_PATH}
		fi
    else
        export CC="$(which gcc 2> /dev/null)"
        export CXX="$(which g++ 2> /dev/null)"
	fi
else
    if (( use_mpi )); 
    then
        TESTENS=""
        for ENS in $ENSEMBLES
        do
            TESTENS+="GOMC_"$ENS"_MPI_Test "
        done
        ENSEMBLES+=$TESTENS
        CMAKEARGS+="-DGOMC_GTEST_MPI=on "
    else
        TESTENS=""
        for ENS in $ENSEMBLES
        do
            TESTENS+="GOMC_"$ENS"_Test "
        done
        ENSEMBLES+=$TESTENS
        CMAKEARGS+="-DGOMC_GTEST=on "
    fi
    export CC="$(which gcc 2> /dev/null)"
    export CXX="$(which g++ 2> /dev/null)"
fi

echo "Ensembles To Compile: $ENSEMBLES"

if (( use_profiler )); then
    if (( use_cuda )); then
      	echo "Enabling NVTX profiling for CUDA "
	    CMAKEARGS+="-DGOMC_NVTX_ENABLED=on "
    else
      	echo "Warning: Cannot enable NVTX profiling without CUDA enabled."
    fi
fi

if (( use_asan )); then
    use_debug=1
    CMAKEARGS+="-DGOMC_ASAN=on "
fi

if (( use_opt )); then
    CMAKEARGS+="-DGOMC_OPT=on "
fi

if (( use_debug )); then
	echo "Enabling Debug Compilation "
	CMAKEARGS+="-DCMAKE_BUILD_TYPE=Debug "
fi

cmake .. $CMAKEARGS
make -j8 $ENSEMBLES
