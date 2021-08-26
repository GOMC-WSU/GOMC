#!/bin/bash
# Check if nvcc is available
use_cuda=0
use_profiler=0
use_gtest=0
use_mpi=0
MPI="off"
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

			# download the download html page
			wget https://nvlabs.github.io/cub/download_cub.html > /dev/null 2>&1

			# find the lines that have the link
			grep "https://github.com/NVlabs/" download_cub.html > link_lines

			# the last line is the easiest to find the link
			awk '/./{line=$0} END{print line}' link_lines > last_line

			# the substring between two quotes is the link!!!!
			LINK="$(awk -F'"' '{ print $2 }' last_line)"
			echo "==== Link found at ${LINK}"

			# remove any temporary files
			rm link_lines
			rm download_cub.html
			rm last_line

			# download the zip file 
			echo "==== Downloading the CUB library... (Shouldn't take too long)"
			wget "${LINK}" > /dev/null 2>&1

			#unzip
			echo "==== Extracting the CUB library..."
			for z in *.zip; do
				unzip "$z" > /dev/null 2>&1
				rm "$z" > /dev/null 2>&1
			done

			# move the cub directory to and remove the rest
			for d in */ ; do
				mv "$d"/cub ../cub > /dev/null 2>&1
				rm -r "$d" > /dev/null 2>&1
			done
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

while getopts 'mpt' opt; do
    case "$opt" in
        p)
            use_profiler=1;;
        m)
            use_mpi=1
            CMAKEARGS+="-DGOMC_MPI=on ";;
        t)
            use_gtest=1;;
        *)  echo 'Error in command line options' >&2
            echo "Available options are: "
            echo "-p (NVTX tags),"
            echo "-t (disables Intel compiler to allow GTests to compile),"
            echo "-m, enables MPI support (Required for Parallel Tempering)"
            echo "For combined usage: -ptm"
            exit 1
    esac
done

shift "$(( OPTIND - 1 ))"

while [ "$#" -ne 0 ]; do
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

mkdir -p bin_MPI
cd bin_MPI

if (( !$use_gtest )); then
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
    if (( $use_mpi )); 
    then
        ENSEMBLES+="GOMC_NVT_MPI_Test "
		ENSEMBLES+="GOMC_NPT_MPI_Test "
        CMAKEARGS+="-DGOMC_GTEST_MPI=on "
    else
        ENSEMBLES+="GOMC_Test"
        CMAKEARGS+="-DGOMC_GTEST=on "
    fi
    export CC="$(which gcc 2> /dev/null)"
    export CXX="$(which g++ 2> /dev/null)"
fi

echo "Ensembles To Compile: $ENSEMBLES"

if (( $use_profiler )); then
    if (( $use_cuda )); then
      	echo "Enabling NVTX profiling for CUDA "
	    CMAKEARGS+="-DGOMC_NVTX_ENABLED=1 "
    else
      	echo "Warning: Cannot enable NVTX profiling without CUDA enabled."
    fi
fi

cmake .. $CMAKEARGS -DCMAKE_BUILD_TYPE=Debug
make -j8 $ENSEMBLES
#make -j8 GOMC_MPI_Test
