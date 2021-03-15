#!/bin/bash
# Check if nvcc is available
use_cuda=0
use_profiler=0
use_gtest=0
MPI="off"
ENSEMBLES=""

if command -v nvcc &> /dev/null
then
        nvcc_version=($(python scripts/get_cuda_version.py))
	use_cuda=1
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
        fi
fi

while getopts 'mpt' opt; do
    case "$opt" in
        p)
            use_profiler=1;;
        m)
            MPI="on";;
        t)
            use_gtest=1
            ENSEMBLES+="GOMC_Test";;
        *)  echo 'Error in command line options' >&2
            echo "Available options are: -p (NVTX tags),"
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
        GOMC_Test)
            ENSEMBLES+="$1"
            use_gtest=1;;
         [!-]*) echo 'Error in Ensembles' >&2
                echo 'Valid Options: {NVT|NPT|GCMC|GEMC|GPU_NVT|GPU_NPT|GPU_GCMC|GPU_GEMC|GOMC_Test}' >&2
                exit 1;;# non-option detected, break out of loop
        --)    shift; break ;; # explicit end of options detected, break
        *)  echo 'Error in command line options' >&2
            exit 1
    esac
    shift
done

echo "Ensembles To Compile: $ENSEMBLES"

mkdir -p bin_MPI
cd bin_MPI
ICC_PATH="$(which icc)"
ICPC_PATH="$(which icpc)"

if (( !$use_gtest )); then
    export CC=${ICC_PATH}
    export CXX=${ICPC_PATH}
fi

if (( $use_profiler )); then
    if (( $use_cuda )); then
      	echo "Enabling NVTX profiling for CUDA "
	cmake .. -DGOMC_MPI=$MPI -DGOMC_NVTX_ENABLED=1
    else
      	echo "Warning: Cannot enable NVTX profiling without CUDA enabled."
	cmake .. -DGOMC_MPI=$MPI
    fi
else
	cmake .. -DGOMC_MPI=$MPI
fi

#make -j8 $ENSEMBLES
make -j8 GOMC_MPI_Test
