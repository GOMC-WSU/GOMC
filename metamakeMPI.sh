#!/bin/bash
# Check if nvcc is available
use_cuda=0
use_profiler=0

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

# check to see if we passed the profile flag
while (( $# > 0 )); do
	if [[ "$1" != --* ]]; then
		echo "ERROR: Expected an option beginning with -- but found $1"
		echo "Available options are: --with-profiling"
		exit 1
	fi
	case "$1" in
	--with-profiling)
	shift
	use_profiler=1
	;;
	*)
		echo "ERROR: unknown option $1"
		echo "Available options are: --with-profiling"
		exit 1
	;;
    esac

    shift
done

mkdir -p bin_MPI
cd bin_MPI
ICC_PATH="$(which icc)"
ICPC_PATH="$(which icpc)"
export CC=${ICC_PATH}
export CXX=${ICPC_PATH}

if (( $use_profiler )); then
    if (( $use_cuda )); then
      	echo "Enabling NVTX profiling for CUDA "
	cmake .. -DGOMC_MPI=on -DGOMC_NVTX_ENABLED=1
    else
      	echo "Warning: Cannot enable NVTX profiling without CUDA enabled."
	cmake .. -DGOMC_MPI=on
    fi
else
	cmake ..
fi

make -j8 
