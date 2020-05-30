#!/bin/bash

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

mkdir -p bin
cd bin
ICC_PATH="$(which icc)"
ICPC_PATH="$(which icpc)"
export CC=${ICC_PATH}
export CXX=${ICPC_PATH}
cmake ..
make -j8
