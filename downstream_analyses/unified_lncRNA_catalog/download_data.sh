#!/usr/bin/env bash

########################
###### Description
########################
# Helps to donwload all necessary data.

# Handle errors
set -e          # exit on any non-0 exit status
set -o pipefail # exit on any non-0 exit status in pipe

# Download data for /Data/Source/Additional/ directory
echo "~~~~~~~~~ Downloading files for /Data/Source/Additional/ directory"
cd ./Data/Source/Additional/
chmod +x README.sh
./README.sh
echo "+++++ DONE +++++"
cd -

# Download data for /Data/Source/masterTable/ directory
echo "~~~~~~~~~ Downloading files for /Data/Source/masterTable/ directory"
cd ./Data/Source/masterTable/
chmod +x README.sh
./README.sh
echo "+++++ DONE +++++"
cd -

# Notify when completed
echo "######## COMPLETED ########"
