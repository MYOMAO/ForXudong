#!/bin/bash

# Source the sphenix_setup.csh script
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.358

# Remove the 'install' directory
rm -r install

# Create a new 'install' directory
mkdir install

# Set environment variables MYINSTALL and LD_LIBRARY_PATH
export MYINSTALL=$PWD/install/
export LD_LIBRARY_PATH=$MYINSTALL/lib:$LD_LIBRARY_PATH

# Update the PATH environment variable
export PATH=$MYINSTALL/bin:$PATH

# Change directory to 'HFMLTrigger_LANL'
cd coresoftware/simulation/g4simulation/

cd g4decayer/

# Clean the project
make clean
# Configure and install the project
./autogen.sh --prefix=$MYINSTALL
make -j10 install

# Change back to the previous directory
cd ..

# Change directory to 'AntiTrigger'
cd g4main/

# Configure and install the AntiTrigger project
make clean
./autogen.sh --prefix=$MYINSTALL
make -j10 install

# Change back to the previous directory
cd ..

cd ../../../

# Add the specified directory to ROOT_INCLUDE_PATH
export ROOT_INCLUDE_PATH=${PWD}/macros/common:$ROOT_INCLUDE_PATH
