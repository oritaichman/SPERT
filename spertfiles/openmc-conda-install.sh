#!/bin/bash

# This script installs OpenMC, its Python API, and all necessary dependencies
# into a conda environment. Because it uses compilers and CMake from the
# Anaconda repository, it's not necessary to have those installed on your
# system. However, at a minimum, you do need to have 'make' installed as well as
# a linker. To get the OpenMC git repository, you'll need git installed as
# well. On Debian derivatives, you can install all of these with:
#
# sudo apt install -y make binutils git

set -e

#===============================================================================
# INSTALLATION OPTIONS

INSTALL_MINICONDA=no                              # Install miniconda?
INSTALL_IN_CONDA=yes                              # Install OpenMC executable/package directly in conda environment?
INSTALL_EDITABLE=no                               # Install in editable mode (overrides INSTALL_IN_CONDA)
INSTALL_PREFIX=$HOME/opt/openmc/openmc-test       # Where to install (if INSTALL_IN_CONDA and INSTALL_EDITABLE are disabled)
WITH_MPI=no                                       # Install with MPI support
PYTHON_VERSION=3.7                                # What version of Python to use
ENVNAME=openmc-py$PYTHON_VERSION                  # conda environment name

#===============================================================================

# Download and install miniconda
if [[ $INSTALL_MINICONDA == yes ]]; then
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b
    $HOME/miniconda3/bin/conda update -y conda
    $HOME/miniconda3/bin/conda init bash
fi

# Make sure conda is activated
. $HOME/miniconda3/etc/profile.d/conda.sh
conda activate

# Create new Python environment to install everything into
conda create -y -n $ENVNAME python=$PYTHON_VERSION
conda activate $ENVNAME

# Install recent version of gcc and CMake. OpenMC requires gcc 4.9+ and CMake
# 3.3+, so if you are on an older distribution, it's possible that your version
# is not sufficient.
conda install -y gcc_linux-64 gxx_linux-64 cmake

# Install OpenMC's Python dependencies
conda install -y numpy scipy cython pandas lxml matplotlib vtk ipython pytest

# The uncertainties package is not in the default conda channel
conda config --append channels conda-forge/label/gcc7
conda install -y uncertainties

# Install MPICH for building with MPI support. Make sure to use the package
# based on gcc7. For depletion, we also need the parallel version of h5py
if [[ $WITH_MPI == yes ]]; then
    conda install -y mpich h5py=*=*mpich*
    export CXX=mpicxx
else
    conda install -y h5py
fi

# Set up temporary directory
WORK_DIR=$(mktemp -d)
if [[ ! "$WORK_DIR" || ! -d "$WORK_DIR" ]]; then
  echo "Could not create temporary directory"
  exit 1
fi

# Make sure temporary directory is deleted if script exits
function cleanup {
    rm -rf "$WORK_DIR"
    echo "Deleted temporary working directory $WORK_DIR"
}
trap cleanup EXIT

# Clone OpenMC
cd $WORK_DIR
git clone https://github.com/openmc-dev/openmc

# Build OpenMC library and executable
cd openmc
mkdir -p build && pushd build
if [[ $INSTALL_IN_CONDA == yes ]]; then
    cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX ..
else
    cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX ..
fi
popd

# Install Python API
if [[ $INSTALL_EDITABLE == yes ]]; then
    pip install --no-dependencies --no-build-isolation -e .
else
    if [[ $INSTALL_IN_CONDA == yes ]]; then
        pip install --no-dependencies --no-build-isolation .
    else
        pip install --target=$INSTALL_PREFIX --no-dependencies --no-build-isolation .
    fi
fi

# Install executable
if [[ $INSTALL_EDITABLE != yes ]]; then
    pushd build
    make install
    popd
fi
