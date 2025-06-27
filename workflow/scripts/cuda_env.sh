#!/bin/bash
# cuda_env_setup.sh - Contains all CUDA environment variables

# Set CUDA path variables
export CUDA_PATH=$(dirname $(dirname $(which nvcc)))
export CUDA_HOME=$CUDA_PATH
export PATH=$CUDA_PATH/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_PATH/lib64:$LD_LIBRARY_PATH
export TORCH_CUDA_ARCH_LIST="7.0"

# Set compiler variables
export CC=$(which gcc)
export CXX=$(which g++)
export CUDAHOSTCXX=$(which g++)