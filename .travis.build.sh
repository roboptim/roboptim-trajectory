#!/bin/sh
set -e

# Directories.
root_dir=`pwd`
build_dir="$root_dir/_travis/build"
install_dir="$root_dir/_travis/install"
core_dir="$build_dir/roboptim-core"
plugin_dir="$build_dir/roboptim-core-plugin-ipopt"

# Shortcuts.
git_clone="git clone --quiet --recursive"

# Create layout.
rm -rf "$build_dir" "$install_dir"
mkdir -p "$build_dir"
mkdir -p "$install_dir"

# Setup environment variables.
export LD_LIBRARY_PATH="$install_dir/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$install_dir/lib/roboptim-core:$install_dir/lib/x86_64-linux-gnu/roboptim-core:$LD_LIBRARY_PATH"
export PKG_CONFIG_PATH="$install_dir/lib/pkgconfig:$install_dir/lib/x86_64-linux-gnu/pkgconfig:$PKG_CONFIG_PATH"

# Checkout Eigen.
cd "$build_dir"
wget "http://bitbucket.org/eigen/eigen/get/3.1.2.tar.gz"
tar xzvf 3.1.2.tar.gz
cd "$build_dir/eigen-eigen-5097c01bcdc4/"
mkdir -p "$build_dir/eigen-eigen-5097c01bcdc4/_build"
cd "$build_dir/eigen-eigen-5097c01bcdc4/_build"
cmake .. -DCMAKE_INSTALL_PREFIX:STRING="$install_dir" \
          -Dpkg_config_libdir:STRING="$install_dir/lib"
make
make install

# Checkout Ipopt.
cd "$build_dir"
wget "http://www.coin-or.org/download/source/Ipopt/Ipopt-3.10.3.tgz"
tar xzvf Ipopt-3.10.3.tgz
cd "$build_dir/Ipopt-3.10.3"
cd ThirdParty/Mumps
./get.Mumps
cd "$build_dir/Ipopt-3.10.3"
# Force GCC for Ipopt, clang will not work.
CC=gcc CXX=g++ ./configure --prefix="$install_dir"
make
make install

echo "Installing dependencies..."

# Checkout roboptim-core.
cd "$build_dir"
$git_clone "git://github.com/roboptim/roboptim-core.git"
cd "$core_dir"
cmake . -DCMAKE_INSTALL_PREFIX:STRING="$install_dir"
make install

# Checkout roboptim-core-plugin-ipopt.
cd "$build_dir"
$git_clone "git://github.com/roboptim/roboptim-core-plugin-ipopt.git"
cd "$plugin_dir"
cmake . -DCMAKE_INSTALL_PREFIX:STRING="$install_dir"
make install

# Build package
echo "Building package..."
cd "$build_dir"
cmake "$root_dir" -DTESTSUITE_SOLVER=ipopt -DCMAKE_INSTALL_PREFIX="$install_dir"
make
make install
# Print error logs when tests fail
export CTEST_OUTPUT_ON_FAILURE=1
make test
