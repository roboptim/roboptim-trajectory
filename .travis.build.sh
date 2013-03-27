#!/bin/sh
set -e

# Directories.
root_dir=`pwd`
build_dir="$root_dir/_travis/build"
install_dir="$root_dir/_travis/install"
core_dir="$build_dir/roboptim-core"

# Shortcuts.
git_clone="git clone --quiet --recursive"

# Create layout.
rm -rf "$build_dir" "$install_dir"
mkdir -p "$build_dir"
mkdir -p "$install_dir"

# Setup environment variables.
export LD_LIBRARY_PATH="$install_dir/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$install_dir/lib/roboptim-core:$LD_LIBRARY_PATH"
export PKG_CONFIG_PATH="$install_dir/lib/pkgconfig:$PKG_CONFIG_PATH"

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

# Checkout roboptim-core
echo "Installing dependencies..."

cd "$build_dir"
$git_clone "git://github.com/roboptim/roboptim-core.git"
cd "$core_dir"
cmake . -DCMAKE_INSTALL_PREFIX:STRING="$install_dir"
make install

cd "$build_dir"
$git_clone "git://github.com/roboptim/roboptim-core-plugin-ipopt.git"
cd "$core_dir"
cmake . -DCMAKE_INSTALL_PREFIX:STRING="$install_dir"
make install

# Build package
echo "Building package..."
cd "$build_dir"
cmake "$root_dir" -DCMAKE_INSTALL_PREFIX="$install_dir"
make
make install
#make test
