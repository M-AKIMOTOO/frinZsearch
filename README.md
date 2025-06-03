# frinZsearch

# CLI11 is a command line parser for C++11
# https://github.com/CLIUtils/CLI11?tab=readme-ov-file 
# Ubuntu 24.04 case
cd CLI11
mkdir build 
cmake .. # If you remove files in the build.
make 
sudo make install

# If g++ version is GCC 8.x, you must add -lstdc++fs to the last g++ compile.

# You should check the library path of your GCC compiler.
g++ -print-search-dirs


# frinZsearch install method, not "make install"
make 
