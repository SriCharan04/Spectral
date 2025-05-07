# Spectral

## Compilation Instructions
### Dependencies required
- Eigen3 (linear algebra)
- HDF5 (data storage)
- ParaView (post-processing)

- If any library is not available in your own system, you can install using `sudo apt install libeigen3-dev libhdf5-dev` for Ubuntu.

### g++ version 17 is necessary for the compilation of this program. Plese install it using `sudo apt install g++` in Ubuntu.
g++ -std=c++17 -O3 -I/usr/include/eigen3 -lhdf5_cpp -lhdf5 -lpthread spectral.cpp -o spectral

### Note: To have high precision in your calculations, make sure to add the flags -mlong-double-128 and -lquadmath. Also, make sure to have the Quadmath library available.
Executing the spectral executive will generate two additional files:
- `spectral_solution.h5` (HDF5 data)
- `solution.xmf` (XDMF metadata)

The rough view of the arrays can be done from command line also, using the commands 
`h5ls -d spectral_solution.h5/error`
`h5ls -d spectral_solution.h5/solution`
`h5ls -d spectral_solution.h5/exact`
## ParaView Visualization

### Basic Workflow
1. **Open ParaView**
2. **File > Open** -> Select `solution.xmf`
3. In Pipeline Browser:
   - Click 👁️ icon next to "ChebyshevGrid"
4. In Dropdown mentu in top left, choose Solution, Error or Exact to view the arrays of your choice
