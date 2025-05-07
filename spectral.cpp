/*
 * Spectral Method Solver for PDEs using Chebyshev Collocation
 * 
 * Implements:
 * - Chebyshev-Gauss-Lobatto grid generation
 * - Spectral differentiation matrices
 * - Direct solution of Helmholtz-type equation
 * - HDF5/XDMF output for ParaView visualization
 * 
 * References:
 * [1] Trefethen, L.N. (2000). Spectral Methods in MATLAB. SIAM.
 * [2] Boyd, J.P. (2001). Chebyshev and Fourier Spectral Methods. Dover.
 * [3] Fornberg, B. (1987). The Pseudospectral Method: Comparisons with Finite Differences
 *     for the Elastic Wave Equation. Geophysics, 52(4), 483-501.
 * [4] HDF5/XDMF Standards: https://www.hdfgroup.org & https://xdmf.org
 */

 #include <iostream>
 #include <fstream>
 #include <vector>
 #include <cmath>
 #include <eigen3/Eigen/Dense>
 #include <H5Cpp.h>
 
 using namespace Eigen;
 typedef long double Real;
 const int N = 32;  // Number of collocation points (Trefethen Chap. 6)
 const Real PI = 3.14159265358979323846264338327950288419716939937510L;
 
 /// Generates Chebyshev-Gauss-Lobatto nodes (Trefethen pg. 53)
 class ChebyshevGrid {
 public:
     std::vector<Real> nodes;
     
     ChebyshevGrid(int n) {
         nodes.resize(n+1);
         for(int j = 0; j <= n; ++j)
             nodes[j] = -cos(PI * j / n);  // CGL node formula
     }
 };
 
 /// Computes 1st and 2nd derivative matrices using spectral collocation
 class ChebyshevDifferentiation {
 public:
     Matrix<Real, Dynamic, Dynamic> D1, D2;
     
     ChebyshevDifferentiation(const ChebyshevGrid& grid) {
         computeD1(grid);  // First derivative matrix
         computeD2(grid);  // Second derivative matrix (Fornberg 1987)
     }
 
 private:
     /// Constructs 1st derivative matrix (Trefethen Program 11)
     void computeD1(const ChebyshevGrid& grid) {
         int n = grid.nodes.size()-1;
         const auto& x = grid.nodes;
         D1 = Matrix<Real, Dynamic, Dynamic>::Zero(n+1, n+1);
         
         for(int i = 0; i <= n; ++i) {
             for(int j = 0; j <= n; ++j) {
                 if(i != j) {
                     Real ci = (i == 0 || i == n) ? 2.0L : 1.0L;  // Boundary weights
                     Real cj = (j == 0 || j == n) ? 2.0L : 1.0L;
                     Real sign = ((i+j) % 2 == 0) ? 1.0L : -1.0L;  // Alternating signs
                     Real diff = x[i] - x[j];
                     
                     if(std::abs(diff) > 1e-12L) {
                         D1(i,j) = (ci/cj) * sign / diff;  // Core differentiation formula
                     }
                 }
             }
         }
         
         // Diagonal entries (Boyd Eqn. 4.19)
         D1(0,0) = (2L*n*n + 1L)/6.0L;
         D1(n,n) = -D1(0,0);
         
         // Interior diagonal terms (Trefethen Lemma 4.1)
         for(int i = 1; i < n; ++i) {
             Real xi_sq = x[i]*x[i];
             Real denom = 1L - xi_sq;
             if(std::abs(denom) > 1e-12L) {
                 D1(i,i) = -x[i]/(2L * denom);
             }
         }
     }
     
     /// Constructs 2nd derivative matrix via negative sum trick (Fornberg 1987)
     void computeD2(const ChebyshevGrid& grid) {
         int n = grid.nodes.size()-1;
         const auto& x = grid.nodes;
         D2 = Matrix<Real, Dynamic, Dynamic>::Zero(n+1, n+1);
         
         for(int i = 0; i <= n; ++i) {
             Real sum = 0;
             for(int j = 0; j <= n; ++j) {
                 if(i != j && std::abs(x[i]-x[j]) > 1e-12L) {
                     D2(i,j) = D1(i,j) * (D1(i,i) - 1L/(x[i]-x[j]));
                     sum += D2(i,j);
                 }
             }
             D2(i,i) = -sum;  // Ensure zero row sums for conservation
         }
     }
 };
 
 /// Solves Helmholtz equation: u'' + k²u = f with Dirichlet BCs
 class SpectralSolver {
 public:
     VectorXd solution, exactSolution, error;
     
     /// Implements tau-method for boundary conditions (Boyd Sec. 6.3)
     void solve(const ChebyshevGrid& grid, const ChebyshevDifferentiation& diff) {
         int n = grid.nodes.size();
         int interior = n - 2;  // Exclude boundary points
         
         // Construct interior operator (Trefethen Chap. 7)
         MatrixXd A = diff.D2.block(1,1,interior,interior).cast<double>() +
                      4.0*PI*PI*MatrixXd::Identity(interior,interior);
         
         // Source term for manufactured solution
         VectorXd b(interior);
         for(int i=0; i<interior; ++i) {
             double x = static_cast<double>(grid.nodes[i+1]);
             b[i] = (4.0*PI*PI*(1.0-x*x) - 2.0)*sin(PI*x) - 4.0*PI*x*cos(PI*x);
         }
 
         // Solve linear system (Trefethen Program 13)
         VectorXd u_interior = A.fullPivLu().solve(b);
 
         // Analytic solution: u(x) = (1-x²)sin(πx)
         exactSolution = VectorXd::Zero(n);
         for(int i=0; i<n; ++i) {
             double x = static_cast<double>(grid.nodes[i]);
             exactSolution[i] = (1.0-x*x)*sin(PI*x);
         }
 
         // Reconstruct full solution with boundary conditions
         solution = VectorXd::Zero(n);
         solution.segment(1, interior) = u_interior;
         error = (solution - exactSolution).cwiseAbs();
     }
     
     /// Writes results to HDF5 format (XDMF standard)
     void write_hdf5(const std::vector<Real>& nodes, const std::string& filename) {
         H5::H5File file(filename, H5F_ACC_TRUNC);
         
         // Store 3D coordinates for ParaView compatibility
         std::vector<double> dnodes(nodes.begin(), nodes.end());
         hsize_t dims[2] = {dnodes.size(), 3};
         std::vector<double> coords(dnodes.size() * 3, 0.0);
         for(size_t i=0; i<dnodes.size(); ++i) {
             coords[3*i] = dnodes[i];  // X-coordinate
         }
         
         H5::DataSet coord_ds = file.createDataSet("coordinates",
             H5::PredType::NATIVE_DOUBLE,
             H5::DataSpace(2, dims));
         coord_ds.write(coords.data(), H5::PredType::NATIVE_DOUBLE);
         
         // Lambda for writing 1D datasets
         auto write_1d = [&](const std::string& name, const VectorXd& data) {
             hsize_t dim = data.size();
             H5::DataSet ds = file.createDataSet(name,
                 H5::PredType::NATIVE_DOUBLE,
                 H5::DataSpace(1, &dim));
             ds.write(data.data(), H5::PredType::NATIVE_DOUBLE);
         };
         
         write_1d("solution", solution);
         write_1d("exact", exactSolution);
         write_1d("error", error);
     }
     
     /// Generates XDMF metadata file for HDF5 visualization
     void write_xdmf(const std::string& h5_filename) {
         std::ofstream xdmf("solution.xmf");
         xdmf << R"(<?xml version="1.0"?>
 <Xdmf Version="3.0">
   <Domain>
     <Grid Name="ChebyshevGrid" GridType="Uniform">
       <Topology TopologyType="Polyline" Dimensions=")" 
           << solution.size() << R"("/>
       <Geometry GeometryType="XYZ">
         <DataItem Dimensions=")"
           << solution.size() << " 3\" NumberType=\"Float\" Format=\"HDF\">\n"
           << "          " << h5_filename << ":/coordinates\n"
           << R"(        </DataItem>
       </Geometry>
       <Attribute Name="Solution" AttributeType="Scalar" Center="Node">
         <DataItem Dimensions=")"
           << solution.size() << R"(" NumberType="Float" Format="HDF">
           )" << h5_filename << R"(:/solution
         </DataItem>
       </Attribute>
       <Attribute Name="Exact" AttributeType="Scalar" Center="Node">
         <DataItem Dimensions=")"
           << exactSolution.size() << R"(" NumberType="Float" Format="HDF">
           )" << h5_filename << R"(:/exact
         </DataItem>
       </Attribute>
       <Attribute Name="Error" AttributeType="Scalar" Center="Node">
         <DataItem Dimensions=")"
           << error.size() << R"(" NumberType="Float" Format="HDF">
           )" << h5_filename << R"(:/error
         </DataItem>
       </Attribute>
     </Grid>
   </Domain>
 </Xdmf>)";
     }
 };
 
 int main() {
     try {
         ChebyshevGrid grid(N);  // Initialize spatial discretization
         ChebyshevDifferentiation diff(grid);  // Compute differentiation operators
         SpectralSolver solver;  // Create PDE solver
         solver.solve(grid, diff);  // Solve boundary value problem
         
         // Write results for visualization
         const std::string h5_filename = "spectral_solution.h5";
         solver.write_hdf5(grid.nodes, h5_filename);
         solver.write_xdmf(h5_filename);
         
         std::cout << "Maximum error: " << solver.error.maxCoeff() << std::endl;
         std::cout << "HDF5/XDMF files written successfully" << std::endl;
         
     } catch (const H5::Exception& error) {
         std::cerr << "HDF5 Error: " << error.getCDetailMsg() << std::endl;
         return 1;
     } catch (const std::exception& e) {
         std::cerr << "Error: " << e.what() << std::endl;
         return 1;
     }
     return 0;
 }