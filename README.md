# CFD-solver-MATLAB
* A 2D Navier-Stokes solver for solving laminar, incompressible flows using finite-volume method and collocated grid arrangement coded in MATLAB
* Capable of solving both steady state and unsteady problems
* Pressure-velocity coupling implemented using SIMPLE algorithm
* Spatial discretization for divergence schemes - Available choices include Upwind, Central differencing, Second order upwind, QUICK and FROMM schemes
* Temporal discretization for unsteady simulation - Implicit Crank-Nicholson
* Cell centred gradient algorithms : Available choices are Gauss cell-based, Gauss node-based and Least squares gradient scheme
* Matrix solvers available : Gauss Siedel, Gauss Jacobi and Incomplete LU decomposition (Fell free to edit the code to implement MATLAB built-in solvers)
* Accepts all-tri mesh as well as all-quad mesh in 2D ASCII Ansys-Fluent mesh file format (.msh)
* You have the option to output the files in Tecplot binary file format

## Instructions:

* Run the file NS_solve.m file to run the solver.
* Some example mesh files and their boundary condition files are provided.
* Set the boundary conditions using the files 'U.bc', 'V.bc', 'P.bc' inside the folder named BC. Check the example boundary condition files. Currently supports fixed value, wall and zero gradient boundary condition
* You can edit the solver setup parameters using the file 'solver_setup' inside the solver_setup folder
* Schemes and SIMPLE algorithm parameters can alo be edited using the files inside the solver_setup folder
  
* Once the simulation is done, the results are available in the variables u,v,p. You have the option to visualize them in Tecplot from the generated Tecplot binary files. Otherwise you can use MATLAB to visualize the flow field using the following function (type them in the MATLAB console):
  * plot_mesh() - To visualize the mesh
  * plot_contour(var) - To visualize the filled contours of a parameter. 'var' can be u,v,p 
  * plot_quiver(u,v) - To visualize the flow field using a quiver plot
* Feel free to use your own post-processing functions to visualize data

## Acknowledgement
* Thanks to Wen Long (https://www.cfd-online.com/Forums/tecplot/103860-mat2tecplot.html) for the MATLAB function mat2tecplot.m that is used to convert data to Tecplot binary format.

## References:
* Error Analysis and Estimation for the Finite Volume Method with Applications to Fluid Flows - PhD thesis by Dr. Hrvoje Jasak (https://foam-extend.fsb.hr/wp-content/uploads/2016/12/Jasak_PhD_1996.pdf)

* The Finite Volume Method in Computational Fluid Dynamics - F. Moukalled, L. Mangani, M. Darwish (https://link.springer.com/book/10.1007/978-3-319-16874-6)

* Computational Methods for Fluid Dynamics - Joel H. Ferziger, Milovan Perić (https://link.springer.com/book/10.1007/978-3-642-56026-2)
