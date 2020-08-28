# CFD-solver-MATLAB
* A 2D Navier-Stokes solver for unsteady,laminar, incompressible flows using finite-volume method and collocated grid arrangement coded in MATLAB
* Pressure-velocity coupling implemented using SIMPLE algorithm.
* Spatial discretization - Central differencing
* Temporal discretization - Implicit Crank-Nicholson
* Accepts all-tri mesh as well as all-quad mesh in 2D ASCII Ansys-Fluent mesh file format (.msh).

## Instructions:

* Run the file NS_solve.m file to run the solver.
* Some example mesh files and their boundary condition files are provided. Try them out with the appropriate flow parameters.
* Set the boundary conditions using the files 'U.bc', 'V.bc', 'P.bc' inside the folder named BC. Check the example boundary condition files. Currently supports fixed value and zero gradient boundary condition
* Don't forget to change the below mentioned parameters according to each specififc case. Specifically change the Reynolds number (Re) and Length scale (L)

* Parameters to set in the code:
  1. mesh_file : Mesh file location
  2. Re : Reynolds number based on the parameter 'L'. Used to compute U_ref
  3. mu : Fluid dynamic viscosity 
  4. rho : Fluid density
  5. L  : Geometric length scale. Usually diameter of the cylinder, airfoil chord length,etc
  6. CFL_max : Maximum CFL number. Time step size varies based on this
  7. alpha_p : Under relaxation factor for pressure (usually 0.1-0.3)
  8. alpha_u : Under relaxation factor for momentum eq (usually 0.3-0.7)
  9. T : Flow end time
  10. flowvis : Flag set to visualize flow field every time step as a quiver plot. Set it to 'T' to activate, 'F' to deactivate
  11. residuals : flag set to visualize residuals every time step
  
* Once the unsteady simulation is done, the results are available in the variables u,v,p at the cell centres x,y. Feel free to visualize the data however you choose

## References:
* Error Analysis and Estimation for the Finite Volume Method with Applications to Fluid Flows - PhD thesis by Dr. Hrvoje Jasak (https://foam-extend.fsb.hr/wp-content/uploads/2016/12/Jasak_PhD_1996.pdf)

* The Finite Volume Method in Computational Fluid Dynamics - F. Moukalled, L. Mangani, M. Darwish (https://link.springer.com/book/10.1007/978-3-319-16874-6)

* Computational Methods for Fluid Dynamics - Joel H. Ferziger, Milovan Perić (https://link.springer.com/book/10.1007/978-3-642-56026-2)
