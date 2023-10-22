Files:
 - Scripts: Chebyshev_Galerkin_method.m, Fundamental_matrix_method.m, Direct_numerical_simulation.m
 - Classes: Chebyshev_basis.m, Fundamental_matrix.m, Ground_state.m
 - Functions: index.m, transform.m, plot_dns.m
 - Data: Chebyshev_integrals.mat, Gaussian_quadrature.mat

This Readme file explains the main features of the eigenvalue problem solvers and direct numerical solver of the problems described in the paper by Leon Kloker and Carina Bringedal titled "Solution approaches for evaporation-driven density instabilities in a slab of saturated
porous media". The equation and figure references below refer to the equations and figures in that paper.

The script Chebyshev_Galerkin_method can be used to solve the eigenvalue problem (40) of the perturbation equations via the Chebyshev-Galerkin method. In this script, the dimensionless height alpha, the time of the ground state t and the wavenumber a of the perturbation can be adjusted and the corresponding value Ra(alpha, t, a) will be computed.

The script Fundamental_matrix_method can be used to solve the eigenvalue problem (40) of the perturbation equations via the fundamental matrix method. In this script, the dimensionless height alpha, the time of the ground state t and the wavenumber a of the perturbation can be adjusted and the corresponding value Ra(alpha, t, a) will be computed. Moreover, the grid of Rayleigh numbers that is used to search the eigenvalue Ra(alpha, t, a) should be adjusted.

The script Direct_numerical_simulation can be used to run a simulation of the full equation system consisting of equation (9),(10),(11) and boundary conditions (13). Here, the dimensionless height alpha, the perturbation wavenumber a, the Rayleigh number Ra, the amount of finite volume cells in the vertical and lateral direction, the end time of the simulation and the timestep can be adjusted. Moreover, the kind of perturbation seed and its amplitude are changeable.

A object of the class Ground_state can be created in order to get the ground-state salt concentration as defined in equation (27) at arbitrary heights z between 0 and alpha and times t between 0 and infinity using the member function get_solution.

The function plot_dns can be used to plot the state vector calculated by the direct numerical simulation similar as done in figures 10, 11 and 12.

The class, function and data files should not be changed. 