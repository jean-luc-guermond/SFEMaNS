// ---------------------------------------------------------------------
// $Id$
//
//    This file is part of SFEMaNS.
//
//    SFEMaNS is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    SFEMaNS is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with SFEMaNS.  If not, see <http://www.gnu.org/licenses/>.
// ---------------------------------------------------------------------

/**
 * @page Getting_Started Getting started 

<h3>Introduction</h3>
SFEMaNS can solve three different types of problems: 
<ol>
<li> Navier-Stokes equations.</li>
<li> Maxwell equations.</li>
<li> MHD equations (Navier-Stokes and Maxwell equations coupled).</li>
</ol>
The user controls the code through four files: 
\c data, \c condlim.f90, \c main.f90, \c read_user_data.f90.
Template files are available in the TEMPLATES directory.
<ol>
<li> \c data contains all the data that are required by the code to run.</li>
<li> \c main.f90 contains the main and user defined instruction to posprocess the results.</li>
<li> \c condlim.f90 enables the user to control the boundary conditions, the initial data
and the various sourcve terms. </li>
<li> \c read_user_data.f90 is a file made available to the user to read extra data.
</ol>

<h3>\c data file</h3>
The file \c data contains all the data needed for the
computation. The user has to specify the geometry of the domain, the
list of conductive and insulating parts, as well as general
information for the parallel runs.
The list of all the data that can be read is in \c TEMPLATES/list_data_queries.

Here is the full list:
!===================================================================================<BR>
!          Data for the mesh<BR>
!===================================================================================<BR>

===Is mesh file formatted (true/false)? <BR>
Mandatory. Logical (.t. or .f.). <BR>
This tells whether the mesh file format is binary or ascii. 

===Directory and name of mesh file <BR>
Mandatory. Character ('directory' 'meshfile'). <BR>
Names of directory and mesh file. 

===Do we read metis partition? (true/false )<BR>
Optional. Logical (.t. or .f.). <BR>
May be useful to restart on a different computer or if the
restart files have been manipulated and interpolated. Default value is .f.

===Which partitioning: "old" or "new" ? <BR>
Optional. Character ('old' or 'new'). <BR>
Only for specialists. Default value is 'new'. 

!===================================================================================<BR>
!          Processor distribution<BR>
!===================================================================================<BR>
===Number of processors in meridian section<BR>
===Number of processors in Fourier space<BR>

!===================================================================================<BR>
!          Fourier modes<BR>
!===================================================================================<BR>
===Number of Fourier modes<BR>
===Select Fourier modes? (true/false)<BR>
===List of Fourier modes (if select_mode=.TRUE.)<BR>

!===================================================================================<BR>
!          Type of problem to be solved<BR>
!===================================================================================<BR>
===Problem type: (nst, mxw, mhd)<BR>

!===================================================================================<BR>
!          Restarts<BR>
!===================================================================================<BR>
===Restart on velocity (true/false)<BR>
===Restart on magnetic field (true/false)<BR>

!===================================================================================<BR>
!          Time integration<BR>
!===================================================================================<BR>
===Time step and number of time iterations<BR>

!===================================================================================<BR>
!          Check numerical stability<BR>
!===================================================================================<BR>
===Check numerical stability (true/false)<BR>

!===================================================================================<BR>
!          Spacial Periodicity<BR>
!===================================================================================<BR>
===How many pieces of periodic boundary?<BR>
===Indices of periodic boundaries and corresponding vectors<BR>

!===================================================================================<BR>
!          Mesh symmetry<BR>
!===================================================================================<BR>
===Is the mesh symmetric (true/false)?<BR>

!===================================================================================<BR>
!          Navier Stokes data<BR>
!===================================================================================<BR>
===Solve Navier-Stokes with u (true) or m (false)?<BR>
===Number of subdomains in Navier-Stokes mesh<BR>
===List of subdomains for Navier-Stokes mesh<BR>

!===================================================================================<BR>
!          Dirichlet BCs for velocity<BR>
!===================================================================================<BR>
===How many boundary pieces for Dirichlet BCs on uradial?<BR>
===List of boundary pieces for Dirichlet BCs on uradial<BR>
===How many boundary pieces for Dirichlet BCs on utheta?<BR>
===List of boundary pieces for Dirichlet BCs on utheta<BR>
===How many boundary pieces for Dirichlet BCs on uzaxis?<BR>
===List of boundary pieces for Dirichlet BCs on uzaxis<BR>

!===================================================================================<BR>
!          Dirichlet BCs for pressure<BR>
!===================================================================================<BR>
===How many boundary pieces for Dirichlet BCs on pressure?<BR>
===List of boundary pieces for Dirichlet BCs on pressure<BR>

!===================================================================================<BR>
!          Reynolds number<BR>
!===================================================================================<BR>
===Reynolds number<BR>

!===================================================================================<BR>
!          Precession<BR>
!===================================================================================<BR>
===Is there a precession term (true/false)?<BR>
===Precession rate<BR>
===Precession angle over pi<BR>

!===================================================================================<BR>
!          Solver parameters for velocity<BR>
!===================================================================================<BR>
===Maximum number of iterations for velocity solver<BR>
===Relative tolerance for velocity solver<BR>
===Absolute tolerance for velocity solver<BR>
===Velocity solver verbose? (true/false)<BR>
===Solver type for velocity (FGMRES, CG, ...)<BR>
===Preconditionner type for velocity solver (HYPRE, JACOBI, MUMPS...)<BR>

!===================================================================================<BR>
!          Solver parameters for pressure<BR>
!===================================================================================<BR>
===Maximum number of iterations for pressure solver<BR>
===Relative tolerance for pressure solver<BR>
===Absolute tolerance for pressure solver<BR>
===Pressure solver verbose? (true/false)<BR>
===Solver type for pressure (FGMRES, CG, ...)<BR>
===Preconditionner type for pressure solver (HYPRE, JACOBI, MUMPS...)<BR>

!===================================================================================<BR>
!          Solver parameters for mass matrix<BR>
!===================================================================================<BR>
===Maximum number of iterations for mass matrix solver<BR>
===Relative tolerance for mass matrix solver<BR>
===Absolute tolerance for mass matrix solver<BR>
===Mass matrix solver verbose? (true/false)<BR>
===Solver type for mass matrix (FGMRES, CG, ...)<BR>
===Preconditionner type for mass matrix solver (HYPRE, JACOBI, MUMPS...)<BR>

!===================================================================================<BR>
!          LES coefficients<BR>
!===================================================================================<BR>
===Use LES? (true/false)<BR>
===Coefficients for LES<BR>

!===================================================================================<BR>
!          Maxwell data. Definition of the domain<BR>
!===================================================================================<BR>
===Solve Maxwell with H (true) or B (false)?<BR>
===Number of subdomains in magnetic field (H) mesh<BR>
===List of subdomains for magnetic field (H) mesh<BR>

!===================================================================================<BR>
!          Interfaces H/H<BR>
!===================================================================================<BR>
===Number of interfaces in H mesh<BR>
===List of interfaces in H mesh<BR>

!===================================================================================<BR>
!          Dirichlet BCs for H (Hxn = )<BR>
!===================================================================================<BR>
===Number of Dirichlet sides for Hxn<BR>
===List of Dirichlet sides for Hxn<BR>

!===================================================================================<BR>
!          Permeability and conductivity for H<BR>
!===================================================================================<BR>
===Is permeability defined analytically (true/false)?<BR>
===Permeability in the conductive part (1:nb_dom_H)<BR>
===Is permeability variable in theta (true/false)?<BR>
===Conductivity in the conductive part (1:nb_dom_H)<BR>

!===================================================================================<BR>
!          Finite element type for H (P1 or P2)<BR>
!===================================================================================<BR>
===Type of finite element for magnetic field<BR>

!===================================================================================<BR>
!          Magnetic Reynolds number<BR>
!===================================================================================<BR>
===Magnetic Reynolds number<BR>

!===================================================================================<BR>
!          Stabilization parameters<BR>
!===================================================================================<BR>
===Stabilization coefficient (divergence)<BR>
===Stabilization coefficient for Dirichlet H and/or interface H/H<BR>

!===================================================================================<BR>
!          Subdomains for phi (if any)<BR>
!===================================================================================<BR>
===Number of subdomains in magnetic potential (phi) mesh<BR>

!===================================================================================<BR>
!          List of subdomains for phi<BR>
!===================================================================================<BR>
===List of subdomains for magnetic potential (phi) mesh<BR>

!===================================================================================<BR>
!          Dirichlet BCs on phi<BR>
!===================================================================================<BR>
===How many boundary pieces for Dirichlet BCs on phi?<BR>
===List of boundary pieces for Dirichlet BCs on phi<BR>

!===================================================================================<BR>
!          H/phi interface<BR>
!===================================================================================<BR>
===Number of interfaces between H and phi<BR>
===List of interfaces between H and phi<BR>

!===================================================================================<BR>
!          Permeability in vacuum<BR>
!===================================================================================
===Permeability in vacuum<BR>

!===================================================================================<BR>
!          Finite element type for phi (scalar potential, use P2)<BR>
!===================================================================================<BR>
===Type of finite element for scalar potential<BR>

!===================================================================================<BR>
!          Stabilization parameters for interface H/phi<BR>
!===================================================================================<BR>
===Stabilization coefficient (interface H/phi)<BR>

!===================================================================================<BR>
!          Solver parameters for phi<BR>
!===================================================================================<BR>
===Maximum number of iterations for Maxwell solver<BR>
===Relative tolerance for Maxwell solver<BR>
===Absolute tolerance for Maxwell solver<BR>
===Maxwell solver verbose? (true/false)<BR>
===Solver type for Maxwell (FGMRES, CG, ...)<BR>
===Preconditionner type for Maxwell solver (HYPRE, JACOBI, MUMPS...)<BR>

!===================================================================================<BR>
!          Data phase and temperature<BR>
!===================================================================================<BR>
===Is there a temperature field?<BR>

!===================================================================================<BR>
!          Gravity coefficient for temperature<BR>
!===================================================================================<BR>
===Non-dimensional gravity coefficient<BR>

!===================================================================================<BR>
!          Temperature diffusivity<BR>
!===================================================================================<BR>
===Diffusivity coefficient for temperature<BR>

!===================================================================================<BR>
!          Dirichlet BCs on temperature<BR>
!===================================================================================<BR>
===How many boundary pieces for Dirichlet BCs on temperature?<BR>
===List of boundary pieces for Dirichlet BCs on temperature<BR>

!===================================================================================<BR>
!          Solver parameters for Temperature<BR>
!===================================================================================<BR>
===Maximum number of iterations for temperature solver<BR>
===Relative tolerance for temperature solver<BR>
===Absolute tolerance for temperature solver<BR>
===Temperature solver verbose? (true/false)<BR>
===Solver type for temperature (FGMRES, CG, ...)<BR>
===Preconditionner type for temperature solver (HYPRE, JACOBI, MUMPS...)<BR>

!===================================================================================<BR>
!          Level set<BR>
!===================================================================================<BR>
===Is there a level set<BR>
===Cmax coefficient for level set<BR>
===Compression factor for level set<BR>
===Density of fluid 0 and fluid 1<BR>
===Dynamic viscosity of fluid 0 and fluid 1<BR>
===How many boundary pieces for Dirichlet BCs on level set?<BR>
===List of boundary pieces for Dirichlet BCs on level set<BR>

!===================================================================================<BR>
!          Solver parameters for level_set<BR>
!===================================================================================<BR>
===Maximum number of iterations for level set solver<BR>
===Relative tolerance for level set solver<BR>
===Absolute tolerance for level set solver<BR>
===Level set solver verbose? (true/false)<BR>
===Solver type for level set (FGMRES, CG, ...)<BR>
===Preconditionner type for level set solver (HYPRE, JACOBI, MUMPS...)<BR>

!===================================================================================<BR>
!          Data for arpack (eigenvalue problems)<BR>
!===================================================================================<BR>
===Do we use Arpack?<BR>
===Number of eigenvalues to compute<BR>
===Maximum number of Arpack iterations<BR>
===Tolerance for Arpack<BR>
===Which eigenvalues (''LM'', ''SM'', ''SR'', ''LR'' ''LI'', ''SI'')<BR>
===Create 2D vtu files for Arpack? (true/false)<BR>

!===================================================================================<BR>
!          Data for post processing<BR>
!===================================================================================<BR>
===Number of planes in real space for Visualization<BR>
===Frequency to write restart file<BR>
===Frequency to write energies<BR>
===Frequency to create plots<BR>
===Just postprocessing without computing? (true/false)<BR>

!===================================================================================
!          Modes to be zeroed out
!===================================================================================
===Should some modes be zeroed out?<BR>
===How many Navier-Stokes modes to zero out?<BR>
===List of Navier-Stokes modes to zero out?<BR>
===How Maxwell modes to zero out?<BR>
===List of Maxwell modes to zero out?<BR>

!===================================================================================<BR>
!          Verbose<BR>
!===================================================================================<BR>
===Verbose timing? (true/false)<BR>
===Verbose divergence? (true/false)<BR>
===Verbose CFL? (true/false)<BR>

!===================================================================================<BR>
!          Stress-free boundary condition for velocity<BR>
!===================================================================================<BR>
===Stress boundary conditions? (true/false)<BR>
===stab_bdy_ns<BR>

 */
