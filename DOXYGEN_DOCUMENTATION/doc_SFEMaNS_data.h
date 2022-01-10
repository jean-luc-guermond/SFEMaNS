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
 * @page doc_SFEMaNS_data Data file

 
The <tt>data</tt> file contains information on the problem approximated
 with the code SFEMaNS. These informations are stocked in a derived data
 type called inputs. For example, the kinetic Reynolds number is saved
 in inputs\%Re. It is possible to add variable, not already included in
 SFEMaNS, with the file <tt>read_user_data.f90</tt>. Such new variables
 are stocked in the derived data type called user. We refer to the
 section \ref doc_SFEMaNS_read_user_data for more information.
The information read by the code SFEMaNS can be
 splitted in the following groups.
<ol>
<li>General settings (mesh, number of Fourier modes, equations approximated, etc.).
<li>Setting for the Navier-Stokes equations.
<li>Setting for the level set equation.
<li>Setting for the temperature equation.
<li>Setting for the Maxwell equations.
<li>Setting for eigenvalue problems with Arpack.
<li>Setting for the outputs computations.
<li>Setting to force Fourier components of the velocity field
 and the magnetic field to zero.
</ol>

The following described the datas associated to each of the
 above groups. In particular, we show how to set these datas
 and specify if it is required to set them. Indeed, some
 variables have default value while other need to be set
 in the data for any problem considered. We note that many
 examples of <tt>data</tt> files are described in the sections
 \ref doc_debug_test
 and 
 \ref doc_example_physical_pb.


Remarks:
<ol>
<li>The set of each data is done in the line under its
 associated key line. These key lines always start with "===".
<li>The datas associated to an equation are not required if this equation
 is not solved.
</ol>


@section data_generat_setting General settings
<ol>
<li>The format of the mesh (formatted or unformatted) is set as follows.
\code
===Is mesh file formatted (true/false)?
\endcode
<li>The directory and the name of the mesh file is set as follows.
\code
===Directory and name of mesh file
\endcode
<li>It is possible to point out that the mesh is symmetric with the following option.
\code
===Is the mesh symmetric (true/false)?
\endcode
This parameter is set to false by default.
<li>The number of processors in the meridian section has to be set.
\code
===Number of processors in meridian section
\endcode
The finite element meridian section is then splitted in as many subsections.
<li>The number of processors in the Fourier space has to be set.
\code
===Number of processors in Fourier space
\endcode
<li>The number of Fourier modes used to approximate the problem is set as follows.
\code
===Number of Fourier modes
\endcode
Remark: the number of Fourier mode is a multiple of the
 number of processors in Fourier space. For instance if
 one use two processors in Fourier space and approximate
 four Fourier modes, each processor in Fourier space is
 approximating the equations for two Fourier modes. The
 total number of processors used is equal to the number
 of precessors in meridian section times the number of
 processors in Fourier space.
<li>It is possible to select a given list of Fourier mode as follows.
\code
===Select Fourier modes? (true/false)
\endcode
This parameter is set to false by default.
<li>If the above parameter is true, the list of Fourier
 modes approximated has to be given.
\code
===List of Fourier modes (if select_mode=.TRUE.)
\endcode
<li>Four type of problems can be considered: nst (hydrodynamic),
 mxw (magnetic), mhd (magnetohydrodynamic) and
 fhd (ferrohydrodynamic). The problem type is set as follows.
\code
===Problem type: (nst, mxw, mhd, fhd)
\endcode
The answer is either 'nst', 'mxw' , 'mhd' or 'fhd'.
<li>The user has to specify if a restart is done on the
 velocity field.
\code
===Restart on velocity (true/false)
\endcode
<li>It is also possible to do a restart with data of
 the magnetic field or the temperature.
\code
===Restart on magnetic field (true/false)
\endcode
\code
===Restart on temperature (true/false)
\endcode
These parameters are set to false by default.
<li>It is possible to use the metis partition of a previous computation.
\code
===Do we read metis partition? (true/false)
\endcode
This parameter is set to false by default. If a restart is done,
 set this parameter to true. Indeed, the metis partition
 contains the information on the splitting of
 the finite element section in subsections.
<li>The number of time step and the number of time iterations
 is set as follows.
\code
===Time step and number of time iterations
\endcode
</ol>

Remark: when doing a restart, the name of the restart files used
 should not display information on the time iteration they had been
 generated. It means the part _I001 (_I002, etc.) should be removed of their name.
 See <a href='doc_SFEMaNS_data.html#data_output'><code>below</code></a>
 for more information on the generation of the restart files.


@section data_period Setting for periodic problem (optional)

<ol>
<li>The number of pairs of boundaries with periodic condition is set as follows.
\code
===How many pieces of periodic boundary?
\endcode
This parameter is set to zero by default.
<li>If pieces of periodic boundary are present, the label
 of the boundaries and the vector that lead to the first
 boundary to the second one are given after the following line.
\code
===Indices of periodic boundaries and corresponding vectors
\endcode
We note that the code needs as much as lines as the number
 of pairs of boundaries with periodic condition.
</ol>




@section data_NS Setting for the Navier-Stokes equations

We refer to the section \ref doc_intro_SFEMaNS_possibilities_nst
 for details on the weak formulation of the Navier-Stokes equations.

@subsection data_NS_general General settings

To approximate solutions of the Navier-Stokes equations, set the following parameters.
<ol>
<li>The dependent variable is either the velocity field or the momentum.
\code
===Solve Navier-Stokes with u (true) or m (false)?
\endcode
It is set to true by default. We note that the momentum
 has to be used for multiphase flow problem.
<li>Number of subdomains where the Navier-Stokes equations are approximated.
\code
===Number of subdomains in Navier-Stokes mesh
\endcode
<li>List of the above subdomains.
\code
===List of subdomains for Navier-Stokes mesh
\endcode
<li>Number of boundary pieces with Dirichlet condition.
\code
===How many boundary pieces for full Dirichlet BCs on velocity?
\endcode
<li>List of the boundary pieces with Dirichlet condition.
\code
===List of boundary pieces for full Dirichlet BCs on velocity
\endcode
<li>Number of boundary pieces with homogeneous Neumann condition.
\code
===How many boundary pieces for homogeneous normal velocity?
\endcode
<li>List of the boundary pieces with homogeneous Neumann condition.
\code
===List of boundary pieces for homogeneous normal velocity
\endcode
<li>Penalty coefficient for enforcing homogeneous Neumann boundary velocity.
\code
===stab_bdy_ns
\endcode
It is set to one by default.
<li>Kinetic Reynolds number \f$\Re\f$.
\code
===Reynolds number
\endcode
<li>Penalty coefficient \f$c_\text{div}\f$, see the section
 \ref doc_intro_SFEMaNS_possibilities_nst_1.
\code
===Coefficient for penalty of divergence in NS?
\endcode
It is set to zero by default.
 It can not be used with the momentum as dependent variable.
</ol>

@subsection data_NS_precession Precession

SFEMaNS can consider precession set up in the precession
 frame when the rotation axis is the vertical axis \f$\textbf{e}_z\f$.
 The method consists of adding the term \f$ 2 \epsilon \textbf{e}_p \times \bu\f$
 in the left handside of the Navier-Stokes equations. The parameter \f$\epsilon\f$ is
 the precession rate. The vector \f$\textbf{e}_p\f$ is equal to 
 \f$\sin(\alpha\pi)\cos(\theta) \textbf{e}_r 
 -\sin(\alpha\pi)\sin(\theta)\textbf{e}_\theta+\cos(\alpha\pi)\textbf{e}_z\f$
 with \f$\alpha\f$ the angle of precession, i.e. the angle between the vertical
 axis and the precession axis. To use this feature, set the following parameters.
<ol>
<li>Is the precession term computed?
\code
===Is there a precession term (true/false)?
\endcode
It is set to false by default.
<li>Precession rate \f$\epsilon\f$.
\code
===Precession rate
\endcode
<li>Precession angle \f$\alpha\f$ in degree.
\code
===Precession angle over pi
\endcode
</ol>

@subsection data_NS_penalty_solid Non axisymmetric solid obstacle 

SFEMaNS can consider non axisymmetric geometry via the use of a penalty method,
 see the section \ref doc_intro_SFEMaNS_possibilities_nst_3
 for more details. To use this feature, set the following parameters.
<ol>
<li>Is the penalty method used?
\code
===Use penalty in NS domain (true/false)?
\endcode
It is set to false by default. If set to true,
 the penalty function \f$\chi\f$ is defined in the subroutine
 <tt>penal_in_real_space</tt> of the file <tt>condlim.f90</tt>.
<li>Does the solid has a nonzero velocity?
\code
===Use nonzero velocity in solids (true/false)?
\endcode
It is set to false by default. If set to true, the velocity \f$\bu_\text{obs}\f$
 is defined in the subroutine
<tt>imposed_velocity_by_penalty</tt> of the file <tt>condlim.f90</tt>.
<li>Should the code compute the following quantity
 \f$ \displaystyle -\frac{2}{\pi}\int_{\Omega_{c,f}} \frac{(1-\chi)(\bu-\bu_\text{obs})
\cdot \textbf{e}_z}{\tau} \f$ with \f$\tau\f$ being the time step.
\code
===Compute z momentum (true/false)?
\endcode
It is set to false by default. The output is written in the file <tt>fort.12</tt> each
 inputs\%freq_en time iterations with inputs\%freq_en the frequency 
 to write energy outputs (defined later).
</ol>

@subsection data_NS_LES LES stabilization with entropy viscosity

We refer to the section \ref doc_intro_SFEMaNS_possibilities_nst_2
 for a description of this method.
 To use this stabilization method, set the following parameters.
<ol>
<li>Is the entropy viscosity (LES) used?
\code
===Use LES? (true/false)
\endcode
It is set to false by default.
<li>Coefficient \f$c_\text{e}\f$.
\code
===Coefficient multiplying residual
\endcode
<li>Optional coefficient \f$c_1\f$.
\code
===Coefficient for explicit LES
\endcode
It is set to zero by default.
</ol>

@subsection data_NS_multiphase Multiphase computations

We refer to the section \ref doc_intro_SFEMaNS_possibilities_nst_4
 for more details on this set up.

To use this feature, set the following parameters.
<ol>
<li>Do we use a second order in time algorithm (BDF2)?
\code
===Do we solve momentum with bdf2 (true/false)?
\endcode
It is set to false by default. Theoretical results about the algorithm's
 stability are only known for BDF1.
<li>Is the momentum equation stabilized with the entropy viscosity?
\code
===Use LES in momentum? (true/false)
\endcode
It is set to true by default.
<li>Do we solve the momentum equation using an artificial compression method?
\code
===Solve Navier-Stokes with art comp scheme (true) or (false)?
\endcode
It is set to false by default, meaning that we use the projection method described in @subpage doc_intro_SFEMaNS_possibilities_nst_4
<li>Set the penalty coefficient for the artificial compression method.
\code
===Penalty coefficient for artifical compression
\endcode
It is set to one by default.
</ol>

@subsection data_NS_solver Informations for the Navier-Stokes equations solvers

@subsubsection data_NS_solver_vel Velocity solver

Set the following information on the solver for the velocity field.
<ol>
<li>Maximum of iterations.
\code
===Maximum number of iterations for velocity solver
\endcode
<li>Relative tolerance.
\code
===Relative tolerance for velocity solver
\endcode
<li>Absolute tolerance.
\code
===Absolute tolerance for velocity solver
\endcode
<li>Type of solver.
\code
===Solver type for velocity (FGMRES, CG, ...)
\endcode
It is possible to set it to GMRES (Generalized Minimal Residual method),
 FGMRES (Flexible Generalized Minimal Residual method)
 and CG (Conjugate Gradient method).
<li>Type of preconditionner.
\code
===Preconditionner type for velocity solver (HYPRE, JACOBI, MUMPS...)
\endcode
</ol>

@subsubsection data_NS_solver_pre Pressure solver

Set the following information on the solver for the pressure.
<ol>
<li>Maximum of iterations.
\code
===Maximum number of iterations for pressure solver
\endcode
<li>Relative tolerance.
\code
===Relative tolerance for pressure solver
\endcode
<li>Absolute tolerance.
\code
===Absolute tolerance for pressure solver
\endcode
<li>Type of solver.
\code
===Solver type for pressure (FGMRES, CG, ...)
\endcode
It is possible to set it to GMRES (Generalized Minimal Residual method),
 FGMRES (Flexible Generalized Minimal Residual method)
 and CG (Conjugate Gradient method).
<li>Type of preconditionner.
\code
===Preconditionner type for pressure solver (HYPRE, JACOBI, MUMPS...)
\endcode
</ol>

@subsubsection data_NS_solver_mass_matrix Mass matrix solver

Set the following information on the solver for the mass matrix.
<ol>
<li>Maximum of iterations.
\code
===Maximum number of iterations for mass matrix solver
\endcode
<li>Relative tolerance.
\code
===Relative tolerance for mass matrix solver
\endcode
<li>Absolute tolerance.
\code
===Absolute tolerance for mass matrix solver
\endcode
<li>Type of solver.
\code
===Solver type for mass matrix (FGMRES, CG, ...)
\endcode
It is possible to set it to GMRES (Generalized Minimal Residual method),
 FGMRES (Flexible Generalized Minimal Residual method)
 and CG (Conjugate Gradient method).
<li>Type of preconditionner.
\code
===Preconditionner type for mass matrix solver (HYPRE, JACOBI, MUMPS...)
\endcode
</ol>




@section data_NS_level_set Setting for the level set equations (multiphase problem)

We refer to the section \ref doc_intro_SFEMaNS_possibilities_nst_4
 for a description of such set ups.


@subsection data_level_set_general General settings

To consider multiphase flow problems, set the following parameters.
<ol>
<li>Is a level set used?
\code
===Is there a level set?
\endcode
It is set to false by default. The following information are only read if it is set to true.
<li>What is the number of fluids?
\code
===How many fluids?
\endcode
As the level set functions represent interfaces between fluids, this number is equal
 to the number of level set plus one.
<li>Define a multiplier that multiplies the minimum mesh size.
 The result is saved in inputs\%h_min_distance. It can be
 used to initialize the level set.
\code
===multiplier for h_min for level set
\endcode
<li>Coefficient \f$c_\text{comp}\in[0,1]\f$.
\code
===Compression factor for level set
\endcode
A typical value is one.
<li>Densities of the different fluids.
\code
===Density of fluid 0, fluid 1, ...
\endcode
<li>Dynamical viscosities of the different fluids.
\code
===Dynamic viscosity of fluid 0, fluid 1, ...
\endcode
<li>If needed, set the electrical conductivity of the different fluids.
\code
===Conductivity of fluid 0, fluid 1, ...
\endcode
<li>Is there surface tension effect?
\code
===Is there a surface tension?
\endcode
The default value is false.
<li>Coefficients of surface tension of each level set.
 These parameters are the inverse of the Weber number.
\code
===Coefficients of surface tension for level set 0, level set 1, ...
\endcode
<li>Is a mass correction applied?
\code
===Do we apply mass correction? (true/false)
\endcode
It is set to true by default.
<li>Number of boundary pieces with Dirichlet condition.
\code
===How many boundary pieces for Dirichlet BCs on level set?
\endcode
<li>List of the boundary pieces with Dirichlet boundaries conditions.
\code
===List of boundary pieces for Dirichlet BCs on level set
\endcode
<li>What function of reconstruction \f$F\f$ is used?
\code
===How are the variables reconstructed from the level set function? (lin, reg)
\endcode
If \f$F\f$ is the identity set it to 'lin', set it to 'reg' otherwise.
<li>Coefficient \f$c_\text{reg}\f$ for reg reconstruction.
\code
===Value of the regularization coefficient in (0,0.5]
\endcode
If used, it is generally set to \f$0.5\f$.
<li>When reconstructing the fluid properties with the function \f$F\f$,
 do we use \f$F(\tilde{\phi})\f$ with \f$\tilde{\phi}=\min(1,\max(0,\phi))\f$?
\code
===Do we kill level set overshoot? (true/false)
\endcode
This parameter is set to false by default.
It allows to ensure that the density does not present overshoot.
<li>The order of the polynomial function used to approximated the level set.
\code
===Do we use P2 finite element for level_set? (true/false)
\endcode
This parameter is set to true by default.
</ol>


@subsection data_level_set_solver Information for the level set equation solver

Set the following information on the solver for the level set equation.
<ol>
<li>Maximum of iterations.
\code
===Maximum number of iterations for level set solver
\endcode
<li>Relative tolerance.
\code
===Relative tolerance for level set solver
\endcode
<li>Absolute tolerance.
\code
===Absolute tolerance for level set solver
\endcode
<li>Type of solver.
\code
===Solver type for level set (FGMRES, CG, ...)
\endcode
It is possible to set it to GMRES (Generalized Minimal Residual method),
 FGMRES (Flexible Generalized Minimal Residual method)
 and CG (Conjugate Gradient method).
<li>Type of preconditionner.
\code
===Preconditionner type for level set solver (HYPRE, JACOBI, MUMPS...)
\endcode
</ol>




@section data_temperature Setting for the temperature equation

We refer to the section \ref doc_intro_SFEMaNS_possibilities_temp
 for details on the weak formulation of the heat equation.

@subsection data_temperature_general General settings
To consider thermal effect, set the following parameters.
<ol>
<li>Is there thermal effect?
\code
===Is there a temperature field?
\endcode
The default value is false. The following information are only read if it is set to true.
<li>The number of subdomains where the heat equation is solved.
\code
===Number of subdomains in temperature mesh
\endcode
<li>The list of the subdomains where the heat equation is solved.
\code
===List of subdomains for temperature mesh
\endcode
<li>The volumetric heat capacity \f$C\f$ in each domains.
\code
===Volumetric heat capacity (1:nb_dom_temp)
\endcode
<li>The thermal conductivity \f$\lambda\f$ in each domains.
\code
===Thermal conductivity (1:nb_dom_temp)
\endcode
<li>It is possible to approximate the heat equation
 with the diffusivity coefficient \f$\kappa=\frac{\lambda}{C}\f$.
\code
===Diffusivity coefficient for temperature (1:nb_dom_temp)
\endcode
In that case, the equations solved is
 \f$\partial_t T + \DIV(T\bu)-\DIV(\kappa\GRAD T)=0\f$.
<li>Nondimensional gravity coefficient that be used to define the
 action of the temperature of the velocity field.
\code
===Non-dimensional gravity coefficient
\endcode
This parameter is saved in inputs\%gravity_coefficient.
<li>Number of boundary pieces with Dirichlet boundary condition.
\code
===How many boundary pieces for Dirichlet BCs on temperature?
\endcode
<li>List of the boundary pieces with Dirichlet boundary condition.
\code
===List of boundary pieces for Dirichlet BCs on temperature
\endcode
<li>Number of interfaces between the velocity field and the temperature field domains.
\code
===Number of interfaces between velocity and temperature only domains (for nst applications)
\endcode
<li>List of interfaces between the velocity field and the temperature field domains.
\code
===List of interfaces between velocity and temperature only domains (for nst applications)
\endcode
</ol>

@subsection data_temperature_solver Information for the temperature equation solver

Set the following information on the solver for the temperature equation.
<ol>
<li>Maximum of iterations.
\code
===Maximum number of iterations for temperature solver
\endcode
<li>Relative tolerance.
\code
===Relative tolerance for temperature solver
\endcode
<li>Absolute tolerance.
\code
===Absolute tolerance for temperature solver
\endcode
<li>Type of solver.
\code
===Solver type for temperature (FGMRES, CG, ...)
\endcode
It is possible to set it to GMRES (Generalized Minimal Residual method),
 FGMRES (Flexible Generalized Minimal Residual method)
 and CG (Conjugate Gradient method).
<li>Type of preconditionner.
\code
===Preconditionner type for temperature solver (HYPRE, JACOBI, MUMPS...)
\endcode
</ol>



@section data_maxwell Setting for the Maxwell equations

We refer to the section \ref doc_intro_SFEMaNS_possibilities_mxw
 for details on the weak formulation of the Maxwell equations.

@subsection data_maxwell_general General settings

The parameters required to approximate solutions of the Maxwell
 equations can be splitted in two main groups.
<ol>
<li>Information on the approximation of the magnetic field in the conducting domain.
<ol>
<li>The dependent variable in the conducting region can be \f$\bH\f$ or \f$\bB\f$.
\code
===Solve Maxwell with H (true) or B (false)?
\endcode
It is set to false by default. When using a magnetic permeability
 that depends of the time and the azimuthal direction, it has to be set to false.
<li>Number of conducting subdomains where the Maxwell equation
 is approximated with the magnetic field \f$\bH\f$ or \f$\bB\f$.
\code
===Number of subdomains in magnetic field (H) mesh
\endcode
<li>List of the conducting subdomains.
\code
===List of subdomains for magnetic field (H) mesh
\endcode
<li>Number of interfaces in the conducting domains.
\code
===Number of interfaces in H mesh
\endcode
 These interfaces represents surfaces with jump in the
 magnetic permeability or surfaces where boundary conditions
 are applied on the velocity field or the temperature field.
<li>List of the above interfaces.
\code
===List of interfaces in H mesh
\endcode
<li>Number of boundary pieces with Dirichlet conditions on \f$\bH\times\textbf{n}\f$.
\code
===Number of Dirichlet sides for Hxn
\endcode
<li>List of the boundary pieces with Dirichlet condition.
\code
===List of Dirichlet sides for Hxn
\endcode
<li>Is the magnetic permeability defined analytically?
\code
===Is permeability defined analytically (true/false)?
\endcode
It is set to false by default. If set to true, the magnetic permeability
 is defined in a function of the file <tt>condlim.f90</tt>.
<li>If the magnetic permeability is defined analytically, does it depends of the time
 and the azimuthal direction?
\code
===Is permeability variable in theta (true/false)?
\endcode
It is set to false by default and the magnetic permeability is defined in the function
 <tt>mu_bar_in_fourier_space</tt> of the file <tt>condlim.f90</tt>.
 If set to true, this function defined the stabilization term
 \f$\overline{\mu}\f$ and the magnetic permeability is defined
 in the function <tt>mu_in_real_space</tt>.
<li>Is the magnetic permeability is computed on gauss points with Gaussian quadratures?
\code
===Use FEM Interpolation for magnetic permeability (true/false)?
\endcode
It is set to true by default.
<li>If the permeability is not defined analytically,
 set the magnetic permeability in each subdomains.
\code
===Permeability in the conductive part (1:nb_dom_H)
\endcode
<li>Electrical conductivity in each subdomains.
\code
===Conductivity in the conductive part (1:nb_dom_H)
\endcode
<li>The order of the polynomial function used to approximated the magnetic field.
\code
===Type of finite element for magnetic field
\endcode
It can be set to one or two.
<li>Magnetic Reynolds number \f$\Rm\f$.
\code
===Magnetic Reynolds number
\endcode
<li>Stabilization coefficient \f$\beta_1\f$, 
 see the section \ref doc_intro_SFEMaNS_possibilities_mxw_1. 
\code
===Stabilization coefficient (divergence)
\endcode
<li>Stabilization coefficient \f$\beta_3\f$,
 see the section \ref doc_intro_SFEMaNS_possibilities_mxw_1.
\code
===Stabilization coefficient for Dirichlet H and/or interface H/H
\endcode
</ol>
<li>Information on the approximation of the scalar potential \f$\phi\f$ in the insulating domain.
<ol>
<li>Number of insulating subdomains where the scalar potential \f$\phi\f$
 is used to approximate the Maxwell equations.
\code
===Number of subdomains in magnetic potential (phi) mesh
\endcode
It is set to zero by default. The following information are only read if it is set to true.
<li>List of the insulating subdomains.
\code
===List of subdomains for magnetic potential (phi) mesh
\endcode
<li>Number of boundary pieces with Dirichlet condition.
\code
===How many boundary pieces for Dirichlet BCs on phi?
\endcode
<li>List of the boundary pieces with Dirichlet condition.
\code
===List of boundary pieces for Dirichlet BCs on phi
\endcode
<li>Number of interfaces between the conducting and insulating subdomains.
\code
===Number of interfaces between H and phi
\endcode
<li>List of the interface between the conducting and insulating subdomains.
\code
===List of interfaces between H and phi
\endcode
<li>Permeability in the insulating domain. It does not depends of the insulating subdomains.
\code
===Permeability in vacuum
\endcode
<li>The order of the polynomial function used to approximated the scalar potential.
\code
===Type of finite element for scalar potential
\endcode
It can be set to one or two. It needs to be larger than the one used for the magnetic field.
<li>Stabilization coefficient \f$\beta_2\f$,
 see the section \ref doc_intro_SFEMaNS_possibilities_mxw_1.
\code
===Stabilization coefficient (interface H/phi)
\endcode
</ol>
<li>In addition, it is possible to use a quasi-static approximation.
\code
===Quasi-static approximation (true) or (false)?
\endcode
It is set to false by default. We refer to the section
 \ref doc_SFEMaNS_condlim for more details on this set up.
</ol>



@subsection data_maxwell_solver Information for the Maxwell equations solver

Set the following information on the solver for the Maxwell equations.
<ol>
<li>Maximum of iterations.
\code
===Maximum number of iterations for Maxwell solver
\endcode
<li>Relative tolerance.
\code
===Relative tolerance for Maxwell solver
\endcode
<li>Absolute tolerance.
\code
===Absolute tolerance for Maxwell solver
\endcode
<li>Type of solver.
\code
===Solver type for Maxwell (FGMRES, CG, ...)
\endcode
It is possible to set it to GMRES (Generalized Minimal Residual method),
 FGMRES (Flexible Generalized Minimal Residual method)
 and CG (Conjugate Gradient method).
<li>Type of preconditionner.
\code
===Preconditionner type for Maxwell solver (HYPRE, JACOBI, MUMPS...)
\endcode
</ol>




@section data_arpack Setting for eigenvalue problem with Arpack

Eigenvalue problems for the Maxwell equations can be considered.
 The following parameters needs to be set.
<ol>
<li>Is Arpack used?
\code
===Do we use Arpack?
\endcode
It is set to false by default. The following information are only read if it is set to true.
<li>Number of eigenvalues to compute.
\code
===Number of eigenvalues to compute
\endcode
<li>Maximum number of Arpack iteration.
\code
===Maximum number of Arpack iteration
\endcode
<li>Tolereance for Arpack solver.
\code
===Tolerance for Arpack
\endcode
<li>Which eigen values are approximated (largest/smallest magnitude, real part or imaginary part).
\code
===Which eigenvalues (''LM'', ''SM'', ''SR'', ''LR'' ''LI'', ''SI'')
\endcode
<li>Generation of visualization files for Paraview.
\code
===Create 2D vtu files for Arpack? (true/false)
\endcode
It is set to false by default.
</ol>



@section data_output Output settings

Information on the outputs computed are provided as follows.
<ol>
<li>Frequency to write restart files of the variables approximated.
\code
===Frequency to write restart file
\endcode
This value is saved in inputs\%freq_restart. The suite files for the
 Navier-Stokes variables have the following name suite_ns_Sxxx_Iyyy.mesh_name where:
<ol>
<li>xxx is the number of the subsection of the meridian plane (000, 001, etc.)
 associated to the restart file.
<li>yyy is the time iteration divided by inputs\%freq_restart when the file was generated. 
<li>mesh_name is the name of the mesh.
</ol>
For the temperature, respectively the magnetic field, the nst
 is replaced by temp, respectively by maxwell.
<li>Frequency to write the energies or other ouputs computed in
 the subroutine <tt>my_post_processing</tt> of the file <tt>main.f90</tt>
\code
===Frequency to write energies
\endcode
This value is saved in inputs\%freq_en.
<li>Frequency to generate visualization files for Paraview.
\code
===Frequency to create plots
\endcode
This value is saved in inputs\%freq_plot
<li>Number of planes used to generate the visualization file.
\code
===Number of planes in real space for Visualization
\endcode
It is set to ten by default.
<li>It is possible to only do postprocessing of restarts files.
\code
===Just postprocessing without computing? (true/false)
\endcode
</ol>

The section \ref doc_SFEMaNS_main_my_post_processing details 
 how to compute these outputs. Additional information can be computed
 by the code in the output file. These informations are written
 each inputs\%freq_en time iterations.
<ol>
<li>Do the code checks the numerical stability?
\code
===Check numerical stability (true/false)
\endcode
The code check that the \f$\bL^2\f$-norm of the velocity field or
 the magnetic field is not larger than a hundred. This value can
 be changed in the file <tt>main.f90</tt>, look for
 inputs\%check_numerical_stability.
<li>Compute the average computational time per iteration.
\code
===Verbose timing? (true/false)
\endcode
<li>Compute the divergence of the velocity field and the magnetic field.
\code
===Verbose divergence? (true/false)
\endcode
<li>For problem involving the Navier-Stokes equations, the CFL can be computed.
\code
===Verbose CFL? (true/false)
\endcode
</ol>


@section data_zero_mode Setting to force Fourier components of velocity field and magnetic field to zero

It is possible to set a list of Fourier components of all variables of
 the Navier-Stokes and Maxwell equations to zero. It is done as follows.
<ol>
<li>Are some Fourier modes set to zero?
\code
===Should some modes be zeroed out?
\endcode
The default value if false. The following information are only read if it is set to true.
<li>Number of Fourier modes to set to zero for the Navier-Stokes equations approximations.
\code
===How many Navier-Stokes modes to zero out?
\endcode
<li>List of above Fourier modes.
\code
===List of Navier-Stokes modes to zero out?
\endcode
<li>Number of Fourier modes to set to zero for the Maxwell equations approximations.
\code
===How Maxwell modes to zero out?
\endcode
<li>List of above Fourier modes.
\code
===List of Maxwell modes to zero out?
\endcode
</ol>

 */
