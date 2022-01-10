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
 * @page doc_debug_test_15 Test 15: Navier-Stokes with Temperature and LES

<h3>Introduction</h3>
In this example, we check the correctness of SFEMaNS for a thermohydrodynamic problem. We use Dirichlet boundary conditions and a stabilization method called entropy viscosity that we refer as LES.  We note this test does not involve manufactured solution and consist to check four quantities, like the kinetic energy of specific Fourier modes, are the same than the values of reference.

We solve the temperature equation:
@f{align*}
\partial_t T+ \bu \cdot \GRAD T - \kappa \LAP T &= 0, \\
T_{|\Gamma} &= T_\text{bdy} , \\
T_{|t=0} &= T_0,
@f}
and the Navier-Stokes equations:
@f{align*}
\partial_t\bu+\left(\ROT\bu\right)\CROSS\bu - \frac{1}{\Re}\LAP \bu +\GRAD p     
&= \DIV (\nu_E \GRAD \bu  ) + \alpha T \textbf{e}_z,
\\ \DIV \bu &= 0, \\
\bu_{|\Gamma} &= \bu_{\text{bdy}} , \\
\bu_{|t=0} &= \bu_0, \\
p_{|t=0} &= p_0,
@f}
in the domain \f$\Omega= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,1] \times [0,2\pi) \times [-1/2,1/2]\} \f$ with  \f$\Gamma= \partial \Omega\f$. We denote by \f$\textbf{e}_z\f$ the unit vector in the vertical direction.
 The term \f$\DIV (\nu_E \GRAD \bu  )\f$ is a stabilization term that involve an artificial viscosity called the entropy viscosity. We refer to the section \ref doc_intro_SFEMaNS_possibilities_nst_2 for the definition of this term.
The data are the boundary datas \f$T_\text{bdy}\f$ and \f$\bu_{\text{bdy}}\f$,
 the initial data \f$T_0\f$, \f$\bu_0\f$ and \f$p_0\f$.
The parameters are the thermal diffusivity \f$\kappa\f$, the kinetic Reynolds number \f$\Re\f$, the thermal gravity number \f$\alpha\f$ and the real \f$c_\text{e}\f$ for the entropy viscosity. We remind that these parameters are dimensionless.

<h3>Manufactured solutions</h3>
As mentionned earlier this test does not involve manufactured solutions. As consequence we do not consider specific source term and only initialize the variables to approximate.
@f{align*}
T(r,\theta,z,t=0) & =  -z - (z-0.5)(z+0.5) \left(1+\cos(\theta)+\sin(\theta)+\cos(2\theta)+\sin(2\theta) \right) ,
\\ u_r(r,\theta,z,t=0) &= 0, 
\\ u_{\theta}(r,\theta,z,t=0) &= 0,
\\ u_z(r,\theta,z,t=0) &=0,
\\ p(r,\theta,z,t=0) &= 0,
@f}
The boundary datas \f$T_\text{bdy}, \bu_{\text{bdy}}\f$ are computed accordingly.

<h3>Generation of the mesh</h3>

The finite element mesh used for this test is named <tt>RECT10_BENCHMARK_CONVECTION_LES.FEM</tt> and 
 has a mesh size of \f$0.1\f$ for the P1 approximation. You can generate this mesh with the files in the following directory:
 ($SFEMaNS_MESH_GEN_DIR)/EXAMPLES/EXAMPLES_MANUFACTURED_SOLUTIONS/RECT10_BENCHMARK_CONVECTION_LES.
The following image shows the mesh for P1 finite elements.
<table width="70%" align="center">
<tr>
    <td align="center">
    @image html  fig_RECT10_BENCHMARK_CONVECTION_LES.png  "Finite element mesh (P1)."
    </td>
</tr>
</table>


<h3>Information on the file <tt>condlim.f90</tt></h3>
The initial conditions, boundary conditions and the forcing term in the Navier-Stokes
 equations are set in the file <tt>condlim_test_15.f90</tt>.
Here is a description of the subroutines and functions of interest.
<ol>
<li>The subroutine <tt>init_velocity_pressure</tt> initializes the velocity field
 and the pressure at the time\f$-dt\f$ and \f$0\f$ with \f$dt\f$ being the time step.
 It is done by using the functions vv_exact and pp_exact as follows:
\code
    time = 0.d0
    DO i= 1, SIZE(list_mode)
       mode = list_mode(i) 
       DO j = 1, 6 
          !===velocity
          un_m1(:,j,i) = vv_exact(j,mesh_f%rr,mode,time-dt)  
          un   (:,j,i) = vv_exact(j,mesh_f%rr,mode,time)
       END DO
       DO j = 1, 2
          !===pressure
          pn_m2(:)       = pp_exact(j,mesh_c%rr,mode,time-2*dt)
          pn_m1  (:,j,i) = pp_exact(j,mesh_c%rr,mode,time-dt)
          pn     (:,j,i) = pp_exact(j,mesh_c%rr,mode,time)
          phin_m1(:,j,i) = pn_m1(:,j,i) - pn_m2(:)
          phin   (:,j,i) = Pn   (:,j,i) - pn_m1(:,j,i)
       ENDDO
    ENDDO
\endcode
<li>The subroutine <tt>init_temperature</tt> initializes the temperature at the time\f$-dt\f$ and \f$0\f$ with \f$dt\f$ being the time step. It is done by using the function temperature_exact as follows:
\code
    time = 0.d0
    DO i= 1, SIZE(list_mode)
       mode = list_mode(i) 
       DO j = 1, 2 
          tempn_m1(:,j,i) = temperature_exact(j, mesh%rr, mode, time-dt)
          tempn   (:,j,i) = temperature_exact(j, mesh%rr, mode, time)
       ENDDO
    ENDDO
\endcode
<li>The function <tt>vv_exact</tt> contains the analytical velocity field.
 It is used to initialize the velocity field and to impose Dirichlet boundary
 conditions on the velocity field. It is set to zero.
\code
    vv(:) = 0.d0
    RETURN
\endcode
<li>The function <tt>pp_exact</tt> contains the analytical pressure. It is used to initialize the pressure and is set to zero.
\code
    vv=0.d0
    RETURN
\endcode
<li>The function <tt>temperature_exact</tt> contains the analytical temperature. It is used to initialize the temperature and to impose Dirichlet boundary condition on the temperature.
<ol>
<li>We construct the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li> We set the temperature to zero.
\code
    vv=0.d0
\endcode
<li>For the Fourier mode \f$m=0\f$  the temperature only depends of the TYPE 1 (cosine).
\code
    IF (m==0 .AND. TYPE==1) THEN
       vv(:)= - z - (z-5d-1)*(z+5d-1)
\endcode
<li>For the Fourier mode \f$m\geq1\f$, the temperature does not depend of the TYPE (1 for cosine and 2 for sine) and is defined as follows:
\code
    ELSE IF (m.GE.1) THEN
       vv= -(z-5d-1)*(z+5d-1)
    END IF
    RETURN
\endcode
</ol>
<li>The function <tt>source_in_temperature</tt> computes the source term denoted \f$f_T\f$ in previous tests, of the temperature equation. As it is not used in this test, we set it to zero.
\code
    vv = 0.d0
    RETURN
\endcode
<li>The function <tt>source_in_NS_momentum</tt> computes the source term \f$\alpha T \textbf{e}_z\f$ of the Navier-Stokes equations depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) as follows:
\code
    IF (TYPE==5) THEN 
       vv = inputs%gravity_coefficient*opt_tempn(:,1,i)
    ELSE IF (TYPE==6) THEN
       vv = inputs%gravity_coefficient*opt_tempn(:,2,i)
    ELSE 
       vv = 0.d0
    END IF
    RETURN
\endcode
</ol>
All the other subroutines present in the file <tt>condlim_test_15.f90</tt> are not used in this test.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file.




<h3>Setting in the data file</h3>
We describe the data file of this test.
 It is called <tt>debug_data_test_15</tt> and can be found in the following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<ol>
<li>We use a formatted mesh by setting:
\code
===Is mesh file formatted (true/false)?
.t.
\endcode
<li>The path and the name of the mesh are specified with the two following lines:
\code
===Directory and name of mesh file
'.' 'RECT10_BENCHMARK_CONVECTION_LES.FEM'
\endcode
where '.' refers to the directory where the data file is, meaning ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<li>We use one processor in the meridian section. It means the finite element mesh is not subdivised. To do so, we write:
\code
===Number of processors in meridian section
1
\endcode
<li>We solve the problem for \f$3\f$ Fourier modes.
\code
===Number of Fourier modes
3
\endcode
<li>We use \f$3\f$ processors in Fourier.
\code
===Number of processors in Fourier space
3
\endcode
It means that each processors is solving the problem for \f$3/3=1\f$ Fourier modes.
<li>We do not select specific Fourier modes to solve.
\code
===Select Fourier modes? (true/false)
.f.
\endcode
<li>We approximate the Navier-Stokes equations by setting:
\code
===Problem type: (nst, mxw, mhd, fhd)
'nst'
\endcode
<li>We do not restart the computations from previous results.
\code
===Restart on velocity (true/false)
.f.
===Restart on magnetic field (true/false)
.f.
\endcode
It means the computation starts from the time \f$t=0\f$.
<li>We use a time step of \f$0.05\f$ and solve the problem over \f$10\f$ time iterations.
\code
===Time step and number of time iterations
5.d-2, 10
\endcode
<li>We set the number of domains and their label, see the files associated to the generation of the mesh,
 where the code approximates the Navier-Stokes equations,
\code
===Number of subdomains in Navier-Stokes mesh
1
===List of subdomains for Navier-Stokes mesh
1
\endcode
<li>We set the number of boundaries with Dirichlet conditions on the velocity field and give their respective labels.
\code
===How many boundary pieces for full Dirichlet BCs on velocity?
2
===List of boundary pieces for full Dirichlet BCs on velocity
2 3
\endcode
<li>We set the kinetic Reynolds number \f$\Re\f$.
\code
===Reynolds number
50.d0
\endcode
<li>We use the entropy viscosity method to stabilize the equation.
\code
===Use LES? (true/false)
.t.
\endcode
<li>We define the coefficient \f$c_\text{e}\f$ of the entropy viscosity.
\code
===Coefficient multiplying residual
1.d0 
\endcode
<li>We give information on how to solve the matrix associated to the time marching of the velocity.
<ol>
<li>
\code
===Maximum number of iterations for velocity solver
100
\endcode
<li>
\code
===Relative tolerance for velocity solver
1.d-6
===Absolute tolerance for velocity solver
1.d-10
\endcode
<li>
\code
===Solver type for velocity (FGMRES, CG, ...)
GMRES
===Preconditionner type for velocity solver (HYPRE, JACOBI, MUMPS...)
MUMPS
\endcode
</ol>
<li>We give information on how to solve the matrix associated to the time marching of the pressure.
<ol>
<li>
\code
===Maximum number of iterations for pressure solver
100
\endcode
<li>
\code
===Relative tolerance for pressure solver
1.d-6
===Absolute tolerance for pressure solver
1.d-10
\endcode
<li>
\code
===Solver type for pressure (FGMRES, CG, ...)
GMRES
===Preconditionner type for pressure solver (HYPRE, JACOBI, MUMPS...)
MUMPS
\endcode
</ol>
<li>We give information on how to solve the mass matrix.
<ol>
<li>
\code
===Maximum number of iterations for mass matrix solver
100
\endcode
<li>
\code
===Relative tolerance for mass matrix solver
1.d-6
===Absolute tolerance for mass matrix solver
1.d-10
\endcode
<li>
\code
===Solver type for mass matrix (FGMRES, CG, ...)
CG
===Preconditionner type for mass matrix solver (HYPRE, JACOBI, MUMPS...)
MUMPS
\endcode
</ol>
<li>We solve the temperature equation (in the same domain than the Navier-Stokes equations).
\code
===Is there a temperature field?
.t.
\endcode
<li>We set the number of domains and their label, see the files associated to the generation of the mesh,
 where the code approximated the temperature equation.
\code
===Number of subdomains in temperature mesh
1
===List of subdomains for temperature mesh
1
\endcode
<li>We set the thermal diffusivity \f$\kappa\f$.
\code
===Diffusivity coefficient for temperature (1:nb_dom_temp)
1.d0
\endcode
<li>We set the thermal gravity number \f$\alpha\f$.
\code
===Non-dimensional gravity coefficient
50.d0
\endcode
<li>We set the number of boundaries with Dirichlet conditions on the temperature and give their respective labels.
\code
===How many boundary pieces for Dirichlet BCs on temperature?
1
===List of boundary pieces for Dirichlet BCs on temperature
2
\endcode
<li>We give information on how to solve the matrix associated to the time marching of the temperature.
<ol>
<li>
\code
===Maximum number of iterations for temperature solver
100
\endcode
<li>
\code
===Relative tolerance for temperature solver
1.d-6
===Absolute tolerance for temperature solver
1.d-10
\endcode
<li>
\code
===Solver type for temperature (FGMRES, CG, ...)
GMRES
===Preconditionner type for temperature solver (HYPRE, JACOBI, MUMPS...)
MUMPS
\endcode
</ol>
<li>To get the total elapse time and the average time in loop minus initialization, we write:
\code
===Verbose timing? (true/false)
.t.
\endcode
These informations are written in the file <tt>lis</tt> when you run the shell <tt>debug_SFEMaNS_template</tt>.
</ol>




<h3> Outputs and value of reference </h3>

The outputs of this test are computed with the file <tt>post_processing_debug.f90</tt> 
that can be found in the following: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. 

To check the well behavior of the code, we compute four quantities:
<ol>
<li>The kinetic energy of the Fourier mode \f$m=0\f$.
<li>The kinetic energy of the Fourier mode \f$m=1\f$.
<li>The kinetic energy of the Fourier mode \f$m=2\f$.
<li>The L2 norm of the velocity divergence divided by the L2 norm of the gradient of the velocity.
</ol>
These quantities are computed at the final time \f$t=1\f$.
 They are compared to reference values to attest of the correctness of the code.  
 These values of reference are in the last lines of the file <tt>debug_data_test_15</tt> in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
============================================
RECT10_BENCHMARK_CONVECTION_LES.FEM, dt=5.d-2, it_max=10
===Reference results
6.33685640432350423E-004    !e_c_u_0
0.12910286398104689         !e_c_u_1
0.10838939329896366         !e_c_u_2
4.93403150219499723E-002    !||div(un)||_L2/|un|_sH1
\endcode


To conclude this test, we show the profile of the approximated, pressure, velocity magnitude
 and temperature at the final time.
 These figures are done in the plane \f$y=0\f$ which
 is the union of the half plane \f$\theta=0\f$ and \f$\theta=\pi\f$.
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_15_pre_tfin.png "Pressure in the plane plane y=0."
    </td>
    <td align="center">
    @image html  fig_test_15_vel_tfin.png  "Velocity magnitude in the plane plane y=0."
    </td>
</tr>
</table>
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_15_temp_tfin.png "Temperature in the plane plane y=0."
    </td>
</tr>
</table>
 */
