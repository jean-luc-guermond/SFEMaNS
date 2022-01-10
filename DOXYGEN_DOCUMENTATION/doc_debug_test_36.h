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
 * @page doc_debug_test_36 Test 36: Viscosity function of temperature

<h3>Introduction</h3>
In this example, we check the correctness of SFEMaNS in cases of problems where the variation of the kinematic viscosity with respect to temperature is taken into account. The viscosity parameter is a function of the temperature instead of being a constant. As the temperature evolves in the system, so does this parameter, which then presents time and space variations. This leads to an additional term in the RHS of the formulation, corresponding to the variable part of the viscous stresses.

The domain of computation is \f$\Omega= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,1) \times [0,2\pi) \times (-1,1)\} \f$. We enforce Dirichlet BCs for the temperature and the velocity on all the boundaries.

We solve the temperature equations:
\f{align*}{
\begin{cases}
\partial_t T + \bu \cdot \GRAD T - \kappa \LAP T &= f_T, \\
T_{|\Gamma} &= T_\text{bdy} , \\
T_{|t=0} &= T_0,
\end{cases}
\f}
in the domain \f$\Omega\f$. 

We solve the Navier-Stokes equations:
\f{align*}{
\begin{cases}
\partial_t\bu + \left(\ROT\bu\right)\CROSS\bu  + \GRAD p - div (2 \nu(T) \GRAD^s \bu)  
&= \bef,
\\ \DIV \bu &= 0, \\
\bu_{|\Gamma} &= \bu_{\text{bdy}} , \\
\bu_{|t=0} &= \bu_0, \\
p_{|t=0} &= p_0,
\end{cases}
\f}
in the domain \f$\Omega\f$ as well, where \f$\GRAD^s \bu = \frac{1}{2} ( \GRAD \bu + (\GRAD \bu)^T )\f$ is the symmetric part of the gradient.

The kinematic viscosity is a function of the temperature
\f{align*}{
\nu(T) = \bar{\nu} + \tilde{\nu}(T) ,
\f}
where \f$\bar{\nu}\f$ is a constant, the maximum value of \f$\nu(T)\f$, and \f$\tilde{\nu}\f$ is the variable part.

The data are the source terms \f$f_T\f$ and \f$\bef\f$, the boundary data \f$T_\text{bdy}\f$ and \f$\bu_{\text{bdy}}\f$, and the initial data \f$T_0\f$, \f$\bu_0\f$ and \f$p_0\f$. The parameters are the thermal diffusivity \f$\kappa\f$ and the maximum kinematic viscosity \f$\bar{\nu}\f$.

<h3>Manufactured solutions</h3>

We approximate the following analytical solution:
@f{align*}{
T(r,\theta,z,t) & = \frac{r^3+e^z}{1+e}\cos(t)^2 ,
\\ u_r(r,\theta,z,t) &= 0,
\\ u_{\theta}(r,\theta,z,t) &= r^2 \sin(t-z),
\\ u_z(r,\theta,z,t) &=  0,
\\ p(r,\theta,z,t) &= 0,
@f}
where the temperature is chosen to stay between 0 and 1 so that \f$\tilde{\nu}\f$ stays negative. The variable part of the kinematic viscosity is defined by 
\f{align*}{
\tilde{\nu}(T) = - \frac{\bar{\nu}}{2} T .
\f}

The source terms \f$f_T, \bef\f$ and the boundary data \f$T_\text{bdy}, \bu_{\text{bdy}}\f$ are computed accordingly.

<h3>Generation of the mesh</h3>

The finite element mesh used for this test is named <tt>RECT_10.FEM</tt> and 
 has a mesh size of \f$0.1\f$ for the P1 approximation. 
The following image shows the mesh for P1 finite elements.
<table width="100%" align="center">
<tr>
    <td align="center">
    @image html  fig_RECT_10.png "Finite element mesh (P1)."
    </td>
</tr>
</table>

<h3>Information on the file <tt>condlim.f90</tt></h3>
The initial conditions, boundary conditions and the forcing terms \f$\textbf{f}\f$ in the Navier-Stokes
 equations and \f$f_T\f$ in the temperature equations are set in the file <tt>condlim_test_36.f90</tt>.
Here is a description of the subroutines and functions of interest.
<ol>
<li>The subroutine <tt>init_velocity_pressure</tt> initializes the velocity field
 and the pressure at the times \f$-dt\f$ and \f$0\f$ with \f$dt\f$ being the time step.
 This is done by using the functions vv_exact and pp_exact as follows:
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
<li>The subroutine <tt>init_temperature</tt> initializes the temperature at the times \f$-dt\f$ and \f$0\f$ with \f$dt\f$ the time step. 
This is done by using the function temperature_exact as follows:
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
 conditions on the velocity field. The \f$r\f$ and \f$z\f$ coordinates of all nodes are defined. The velocity is azimuthal and carried by the mode 0.
\code
    r = rr(1,:)
    z = rr(2,:)

    IF ((TYPE == 3) .AND. (m == 0)) THEN
       vv = r**2 * sin(t - z)
    ELSE
       vv = 0.d0
    END IF
\endcode
<li>The function <tt>pp_exact</tt> contains the analytical pressure. It is used to initialize the pressure. The pressure is here equal to 0.
\code
    vv = 0.d0
\endcode
<li>The function <tt>temperature_exact</tt> contains the analytical temperature. It is used to initialize the temperature and to impose Dirichlet boundary condition on the temperature. Again, the coordinates array are defined for simplicity. The temperature is carried by the mode 0. 
\code
    r = rr(1,:)
    z = rr(2,:)

    IF ((TYPE == 1) .AND. (m == 0)) THEN
       vv = (r**3 + exp(z)) / (1.d0 + exp(1.d0)) * cos(t)**2
    ELSE
       vv = 0.d0
    END IF
\endcode
<li>The function <tt>source_in_temperature</tt> computes the source term \f$f_T = \partial_t T + \bu \cdot \GRAD T - \kappa \LAP T\f$ of the temperature equations. The second term is null. The diffusivity parameter \f$\kappa\f$ is defined in the data file. The source term is carried by the mode 0.
\code
    IF ((TYPE == 1) .AND. (m == 0)) THEN
       vv = - (r**3 + exp(z)) / (1.d0 + exp(1.d0)) * sin(2 * t) &
            - inputs%temperature_diffusivity(1) * (9 * r + exp(z)) / (1.d0 + exp(1.d0)) * cos(t)**2
    ELSE
       vv = 0.d0
    END IF
\endcode
<li>The function <tt>source_in_NS_momentum</tt> computes the source term \f$\bef = \partial_t\bu + \left(\ROT\bu\right)\CROSS\bu  + \GRAD p - div (2 \nu(T) \GRAD^s \bu)\f$ of the Navier-Stokes equations. The implementation is done term by term (3 terms, the pressure term being null).
\code
    r = rr(1,:)
    z = rr(2,:)
    
    nu_bar = inputs%Re**(-1)
    nu = nu_bar - (nu_bar/2) * (r**3 + exp(z)) / (1.d0 + exp(1.d0)) * cos(time)**2

    ! du/dt

    IF ((TYPE == 3) .AND. (mode == 0)) THEN
       vv = r**2 * cos(time - z)
    ELSE 
       vv = 0.d0
    END IF
    
    ! (curl u) cross u

    IF ((TYPE == 1) .AND. (mode == 0)) THEN
       vv = vv - 3 * r**3 * sin(time - z)**2
    ELSE IF ((TYPE == 5) .AND. (mode == 0)) THEN
       vv = vv + r**4 * cos(time - z) * sin(time - z)
    END IF

    ! - 2 * div (nu * grads u)

    IF ((TYPE == 3) .AND. (mode == 0)) THEN
       vv = vv &
            + (nu_bar/2)/(1.d0 + exp(1.d0)) * r**2*cos(time)**2 * (3*r*sin(time - z) - exp(z)*cos(time - z)) & 
            - nu * (3.d0 - r**2) * sin(time - z)
    END IF  
\endcode
</ol>
All the other subroutines present in the file <tt>condlim_test_36.f90</tt> are not used in this test.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file.


<h3>Setting in the data file</h3>
We describe the data file of this test. It is called <tt>debug_data_test_36</tt> and can be found in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<ol>
<li>We use a formatted mesh by setting:
\code
===Is mesh file formatted (true/false)?
.t.
\endcode
<li>The path and the name of the mesh are specified with the two following lines:
\code
===Directory and name of mesh file
'.' 'RECT_10.FEM'
\endcode
where '.' refers to the directory where the data file is, meaning ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<li>We use \f$4\f$ processors in the meridian section.
\code
===Number of processors in meridian section
4
\endcode
<li>We solved the problem for \f$2\f$ Fourier modes.
\code
===Number of Fourier modes
2
\endcode
The Fourier modes are not detailed so the first 2 modes \f$0,1\f$ are solved. We verify that the second mode (\f$1\f$) is null.
<li>We use \f$2\f$ processors in Fourier space.
\code
===Number of processors in Fourier space
2
\endcode
It means that each processors is solving the problem for \f$2/2=1\f$ Fourier modes.
<li>We approximate the Navier-Stokes equations by setting:
\code
===Problem type: (nst, mxw, mhd, fhd)
'nst'
\endcode
<li>We do not restart the computations from previous results.
\code
===Restart on velocity (true/false)
.f.
\endcode
It means the computation starts from the time \f$t=0\f$.
<li>We use a time step of \f$0.01\f$ and solve the problem over \f$10\f$ time iterations.
\code
===Time step and number of time iterations
1.d-2 10
\endcode
<li>We set the number of domains and their label, see the files associated to the generation of the mesh, where the code approximates the Navier-Stokes equations.
\code
===Number of subdomains in Navier-Stokes mesh
1
===List of subdomains for Navier-Stokes mesh
1
\endcode
<li>We set the number of boundaries with Dirichlet conditions on the velocity field and give their respective labels.
\code
===How many boundary pieces for full Dirichlet BCs on velocity?
3
===List of boundary pieces for full Dirichlet BCs on velocity
2 3 4
\endcode
<li>We set the constant part of the kinetic viscosity \f$\bar{\nu}\f$ and indicate that the kinematic viscosity is variable.
\code
===Kinematic viscosity
10.d0
===Variable viscosity (true/false)?
.t.
\endcode
<li>The penalty term of the divergence of \f$\bu\f$ has a coefficient of 1.
\code
===Coefficient for penalty of divergence in NS?
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
<li>We solve the temperature equation.
\code
===Is there a temperature field?
.t.
\endcode
<li>We set the number of domains and their label, see the files associated to the generation of the mesh, where the code approximates the temperature equation.
\code
===Number of subdomains in temperature mesh
1
===List of subdomains for temperature mesh
1
\endcode
<li>We set the thermal diffusivity parameter \f$\kappa\f$.
\code
===Diffusivity coefficient for temperature (1:nb_dom_temp)
3.d0
\endcode
<li>We set the number of boundaries with Dirichlet conditions on the velocity and give their respective labels.
\code
===How many boundary pieces for Dirichlet BCs on temperature?
3
===List of boundary pieces for Dirichlet BCs on temperature
2 3 4
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
that can be found in: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. 

To check the well behavior of the code, we compute four quantities:
<ol>
<li>The L2-norm of error on u divided by L2-norm of u exact.
<li>The H1-norm of error on u divided by H1-norm of u exact.
<li>The L2-norm of error on p.
<li>The H1-norm of error on T divided by H1-norm of T exact.
</ol>
These quantities are computed at the final time \f$t=0.1\f$.
 They are compared to reference values to attest of the correct behavior of the code.  
 These values of reference are in the last lines of the file <tt>debug_data_test_36</tt> in ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
============================================
(RECT_10.FEM)
===Reference results
5.858165337128355E-006  L2-norm of error on u / L2-norm of u exact
6.849107330069875E-005  H1-norm of error on u / H1-norm of u exact
4.361164116502296E-005  L2-norm of error on p
1.533231503293184E-006  L2-norm of error on T / L2-norm of T exact
\endcode

To conclude this test, we display the profile of the approximated pressure, velocity and temperature at the final time.
 The figure is done in the plane \f$y=0\f$ which is the union of the half plane \f$\theta=0\f$ and \f$\theta=\pi\f$.
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_36_pre_tfin.png "Pressure in the plane y=0."
    </td>
</tr>
<tr>
    <td align="center">
     @image html  fig_test_36_vel_tfin.png "Azimuthal compound of the velocity in the plane y=0."
    </td>
</tr>
<tr>
    <td align="center">
     @image html  fig_test_36_temp_tfin.png "Temperature in the plane y=0."
    </td>
</tr>
</table>
 */
