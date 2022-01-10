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
 * @page doc_debug_test_35 Test 35: Robin boundary conditions for Temperature

<h3>Introduction</h3>
In this example, we check the correctness of SFEMaNS for the use of Robin boundary conditions on the temperature. This boundary condition is usefull when the cooling or the heating of a body by convection at the boundaries is considered. A mix of Dirichlet and Robin boundary conditions is used in this test.

The domain of computation is \f$\Omega= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in (1/2,1) \times [0,2\pi) \times (0,1)\} \f$. We enforce Dirichlet BCs on the union of the boundaries \f$r=1/2\f$ and \f$z=0\f$, called \f$\Gamma_D\f$. We enforce Robin BCs on the boundaries \f$r=1\f$, called \f$\Gamma_{R,1}\f$, and \f$z=1\f$, called \f$\Gamma_{R,2}\f$. We call \f$\Gamma = \partial \Omega\f$.

We solve the temperature equations:
\f{align*}{
\begin{cases}
\rho c\partial_t T + \rho c \bu \cdot \GRAD T - \DIV (\lambda \GRAD T) &= f_T, \\
T_{|\Gamma_D} &= T_\text{bdy} , \\
(- h_1 T - \lambda \GRAD T \SCAL \bn)_{|\Gamma_{R,1}} &= - h_1 T_1 , \\
(- h_2 T - \lambda \GRAD T \SCAL \bn)_{|\Gamma_{R,2}} &= - h_2 T_1 , \\
T_{|t=0} &= T_0,
\end{cases}
\f}
in the domain \f$\Omega\f$. The parameters \f$h_1\f$, \f$h_2\f$ and \f$T_1\f$ are constants that represent the convection ceofficients and the exterior temperature respectively.

We solve the Navier-Stokes equations:
\f{align*}{
\begin{cases}
\partial_t\bu+\left(\ROT\bu\right)\CROSS\bu - \frac{1}{\Re}\LAP \bu +\GRAD p     
&= \bef,
\\ \DIV \bu &= 0, \\
\bu_{|\Gamma} &= \bu_{\text{bdy}} , \\
\bu_{|t=0} &= \bu_0, \\
p_{|t=0} &= p_0,
\end{cases}
\f}
in the domain \f$\Omega\f$.

The data are the source terms \f$f_T\f$ and \f$\bef\f$, the boundary data \f$T_1\f$, \f$T_\text{bdy}\f$ and \f$\bu_{\text{bdy}}\f$, and the initial data \f$T_0\f$, \f$\bu_0\f$ and \f$p_0\f$. The parameters are the density \f$\rho\f$, the heat capacity \f$c\f$, the thermal conductivity \f$\lambda\f$, the convection coefficients \f$h_1\f$ and \f$h_2\f$ and the kinetic Reynolds number \f$\Re\f$.

<h3>Manufactured solutions</h3>

We approximate the following analytical solution:
@f{align*}{
T(r,\theta,z,t) & = \left( \exp \left( - \frac{h_1 r+h_2 z}{\lambda} \right) + (r-1)^2 + (z-1)^2 \right)(1+\cos(\theta))\cos(t) + T_1 ,
\\ u_r(r,\theta,z,t) &= 0,
\\ u_{\theta}(r,\theta,z,t) &= 0,
\\ u_z(r,\theta,z,t) &=  0,
\\ p(r,\theta,z,t) &= 0,
@f}
where the temperature is chosen to satisfy the Robin BCs.

The source terms \f$f_T, \bef\f$ and the boundary data \f$T_\text{bdy}, \bu_{\text{bdy}}\f$ are computed accordingly.

<h3>Generation of the mesh</h3>

The finite element mesh used for this test is named <tt>SOLID_FLUID_10.FEM</tt> and 
 has a mesh size of \f$0.1\f$ for the P1 approximation. 
You can generate this mesh with the files in the following directory:
 ($SFEMaNS_MESH_GEN_DIR)/EXAMPLES/EXAMPLES_MANUFACTURED_SOLUTIONS/SOLID_FLUID_10.
The following image shows the mesh for P1 finite elements.
<table width="100%" align="center">
<tr>
    <td align="center">
    @image html  fig_SOLID_FLUID_10.png  "Finite element mesh (P1)."
    </td>
</tr>
</table>

Only the right part of the mesh is used for the computation in this test: \f$ r \ge 1/2 \f$.

<h3>Information on the file <tt>condlim.f90</tt></h3>
The initial conditions, boundary conditions and the forcing terms \f$\textbf{f}\f$ in the Navier-Stokes
 equations and \f$f_T\f$ in the temperature equations are set in the file <tt>condlim_test_35.f90</tt>.
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
 conditions on the velocity field. The velocity is equal to 0 (for every type or mode).
\code
    vv = 0.d0
\endcode
<li>The function <tt>pp_exact</tt> contains the analytical pressure. It is used to initialize the pressure. The pressure is equal to 0 (for every type or mode).
\code
    vv = 0.d0
\endcode
<li>The function <tt>temperature_exact</tt> contains the analytical temperature. It is used to initialize the temperature and to impose Dirichlet boundary condition on the temperature.
<ol>
<li>We construct the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>We define the thermal conductivity, the convection coefficients and the exterior temperature.
\code
    lambda = inputs%temperature_diffusivity(1)
    h1 = inputs%convection_coeff(1)
    h2 = inputs%convection_coeff(2)
    T1 = inputs%exterior_temperature(1)
\endcode
<li>We define the temperature depending on its TYPE (1 and 2 for cosine and sine) and on its mode as follows:
\code
    x = (h1 * r + h2 * z) / lambda

    IF (TYPE==1) THEN
       IF (m==0) THEN
          vv = ( exp(-x) + (r - 1d0)**2 * (z - 1d0)**2 ) * cos(t) + T1
       ELSE IF (m==1) THEN
          vv = ( exp(-x) + (r - 1d0)**2 * (z - 1d0)**2 ) * cos(t)
       ELSE
          vv = 0d0
       END IF
    ELSE
       vv = 0.d0
    END IF
\endcode
</ol>
<li>The function <tt>source_in_temperature</tt> computes the source term \f$f_T\f$ of the temperature equations.
<ol>
<li>We construct the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>We define the volumetric heat capacity (the product \f$\rho c\f$), the thermal conductivity and the convection coefficients.
\code
    c = inputs%vol_heat_capacity(1)
    lambda = inputs%temperature_diffusivity(1)
    h1 = inputs%convection_coeff(1)
    h2 = inputs%convection_coeff(2)
\endcode
<li>The source term \f$ f_T = \rho c \partial_t T  + \rho c \bu \cdot \GRAD T - \DIV(\lambda\GRAD T) \f$ is defined in two parts since the velocity is null:
\code
    x = (h1 * r + h2 * z) / lambda

    ! source = c pdt T - lambda laplacien T

    ! c pdt T

    IF (TYPE==1) THEN
       IF ((m==0).OR.(m==1)) THEN
          vv = - c * ( exp(-x) + (r - 1d0)**2 * (z - 1d0)**2 ) * sin(t)
       ELSE
          vv = 0d0
       END IF
    ELSE
       vv = 0d0
    END IF

    ! - lambda laplacien T

    IF (TYPE==1) THEN
       IF (m==0) THEN
          vv = vv - lambda * ( &
               1d0/r * ( -h1/lambda * exp(-x) + 2 * (r - 1d0) * (z - 1d0)**2 ) * cos(t) &
               + ( (h1/lambda)**2 * exp(-x) + 2 * (z - 1d0)**2 ) * cos(t) &
               + ( (h2/lambda)**2 * exp(-x) + 2 * (r - 1d0)**2 ) * cos(t) )
       ELSE IF (m==1) THEN
          vv = vv - lambda * ( &
               1d0/r * ( -h1/lambda * exp(-x) + 2 * (r - 1d0) * (z - 1d0)**2 ) * cos(t) &
               + ( (h1/lambda)**2 * exp(-x) + 2 * (z - 1d0)**2 ) * cos(t) &
               - 1d0/r**2 * ( exp(-x) + (r - 1d0)**2 * (z - 1d0)**2 ) * cos(t) &
               + ( (h2/lambda)**2 * exp(-x) + 2 * (r - 1d0)**2 ) * cos(t) )
       END IF
    END IF
\endcode
</ol>
<li>The function <tt>source_in_NS_momentum</tt> computes the source term \f$\bef\f$ of the Navier-Stokes equations. It is null.
\code
    vv = 0d0
\endcode
</ol>
All the other subroutines present in the file <tt>condlim_test_35.f90</tt> are not used in this test.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file.




<h3>Setting in the data file</h3>
We describe the data file of this test.
 It is called <tt>debug_data_test_35</tt> and can be found in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<ol>
<li>We use a formatted mesh by setting:
\code
===Is mesh file formatted (true/false)?
.t.
\endcode
<li>The path and the name of the mesh are specified with the two following lines:
\code
===Directory and name of mesh file
'.' 'SOLID_FLUID_10.FEM'
\endcode
where '.' refers to the directory where the data file is, meaning ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<li>We use \f$2\f$ processors in the meridian section.
\code
===Number of processors in meridian section
2
\endcode
<li>We solved the problem for \f$2\f$ Fourier modes.
\code
===Number of Fourier modes
2
\endcode
The Fourier modes are not detailed so the first 2 modes \f$0,1\f$ are solved.
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
<li>We use a time step of \f$0.01\f$ and solve the problem over \f$100\f$ time iterations.
\code
===Time step and number of time iterations
1.d-2 100
\endcode
<li>We set the number of domains and their label, see the files associated to the generation of the mesh,
 where the code approximates the Navier-Stokes equations.
\code
===Number of subdomains in Navier-Stokes mesh
1
===List of subdomains for Navier-Stokes mesh
2
\endcode
<li>We set the number of boundaries with Dirichlet conditions on the velocity field and give their respective labels.
\code
===How many boundary pieces for full Dirichlet BCs on velocity?
4
===List of boundary pieces for full Dirichlet BCs on velocity
3 5 2 4
\endcode
<li>We set the kinetic Reynolds number \f$\Re\f$.
\code
===Reynolds number
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
<li>We set the number of domains and their label, see the files associated to the generation of the mesh,
 where the code approximates the temperature equation.
\code
===Number of subdomains in temperature mesh
1
===List of subdomains for temperature mesh
2
\endcode
<li>We set the density \f$\rho\f$, the heat capacity \f$c\f$ and the thermal conductivity \f$\lambda\f$.
\code
===Density (1:nb_dom_temp)
1.d0
===Heat capacity (1:nb_dom_temp)
2.d0
===Thermal conductivity (1:nb_dom_temp)
10.d0
\endcode
<li>We set the number of boundaries with Dirichlet conditions on the velocity and give their respective labels.
\code
===How many boundary pieces for Dirichlet BCs on temperature?
2
===List of boundary pieces for Dirichlet BCs on temperature
3 4
\endcode
<li>We set the number of boundaries with Robin conditions on the velocity and give their respective labels.
\code
===How many boundary pieces for Dirichlet BCs on temperature?
2
===List of boundary pieces for Dirichlet BCs on temperature
5 2
\endcode
<li>We set the convection coefficients \f$h_1\f$ and and \f$h_2\f$, and the exterior temperature \f$T_1\f$.
\code
===Convection heat transfert coefficient (1:temperature_nb_robin_sides)
5d0 2d0
===Exterior temperature (1:temperature_nb_robin_sides)
3d0 3d0
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
<li>The L2-norm of error on u.
<li>The L2-norm of error on p.
<li>The L2-norm of error on T divided by L2-norm of T exact.
<li>The H1-norm of error on T divided by H1-norm of T exact.
</ol>
These quantities are computed at the final time \f$t=1\f$.
 They are compared to reference values to attest of the correct behavior of the code.  
 These values of reference are in the last lines of the file <tt>debug_data_test_35</tt> in ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
============================================
(SOLID_FLUID_10.FEM)
===Reference results
0.000000000000000E+000  L2-norm of error on u
0.000000000000000E+000  L2-norm of error on p
3.017387149621566E-007  L2-norm of error on T / L2-norm of T exact
1.936024637254978E-005  H1-norm of error on T / H1-norm of T exact
\endcode

To conclude this test, we display the profile of the approximated temperature at the final time.
 The figure is done in the plane \f$y=0\f$ which is the union of the half plane \f$\theta=0\f$ and \f$\theta=\pi\f$.
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_35_temp_tfin.png "Temperature in the plane y=0."
    </td>
</tr>
</table>
 */
