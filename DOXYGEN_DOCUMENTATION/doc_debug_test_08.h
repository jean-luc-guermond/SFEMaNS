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
 * @page doc_debug_test_08 Test 8: Navier-Stokes with Temperature

<h3>Introduction</h3>
In this example, we check the correctness of SFEMaNS for a thermohydrodynamic problem involving Dirichlet boundary conditions.

We solve the temperature equation:
@f{align*}
\partial_t T+ \bu \cdot \GRAD T - \kappa \LAP T &= f_T, \\
T_{|\Gamma} &= T_\text{bdy} , \\
T_{|t=0} &= T_0,
@f}
and the Navier-Stokes equations:
@f{align*}
\partial_t\bu+\left(\ROT\bu\right)\CROSS\bu - \frac{1}{\Re}\LAP \bu +\GRAD p     
&= \alpha T \textbf{e}_z + \bef,
\\ \DIV \bu &= 0, \\
\bu_{|\Gamma} &= \bu_{\text{bdy}} , \\
\bu_{|t=0} &= \bu_0, \\
p_{|t=0} &= p_0,
@f}
in the domain \f$\Omega= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,1/2] \times [0,2\pi) \times [0,1]\} \f$ with  \f$\Gamma= \partial \Omega\f$. We denote by \f$\textbf{e}_z\f$ the unit vector in the vertical direction.
The data are the source terms \f$f_T\f$ and \f$\bef\f$,
 the boundary datas \f$T_\text{bdy}\f$ and \f$\bu_{\text{bdy}}\f$,
 the initial data \f$T_0\f$, \f$\bu_0\f$ and \f$p_0\f$.
The parameters are the thermal diffusivity \f$\kappa\f$, the kinetic Reynolds number \f$\Re\f$ and a thermal gravity number \f$\alpha\f$. Recall that these parameters are dimensionless.

<h3>Manufactured solutions</h3>
We approximate the following analytical solutions:
@f{align*}
T(r,\theta,z,t) & =  \left(r^2 z + r^2z^2(\cos(\theta)+2\sin(2\theta))\right)\cos(t),
\\ u_r(r,\theta,z,t) &= r^3\cos(2\pi z)\sin(t) , 
\\ u_{\theta}(r,\theta,z,t) &= r^2z\sin(t), 
\\ u_z(r,\theta,z,t) &= -\frac{2r^2}{\pi}\sin(2\pi z)\sin(t), 
\\ p(r,\theta,z,t) &= \frac{1}{2} \left( r^6\cos(2\pi z)^2 + r^4z^2 + \frac{4r^4}{\pi^2}\sin(2\pi z)^2
  \right) \sin(t)^2,
@f}
where the source terms \f$f_T, \bef\f$ and the boundary datas \f$T_\text{bdy}, \bu_{\text{bdy}}\f$ are computed accordingly.

<h3>Generation of the mesh</h3>

The finite element mesh used for this test is named <tt>Mesh_10_form.FEM</tt> and 
 has a mesh size of \f$0.1\f$ for the P1 approximation. 
You can generate this mesh with the files in the following directory:
 ($SFEMaNS_MESH_GEN_DIR)/EXAMPLES/EXAMPLES_MANUFACTURED_SOLUTIONS/Mesh_10_form.
The following image shows the mesh for P1 finite elements.
<table width="70%" align="center">
<tr>
    <td align="center">
    @image html  fig_Mesh_10_form.png  "Finite element mesh (P1)."
    </td>
</tr>
</table>


<h3>Information on the file <tt>condlim.f90</tt></h3>
The initial conditions, boundary conditions and the forcing terms \f$\textbf{f}\f$ in the Navier-Stokes
 equations and \f$f_T\f$ in the temperature equation are set in the file <tt>condlim_test_8.f90</tt>.
Here is a description of the subroutines and functions of interest.
<ol>
<li>First we define the number \f$\pi\f$ at the begining of the module so that 
every subroutines has access to this real number.
\code
  REAL(KIND=8) :: pi=ACOS(-1.d0)
\endcode
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
<li>The subroutine <tt>init_temperature</tt> initializes the temperature at the times \f$-dt\f$ and \f$0\f$ with \f$dt\f$ being the time step. 
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
 conditions on the velocity field.
<ol>
<li>We construct the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>If the Fourier mode m is not equal to 0, the velocity field is set to zero.
\code 
    IF (m/=0) THEN
       vv = 0
       RETURN
    END IF
\endcode
<li>For the Fourier mode \f$m=0\f$, we define the velocity field depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) as follows:
\code
    IF (TYPE == 1) THEN
       vv(:) =  r**3*COS(2*PI*z)
    ELSE IF (TYPE == 2) THEN
       vv(:) = 0
    ELSE IF (TYPE == 3) THEN
       vv(:) = r**2*z
    ELSE IF (TYPE == 4) THEN
       vv(:) = 0.d0
    ELSE IF (TYPE == 5) THEN
       vv(:) = -(2*r**2/PI)*SIN(2*PI*z)
    ELSE IF (TYPE == 6) THEN
       vv(:) = 0.d0
    ENDIF
    vv(:) = vv(:) * SIN(t)
    RETURN
\endcode
where \f$t\f$ is the time.
</ol>
<li>The function <tt>pp_exact</tt> contains the analytical pressure. It is used to initialize the pressure.
<ol>
<li>We define the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>If the Fourier mode m is not equal to 0, the pressure is set to zero.
\code 
    IF (m/=0) THEN
       vv = 0 
       RETURN
    END IF
\endcode
<li>For the Fourier mode \f$m=0\f$, the pressure only depends of the TYPE 1 (cosine) so we write:
\code
    IF (TYPE == 1) THEN
       vv(:) =  0.5d0*((r**3*COS(2*PI*z))**2 &
            + (r**2*z)**2 &
            + ((2*r**2/PI)*SIN(2*PI*z))**2)*SIN(t)**2
    ELSE IF (TYPE == 2) THEN
       vv(:) = 0
    END IF
    RETURN
\endcode
where \f$t\f$ is the time.
 We note that the component sine of a mode \f$m=0\f$ is always zero because \f$\sin(0)=0\f$.
</ol>
<li>The function <tt>temperature_exact</tt> contains the analytical temperature. It is used to initialize the temperature and to impose Dirichlet boundary condition on the temperature.
<ol>
<li>We construct the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>If the Fourier mode m is strictly larger than 2, the temperature is set to zero.
\code
    IF (m>2) THEN
       vv =0
       RETURN
\endcode
<li>For the Fourier mode \f$m\leq2\f$, we set the analytical temperature depending of its TYPE (1 for cosine and 2 for sine) as follows:
\code
    ELSE
       !temperature = (r**2*z + r**2*z**2(cos(theta)+2*sin(2*theta))cos(t)
       IF(m==0 .AND. TYPE==1) THEN
          vv = r**2*z*COS(t)
       ELSE IF (m==0 .AND. TYPE==2) THEN
          vv = 0.d0
       ELSE IF (m==1 .AND. TYPE==1) THEN
          vv = r**2*z**2*COS(t)
       ELSE IF (m==1 .AND. TYPE==2) THEN
          vv = 0
       ELSE IF (m==2 .AND. TYPE==1) THEN
          vv = 0
       ELSE IF (m==2 .AND. TYPE==2) THEN
          vv = 2*r**2*z**2*COS(t)
       ELSE
          CALL error_petsc('BUG in temperature_exact')
       END IF
    END IF
    RETURN
\endcode
where \f$t\f$ is the time.
We note that the subroutine <tt>error_petsc</tt> is here to stop the code in the case the input TYPE is
 different than 1 or 2 (for cosine and sine).
</ol>
<li>The function <tt>source_in_temperature</tt> computes the source term \f$f_T\f$ of the temperature equation.
<li>The function <tt>source_in_NS_momentum</tt> computes the source term \f$\alpha T \textbf{e}_z+\bef\f$ of the Navier-Stokes equations.
</ol>
All the other subroutines present in the file <tt>condlim_test_8.f90</tt> are not used in this test.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file.




<h3>Setting in the data file</h3>
We describe the data file of this test.
 It is called <tt>debug_data_test_8</tt> and can be found in the  following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<ol>
<li>We use a formatted mesh by setting:
\code
===Is mesh file formatted (true/false)?
.t.
\endcode
<li>The path and the name of the mesh are specified with the two following lines:
\code
===Directory and name of mesh file
'.' 'Mesh_10_form.FEM'
\endcode
where '.' refers to the directory where the data file is, meaning ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<li>We use one processor in the meridian section. It means the finite element mesh is not subdivised.
\code
===Number of processors in meridian section
1
\endcode
<li>We solve the problem for \f$3\f$ Fourier modes.
\code
===Number of Fourier modes
3
\endcode
<li>We use \f$3\f$ processors in Fourier space.
\code
===Number of processors in Fourier space
3
\endcode
It means that each processors is solving the problem for \f$3/3=1\f$ Fourier modes.
<li>We select specific Fourier modes to solve.
\code
===Select Fourier modes? (true/false)
.t.
\endcode
<li>We give the list of the Fourier modes to solve.
\code
===List of Fourier modes (if select_mode=.TRUE.)
0 1 2
\endcode
We note that setting select Fourier modes to false would give the same result 
 as we select the first three Fourier modes.
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
<li>We use a time step of \f$0.02\f$ and solve the problem over \f$50\f$ time iterations.
\code
===Time step and number of time iterations
.02d0, 50
\endcode
<li>We set the number of domains and their label, see the files associated to the generation of the mesh,
 where the code approximates the Navier-Stokes equations.
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
5 2 4
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
<li>We solve the temperature equation (in the same domain than the Navier-Stokes equations).
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
1
\endcode
<li>We set the thermal diffusivity \f$\kappa\f$.
\code
===Diffusivity coefficient for temperature
1.d-1
\endcode
<li>We set the thermal gravity number \f$\alpha\f$.
\code
===Non-dimensional gravity coefficient
10.d0
\endcode
<li>We set the number of boundaries with Dirichlet conditions on the temperature and give their respective labels.
\code
===How many boundary pieces for Dirichlet BCs on temperature?
3
===List of boundary pieces for Dirichlet BCs on temperature
5 2 4
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
that can be found in the following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. 

To check the well behavior of the code, we compute four quantities:
<ol>
<li>The L2 norm of the error on the temperature.
<li>The H1 norm of the error on the temperature.
<li>The L2 norm of the divergence of the velocity field.
<li>The L2 norm of the error on the pressure.
</ol>
These quantities are computed at the final time \f$t=1\f$.
 They are compared to reference values to attest of the correctness of the code.  
 These values of reference are in the last lines of the file <tt>debug_data_test_8</tt> in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
============================================
(Mesh_10_form.FEM)
===Reference results
2.411363013813094E-006  L2 error on temperature
1.347463767239025E-004  H1 error on temperature
5.567600159075602E-003  L2 norm of divergence
1.791439464958010E-003  L2 error on pressure
\endcode


To conclude this test, we show the profile of the approximated, pressure, velocity magnitude
 and temperature at the final time.
 These figures are done in the plane \f$y=0\f$ which
 is the union of the half plane \f$\theta=0\f$ and \f$\theta=\pi\f$.
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_08_pre_tfin.png "Pressure in the plane plane y=0."
    </td>
    <td align="center">
    @image html  fig_test_08_vel_tfin.png  "Velocity magnitude in the plane plane y=0."
    </td>
</tr>
</table>
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_08_temp_tfin.png "Temperature in the plane plane y=0."
    </td>
</tr>
</table>
 */
