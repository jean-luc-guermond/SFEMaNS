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
 * @page doc_debug_test_24 Test 24: Navier-Stokes with penalty for non axisymetric domain

<h3>Introduction</h3>
In this example, we check the correctness behavior of SFEMaNS for a hydrodynamic problem with a solid obstacle involving Dirichlet boundary conditions.

We solve the Navier-Stokes equations:
@f{align*}
\partial_t\bu+\left(\ROT\bu\right)\CROSS\bu - \frac{1}{\Re}\LAP \bu +\GRAD p     &=\bef \text{ in } \Omega_1,
\\ \bu & = 0 \text{ in } \Omega_2,
\\ \DIV \bu &= 0, \\
\bu_{|\Gamma} &= \bu_{\text{bdy}} , \\
\bu_{|t=0} &= \bu_0, \\
p_{|t=0} &= p_0,
@f}
in the domain \f$\Omega= \Omega_1 \cup \Omega_2\f$  with \f$ \Omega_1=\{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0.1,1/2] \times [0,2\pi) \times [0,1]\}\f$ and \f$\Omega_2=\{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [1/2,1] \times [0,2\pi) \times [0,1]\}  \f$. We also define \f$\Gamma= \partial \Omega \f$. We note that the condition \f$\bu=0\f$ in \f$\Omega_2\f$ is imposed via a penalty method that involves a penalty function \f$\chi\f$ equal to 1 in \f$\Omega_1\f$ and zero elsewhere.
The data are the source term \f$\bef\f$, the penalty function \f$\chi\f$, the boundary data \f$\bu_{\text{bdy}}\f$, the initial datas \f$\bu_0\f$ and \f$p_0\f$. The parameter \f$\Re\f$ is the kinetic Reynolds number.

<h3>Manufactured solutions</h3>
We approximate the following analytical solutions:
@f{align*}
u_r(r,\theta,z,t) &= (2r-1)^2 \sin(r+t) \mathbb{1}_{r\geq0.5},
\\ u_{\theta}(r,\theta,z,t) &= 0,
\\ u_z(r,\theta,z,t) &= \left( (2-1/r) (6r-1) \cos(r+t) + (r-0.5) \sin(2\theta)  \right) \mathbb{1}_{r\geq0.5},
\\ p(r,\theta,z,t) &= r^2 z^3\cos(t) + r \cos(\theta) ,
@f}
with \f$ \mathbb{1}_{r\geq0.5} \f$ the function equals to \f$r\geq0.5\f$ if \f$ \f$ and \f$0\f$ elsewhere. The source term \f$\bef\f$ and the boundary data \f$ \bu_{\text{bdy}}\f$ are computed accordingly.

<h3>Generation of the mesh</h3>

The finite element mesh used for this test is named <tt>cylinder_0.05.FEM</tt> and 
 has a mesh size of \f$0.05\f$ for the P1 approximation. You can generate this mesh with the files in the following directory:
 ($SFEMaNS_MESH_GEN_DIR)/EXAMPLES/EXAMPLES_MANUFACTURED_SOLUTIONS/cylinder_0.05.
The following image shows the mesh for P1 finite elements.
<table width="70%" align="center">
<tr>
    <td align="center">
    @image html  fig_cylinder_0.05.FEM.png  "Finite element mesh (P1)."
    </td>
</tr>
</table>


<h3>Information on the file <tt>condlim.f90</tt></h3>

The initial conditions, boundary conditions, the forcing term and the penalty function are set in the file <tt>condlim_test_24.f90</tt>.
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
<li>The function <tt>vv_exact</tt> contains the analytical velocity field.
 It is used to initialize the velocity field and to impose Dirichlet boundary
 conditions on the velocity field.
<ol>
<li>First we define the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>We define the velocity field depending of the Fourier mode and its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) as follows:
\code
    IF (TYPE==1.AND.m==0) THEN
       DO n = 1, SIZE(rr,2)
          IF (rr(1,n)>0.5d0) THEN
             vv(n) = (2*rr(1,n)-1)**2*SIN(rr(2,n)+t)
          ELSE
             vv(n) = 0.d0
          END IF
       END DO
    ELSE IF (TYPE==5.AND.m==0) THEN
       DO n = 1, SIZE(rr,2)
          IF (rr(1,n)>0.5d0) THEN
             vv(n) = (2-1.d0/rr(1,n))*(6*rr(1,n)-1)*COS(rr(2,n)+t)
          ELSE
             vv(n) = 0.d0
          END IF
       END DO
    ELSE IF (TYPE==6.AND.m==2) THEN
       DO n = 1, SIZE(rr,2)
          IF (rr(1,n)>0.5d0) THEN
             vv(n) = rr(1,n)-0.5d0
          ELSE
             vv(n) = 0.d0
          END IF
       END DO
    ELSE
       vv = 0.d0
    END IF
    RETURN
\endcode
where \f$t\f$ is the time.
</ol>
<li>The function <tt>pp_exact</tt> contains the analytical pressure.
 It is used to initialize the pressure.
<ol>
<li>First we define the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>We define the pressure depending of the Fourier mode and its TYPE (1 for cosine and 2 for sine) as follows:
\code    
    IF (TYPE==1.AND.m==0) THEN
       vv(:) = r**2*z**3*COS(t)
    ELSE IF (TYPE==1.AND.m==1) THEN
       vv(:) = r
    ELSE
       vv = 0.d0
    END IF
    RETURN
\endcode
where \f$t\f$ is the time.
</ol>
<li>The function <tt>source_in_NS_momentum</tt> computes the source term \f$\bef\f$ of the Navier-Stokes equations.
<li> The function <tt>penal_in_real_space</tt> define the penalty function \f$\chi\f$ in the real space (depending of the node in the meridian plan and its angle n). It is done as follows:
\code
    DO n = nb, ne
       n_loc = n - nb + 1
       IF (rr_gauss(1,n_loc).LE.0.5d0) THEN
          vv(:,n_loc) = 0.d0
       ELSE
          vv(:,n_loc) = 1.d0
       END IF
    END DO
    RETURN
\endcode
As defined earlier, this function is equal to one when the cylindrical coordinate r is smaller than 0.5 and else is equal to 1.
<li>The function <tt>imposed_velocity_by_penalty</tt> defines the velocity in the solid domain \f$\Omega_2\f$. It is set to zero as follows:
\code
    vv=0.d0
    RETURN
\endcode
</ol>
All the other subroutines present in the file <tt>condlim_test_24.f90</tt> are not used in this test.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file.

<h3>Setting in the data file</h3>
We describe the data file of this test.
 It is called <tt>debug_data_test_24</tt> and can be found in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<ol>
<li>We use a formatted mesh by setting:
\code
===Is mesh file formatted (true/false)?
.t.
\endcode
<li>The path and the name of the mesh are specified with the two following lines:
\code
===Directory and name of mesh file
'.' 'cylinder_0.05.FEM'
\endcode
where '.' refers to the directory where the data file is, meaning ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<li>We use two processors in the meridian section. It means the finite element mesh is subdivised in two.
\code
===Number of processors in meridian section
2
\endcode
<li>We solve the problem for \f$6\f$ Fourier modes.
\code
===Number of Fourier modes
6
\endcode
<li>We use \f$6\f$ processors in Fourier space.
\code
===Number of processors in Fourier space
6
\endcode
It means that each processors is solving the problem for \f$6/6=1\f$ Fourier modes.
<li>We do not select specific Fourier modes to solve.
\code
===Select Fourier modes? (true/false)
.f.
\endcode
As a consequence, the code approximates the problem on the first \f$6\f$ Fourier modes.

<li>We approximate the Navier-Stokes equations by setting:
\code
===Problem type: (nst, mxw, mhd, fhd)
'nst'
\endcode
<li>We approximate the Navier-Stokes equations with the velocity field as dependent variable.
\code
===Solve Navier-Stokes with u (true) or m (false)?
.t.
\endcode
We note this data is set to true by default. The momentum \f$m\f$ is only used for multiphase flow problem.
<li>We do not restart the computations from previous results.
\code
===Restart on velocity (true/false)
.f.
\endcode
 It means the computation starts from the time \f$t=0\f$.
<li>We use a time step of \f$0.0005\f$ and solve the problem over \f$100\f$ time iterations.
\code
===Time step and number of time iterations
0.0005d0 100
\endcode
<li>We use a penalty function function to take into account the presence of a solid obstacle.
\code
===Use penalty in NS domain (true/false)?
.t.
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
2
===List of boundary pieces for full Dirichlet BCs on velocity
1 2
\endcode
<li>We set the kinetic Reynolds number \f$\Re\f$.
\code
===Reynolds number
100.d0
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
<li>The L2 norm of the error on the velocity field.
<li>The H1 norm of the error on the velocity field.
<li>The L2 norm of the divergence of the velocity field.
<li>The L2 norm of the error on the pressure outside the obstacle.
</ol>
These quantities are computed at the final time \f$t=0.05\f$.
 They are compared to reference values to attest of the correctness of the code.
  
 These values of reference are in the last lines of the file <tt>debug_data_test_24</tt> in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
============================================
(cylinder_0.05.FEM)
===Reference results
5.35272511967831415E-003  L2 error on velocity
0.41315380860605949       H1 error on velocity
0.19511906134562279       L2 norm of divergence
5.07028095271459204E-003  L2 error on pressure outter obstacle
\endcode


To conclude this test, we show the profile of the approximated pressure and velocity magnitude at the final time. These figures are done in the plane \f$y=0\f$ which is the union of the half plane \f$\theta=0\f$ and \f$\theta=\pi\f$.
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_24_pre_tfin.png "Pressure in the plane plane y=0."
    </td>
    <td align="center">
    @image html  fig_test_24_vel_tfin.png  "Velocity magnitude in the plane plane y=0."
    </td>
</tr>
</table>
 */
