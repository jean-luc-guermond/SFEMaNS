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
 * @page doc_debug_test_16 Test 16: Navier-Stokes Neumann bdy with precession

<h3>Introduction</h3>
In this example, we check the correctness of SFEMaNS for a hydrodynamic problem of a precession set up involving Neumann boundary conditions. The main rotation axis is the vertical axis while the precession axis is along the unit normal vector \f$\textbf{e}_x\f$ associated to the x cartesian coordinate. The computation is done in the precession frame, meaning the walls only see the main rotation along the vertical axis.
 We note this test does not involve manufactured solution and consist to check four quantities, like the total kinetic energy, are the same than the values of reference.

We solve the Navier-Stokes equations:
@f{align*}
\partial_t\bu+\left(\ROT\bu + 2 \epsilon \textbf{e}_x \right)\CROSS\bu - \frac{1}{\Re}\LAP \bu +\GRAD p   &=0,
\\ \DIV \bu &= 0, \\
\bu \cdot \textbf{n}_{|\Gamma} &= 0, \\
\bu_{|t=0} &= \bu_0, \\
p_{|t=0} &= p_0,
@f} 
in the domain \f$\Omega= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,1] \times [0,2\pi) \times [0,.8]\; | \; r^2 + \frac{z^2}{0.8^2} =1\} \f$ with  \f$\Gamma= \partial \Omega \f$.
The data are the initial datas \f$\bu_0\f$ and \f$p_0\f$. The parameter \f$\Re\f$ is the kinetic Reynolds number, \f$\epsilon\f$ is the precession rate and \f$\alpha\f$ is the precession angle.

<h3>Manufactured solutions</h3>
As mentionned earlier this test does not involve manufactured solutions. As consequence we do not consider specific source term and only initialize the variables to approximate.
@f{align*}
u_r(r,\theta,z,t=0) &= 0,
\\ u_{\theta}(r,\theta,z,t=0) &= 0.1r,
\\ u_z(r,\theta,z,t=0) &=0,
\\ p(r,\theta,z,t=0) &= \frac{(0.1 r)^2}{2} .
@f}

<h3>Generation of the mesh</h3>

The finite element mesh used for this test is named <tt>ELL_b0p8_10_form.FEM</tt>.
 The mesh size for the P1 approximation \f$0.1\f$ on the boundary of the domain and \f$0.033\f$ at \f$(r,z)=(0,0)\f$. You can generate this mesh with the files in the following directory:
 ($SFEMaNS_MESH_GEN_DIR)/EXAMPLES/EXAMPLES_MANUFACTURED_SOLUTIONS/ELL_b0p8_10_form.
The following image shows the mesh for P1 finite elements.
<table width="70%" align="center">
<tr>
    <td align="center">
    @image html  fig_ELL_b0p8_10_form.png  "Finite element mesh (P1)."
    </td>
</tr>
</table>


<h3>Information on the file <tt>condlim.f90</tt></h3>

The initial conditions, boundary conditions and the forcing term \f$\textbf{f}\f$ in the Navier-Stokes
 equations are set in the file <tt>condlim_test_16.f90</tt>.
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
 It is used to initialize the velocity field.
<ol>
<li>First we define the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>For a time strictly larger than 0, we set the valocity field to zero.
\code    
    IF (t>1.d-14) THEN
          vv = 0.d0
\endcode
We note that the above line are not required since this test does not involved Dirichlet boundary conditions. As a consequence this function is only used to initialize the velocity field.
<li>If the Fourier mode m is not equal to 0, the velocity field is set to zero.
\code
    ELSE
       IF (m/=0) THEN
          vv = 0.d0
\endcode
<li>For the Fourier mode \f$m=0\f$, we define the velocity field depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) as follows:
\code
       ELSE
          IF (TYPE==3) THEN
             vv = 0.1d0*r
          ELSE
             vv = 0.d0
          END IF
       END IF
    END IF
    RETURN
\endcode
</ol>
<li>The function <tt>pp_exact</tt> contains the analytical pressure.
 It is used to initialize the pressure.
<ol>
<li>First we define the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>For a time strictly larger than 0, we set the pressure to zero.
\code
    IF (t>1.d-14) THEN
       vv = 0.d0
    ELSE
\endcode
We note that the above line are not required since this function is only used to initialize the pressure.
<li>If the Fourier mode m is not equal to 0, the pressure is set to zero.
\code
       IF (m/=0) THEN
          vv = 0.d0
\endcode
<li>For the Fourier mode \f$m=0\f$, we define the pressure depending of its TYPE (1 for cosine and 2 for sine) as follows:
\code
       ELSE
          IF (TYPE==1) THEN
             vv = (0.1d0*r)**2/2.d0
          ELSE
             vv = 0.d0
          END IF
       END IF
    END IF
    RETURN
\endcode
We note that the sine part of a Fourier mode \f$m=0\f$ is always zero.
</ol>
<li>The function <tt>source_in_NS_momentum</tt>is used to define the source term \f$\bef\f$ of the Navier-Stokes equations. As this term is not used for this test, it is set to zero.
\code
    vv=0.d0
    RETURN
\endcode
</ol>
All the other subroutines present in the file <tt>condlim_test_16.f90</tt> are not used in this test.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file.




<h3>Setting in the data file</h3>
We describe the data file of this test.
 It is called <tt>debug_data_test_16</tt> and can be found in the following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<ol>
<li>We use a formatted mesh by setting:
\code
===Is mesh file formatted (true/false)?
.t.
\endcode
<li>The path and the name of the mesh are specified with the two following lines:
\code
===Directory and name of mesh file
'.' 'ELL_b0p8_10_form.FEM'
\endcode
where '.' refers to the directory where the data file is, meaning ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<li>We use two processors in the meridian section. It means the finite element mesh is subdivised in two. To do so, we write:
\code
===Number of processors in meridian section
2
\endcode
<li>We solve the problem for \f$8\f$ Fourier modes.
\code
===Number of Fourier modes
8
\endcode
<li>We use \f$4\f$ processors in Fourier.
\code
===Number of processors in Fourier space
4
\endcode
It means that each processors is solving the problem for \f$8/4=2\f$ Fourier modes.
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
\endcode
It means the computation starts from the time \f$t=0\f$.
<li>We use a time step of \f$0.1\f$ and solve the problem over \f$20\f$ time iterations.
\code
===Time step and number of time iterations
1d-1, 20 
\endcode
<li>We set the number of domains and their label, see the files associated to the generation of the mesh,
 where the code approximates the Navier-Stokes equations,
\code
===Number of subdomains in Navier-Stokes mesh
1
===List of subdomains for Navier-Stokes mesh
1
\endcode
<li>We set the number of boundaries with Dirichlet conditions on the velocity field.
\code
===How many boundary pieces for full Dirichlet BCs on velocity?
0
\endcode
<li>We set the number of boundaries with homogeneous Neumann conditions on the velocity field and give their respective labels.
\code
===How many boundary pieces for homogeneous normal velocity?
1
===List of boundary pieces for homogeneous normal velocity
2
\endcode
<li>We set the kinetic Reynolds number \f$\Re\f$.
\code
===Reynolds number
42.d0
\endcode
<li>We want to add the term \f$ 2\epsilon \textbf{e}_x \times \bu \f$ in the left hand side of the Navier-Stokes equations. Such features is already programmed in SFEMaNS via the following lines in your data file.
<ol>
<li>We set the option precession to true.
\code
===Is there a precession term (true/false)?
.t.
\endcode
It adds the term \f$ 2\epsilon \left(\sin(\alpha\pi)\textbf{e}_x+\cos(\alpha\pi)\textbf{e}_z\right) \times \bu \f$ in the left hand side of the Navier-Stokes equations.
<li>We set the precession rate \f$\epsilon\f$ and the precesion angle \f$\alpha\f$.
\code
===Precession rate
0.25d0
===Precession angle over pi
0.5d0
\endcode
</ol>
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
<li>The total kinetic energy \f$ \displaystyle 0.5 {\Vert \bu \Vert}_{\bL^2(\Omega)}^2\f$.
<li>The quantity \f$M_x= \displaystyle \int_\Omega
 \left(-z(\sin(\theta)u_r + \cos(\theta) u_\theta) +r \sin(\theta)u_z   \right) r   dr d\theta dz\f$.
<li>The quantity \f$M_y= \displaystyle \int_\Omega 
\left(z(\cos(\theta)u_r - \sin(\theta) u_\theta ) -r \cos(\theta)u_z  \right) r   dr d\theta dz\f$.
<li>The quantity \f$M_z= \displaystyle \int_\Omega  u_\theta r^2 dr d\theta dz \f$.
</ol>
These quantities are computed at the final time \f$t=2\f$. We note the last three quantities are computed via the subroutine <tt>angular_momentum </tt> of the file <tt>tn_axi.f90</tt>.
 They are compared to reference values to attest of the correctness of the code.  
 These values of reference are in the last lines of the file <tt>debug_data_test_16</tt> in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
============================================
ELL_b0p8_10_form.FEM, dt=1d-1, it_max=20
===Reference results
6.67555315567430665E-003  !Total kinetic energy at t=2
9.61565539080621234E-004  !Mx
4.87365427729861689E-002  !My
0.12184513917556984       !Mz
\endcode


To conclude this test, we show the profile of the approximated pressure and velocity magnitude at the final time. These figures are done in the plane \f$y=0\f$ which is the union of the half plane \f$\theta=0\f$ and \f$\theta=\pi\f$.
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_16_pre_tfin.png "Pressure in the plane plane y=0."
    </td>
    <td align="center">
    @image html  fig_test_16_vel_tfin.png  "Velocity magnitude in the plane plane y=0."
    </td>
</tr>
</table>
 */
