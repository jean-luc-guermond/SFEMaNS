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
 * @page doc_debug_test_28 Test 28 Hydrodynamic Simulation of the VKS Experiment


<h3>Introduction</h3>
In this example, we check the correctness of SFEMaNS for a hydrodynamic problem with a moving solid obstacle involving Dirichlet boundary conditions. The set up consists of a fluid driven by contra rotating impellers in a cylindrical container. In the litterature, it is referred as Von Karman Sodium. Here we study the case with impellers called TM73, we refer to the paper <a href='http://www.math.tamu.edu/~guermond/PUBLICATIONS/Nore_castanon_cappanera_guermond_EPL_2016.pdf'><code>Direct numerical simulation of the axial dipolar dynamo in the Von Karman Sodium experiment</code> </a> (Nore et al. 2016) for more information on this set up.

This test does not involve manufactured solutions and consists of checking four quantities, like the \f$\bL^2\f$ norm of the velocity, are the same as the reference values.

We solve the Navier-Stokes equations:
@f{align*}
\partial_t\bu+\left(\ROT\bu\right)\CROSS\bu - \frac{1}{\Re}\LAP \bu +\GRAD p 
    &=\bef  &\text{ in } \Omega_\text{fluid},
\\ \bu & = r \omega  \be_\theta &\text{ in } \Omega_\text{imp_bot},
\\ \bu & = -r \omega \be_\theta &\text{ in } \Omega_\text{imp_top},
\\ \DIV \bu &= 0, &\\
\bu_{|\Gamma} &= \bu_{\text{bdy}} ,& \\
\bu_{|t=0} &= \bu_0, &\\
p_{|t=0} &= p_0,&
@f}
in the domain \f$\Omega= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,1] \times [0,2\pi) \times [-1,1]\}\f$. This domain is the union of a fluid domain \f$\Omega_\text{fluid}\f$ and two solid domains, \f$\Omega_\text{imp_top}\f$ and \f$\Omega_\text{imp_bot}\f$ that represent the impellers. These subdomains depend of time as the impellers are contra-rotating with the angular velocity \f$\omega\f$.
 We also define \f$\Gamma= \partial \Omega \f$. 
 The data are the source term \f$\bef\f$, the angular velocity \f$\omega\f$, the penalty function \f$\chi\f$, the boundary data \f$\bu_{\text{bdy}}\f$, the initial datas \f$\bu_0\f$ and \f$p_0\f$. The parameter \f$\Re\f$ is the kinetic Reynolds number.

Remark: The velocity field is forced to match the velocity of the impellers in the solid subdomains with a penalty method. This method involves a penalty function \f$\chi\f$ equal to 1 in \f$\Omega_\text{fluid}\f$ and zero elsewhere.


<h3>Manufactured solutions</h3>
As mentionned earlier this test does not involve manufactured solutions. As a consequence, we do not consider specific source term and only initialize the dependent variables.

The initial velocity field and pressure are initialized as follows:
@f{align*}
u_r(r,\theta,z,t) &= 0,
\\ u_{\theta}(r,\theta,z,t) &= 
\begin{cases}
 -r \omega \be_\theta & \text{ in } \Omega_\text{imp_top}, \\
0  &\text{ in } \Omega_\text{fluid}, \\
  r \omega \be_\theta & \text{ in } \Omega_\text{imp_bot},
\end{cases}
\\ u_z(r,\theta,z,t) &= 0,
\\ p(r,\theta,z,t) &= 0,
@f}
The penalty function \f$\chi\f$ is defined such that it is equal to 1 in the fluid domain \f$\Omega_1\f$ and \f$0\f$ elsewhere. We note that we use a smooth penalty function.


<h3>Generation of the mesh</h3>

The finite element mesh used for this test is named <tt>mesh_T28_0_04_04_ext3.FEM</tt> and 
 has a mesh size of \f$0.04\f$ for the P1 approximation. You can generate this mesh with the files in the following directory:
 ($SFEMaNS_MESH_GEN_DIR)/EXAMPLES/EXAMPLES_MANUFACTURED_SOLUTIONS/mesh_T28_0_04_04_ext3.
The following image shows the mesh for P1 finite elements.
<table width="70%" align="center">
<tr>
    <td align="center">
    @image html  fig_mesh_T28_0_04_04_ext3.png  "Finite element mesh (P1)."
    </td>
</tr>
</table>

The following images show a 3D representation of the VKS set up and the shape of the impellers that drive the fluid.
<table width="70%" align="center">
<tr>
    <td align="center">
    @image html  vks_hydro_setting_small.png "VKS Setting."
    </td>

    <td align="center">
    @image html  tm73_small.png  "Impeller TM73."
    </td>
</tr>
</table>


<h3>Information on the file <tt>condlim.f90</tt></h3>

The initial conditions, boundary conditions, the forcing term and the penalty function are set in the file <tt>condlim_test_28.f90</tt>.
Here is a description of the subroutines and functions of interest.
<ol>
<li>First we define a set of shape parameters at the begining of the module so that every subroutines has access to these parameters. These parameters are used to define the penalty function, meaning the shape of the impellers, or to set the velocity of the impellers.
<li>The subroutine <tt>init_velocity_pressure</tt> initializes the velocity field
 and the pressure at the time\f$-dt\f$ and \f$0\f$ with \f$dt\f$ being the time step.
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
<li>The function <tt>vv_exact</tt> contains the analytical velocity field.
 It is used to initialize the velocity field and to impose Dirichlet boundary
 conditions on the velocity field.
<ol>
<li>First we set the velocity to zero.
\code
    vv=0.0  
\endcode
<li>We define the velocity field depending of the Fourier mode and its TYPE (1 and 2 for the radial component, cosine and sine, 3 and 4 for the azimuthal component, cosine and sine, 5 and 6 for the vertical component, cosine and sine) as follows:
\code
    IF (type==3 .AND. m==0) THEN
       DO n = 1, SIZE(rr,2)
          r= rr(1,n)
          z= rr(2,n)
          ! Are we in the Bottom propeller?
          IF ( if_bottom_prop .AND. r <disk_radius .AND. z <  top_of_blade_bot ) then 
             vv(n)=solid_vel*r
          END IF
          
          !are we in the top Propeller?
          IF ( if_top_prop  .AND. r <disk_radius .AND. z  >   bot_of_blade_top) then 
             vv(n)=-solid_vel*r
             
          END IF
       END DO
    END IF    
    RETURN
\endcode
where solid_vel is the angular velocity of the impeller and is set to 1.
</ol>
<li>The function <tt>pp_exact</tt> contains the analytical pressure.
 It is used to initialize the pressure to zero.
\code
    vv=0.d0
    RETURN
\endcode
<li>The function <tt>source_in_NS_momentum</tt> computes the source term \f$\bef\f$ of the Navier-Stokes equations. It is set to zero.
<li> The function <tt>penal_in_real_space</tt> defines the penalty function \f$\chi\f$ in the real space (depending of the node in the meridian plan and its angle n). This is done by calling the function <tt>smooth_penal_in_real_space</tt> as follows:
\code
      vv=smooth_penal_in_real_space(mesh,rr_gauss,angles,nb_angles,nb,ne,time)
\endcode
This function is equal to zero in the impeller and 1 elsewhere. It is defined with the use other functions such as <tt>smooth_bottom_propeller</tt> or <tt>smooth_top_propeller</tt>.
<li>The function <tt>imposed_velocity_by_penalty</tt> is used to set a non zero velocity in the impellers.
<ol>
<li>First we set the empellers velocity to zero.
\code
    vv=0.d0
\endcode
<li>We do a loop on the node of the mesh and define the radial and vertical coordinates r, z.
\code
    DO n = 1, SIZE(rr,2)
       r= rr(1,n)
       z= rr(2,n)
\endcode
<li>If the node considered has a vertical coordinate smaller than \f$-0.5\f$, it is in the bottom impeller of angular velocity solid_vel.
\code
       IF (z<-0.5d0) THEN
          vv(n,3) =  solid_vel*rr(1,n)
\endcode
<li>Else the node is in the top impeller of angular velocity -solid_vel.
\code
       ELSE
          vv(n,3) =  -solid_vel*rr(1,n)
       ENDIF
    END DO
    RETURN
\endcode
</ol>
</ol>
All the other subroutines present in the file <tt>condlim_test_28.f90</tt> are not used in this test.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file.

<h3>Setting in the data file</h3>
We describe the data file of this test.
 It is called <tt>debug_data_test_28</tt> and can be found in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<ol>
<li>We use a formatted mesh by setting:
\code
===Is mesh file formatted (true/false)?
.t.
\endcode
<li>The path and the name of the mesh are specified with the two following lines:
\code
===Directory and name of mesh file
'.'  mesh_T28_0_04_04_ext3.FEM 
\endcode
where '.' refers to the directory where the data file is, meaning ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<li>We use one processor in the meridian section. It means the finite element mesh is not subdivised.
\code
===Number of processors in meridian section
1
\endcode
<li>We solve the problem for \f$64\f$ Fourier modes.
\code
===Number of Fourier modes
64
\endcode
<li>We use \f$8\f$ processors in Fourier space.
\code
===Number of processors in Fourier space
8
\endcode
It means that each processors is solving the problem for \f$64/8=8\f$ Fourier modes.
<li>We do not select specific Fourier modes to solve.
\code
===Select Fourier modes? (true/false)
.f.
\endcode
As a consequence, the code approximates the problem on the first \f$64\f$ Fourier modes.
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
0.01  10 !628 iterations = one turn since omega=1.0
\endcode
As the angular speed is set to one, we need 628 time steps to have a full rotation of the impellers.
<li>We set the number of domains and their label, see the files associated to the generation of the mesh, where the code approximates the Navier-Stokes equations.
\code
===Number of subdomains in Navier-Stokes mesh
7
===List of subdomains for Navier-Stokes mesh
1 2 3 4 5 6 7
\endcode
<li>We set the number of boundaries with Dirichlet conditions on the velocity field and give their respective labels.
\code
===How many boundary pieces for full Dirichlet BCs on velocity?
3
===List of boundary pieces for full Dirichlet BCs on velocity
2 10 4 
\endcode
<li>We set the kinetic Reynolds number \f$\Re\f$.
\code
===Reynolds number
20.d0
\endcode
<li>We use a penalty function function to take into account the presence of impellers.
\code
===Use penalty in NS domain (true/false)?
.t.
\endcode
<li>The solid is moving (contra rotating impellers), so we need to set:
\code
===Use nonzero velocity in solids (true/false)?
.t.
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

To check the correctness of the code, we compute four quantities:
<ol>
<li>The L2 norm of the divergence of the velocity field.
<li>The L2 norm of the divergence of the velocity field divided by the H1 norm of the veloctity field.
<li>The L2 norm of the velocity field.
<li>The H1 norm of the pressure.
</ol>
These quantities are computed at the final time \f$t=0.1\f$.
 They are compared to reference values to attest of the correctness of the code.
  
 These values of reference are in the last lines of the file <tt>debug_data_test_28</tt> in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
============================================
mesh_T28_0_04_04_ext3.FEM
===Reference results
8.59642068557655986E-002 !L2-norm on div of  u                   
2.32688350427632945E-002 !L2-norm of div of  u err/norm         
0.47840385650238321      !L2-norm  of u                               
1.7279248727889553       !H1-norm  of p
\endcode


To conclude this test, we show the profile of the approximated pressure and velocity magnitude at the final time. These figures are done in the plane \f$y=0\f$ which is the union of the half plane \f$\theta=0\f$ and \f$\theta=\pi\f$.
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_28_pre_tfin.png "Pressure in the plane plane y=0."
    </td>
    <td align="center">
    @image html  fig_test_28_vel_tfin.png  "Velocity magnitude in the plane plane y=0."
    </td>
</tr>
</table>

*/
