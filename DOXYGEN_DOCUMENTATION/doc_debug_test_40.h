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
 * @page doc_debug_test_40 Test 40: Navier-Stokes nonaxisymetric with variable density and viscosity

<h3>Introduction</h3>
In this example, we check the correctness of SFEMaNS for a hydrodynamic problem with variable density and dynamical viscosity involving Dirichlet boundary conditions. The variables are all made dependant of the azimuthal component (nonaxisymmetric setting). The mass equation is not approximated. We consider a level set, solution of the same advection equation, that is used to reconstruct the density and the other fluid's properties.

We solve the level set equation:
@f{align*}
\partial_t \varphi + \bu \cdot \GRAD \varphi = f_\varphi.
@f}

We recontruct the density and dynamical viscosity as follows:
@f{align*}
\rho=\rho_1+(\rho_2-\rho_1) F(\varphi),
\\ \eta=\eta_1+(\eta_2-\rho_1) F(\varphi),
@f}
with \f$\eta_i\f$ and \f$\rho_i\f$ data to define. The function \f$F(\varphi)\f$ is either the identity (called linear reconstruction) or a piece-wise polynomial function (called reg reconstruction).

We solve the Navier-Stokes equations with the momentum \f$\bm=\rho\bu\f$ as dependent variable:
@f{align*}
\partial_t\bm + \DIV(\bm\otimes\bu) - \frac{2}{\Re} \DIV(\eta \epsilon(\bu)) +\GRAD p     &=\bef,
\\ \DIV \bu &= 0, \\
\bu_{|\Gamma} &= \bu_{\text{bdy}} , \\
\bu_{|t=0} &= \bu_0, \\
p_{|t=0} &= p_0,
@f}
where \f$\epsilon(\bu)=\frac{1}{2}(\GRAD \bu + \GRAD^s \bu)\f$ is the strain rate tensor.

These equations are solved in the domain \f$\Omega= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,1] \times [0,2\pi) \times [-1,1]\} \f$ with  \f$\Gamma= \partial \Omega \f$. 
The data are the source terms \f$f_\varphi\f$ and \f$\bef\f$, the boundary data \f$\bu_{\text{bdy}}\f$,
the initial datas \f$\bu_0\f$ and \f$p_0\f$. The parameter \f$\Re\f$ is the kinetic Reynolds number, the densities \f$\rho_i\f$ and the dynamical viscosities \f$\eta_i\f$.

Remarks:
<ol>
<li>The level set and the momentum equations can be stabilized with the entropy viscosity method used in test 15 (called LES).
<li>For physical problem, the level set has to take values in [0,1] such that the interface is represented by \f$\varphi^{-1}(\{1/2\})\f$. The fluids area are respectively represented by \f$\varphi^{-1}(\{0\})\f$ and \f$\varphi^{-1}(\{1\})\f$. Notice that this test does not consider immiscible fluids. It involves manufactured solution with a smooth variable density.
<li>A compression term can also be added in the level set equation. This term allows the level set to remains sharp near the fluids interface represented by \f$ \varphi^{-1}(\{1/2\}) \f$.
</ol>
 We refer to the section \ref doc_intro_SFEMaNS_possibilities_nst_4 for more details on the algorithms implemented in SFEMaNS for multiphase problem.

<h3>Manufactured solutions</h3>
We approximate the following analytical solutions:
@f{align*}
u_r(r,\theta,z,t) &=  r \cos(t)\sin(2\theta)/2,
\\ u_{\theta}(r,\theta,z,t) &= 3 r \cos(t)/2 + r \cos(t) \cos(2\theta)/2,
\\ u_z(r,\theta,z,t) &=1 ,
\\ p(r,\theta,z,t) &=  r^2 z^3 \cos(t) + r^2 \cos(t-z)\sin(\theta) + r z \sin(t-r)\cos(2\theta),
\\ \varphi(r,\theta,z,t) &= 0.5 + (\cos(t-z) + \sqrt{2}\cos(\sqrt{2}\sin(t)) r \cos(\theta) 
                                   +\sin(\sqrt{2}\sin(t)) r \sin(\theta))/8,
\\ \rho(r,\theta,z,t) &=1 + 99\varphi(r,\theta,z,t),
\\ \eta(r,\theta,z,t) &=1 + 99\varphi(r,\theta,z,t).
@f}
where the source terms \f$f_\varphi\f$, \f$\bef\f$ and the boundary data \f$ \bu_{\text{bdy}}\f$ are computed accordingly.

<h3>Generation of the mesh</h3>

The finite element mesh used for this test is named <tt>RECT20.FEM</tt> and 
 has a mesh size of \f$0.05\f$ for the P1 approximation. You can generate this mesh with the files in the following directory:
 ($SFEMaNS_MESH_GEN_DIR)/EXAMPLES/EXAMPLES_MANUFACTURED_SOLUTIONS/RECT10_form after modifying the topology file.
The following image shows the mesh for P1 finite elements.
<table width="70%" align="center">
<tr>
    <td align="center">
    @image html  fig_RECT_20.png  "Finite element mesh (P1)."
    </td>
</tr>
</table>


<h3>Information on the file <tt>condlim.f90</tt></h3>

The initial conditions, boundary conditions and the forcing terms are set in the file <tt>condlim_test_40.f90</tt>.
Here is a description of the subroutines and functions of interest.
<ol>
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
<li>The function <tt>init_level_set</tt> initializes the level set at the time\f$-dt\f$ and \f$0\f$ with \f$dt\f$ being the time step. This is done by using the function level_set_exact as follows:
\code
    time = 0.d0
    DO i= 1, SIZE(list_mode)
       mode = list_mode(i) 
       DO j = 1, 2
          !===level_set
          DO n = 1, inputs%nb_fluid -1
             level_set_m1(n,:,j,i) = level_set_exact(n,j,pp_mesh%rr,mode,time-dt)  
             level_set   (n,:,j,i) = level_set_exact(n,j,pp_mesh%rr,mode,time)
          END DO
       END DO
    END DO
\endcode
We note there is one level set per interface and the different phases are stratified. So for nb_fluid fluids, we have inputs nb_fluid-1 interfaces.
<li>The function <tt>vv_exact</tt> contains the analytical velocity field.
 It is used to initialize the velocity field and to impose Dirichlet boundary
 conditions on the velocity field.
<ol>
<li>First we define the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>We define the velocity field depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) as follows:
\code
    vv = 0.d0
    IF (m==0) THEN
       IF(TYPE==3) THEN
          vv = 1.5d0*r*COS(t)
       ELSE IF (TYPE==5) THEN
          vv = 1.d0
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (m==2) THEN
       IF (TYPE==2) THEN
          vv = 0.5d0*r*COS(t)
       ELSE IF (TYPE==3) THEN
          vv = 0.5d0*r*COS(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE
       vv = 0.d0
    END IF
    RETURN
\endcode
where \f$t\f$ is the time.
</ol>
<li>The function <tt>pp_exact</tt> contains the analytical pressure.
It is defined as follows.
\code
    IF (m==0.AND.TYPE==1) THEN
       vv(:) = rr(1,:)**2*rr(2,:)**3*COS(t)
    ELSE IF (m==1 .AND. TYPE==2) THEN
       vv = rr(1,:)**2*COS(t-rr(2,:))
    ELSE IF (m==2 .AND. TYPE==1) THEN
       vv = rr(1,:)*rr(2,:)*SIN(t-rr(1,:))
    ELSE
       vv = 0.d0
    END IF
    RETURN
\endcode
<li>The function <tt>level_set_exact</tt> is used to initialize the level set.
<ol>
<li>We define the level set of the mode and its TYPE (1 for cosine and 2 for sine) as follows:
\code
    IF (interface_nb==1) THEN
       IF (m==0 .AND. TYPE==1) THEN
          vv = 0.5d0 + 0.125d0*COS(t-z)
       ELSE IF (m==1 .AND. TYPE==1) THEN
          vv = 0.125d0*SQRT(2.d0)*COS(SQRT(2.d0)*SIN(t))*r
       ELSE IF (m==1 .AND. TYPE==2) THEN
          vv = 0.125d0*SIN(SQRT(2.d0)*SIN(t))*r
       ELSE
          vv = 0.d0
       END IF
\endcode
<li>If more than one level set is considered, the computation is stopped.
\code
    ELSE 
       CALL error_petsc(' BUG in level_set_exact, we should compute only 1 level set')
    END IF
    RETURN
\endcode
Indeed with two fluids there is one interface. So we should compute only one level set.
<li>The function <tt>source_in_NS_momentum</tt> computes the source term \f$\bef\f$ of the Navier-Stokes equations.
<li>The function <tt>source_in_level_set</tt> computes the source term \f$f_\varphi\f$ of the level set equations. It is equal to zero.
</ol>
All the other subroutines present in the file <tt>condlim_test_40.f90</tt> are not used in this test.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file.


<h3>Setting in the data file</h3>
We describe the data file of this test.
 It is called <tt>debug_data_test_40</tt> and can be found in the following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<ol>
<li>We use a formatted mesh by setting:
\code
===Is mesh file formatted (true/false)?
.t.
\endcode
<li>The path and the name of the mesh are specified with the two following lines:
\code
===Directory and name of mesh file
'.' 'RECT20.FEM'
\endcode
where '.' refers to the directory where the data file is, meaning ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<li>We use four processor in the meridian section. It means the finite element mesh is subdivised in four.
\code
===Number of processors in meridian section
4
\endcode
<li>We solve the problem for \f$4\f$ Fourier modes.
\code
===Number of Fourier modes
4
\endcode
<li>We use \f$4\f$ processors in Fourier space.
\code
===Number of processors in Fourier space
4
\endcode
It means that each processors is solving the problem for \f$4/4=1\f$ Fourier modes.
<li>We do not select specific Fourier modes to solve.
\code
===Select Fourier modes? (true/false)
.f.
\endcode
As a consequence, the code approximates the problem on the first \f$4\f$ Fourier modes.
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
<li>We use a time step of \f$0.00125\f$ and solve the problem over \f$400\f$ time iterations.
\code
===Time step and number of time iterations
0.00125d0, 400
\endcode
<li>We do not apply mass correction on the level set.
\code
===Do we apply mass correction? (true/false)
.f.
\endcode
The default value is true.
<li>We do not kill the overshoot of the level set with respect of the interval \f$[0,1]\f$.
\code
===Do we kill level set overshoot? (true/false)
.f.
\endcode
The default value if false so these two lines are not required.
<li> We use P2 finite element to approximate the level set.
\code
===Do we use P2 finite element for level_set? (true/false)
.t.
\endcode
The default value is true.
<li>We use a BDF1 approximation of the time derivative in the momentum equation.
\code
===Do we solve momentum with bdf2 (true/false)?
.f.
\endcode
The default value is false.
<li>We use a BDF2 approximation of the time derivative in the level set equation.
\code
===Do we solve level set with bdf2 (true/false)?
.t.
\endcode
The default value is false.
<li>We use an artificial compression technique to approximate the Navier-Stokes equations.
\code
===Solve Navier-Stokes with art comp scheme (true) or (false)?
.t.
\endcode
The default value is false. This method can only be used when the momentum is used
to approximate the Navier-Stokes equations.
<li>We set the penalty coeffcient for the artifical compression method.
\code
===Penalty coefficient for artifical compression
1.0d0
\endcode
The default value is one.
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
2 4 5
\endcode
<li>We use the momentum as dependent variable for the Navier-Stokes equations.
\code
===Solve Navier-Stokes with u (true) or m (false)?
.f.
\endcode
If the density or the viscosity are variable, this parameter needs to be false. 
The default value is true (constant density and viscosity).
<li>We set the kinetic Reynolds number \f$\Re\f$.
\code
===Reynolds number
250.d0
\endcode
<li>We use the entropy viscosity method to stabilize the level set equation.
\code
===Use LES? (true/false)
.t.
\endcode
<li>We use the entropy viscosity method to stabilize the momentum equation.
\code
==Use LES in momentum? (true/false)
.t.
\endcode
<li>We define the coefficient \f$c_\text{e}\f$ of the entropy viscosity.
\code
===Coefficient multiplying residual
1.0d0
\endcode
<li>We give information on how to solve the matrix associated to the time marching of the velocity (or momentum in this case).
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
<li>We solve the level set equation.
\code
===Is there a level set?
.t.
\endcode
<li>We want to consider one level set \f$\varphi\f$, so we set:
\code
===How many fluids?
2
\endcode
We note this test does not consider two immiscible fluids.
<li>We use compression the compression technique to get a sharp interface between immiscible fluids.
\code
===Compression factor for level set
1.d0
\endcode
<li>We define the parameters \f$(\rho_1,\rho_2)\f$ used to reconstruct the density.
\code
===Density of fluid 0, fluid 1, ...
1.d0 100.d0
\endcode
<li>We define the parameters \f$(\eta_1,\eta_2)\f$ used to reconstruct the dynamical viscosity.
\code
===Dynamic viscosity of fluid 0, fluid 1, ...
1.d0 100.d0
\endcode
<li>We define a multiplier coefficient.
\code
===multiplier for h_min for level set
1.d0
\endcode
This multiplier times the smallest mesh size is stored in the variable <tt>inputs\%h_min_distance</tt>. It can be used in the condlim file to set the wideness of the initial interface. It is not used in this case as the level set does not represent an interface between two immiscible fluids.
<li>We use a linear reconstruction, meaning \f$F(\varphi)=\varphi\f$.
\code
===How are the variables reconstructed from the level set function? (lin, reg)
'lin'
\endcode
<li>We  set the number of boundaries with Dirichlet conditions on the level set and give their respective labels.
\code
===How many boundary pieces for Dirichlet BCs on level set?
1
===List of boundary pieces for Dirichlet BCs on level set
4
\endcode
<li>We give information on how to solve the matrix associated to the time marching of the level set.
<ol>
<li>
\code
===Maximum number of iterations for level set solver
100
\endcode
<li>
\code
===Relative tolerance for level set solver
1.d-6
===Absolute tolerance for level set solver
1.d-10
\endcode
<li>
\code
===Solver type for level set (FGMRES, CG, ...)
GMRES
===Preconditionner type for level set solver (HYPRE, JACOBI, MUMPS...)
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
</ol>


<h3> Outputs and value of reference </h3>

The outputs of this test are computed with the file <tt>post_processing_debug.f90</tt> 
that can be found in the following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. 

To check the well behavior of the code, we compute four quantities:
<ol>
<li>The L2 norm of the error on the velocity field.
<li>The H1 norm of the error on the velocity field.
<li>The L2 norm of the error on the level set.
<li>The L2 norm of the error on the pressure.
</ol>
These quantities are computed at the final time \f$t=0.5\f$.
 They are compared to reference values to attest of the correctness of the code.  
 These values of reference are in the last lines of the file <tt>debug_data_test_40</tt> in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
============================================
(RECT20.FEM)
===Reference results
8.7749796248369450E-004  L2 error on velocity
1.0482305912429532E-002  H1 error on velocity
2.6353807435704957E-005  L2 error on level set
8.1493730051115043E-002  L2 error on pressure
\endcode


To conclude this test, we show the profile of the approximated level set, density, pressure and velocity magnitude at the final time. These figures are done in the plane \f$y=0\f$ which is the union of the half plane \f$\theta=0\f$ and \f$\theta=\pi\f$.
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_40_level_set_tfin.png "Level set in the plane plane y=0."
    </td>
    <td align="center">
    @image html  fig_test_40_density_tfin.png  "Density in the plane plane y=0."
    </td>
</tr>
<tr>
    <td align="center">
     @image html  fig_test_40_pre_tfin.png "Pressure in the plane plane y=0."
    </td>
    <td align="center">
    @image html  fig_test_40_vel_tfin.png  "Velocity magnitude in the plane plane y=0."
    </td>
</tr>
</table>
 */
