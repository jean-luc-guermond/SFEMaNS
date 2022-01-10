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
 * @page doc_debug_test_13 Test 13: MHD periodic and Neumann bdy

<h3>Introduction</h3>
In this example, we check the correctness of SFEMaNS for a magnetohydrodynamic problem. This test uses Dirichlet, Neumann and periodic boundary conditions. We note this test does not involve manufactured solutions and consist of checking four quantities, like the \f$\bL^2\f$ norm of the velocity, are the same than the values of reference.

We solve the Navier-Stokes equations:
@f{align*}
\partial_t\bu+\left(\ROT\bu + 2 \epsilon \textbf{e}_z \right)\CROSS\bu - \frac{1}{\Re}\LAP \bu +\GRAD p     
&= (\ROT \textbf{H}) \times (\mu^c \textbf{H}),
\\ \DIV \bu &= 0, \\
\bu_{|\{z=0\}} &= \bu_{|\{z=1\}} , \\
\bu_{|\Gamma} &= \bu_{\text{bdy}} , \\
\bu_{|t=0} &= \bu_0, \\
p_{|t=0} &= p_0,
@f}
and the Maxwell equations:
@f{align*}
\partial_t (\mu^c \mathbf{H}) + \nabla \times \left(\frac{1}{\Rm \sigma}\nabla \times \mathbf{H}  \right) 
& = \nabla\times (\bu \times \mu^c \mathbf{H}), \\
\text{div} (\mu^c \mathbf {H}) &= 0   ,\\
\mathbf{H}_{|\{z=0\}} &= \mathbf{H}_{|\{z=1\}} , \\
 \left( \frac{1}{\Rm \sigma} \left( \ROT (\mathbf{H}) - \mathbf{j}  \right) - \bu \times \mu \mathbf{H}
 \right)  \times \bn_{|\Gamma} & = {\textbf{a} \times \bn}_{|\Gamma},\\
\bH_{|t=0}&= \bH_0.
@f}
Theses equations are solved in the domain \f$\Omega= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,1/2] \times [0,2\pi) \times [0,1]\} \f$ with  \f$\Gamma= \partial \Omega\setminus\{ \{z=0\} \cup \{z=1\} \}  \f$.
The data are the boundary datas \f$\bu_{\text{bdy}}, \textbf{a}\f$,
 the initial data \f$\bu_0, p_0, \bH_0\f$.
 The parameters are the kinetic Reynolds number \f$\Re\f$, the magnetic Reynolds number \f$\Rm\f$, the magnetic permeability \f$\mu^c\f$ and the conductivity  \f$\sigma\f$ of the fluid.


<h3>Manufactured solutions</h3>
As mentionned earlier this test does not involve manufactured solutions. As a consequence, we do not consider specific source term and only initialize the variables to approximate.

The initial velocity field and pressure are set to zero.
@f{align*}
 u_r(r,\theta,z,t=0) &= 0.5 - r, 
\\ u_{\theta}(r,\theta,z,t=0) &= (r-0.5)r\sin(2\pi z), 
\\ u_z(r,\theta,z,t=0) &= 0, 
\\ p(r,\theta,z,t=0) &= 0.
@f}

The magnetic field is iniatialized as follows:
@f{align*}
\\ H_r(r,\theta,z,t=0) &= 0,
\\ H_{\theta}(r,\theta,z,t=0) &= r,
\\ H_z(r,\theta,z,t=0) &= 1 + r(r-0.5)\left(\cos(\theta)+\sin(\theta)+ \cos(2\theta)+\sin(2\theta)\right).
@f}
The boundary datas \f$\bu_{\text{bdy}}, \textbf{a}\f$ are set to zero.


<h3>Generation of the mesh</h3>

The finite element mesh used for this test is named <tt>Mesh_20_form.FEM</tt>. 
 The mesh size is \f$0.05\f$ for the P1 approximation. You can generate this mesh with the files in the following directory:
 ($SFEMaNS_MESH_GEN_DIR)/EXAMPLES/EXAMPLES_MANUFACTURED_SOLUTIONS/Mesh_20_form.
The following image shows the mesh for P1 finite elements.
<table width="70%" align="center">
<tr>
    <td align="center">
    @image html  fig_Mesh_20_form.png  "Finite element mesh (P1)."
    </td>
</tr>
</table>


<h3>Information on the file <tt>condlim.f90</tt></h3>
The initial conditions, boundary conditions and the forcing terms are set in the file <tt>condlim_test_13.f90</tt>.
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
<li>We construct the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>For a time stritly larger than 0, we set the Dirichlet boundary condition to zero.
\code
    IF (t>1.d-14) THEN
       vv = 0.d0
\endcode
<li>For the initialization, if the Fourier mode m is not equal to 0, the velocity field is set to zero.
\code
    ELSE
       IF (m/=0) THEN
          vv = 0.d0
\endcode
<li>For the Fourier mode \f$m=0\f$, we define the velocity field depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) as follows:
\code
       ELSE
          IF (TYPE==1) THEN
             vv = 0.5-r
          ELSE IF (TYPE==3) THEN
             vv = (r-0.5)*r*SIN(2*PI*z)
          ELSE
             vv = 0.d0
          END IF
       END IF
    END IF
    RETURN
\endcode
</ol>
<li>The function <tt>pp_exact</tt> contains the analytical pressure.
 It is used to initialize the pressure and is set to zero.
\code
    vv=0.d0
    RETURN
\endcode
<li>The function <tt>source_in_NS_momentum</tt> is used to define the source term \f$\bef\f$ of the Navier-Stokes equations. As this term is not used for this test, it is set to zero.
\code
    vv=0.d0
    RETURN
\endcode
<li>The subroutine <tt>init_maxwell</tt> initializes the magnetic field at the time\f$-dt\f$ and \f$0\f$ with \f$dt\f$ being the time step.
 It is done by using the functions Hexact as follows:
\code
    time = -dt
    DO k=1,6
       DO i=1, SIZE(list_mode)
          Hn1(:,k,i) = Hexact(H_mesh,k, H_mesh%rr, list_mode(i), mu_H_field, time)
          IF (inputs%nb_dom_phi>0) THEN
             IF (k<3) THEN
                phin1(:,k,i) = Phiexact(k, phi_mesh%rr, list_mode(i) , mu_phi, time)
             ENDIF
          END  IF
       ENDDO
    ENDDO

    time = time + dt 
    DO k=1,6
       DO i=1, SIZE(list_mode)
          Hn(:,k,i) = Hexact(H_mesh,k, H_mesh%rr, list_mode(i), mu_H_field, time)
          IF (inputs%nb_dom_phi>0) THEN
             IF (k<3) THEN
                phin(:,k,i) = Phiexact(k, phi_mesh%rr, list_mode(i), mu_phi, time)
             ENDIF
          END  IF
       ENDDO
\endcode
We note there is no scalar potential in this test so \f$\text{inputs%nb_dom_phi}=0\f$.
<li>The function <tt>Hexact</tt> contains the analytical magnetic field.
 It is used to initialize the magnetic field.
<ol>
<li>For the Fourier mode \f$m\in\{1,2\}\f$, we define the magnetic field depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) as follows:
\code
    IF (m/=0) THEN
       IF (TYPE==5 .OR. TYPE==6) THEN
          vv = rr(1,:)*(rr(1,:)-0.5)
       ELSE
          vv = 0.d0
       END IF
\endcode
<li>For the Fourier mode \f$m=0\f$, we define the magnetic field depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) as follows:
\code
    ELSE
       IF (TYPE==3) THEN
          vv = rr(1,:)
       ELSE IF (TYPE==5) THEN
          vv = 1.d0
       ELSE
          vv = 0.d0
       END IF
    END IF
    RETURN
\endcode
</ol>
<li>The function <tt>Jexact_gauss</tt> is used to define the source term \f$\textbf{j}\f$ that is not used in this test. So we set it to zero.
\code
   vv=0.d0
   RETURN
\endcode
<li>The function <tt>Eexact_gauss</tt> is used to define the boundary data \f$\textbf{a}\f$. It is set to zero.
\code
   vv=0.d0
   RETURN
\endcode
</ol>
All the other subroutines present in the file <tt>condlim_test_13.f90</tt> are not used in this test.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file.




<h3>Setting in the data file</h3>
We describe the data file of this test.
 It is called <tt>debug_data_test_13</tt> and can be found in the following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<ol>
<li>We use a formatted mesh by setting:
\code
===Is mesh file formatted (true/false)?
.t.
\endcode
<li>The path and the name of the mesh are specified with the two following lines:
\code
===Directory and name of mesh file
'.' 'Mesh_20_form.FEM'
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
<li>We approximate the Maxwell and the Navier-Stokes equations by setting:
\code
===Problem type: (nst, mxw, mhd, fhd)
'mhd'
\endcode
<li>We do not restart the computations from previous results.
\code
===Restart on velocity (true/false)
.f.
===Restart on magnetic field (true/false)
.f.
\endcode
It means the computation starts from the time \f$t=0\f$.
<li>We use a time step of \f$0.02\f$ and solve the problem over \f$10\f$ time iterations.
\code
===Time step and number of time iterations
2d-2, 20
\endcode
<li>We set periodic boundary condition.
<ol>
<li>We set the number of pair of boundaries that has to be periodic.
\code
===How many pieces of periodic boundary?
1
\endcode
<li>We give the label of the boundaries and the vector that lead to the first boundary to the second one.
\code
===Indices of periodic boundaries and corresponding vectors
4 2 .0d0 1.d0
\endcode
We note that we need as much as lines as the number of pairs of boundaries with periodic condition.
</ol>
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
1
===List of boundary pieces for full Dirichlet BCs on velocity
5
\endcode
<li>We set the kinetic Reynolds number \f$\Re\f$.
\code
===Reynolds number
1.d1
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
<li>We set the number of domains and their label, see the files associated to the generation of the mesh,
 where the code approximates the magnetic field.
\code
===Number of subdomains in magnetic field (H) mesh
1
===List of subdomains for magnetic field (H) mesh
1
\endcode
<li>We set the number of interface in H_mesh.
\code
===Number of interfaces in H mesh
0
\endcode
Such interfaces represent interfaces with discontinuities in magnetic permeability or interfaces between the magnetic field mesh and the temperature or the velocity field meshes.
<li>We set the number of boundaries with Dirichlet conditions on the magnetic field and give their respective labels.
\code
===Number of Dirichlet sides for Hxn
0
===List of Dirichlet sides for Hxn
0
\endcode
We note that the lines regarding the list of Dirichlet sides are not read as their is no such sides.
<li>We set the magnetic permeability in each domains where the magnetic field is approximated.
\code
===Permeability in the conductive part (1:nb_dom_H)
1.d0
\endcode
<li>We set the conductivity in each domains where the magnetic field is approximated.
\code
===Conductivity in the conductive part (1:nb_dom_H)
1.d0
\endcode
<li>We set the type of finite element used to approximate the magnetic field.
\code
===Type of finite element for magnetic field
2
\endcode
<li>We set the magnetic Reynolds number \f$\Rm\f$.
\code
===Magnetic Reynolds number
1.d0
\endcode
<li>We set stabilization coefficient for the divergence of the magnetic field and the penalty of the Dirichlet and interface terms.
\code
===Stabilization coefficient (divergence)
1.d0
===Stabilization coefficient for Dirichlet H and/or interface H/H
1.d0
\endcode
We note these coefficients are usually set to \f$1\f$.
<li>We set the number of domains and their label, see the files associated to the generation of the mesh,
 where the code approximates the scalar potential.
\code
===Number of subdomains in magnetic potential (phi) mesh
0
\endcode
<li>We give information on how to solve the matrix associated to the time marching of the Maxwell equations.
<ol>
<li>
\code
===Maximum number of iterations for Maxwell solver
100
\endcode
<li>
\code
===Relative tolerance for Maxwell solver
1.d-6 
===Absolute tolerance for Maxwell solver
1.d-10
\endcode
<li>
\code
===Solver type for Maxwell (FGMRES, CG, ...)
GMRES
===Preconditionner type for Maxwell solver (HYPRE, JACOBI, MUMPS...)
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
<li>The H1 norm of the velocity.
<li>The L2 norm of the magnetic field divergence.
<li>The L2 norm of the magnetic field.
<li>The L2 norm of the pressure.
</ol>
These quantities are computed at the final time \f$t=0.2\f$.
 They are compared to reference values to attest of the correctness of the code.  
 These values of reference are in the last lines of the file <tt>debug_data_test_13</tt> in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
============================================
Mesh_20_form.FEM, P2
===Reference results
3.506833380349648E-002      H1 norm on velocity     
3.720369285322975E-006      L2 norm of div(Hn)  
0.886235556266004           L2 norm of Hn
2.313775787175324E-003      L2 norm of pressure
\endcode


To conclude this test, we show the profile of the approximated pressure, velocity magnitude
 and magnetic field magnitude at the final time.
 These figures are done in the plane \f$y=0\f$ which
 is the union of the half plane \f$\theta=0\f$ and \f$\theta=\pi\f$.
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_13_pre_tfin.png "Pressure in the plane plane y=0."
    </td>
    <td align="center">
    @image html  fig_test_13_vel_tfin.png  "Velocity magnitude in the plane plane y=0."
    </td>
</tr>
</table>
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_13_mag_tfin.png "Magnetic field magnitude in the plane plane y=0."
    </td>
</tr>
</table>
 */
