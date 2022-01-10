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
 * @page doc_debug_test_17 Test 17: Maxwell Equations with vacuum and variable magnetic permeability in (r,z)

<h3>Introduction</h3>
In this example, we check the correctness of SFEMaNS for stationary magnetic problem involving a conducting with a variable magnetic permeability in (r,z) and an insulating domain. The magnetic permeability in the conducting region is given. We consider Dirichlet boundary conditions. We use P2 finite elements. This set up is approximated over one time iteration with a large time step.


We solve the Maxwell equations with the magnetic field \f$\bB=\mu \bH\f$ as dependent variable:
\f{align}{
\begin{cases}
\partial_t (\mathbf{B}) + \nabla \times \left(\frac{1}{\Rm \sigma}\nabla \times( \frac{1}{\mu^c}\mathbf{B} ) \right)  =
 \nabla\times (\bu^\text{given} \times  \mathbf{B}) 
+ \nabla \times \left(\frac{1}{\Rm \sigma} \nabla\times \mathbf{j} \right) & \text{in } \Omega_1, \\
\text{div} (\mathbf{B}) = 0  & \text{in } \Omega_1 ,\\
-\partial_t \DIV (\mu^v \GRAD( \phi)) = 0 & \text{in } \Omega_2 ,\\
\mathbf{H}\times  \bn^c + \nabla \phi \times \bn^v = 0 & \text{on } \Sigma, \\
 \bB \cdot  \bn^c + \mu ^v \nabla \phi \cdot \bn^v = 0 & \text{on } \Sigma, \\
\phi = \phi_{\text{bdy}}  & \text{on } \Gamma,\\
\bH_{|t=0}= \bH_0, \\
\phi_{|t=0}= \phi _0,
\end{cases}
\f}
in the domain \f$\Omega= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,10] \times [0,2\pi) \times [-10,10] \; | \; r^2+z^2 \leq 100\} \f$.
 This domain is the union of a conducting domain
\f$\Omega_1= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,1] \times [0,2\pi) \times [-1,1]\} \f$
 and an insulating domain
\f$\Omega_2\f$ that is the complementary of \f$\Omega_1\f$ in \f$\Omega\f$.
We also set \f$\Sigma= \Omega_1 \cap \Omega_2 \f$ and \f$\Gamma= \partial \Omega_2\setminus \Sigma \f$.
 We define the outward normals \f$\bn^c, \bn^v\f$ to the surface \f$\Sigma\f$.
The data are the source term \f$\mathbf{j}\f$, the given velocity \f$\bu^\text{given}\f$,
 the boundary data \f$\phi_{\text{bdy}}\f$,
the initial datas \f$\bH_0\f$ and \f$\phi _0\f$. The parameter \f$\Rm\f$ is the magnetic Reynolds number.
 The parameter \f$\mu^c\f$, respectively \f$\mu^v\f$, is the magnetic permeability of the conducting region \f$\Omega_1\f$, respectively of the vacuum region \f$\Omega_2\f$.
 The parameter \f$\sigma\f$ is the conductivity in the conducting region \f$\Omega_1\f$.

<h3>Manufactured solutions</h3>

We define the magnetic permeability in the conducting region as follows:
@f{align*}
\mu^c(r,z)=\frac{1}{1+f_{17}(r,z)},
@f}
with
@f{align*}
f_{17}(r,z)=2^6\left(1-\frac{1}{\text{ratio}_\mu}\right) \left(r(1-r)(z^2-1)  \right)^3,
@f}
where \f$\text{ratio}_\mu=50\f$ represents the fraction of the minimum over the maximum of the magnetic permeability in the conducting region.

Tha manufactured solutions considered in this test involve two bessel functions denoted \f$BESSJ0\f$ and \f$BESSJ1\f$.They are defined for any positive real argument \f$x\f$ and satisfy the relation \f$ \partial_x BESSJ0 = - BESSJ1\f$. These functions are defined in the file <tt>bessel.f90</tt> in the following directory: ($SFEMaNS_DIR)/FEMSUB.

We approximate the following analytical solutions:
@f{align*}
H_r(r,\theta,z,t) &= -(1+f_{17}(r,z)) \cosh(z) BESSJ(1,r),
\\ H_{\theta}(r,\theta,z,t) &=0,
\\ H_z(r,\theta,z,t) &= (1+f_{17}(r,z))  \sinh(z) BESSJ(0,r) ,
\\ \phi(r,\theta,z,t) &=  \cosh(z) BESSJ(0,r),
@f}
and remind that \f$\bB=\mu\bH\f$.

 We also set the given velocity field to zero.
@f{align*}
u^\text{given}_r(r,\theta,z,t) &= 0, 
\\ u^\text{given}_{\theta}(r,\theta,z,t) &= 0, 
\\ u^\text{given}_z(r,\theta,z,t) &= 0.
@f}
The source term \f$ \mathbf{j}\f$ and the boundary data \f$\phi_{\text{bdy}}\f$ are computed accordingly.


<h3>Generation of the mesh</h3>
The finite element mesh used for this test is named <tt>mesh_17_0.1.FEM</tt>. 
 The mesh size for the P1 approximation is \f$0.1\f$ in \f$\Omega_1\f$ and \f$1\f$ on the boundary \f$\Gamma\f$. You can generate this mesh with the files in the following directory:
 ($SFEMaNS_MESH_GEN_DIR)/EXAMPLES/EXAMPLES_MANUFACTURED_SOLUTIONS/mesh_17_0.1.FEM.
The following images show the mesh for P1 finite elements for the conducting and vacuum regions.

<table width="70%" align="center">
<tr>
    <td align="center">
    @image html  fig_mesh_17_0.1.png  "Finite element mesh (P1) of the conducting region."
    </td>
    <td align="center">
    @image html  fig_mesh_17_0.1_vacuum.png  "Finite element mesh (P1) of the vacuum region."
    </td>
</tr>
</table>




<h3>Information on the file <tt>condlim.f90</tt></h3>
The initial conditions, boundary conditions and the forcing term \f$\textbf{j}\f$ in the Maxwell
 equations are set in the file <tt>condlim_test_17.f90</tt>. This condlim also uses functions of the file <tt>test_17.f90</tt> to define the function \f$f_{17}\f$, the magnetic permeability and its gradient. Here is a description of the subroutines and functions of interest of the file <tt>condlim_test_17.f90</tt>.
<ol>
<li>The subroutine <tt>init_maxwell</tt> initializes the magnetic field and the scalar potential 
 at the time\f$-dt\f$ and \f$0\f$ with \f$dt\f$ being the time step.
 It is done by using the functions Hexact, Phiexact as follows:
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
    ENDDO
\endcode
<li>The function <tt>Hexact</tt> contains the analytical magnetic field.
 It is used to initialize the magnetic field.
<ol>
<li>If the Fourier mode m is not equal to 0, we set the magnetic field to zero.
\code
    IF (m/=0) THEN
       vv(:) = 0.d0
       RETURN
    END IF
\endcode
<li>We define the function \f$f_{17}\f$ with the function <tt>f_test_T17</tt> from the file <tt>test_17.f90</tt>.
\code
    vv =  f_test_T17(rr(1,:),rr(2,:))
\endcode
<li>For the Fourier mode \f$m=0\f$, we define the magnetic field depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) as follows:
\code
    IF (TYPE==1) THEN !Bessel functions defined for reals
       vv = COSH(rr(2,:))*(1+vv)
       DO n = 1, SIZE(rr,2)
          vv(n) = -BESSJ1(rr(1,n))*vv(n)
       END DO
    ELSE IF (TYPE==5) THEN
       vv = SINH(rr(2,:))*(1+vv) 
       DO n = 1, SIZE(rr,2)
          vv(n) = BESSJ0(rr(1,n))*vv(n)
       END DO
    ELSE
       vv = 0.d0
    END IF
    RETURN
\endcode
</ol>
<li>The function <tt>Phiexact </tt> contains the analytical scalar potential.
 It is used to initialize the scalar potential and to impose Dirichlet boundary
 conditions on the scalar potential.
<ol>
<li>If the Fourier mode m is not equal to 0, we set the scalar potential to zero.
\code
   IF (m/=0) THEN
      vv = 0.d0
      RETURN
   END IF
\endcode
<li>We define the radial and vertical coordinates \f$(r, z)\f$.
\code
   r = rr(1,:)
   z = rr(2,:)
\endcode
<li>For the Fourier mode \f$m=0\f$, we define the scalar potential depending of its TYPE (1 for cosine and 2 for sine) as follows:
\code
   IF (TYPE==1) THEN !Bessel functions defined for reals
      DO n = 1, SIZE(rr,2)
         vv(n) = BESSJ0(r(n))*COSH(z(n)) 
      END DO
   ELSE
      vv = 0.d0
   END IF
   RETURN
\endcode
</ol>
<li>The function <tt>Vexact</tt> is used to define the given velocity field \f$\bu^\text{given}\f$. It is set to zero.
\code
   vv = 0.d0
\endcode
<li>The function <tt>Jexact_gauss</tt> is used to define the source term \f$\textbf{j}\f$.
 We remind that \f$ \bu^\text{given}=0 \f$ and \f$\partial_t \textbf{H}=0\f$ so it is defined such that \f$\textbf{j}=\ROT \bH\f$.
<ol>
<li>If the Fourier mode m is not equal to 0, it is set to zero.
\code
   IF (m/=0) THEN
      vv = 0.d0
      RETURN
   END IF
\endcode
<li>We define the radial and vertical coordinates \f$r, z\f$.
\code
   r = rr(1)
   z = rr(2)
\endcode
<li>For the Fourier mode \f$m=0\f$, we define the source term depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) as follows:
\code
   IF (TYPE==3) THEN
      vv = -dfdz_test_T17(r,z)*BESSJ1(r)*COSH(z) &
           -dfdr_test_T17(r,z)*BESSJ0(r)*SINH(z)
   ELSE
      vv = 0.d0
   END IF
   RETURN
\endcode
</ol>
<li>The function <tt>mu_bar_in_fourier_space</tt> defines the magnetic permeability that depends of the radial and vertical coordinates. It is done by using the function <tt>mu_bar_in_fourier_space_anal_T17</tt> of the file <tt>test_17.f90</tt>.
\code
    IF( PRESENT(pts) .AND. PRESENT(pts_ids) ) THEN
       vv=mu_bar_in_fourier_space_anal_T17(H_mesh,nb,ne,pts,pts_ids) 
    ELSE
       vv=mu_bar_in_fourier_space_anal_T17(H_mesh,nb,ne) 
    END IF
    RETURN
\endcode
When the magnetic permeability depends of the time or the azimuthal direction (set "===Is permeability variable in theta (true/false)?" to true in the data file), this function is used to define a stabilization term mu_bar. It needs to be smaller than mu everywhere in the domain. The magnetic permeability is then defined with the function <tt>mu_in_real_space</tt>.
<li>The function <tt>grad_mu_bar_in_fourier_space</tt> defines the gradient of the magnetic permeability. It is done by using the function <tt>grad_mu_bar_in_fourier_space_anal_T17</tt> of the file <tt>test_17.f90</tt>.
\code
    vv=grad_mu_bar_in_fourier_space_anal_T17(pt,pt_id) 
    RETURN
\endcode
When the magnetic permeability depends of the time or the azimuthal direction, this function is used to define the gradient of a stabilization term mu_bar.
</ol>
All the other subroutines present in the file <tt>condlim_test_17.f90</tt> are not used in this test.  We remind that the bessel function \f$BESSJ0\f$ and \f$BESSJ1\f$ are defined in the file <tt>bessel.f90</tt> in the following directory: ($SFEMaNS_DIR)/FEMSUB.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file. 
 The following describes the functions of interest of the file <tt>test_17.f90</tt>.
<ol>
<li>First we define two real numbers \f$\text{ratio}_\mu\f$ and \f$b\f$.
\code
  REAL (KIND=8), PARAMETER, PRIVATE :: ratio_mu_T17 = 50.d0 ! the variation of mu
  REAL (KIND=8), PRIVATE            :: b_factor_T17_anal = (2**6) * (1.d0 - 1/ratio_mu_T17)
\endcode
<li>The function <tt>f_test_T17</tt> defines the function \f$f_{17}\f$.
\code
    vv = b_factor_t17_anal*(r*(1-r)*(z**2-1))**3
    RETURN
\endcode
<li>The function <tt>dfdr_test_T17</tt> defines the r-derivative of \f$f_{17}\f$.
\code
    vv = 3 * b_factor_t17_anal * (z**2-1)**3 * (r*(1-r))**2 * (1-2*r)
    RETURN
\endcode
<li>The function <tt>dfdz_test_T17</tt> defines the z-derivative of \f$f_{17}\f$.
\code
    vv = 3*b_factor_t17_anal*(r*(1-r))**3*(z**2-1)**2*(2*z)
    RETURN
\endcode
<li>The function <tt>mu_bar_in_fourier_space_anal_T17</tt> defines the magnetic permeability.
\code
    IF( PRESENT(pts) .AND. PRESENT(pts_ids) ) THEN !Computing mu at  pts
       r=pts(1,nb:ne)
       z=pts(2,nb:ne)
    ELSE
       r=H_mesh%rr(1,nb:ne) !Computing mu  at nodes
       z=H_mesh%rr(2,nb:ne)
    END IF
    vv=1.d0/(1.d0+f_test_T17(r,z))
    RETURN
\endcode
<li>The function <tt>grad_mu_bar_in_fourier_space_anal_T17</tt> defines the gradient of the magnetic permeability on gauss points.
\code
    r=pt(1)
    z=pt(2)
    tmp=(1.d0 +f_test_T17(r,z))**2
    vv(1)=-dfdr_test_T17(r(1),z(1))/tmp(1)
    vv(2)=-dfdz_test_T17(r(1),z(1))/tmp(1)
    RETURN
\endcode
</ol>








<h3>Setting in the data file</h3>
We describe the data file of this test.
 It is called <tt>debug_data_test_17</tt> and can be found in the following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<ol>
<li>We use a formatted mesh by setting:
\code
===Is mesh file formatted (true/false)?
.t.
\endcode
<li>The path and the name of the mesh are specified with the two following lines:
\code
===Directory and name of mesh file
'.' 'mesh_17_0.1.FEM'
\endcode
where '.' refers to the directory where the data file is, meaning ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<li>We use one processor in the meridian section. It means the finite element mesh is not subdivised.
\code
===Number of processors in meridian section
1
\endcode
<li>We solve the problem for \f$1\f$ Fourier modes.
\code
===Number of Fourier modes
1
\endcode
<li>We use \f$1\f$ processors in Fourier space.
\code
===Number of processors in Fourier space
1
\endcode
It means that each processors is solving the problem for \f$1/1=1\f$ Fourier modes.
<li>We do not select specific Fourier modes to solve.
\code
===Select Fourier modes? (true/false)
.f.
\endcode
As a consequence, the code approximates the problem on the first Fourier mode \f$m=0\f$.
<li>We approximate the Maxwell equations by setting:
\code
===Problem type: (nst, mxw, mhd, fhd)
'mxw'
\endcode
<li>We do not restart the computations from previous results.
\code
===Restart on velocity (true/false)
.f.
===Restart on magnetic field (true/false)
.f.
\endcode
 It means the computation starts from the time \f$t=0\f$.
<li>We use a time step of \f$10^{40}\f$ and solve the problem over \f$1\f$ time iteration.
\code  
===Time step and number of time iterations
1d40, 1 
\endcode
<li>We use the magnetic field \f$\textbf{B}\f$ as dependent variable for the Maxwell equations.
\code
===Solve Maxwell with H (true) or B (false)?
.f.
\endcode
This parameter is set to false by default.
<li>We set the number of domains and their label, see the files associated to the generation of the mesh, where the code approximates the Maxwell equations.
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
<li>We set the number of boundaries with Dirichlet conditions on the magnetic field.
\code
===Number of Dirichlet sides for Hxn
0
\endcode
<li>The magnetic permeability is defined with a function of the file <tt>condlim_test_17.f90</tt>.
\code
===Is permeability defined analytically (true/false)?
.t.
\endcode
<li>We construct the magnetic permeability on the gauss points without using its value on the nodes and a finite element interpolation.
\code
===Use FEM Interpolation for magnetic permeability (true/false)?
.f.
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
<li>We set the number of domains and their label, see the files associated to the generation of the mesh,
 where the code approximates the scalar potential.
\code
===Number of subdomains in magnetic potential (phi) mesh
1
===List of subdomains for magnetic potential (phi) mesh
2
\endcode
<li>We set the number of boundaries with Dirichlet conditions on the scalar potential and give their respective labels.
\code
===How many boundary pieces for Dirichlet BCs on phi?
1
===List of boundary pieces for Dirichlet BCs on phi
3
\endcode
<li>We set the number of interface between the magnetic field and the scalar potential and give their respective labels.
\code
===Number of interfaces between H and phi
1
===List of interfaces between H and phi
2
\endcode
<li>We set the magnetic permeability in the vacuum domain.
\code
===Permeability in vacuum
1
\endcode
<li>We set the type of finite element used to approximate the scalar potential.
\code
===Type of finite element for scalar potential
2
\endcode
<li>We set a stabilization coefficient used to penalize term across the interface between the magnetic field and the scalar potential.
\code
===Stabilization coefficient (interface H/phi)
1.d0
\endcode
We note this coefficient is usually set to \f$1\f$.
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
</ol>


<h3> Outputs and value of reference </h3>

The outputs of this test are computed with the file <tt>post_processing_debug.f90</tt> 
that can be found in the following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. 

To check the well behavior of the code, we compute four quantities:
<ol>
<li>The L2 norm of the error on the divergence of the magnetic field \f$\textbf{B}=\mu\textbf{H}\f$.
<li>The L2 norm of the error on the curl of the magnetic field \f$\textbf{H}\f$.
<li>The L2 norm of the error on the magnetic field \f$\textbf{H}\f$.
<li>The L2 norm of the error on the scalar potential.
</ol>

 They are compared to reference values to attest of the correctness of the code.  
 These values of reference are in the last lines of the file <tt>debug_data_test_17</tt> in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
============================================
mesh_17_0.1.FEM, dt=1d40, it_max=1
===Reference results
3.658447514681710E-003 !L2 norm of div(mu (H-Hexact))
4.500377040641551E-004 !L2 norm of curl(H-Hexact)
5.441900368783823E-003 !L2 norm of H-Hexact
1.586225425289322E-003 !L2 norm of phi-phiexact
\endcode


To conclude this test, we show the profile of the approximated  magnetic field \f$\textbf{H}\f$ and
 scalar potential \f$\phi\f$ at the final time.
 These figures are done in the plane \f$y=0\f$ which is 
 the union of the half plane \f$\theta=0\f$ and \f$\theta=\pi\f$.
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_17_mag_tfin.png "Magnetic field magnitude in the plane plane y=0."
    </td>
    <td align="center">
    @image html  fig_test_17_phi_tfin.png  "Scalar potential in the plane plane y=0."
    </td>
</tr>
</table>


*/
