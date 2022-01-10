1// ---------------------------------------------------------------------
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
 * @page doc_debug_test_04 Test 4: Maxwell with vacuum periodic P1-P2

<h3>Introduction</h3>
In this example, we check the correctness of SFEMaNS for a magnetic problem involving a conducting and an insulating domain. We consider periodic and Dirichlet boundary conditions.

We solve the Maxwell equations:
\f{align}{
\begin{cases}
\partial_t (\mu^c \mathbf{H}) + \nabla \times \left(\frac{1}{\Rm \sigma}\nabla \times \mathbf{H}  \right)  =
 \nabla\times (\bu^\text{given} \times \mu^c \mathbf{H}) 
+ \nabla \times \left(\frac{1}{\Rm \sigma} \nabla\times \mathbf{j} \right) & \text{in } \Omega_1, \\
\text{div} (\mu^c \mathbf {H}) = 0  & \text{in } \Omega_1 ,\\
- \partial_t \DIV (\mu^v \GRAD( \phi)) = 0 & \text{in } \Omega_2 ,\\
\mathbf{H}_{|\{z=0\}} = \mathbf{H}_{|\{z=1\}}, & \\
\phi_{|\{z=0\}} = \phi_{|\{z=1\}}, & \\
\bH \times  \bn^c + \nabla \phi \times \bn^v = 0 & \text{on } \Sigma, \\
\mu^c \bH \cdot  \bn^c + \mu ^v \nabla \phi \cdot \bn^v = 0 & \text{on } \Sigma, \\
\phi = \phi_{\text{bdy}}  & \text{on } \Gamma_2,\\
\bH_{|t=0}= \bH_0, \\
\phi_{|t=0}= \phi _0,
\end{cases}
\f}
in the domain \f$\Omega= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,1] \times [0,2\pi) \times [0,1]\} \f$.
 This domain is the union of a conducting domain
\f$\Omega_1= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,0.5] \times [0,2\pi) \times [0,1]\} \f$
 and an insulating domain
\f$\Omega_2= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0.5,1] \times [0,2\pi) \times [0,1]\} \f$.
We also set \f$\Sigma= \Omega_1 \cap \Omega_2 \f$, \f$\Gamma_1= \partial \Omega_1\setminus \Sigma \f$ 
and \f$\Gamma_2= \partial \Omega_2\setminus \Sigma \f$.
The data are the source term \f$\mathbf{j}\f$, the given velocity \f$\bu^\text{given}\f$,
 the boundary data \f$\phi_{\text{bdy}}\f$,
the initial datas \f$\bH_0\f$ and \f$\phi _0\f$. The parameter \f$\Rm\f$ is the magnetic Reynolds number.
 The parameter \f$\mu^c\f$, respectively \f$\mu^v\f$, is the magnetic permeability of the conducting region \f$\Omega_1\f$, respectively of the vacuum region \f$\Omega_2\f$.
 The parameter \f$\sigma\f$ is the conductivity in the conducting region \f$\Omega_1\f$.



<h3>Manufactured solutions</h3>
Tha manufactured solutions considered in this test involve the modified bessel function of the third kind. These functions are defined for an integer \f$n\f$ and for any positive real argument \f$x\f$. In the following, we denote this function \f$BESSK(n,x)\f$. They are defined in the file <tt>bessel.f90</tt> in the following directory: ($SFEMaNS_DIR)/FEMSUB. We refer to C.W. CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS, MATHEMATICAL TABLES, VOL.5, 1962 for more details on these functions and their definition.

We introduce the following real numbers:
@f{align*}
f_{1,r_0} &= \frac{1}{k_0} \left(\frac{2}{r_0} -\frac{m_0}{r_0} BESSK(m_0,z_0) \right),
\\f_{2,r_0} &= \frac{m_0}{k_0 r_0} - \left( -BESSK(m_0-1,z_0) -
 \frac{m_0}{k_0 r_0} BESSK(m_0,z_0) \right),
\\df_{2,r_0}&= \frac{m_0}{r_0} f_{1,r_0} - f_{2,r_0} -k_0 BESSK(m_0,z_0),
\\A &= 3 f_{2,r_0} - r_0 df_{2,r_0},
\\B &= r_0 df_{2,r_0} - 2 f_{2,r_0},
\\f_1 &=  (\frac{r}{r_0})^2 f_{1,r_0},
\\f_2 &=  (\frac{r}{r_0})^2 \left( A+ B\frac{r}{r_0} \right),
\\f_4 &= \frac{r}{r_0^2} \left( 3 A+ 4 B\frac{r}{r_0} \right) 
- m_0\frac{r}{r_0^2} f_{1,r_0},
\\df_3 &= \frac{2r}{r_0^2},
@f}
with \f$m_0=1\f$, \f$r_0=0.5\f$, \f$k_0=2\pi\f$ and \f$z_0=r_0 k_0\f$.

We approximate the following analytical solutions:
@f{align*}
H_r(r,\theta,z,t) &= \left(\frac{m_0 r}{r_0^2} -k_0 f_2\right) 
\cos(k_0 z) \cos(m_0\theta)\cos(t),
\\ H_{\theta}(r,\theta,z,t) &= \left(k_0 f_1-df_3\right)
 \cos(k_0 z)\sin(m_0\theta) \cos(t),
\\ H_z(r,\theta,z,t) &= f_4 \sin(k_0 z) \cos(m_0\theta)\cos(t),
\\ \phi(r,\theta,z,t) &= BESS(m_0,k_0 r) \cos(k_0 z) \cos(m_0\theta) \cos(t).
@f}

We also set the given velocity field as follows:
@f{align*}
u^\text{given}_r(r,\theta,z,t) &= 0, 
\\ u^\text{given}_{\theta}(r,\theta,z,t) &= 0, 
\\ u^\text{given}_z(r,\theta,z,t) &= 0.
@f}
The source term \f$ \mathbf{j} \f$ and the boundary data \f$\phi_{\text{bdy}}\f$ are computed accordingly.



<h3>Generation of the mesh</h3>
The finite element mesh used for this test is named <tt>Mesh_10_form.FEM</tt> and 
 has a mesh size of \f$0.1\f$ for the P1 approximation. You can generate this mesh with the files in the following directory:
 ($SFEMaNS_MESH_GEN_DIR)/EXAMPLES/EXAMPLES_MANUFACTURED_SOLUTIONS/Mesh_10_form.
The following images shows the mesh for P1 finite elements for the conducting and vacuum regions.

<table width="70%" align="center">
<tr>
    <td align="center">
    @image html  fig_Mesh_10_form.png  "Finite element mesh (P1) of the conducting region."
    </td>
    <td align="center">
    @image html  fig_Mesh_10_form_vacuum.png  "Finite element mesh (P1) of the vacuum region."
    </td>
</tr>
</table>







<h3>Information on the file <tt>condlim.f90</tt></h3>
The initial conditions, boundary conditions and the forcing term \f$\textbf{j}\f$ in the Maxwell
 equations are set in the file <tt>condlim_test_4.f90</tt>.
Here is a description of the subroutines and functions of interest.
<ol>
<li>First we define the numbers \f$r_0\f$, \f$k_0\f$ and \f$m_0\f$ at the begining of the module
 so that every subroutines has access to these real numbers.
\code
  REAL(KIND=8), PRIVATE  :: r0 = 0.5d0
  REAL(KIND=8), PRIVATE  :: k0=6.28318530717958d0
  INTEGER     , PRIVATE  :: m0=1
\endcode
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
 It is used to initialize the magnetic field and to impose Dirichlet boundary
 conditions on \f$\textbf{H}\times\textbf{n}\f$ with \f$\textbf{n}\f$ the outter normal vector.
<ol>
<li>We define the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>If the Fourier mode is not equal to \f$m_0\f$, we set the magnetic field to zero.
\code
    IF (m/=m0) THEN
       vv = 0
       RETURN
    END IF
\endcode
<li>We define various real numbers.
\code
    !value in r0
    z0 = r0*k0
    f3_r0  = 1.d0
    df3_r0 = 2.d0/r0  
    f1_r0  = 1.d0/k0*(df3_r0-m0/r0*BESSK(m0,z0))
    !f2_r0  = m0/(k0*r0)*f3_r0 - (BESSI(m0-1,z0)-m0/(k0*r0)*BESSI(m0,z0))
    f2_r0  = m0/(k0*r0)*f3_r0 - (-BESSK(m0-1,z0)-m0/(k0*r0)*BESSK(m0,z0))
    df2_r0 = (m0/r0*f1_r0-f2_r0/r0-k0*BESSK(m0,z0))
    !to evaluate f2
    A = 3*f2_r0 - r0*df2_r0
    B = r0*df2_r0-2*f2_r0    
    !function for all points
    f1  = (r/r0)**2*f1_r0
    f2  = (r/r0)**2*(A+B*r/r0)
    f3  = (r/r0)**2
    df3 = 2*r/r0**2
    df2 = r/r0**2*(2*A+3*B*r/r0)
    f4  = r/r0**2*(A+B*r/r0) + df2 - m0*r/r0**2*f1_r0 
\endcode
<li>For the Fourier modes \f$m_0\f$, we define the magnetic field depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) as follows:
\code
    IF (TYPE == 1) THEN
       vv = (m0*r/r0**2-k0*f2)*COS(k0*z) 
    ELSEIF (TYPE == 2) THEN
       vv = 0.d0
    ELSEIF (TYPE ==3) THEN
       vv = 0.d0
    ELSEIF (TYPE == 4)  THEN
       vv = (k0*f1-df3)*COS(k0*z)
    ELSEIF (TYPE == 5) THEN
       vv = f4*SIN(k0*z)
    ELSEIF (TYPE == 6) THEN
       vv = 0.d0
    ENDIF
    vv = vv * COS(t)
    RETURN
\endcode
</ol>
<li>The function <tt>Phiexact </tt> contains the analytical scalar potential.
 It is used to initialize the scalar potential and to impose Dirichlet boundary
 conditions on the scalar potential.
<ol>
<li>We define the radial and vertical coordinates r, z.
\code
   r = rr(1,:)
   z = rr(2,:)
\endcode
<li>We set the scalar potential to zero.
\code
   vv(:) = 0.d0
\endcode
<li>If we consider a Fourier mode different than \f$m_0\f$ we are done.
\code
   IF (m/=m0) RETURN
\endcode
<li>For the Fourier modes \f$m_0\f$, we define the scalar potential depending of its TYPE (1 for cosine and 2 for sine) as follows:
\code
   IF  (TYPE == 1) THEN
      DO n= 1, SIZE(rr,2)
         vv(n) = BESSK(m0,k0*r(n))*COS(k0*z(n))
      END DO
   ELSEIF (TYPE == 2) THEN
      vv = 0.d0
   ENDIF
   vv = vv*COS(t)
   RETURN
\endcode
</ol>
<li>The function <tt>Vexact</tt> is used to define the given velocity field \f$\bu^\text{given}\f$. It is set to zero.
\code
    vv = 0.d0
    RETURN
\endcode
<li>The function <tt>Jexact_gauss</tt> is used to define the source term \f$\textbf{j}\f$. It uses the function <tt>Eexact_gauss</tt>.
<li>The function <tt>Eexact_gauss</tt> is used to define the part of Jexact that matches the time derivative of the magnetic field.
</ol>
All the other subroutines present in the file <tt>condlim_test_4.f90</tt> are not used in this test. We remind that the bessel function \f$BESSK\f$ are defined in the file <tt>bessel.f90</tt> in the following directory: ($SFEMaNS_DIR)/FEMSUB.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file.







<h3>Setting in the data file</h3>
We describe the data file of this test.
 It is called <tt>debug_data_test_4</tt> and can be found in the following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
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
<li>We use two processors in the meridian section. It means the finite element mesh is subdivised in two.
\code
===Number of processors in meridian section
2
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
<li>We select specific Fourier modes to solve by setting:
\code
===Select Fourier modes? (true/false)
.t.
\endcode
<li>We give the list of the Fourier modes to solve.
\code
===List of Fourier modes (if select_mode=.TRUE.)
1
\endcode
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
We note that we need to say if we use a restart for the velocity. If yes, the \f$\bu^\text{given}\f$ comes from a suite file and not the function <tt>Vexact</tt>.
<li>We use a time step of \f$0.01\f$ and solve the problem over \f$100\f$ time iterations.
\code
===Time step and number of time iterations
.01d0, 100
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
4 2 .0 1.
\endcode
We note that we need as much as lines as the number of pairs of boundaries with periodic condition.
</ol>
<li>We set the number of domains and their label, see the files associated to the generation of the mesh,
 where the code approximates the magnetic field.
\code
===Number of subdomains in magnetic field (H) mesh
1
===List of subdomains for magnetic field (H) mesh
1
\endcode
<li>We set the number of interface in H_mesh and give their respective labels.
\code
===Number of interfaces in H mesh
0
===List of interfaces in H mesh
xx
\endcode
Such interfaces represent interfaces with discontinuities in magnetic permeability or interfaces between the magnetic field mesh and the temperature or the velocity field meshes.
<li>We set the number of boundaries with Dirichlet conditions on the magnetic field and give their respective labels.
\code
===Number of Dirichlet sides for Hxn
0
===List of Dirichlet sides for Hxn
xx
\endcode
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
1
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
5
\endcode
<li>We set the magnetic permeability in the vacuum domain.
\code
===Permeability in vacuum
1.d0
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
<li>The L2 norm of the error on the magnetic field \f$\textbf{H}\f$.
<li>The L2 norm of the error on the curl of the magnetic field \f$\textbf{H}\f$.
<li>The L2 norm of the divergence of the magnetic field \f$\textbf{B}=\mu\bH\f$.
<li>The H1 norm of the error on the scalar potential.
</ol>

 They are compared to reference values to attest of the correctness of the code.  
 These values of reference are in the last lines of the file <tt>debug_data_test_4</tt> in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
============================================
Mesh_10_form.FEM, P1P2
===Reference results
0.117284337370832 L2 error on Hn     
0.158175579990508 L2 error on Curl(Hn)  
0.276743273545775 L2 norm of Div(mu Hn)
0.236676729274824 H1 error on phin
\endcode


To conclude this test, we show the profile of the approximated magnetic field \f$\textbf{H}\f$ and
  scalar potential \f$\phi\f$ at the final time.
 These figures are done in the plane \f$y=0\f$ which is 
 the union of the half plane \f$\theta=0\f$ and \f$\theta=\pi\f$.
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_04_mag_tfin.png "Magnetic field magnitude in the plane plane y=0."
    </td>
    <td align="center">
    @image html  fig_test_04_phi_tfin.png  "Scalar potential in the plane plane y=0."
    </td>
</tr>
</table>
 */
