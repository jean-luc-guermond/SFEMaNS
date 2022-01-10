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
 * @page doc_debug_test_06 Test 6: Maxwell variable permeability with vacuum P1-P2

<h3>Introduction</h3>
In this example, we check the correctness of SFEMaNS for a magnetic problem involving a conducting and an insulating domain. The conducting domain is splitted into two subdomains with different magnetic permeability. We consider Dirichlet boundary conditions. We use P1 finite elements for the magnetic field and P2 finite element for the scalar potential. This set up can be referred as Durand sphere in the litterature.


We solve the Maxwell equations:
\f{align}{
\begin{cases}
\partial_t (\mu^c \mathbf{H}) + \nabla \times \left(\frac{1}{\Rm \sigma}\nabla \times \mathbf{H}  \right)  =
 \nabla\times (\bu^\text{given} \times \mu^c \mathbf{H}) 
+ \nabla \times \left(\frac{1}{\Rm \sigma} \nabla\times \mathbf{j} \right) & \text{in } \Omega_1, \\
\text{div} (\mu^c \mathbf {H}) = 0  & \text{in } \Omega_1 ,\\
- \partial_t \DIV (\mu^v \GRAD( \phi)) = 0 & \text{in } \Omega_2 ,\\
\bH_1 \times  \bn_1 +  \bH_2 \times  \bn_2 = 0 & \text{on } \Sigma_\mu,\\
\mu^c_1\bH _1 \cdot  \bn_1 +  \mu^c_2 \bH _2 \cdot  \bn_2 = 0 & \text{on } \Sigma_\mu,\\
\bH \times  \bn^c + \nabla \phi \times \bn^v = 0 & \text{on } \Sigma, \\
\mu^c \bH \cdot  \bn^c + \mu ^v \nabla \phi \cdot \bn^v = 0 & \text{on } \Sigma, \\
\phi = \phi_{\text{bdy}}  & \text{on } \Gamma_2,\\
\bH_{|t=0}= \bH_0, \\
\phi_{|t=0}= \phi _0,
\end{cases}
\f}
in the domain \f$\Omega= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,10] \times [0,2\pi) \times [0,10] \; | \; r^2+z^2 \leq 100\} \f$.
 This domain is the union of a conducting domain
\f$\Omega_1= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,1] \times [0,2\pi) \times [0,1] \; | \; r^2+z^2\leq 1\} \f$
 and an insulating domain
\f$\Omega_2= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,10] \times [0,2\pi) \times [0,10] \; | \; 1\leq r^2+z^2\leq 100\} \f$.
We also set \f$\Sigma= \Omega_1 \cap \Omega_2 \f$, \f$\Gamma_1= \partial \Omega_1\setminus \Sigma \f$ 
and \f$\Gamma_2= \partial \Omega_2\setminus \Sigma \f$ and 
\f$\Sigma_\mu = \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,1] \times [0,2\pi) \times [0,1] \; | \; r^2+z^2 = 0.5^2\} \f$.
 The variable \f$\bH_1\f$, respectively \f$\bH_2\f$, represents the magnetic field in the area \f$\{r^2+z^2 \leq 0.5^2\} \f$, respectively \f$\{0.5^2\leq r^2+z^2  \leq 1\} \f$. By analogy, we define outward normals \f$\bn_1, \bn_2\f$ to the surface \f$\Sigma_\mu\f$.
The data are the source term \f$\mathbf{j}\f$, the given velocity \f$\bu^\text{given}\f$,
 the boundary data \f$\phi_{\text{bdy}}\f$,
the initial datas \f$\bH_0\f$ and \f$\phi _0\f$. The parameter \f$\Rm\f$ is the magnetic Reynolds number.
 The parameter \f$\mu^c\f$, respectively \f$\mu^v\f$, is the magnetic permeability of the conducting region \f$\Omega_1\f$, respectively of the vacuum region \f$\Omega_2\f$.
 The parameter \f$\sigma\f$ is the conductivity in the conducting region \f$\Omega_1\f$.

<h3>Manufactured solutions</h3>
We approximate the following analytical solutions:
@f{align*}
H_r(r,\theta,z,t) &= 
\begin{cases} 0 \quad \text{if} \quad \sqrt{r^2+z^2}\leq 0.5,
\\  3 \text{capC} (\frac{a}{\sqrt{r^2+z^2}})^3 \cos(\varphi)\sin(\varphi) \quad \text{if} \quad \sqrt{r^2+z^2} > 0.5,
\end{cases}
\\ H_{\theta}(r,\theta,z,t) &=0,
\\ H_z(r,\theta,z,t) &=
\begin{cases} -\text{capA}   \quad \text{if} \quad \sqrt{r^2+z^2} \leq 0.5,
\\ - \text{capB} +\text{capC} (\frac{a}{\sqrt{r^2+z^2}})^3 ( 3 \cos(\varphi)^2-1)  \quad \text{if} \quad \sqrt{r^2+z^2} > 0.5,
\\
\end{cases},
\\ \phi(r,\theta,z,t) &= \sqrt{r^2+z^2} \cos(\varphi) - \frac{\text{capD} \cos(\varphi) a^3}{r^2+z^2},
@f}
where \f$\varphi\f$ is the polar angle of the spherical coordinate system and where \f$\text{capA}, \text{capB}, \text{capC},\f$ and \f$\text{capD}\f$ are defined as follows:
@f{align*}
\text{capA} &=  \dfrac{-9 \mu_\max^c \mu_0}{(2\mu_\max^c\mu_0)(\mu_\max^c+2\mu_0)-2(\mu_\max^c-\mu_0)(\frac{a}{b})^3},
\\ \text{capB} &= \dfrac{(2\mu_\max^c+\mu_0)(\mu_\max^c-\mu_0)(\frac{b^3}{a^3}-1)}{(2\mu_\max^c\mu_0)(\mu_\max^c+2\mu_0)-2(\mu_\max^c-\mu_0)(\frac{a}{b})^3}  ,
\\ \text{capC} &= \left(1-\frac{\mu_0}{\mu_\max^c} \right) \frac{\text{capA}}{3}  ,
\\ \text{capD} &= \left(2+\frac{\mu_0}{\mu_\max^c} \right) \frac{\text{capA}}{3}  ,
@f}
with \f$\mu_\max^c = \max(\mu^c)\f$, \f$\mu_0=1\f$, \f$a=0.5\f$ and \f$b=1\f$.

 We also set the given velocity field to zero.
@f{align*}
u^\text{given}_r(r,\theta,z,t) &= 0, 
\\ u^\text{given}_{\theta}(r,\theta,z,t) &= 0, 
\\ u^\text{given}_z(r,\theta,z,t) &= 0.
@f}
The source term \f$ \mathbf{j}\f$ is set to zero and the boundary data \f$\phi_{\text{bdy}}\f$ is computed accordingly.


<h3>Generation of the mesh</h3>
The finite element mesh used for this test is named <tt>sphere_0.05_form.FEM</tt>. 
 The mesh size for the P1 approximation is \f$0.05\f$ for \f$r\leq1\f$ and \f$1\f$ for \f$r=10\f$. You can generate this mesh with the files in the following directory:
 ($SFEMaNS_MESH_GEN_DIR)/EXAMPLES/EXAMPLES_MANUFACTURED_SOLUTIONS/sphere_0.05_form.
The following images show the mesh for P1 finite elements for the conducting and vacuum regions.

<table width="70%" align="center">
<tr>
    <td align="center">
    @image html  fig_sphere_0.05_form.png  "Finite element mesh (P1) of the conducting region."
    </td>
    <td align="center">
    @image html  fig_sphere_0.05_form_vacuum.png  "Finite element mesh (P1) of the vacuum region."
    </td>
</tr>
</table>
We note that the conducting region is the union of two domain of respective index 1 and 2. Their interface is \f$r^2+z^2=0.5^2\f$ of index 2. It allows us to set different magnetic permeability (and also electrical conductivity) in each subregions of the conducting region.






<h3>Information on the file <tt>condlim.f90</tt></h3>
The initial conditions, boundary conditions and the forcing term \f$\textbf{j}\f$ in the Maxwell
 equations are set in the file <tt>condlim_test_6.f90</tt>.
Here is a description of the subroutines and functions of interest.
<ol>
<li>The subroutine <tt>init_maxwell</tt> initializes the magnetic field and the scalar potential 
 at the time\f$-dt\f$ and \f$0\f$ with \f$dt\f$ being the time step.
 It is done by using the functions Hexact, Phiexact as follows:
\code
    Hn1 = 0.d0
    Hn = 0.d0
    phin1 = 0.d0
    phin = 0.d0
    RETURN
\endcode
<li>The function <tt>Hexact</tt> contains the analytical magnetic field.
 It is used to initialize the magnetic field and to impose Dirichlet boundary
 conditions on \f$\textbf{H}\times\textbf{n}\f$ with \f$\textbf{n}\f$ the outter normal vector.
<ol>
<li>First we define the parameters \f$aa, bb\f$ and \f$\mu_0\f$.
\code
    REAL(KIND=8)                                      :: aa=0.5d0, bb=1.d0, mu0=1.d0
\endcode
<li>If the Fourier mode m is not equal to 0, we set the magnetic field to zero.
\code
    IF (m/=0) THEN
       vv = 0.d0
       RETURN
    END IF
\endcode
<li>If the the processor has no point associated to the magnetic field domain, it does nothing.
\code
    IF (SIZE(rr,2)==0) RETURN
\endcode
We note such case can not happen as the every processor has nodes of the magnetic field domain and the scalar potential domain.
<li>We check if the number of pair of coordinates (r,z) corresponds to the number of point of the mangetic field mesh.
\code
    IF (SIZE(rr,2)/=H_mesh%np) THEN
       CALL error_petsc(' BUG in Hexact')
    END IF
\endcode
<li>We define \f$\mu_\max^c\f$, the radial and vertical coordinates \f$r, z\f$.
\code
    !CHAMP DURAND, H/H/phi configuration, H multi-valued at interface.
    ! mu1 = mu3 = mu_0 = 1
    ! mu2 = mu
    mu = MAXVAL(mu_H_field)
    r = H_mesh%rr(1,:)
    z = H_mesh%rr(2,:)
\endcode
<li>We define the polar angle \f$\varphi\f$ and the distance to the origin \f$\rho\f$.
\code
    theta = ATAN2(r,z)
    rho = SQRT(r**2+z**2)
\endcode
<li>We define the parameters \f$\text{capA}, \text{capB}, \text{capC}\f$ and \f$\text{capD}\f$.
\code
    capA = -9*mu*mu0/((2*mu+mu0)*(mu+2*mu0) - 2*(mu-mu0)**2*(aa/bb)**3)
    capD = (2*mu+mu0)*(mu-mu0)*((bb/aa)**3-1.d0)/((2*mu+mu0)*(mu+2*mu0) - 2*(mu-mu0)**2*(aa/bb)**3)
    capC = (1.d0 - mu0/mu)*capA/3
    capB = (2.d0 + mu0/mu)*capA/3
\endcode
<li>For the Fourier mode \f$m=0\f$, we define the magnetic field depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) and the distance to the origin as follows:
\code
    DO n = 1, H_mesh%np
       IF (rho(n) .LE. aa) THEN
          IF (TYPE==1) THEN
             vv(n) =0.d0
          ELSE IF(TYPE==5) THEN
             vv(n) = -capA
          ELSE
             vv(n) = 0.d0
          END IF
       ELSE
          IF (TYPE==1) THEN
             vv(n) = 3*capC*(aa/rho(n))**3*COS(theta(n))*SIN(theta(n))
          ELSE IF(TYPE==5) THEN
             vv(n) = -capB + capC*(aa/rho(n))**3*(3.d0*COS(theta(n))**2-1.d0)
          ELSE
             vv(n) = 0.d0
          END IF
       END IF
    END DO
\endcode
<li>The presence of an interface in the magnetic field, due to a discontinuity of the magnetic permeability, 
 induces a multi-valued magnetic field on this interface. 
So we ensure that the magnetic field is correctly define on each sides of the interface.
\code
    DO ms = 1, H_mesh%mes !Do loop on interface, since H multi-valued
       DO ns = 1, H_mesh%gauss%n_ws
          n = H_mesh%jjs(ns,ms)
          IF (H_mesh%i_d(H_mesh%neighs(ms)) == 1) THEN 
             IF (TYPE==1) THEN
                vv(n) = 0.d0
             ELSE IF(TYPE==5) THEN
                vv(n) = -capA
             ELSE
                vv(n) = 0.d0
             END IF
          ELSE
             IF (TYPE==1) THEN
                vv(n) = 3*capC*(aa/rho(n))**3*COS(theta(n))*SIN(theta(n))
             ELSE IF(TYPE==5) THEN
                vv(n) = -capB + capC*(aa/rho(n))**3*(3.d0*COS(theta(n))**2-1.d0)
             ELSE
                vv(n) = 0.d0
             END IF
          END IF
       END DO
    END DO
    RETURN
\endcode
</ol>
<li>The function <tt>Phiexact </tt> contains the analytical scalar potential.
 It is used to initialize the scalar potential and to impose Dirichlet boundary
 conditions on the scalar potential.
<ol>
<li>First we define the parameters \f$a, b\f$ and \f$\mu_0\f$.
\code
   REAL(KIND=8)                                      :: a=0.5d0, b=1.d0, mu0=1.d0
\endcode
<li>If the Fourier mode m is not equal to 0, we set the scalar potential to zero.
\code
   IF (m/=0) THEN
      vv = 0.d0
      RETURN
   END IF
\endcode
<li>We define \f$\mu_\max^c\f$, the radial and vertical coordinates \f$r, z\f$ and the distance to the origin \f$\rho\f$.
\code
   !CHAMP DURAND
   mu = MAXVAL(inputs%mu_H)
   r = rr(1,:)
   z = rr(2,:)   
   rho = SQRT(r**2+z**2)
\endcode
<li>We define the polar angle \f$\varphi\f$.
\code
   DO n = 1, SIZE(rho)
      IF (rho(n).LE.1.d-10) THEN
         theta(n) = 0.d0
      ELSE
         theta(n) = ATAN2(r(n),z(n))
      END IF
   END DO
\endcode
<li>We define the parameters \f$\text{capA}\f$ and \f$\text{capD}\f$.
\code
   capA = -9*mu*mu0/((2*mu+mu0)*(mu+2*mu0) - 2*(mu-mu0)**2*(a/b)**3)
   capD = (2*mu+mu0)*(mu-mu0)*((b/a)**3-1.d0)/((2*mu+mu0)*(mu+2*mu0) - 2*(mu-mu0)**2*(a/b)**3)
\endcode
<li>For the Fourier mode \f$m=0\f$, we define the scalar potential depending of its TYPE (1 for cosine and 2 for sine)
 and the distance to the origin.
\code
   DO n = 1, SIZE(rho)
      IF (TYPE==1 .AND. rho(n).LE. (a+1.d-1)) THEN
         vv(n) =  -capA*rho(n)*COS(theta(n)) 
      ELSE IF (TYPE==1 .AND. rho(n) .GE. (b-1.d-1)) THEN
         vv(n) = (rho(n)*COS(theta(n)) - capD*COS(theta(n))*a**3/rho(n)**2) !*(t/t_al)**3/(1.d0+(t/t_al)**3)
      ELSE
         vv(n) = 0.d0
      END IF
   END DO
   RETURN
\endcode
We note that the condition \f$\rho(n) .LE. (a+1.d-1)\f$ is never satisfied for this test (but it will be for the test 7).
</ol>
<li>The function <tt>Vexact</tt> is used to define the given velocity field \f$\bu^\text{given}\f$. It is set to zero.
\code
   vv = 0.d0
\endcode
<li>The function <tt>Jexact_gauss</tt> is used to define the source term \f$\textbf{j}\f$.
 Since \f$\ROT \textbf{H}=0\f$, \f$ \bu^\text{given}=0 \f$ and \f$\partial_t \textbf{H}=0\f$, we have \f$\textbf{j}=0\f$.
\code
   vv=0.d0
   RETURN
\endcode
<li>The function <tt>Eexact_gauss</tt> is not used in this test (no Neumann boundary conditions). So it is set to zero.
\code
    vv = 0.d0
    RETURN
\endcode
</ol>
All the other subroutines present in the file <tt>condlim_test_6.f90</tt> are not used in this test.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file.







<h3>Setting in the data file</h3>
We describe the data file of this test.
 It is called <tt>debug_data_test_6</tt> and can be found in the following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<ol>
<li>We use a formatted mesh by setting:
\code
===Is mesh file formatted (true/false)?
.t.
\endcode
<li>The path and the name of the mesh are specified with the two following lines:
\code
===Directory and name of mesh file
'.' 'sphere_0.05_form.FEM'
\endcode
where '.' refers to the directory where the data file is, meaning ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<li>We use two processors in the meridian section. It means the finite element mesh is subdivised in two.
\code
===Number of processors in meridian section
2
\endcode
<li>We solve the problem for \f$1\f$ Fourier mode.
\code
===Number of Fourier modes
1
\endcode
<li>We use \f$1\f$ processor in Fourier space.
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
<li>We use a time step of \f$1.d10\f$ and solve the problem over \f$1\f$ time iterations.
\code
===Time step and number of time iterations
1.d10, 1
\endcode
<li>We approximate the Maxwell equations with the variable \f$\textbf{H}\f$.
\code
===Solve Maxwell with H (true) or B (false)?
.t.
\endcode
<li>We set the number of domains and their label, see the files associated to the generation of the mesh,
 where the code approximates the magnetic field.
\code
===Number of subdomains in magnetic field (H) mesh
2
===List of subdomains for magnetic field (H) mesh
1 2
\endcode
<li>We set the number of interface in H_mesh and give their respective labels.
\code
===Number of interfaces in H mesh
1
===List of interfaces in H mesh
2
\endcode
Such interfaces represent interfaces with discontinuities in magnetic permeability or interfaces between the magnetic field mesh and the temperature or the velocity field meshes.
<li>We set the number of boundaries with Dirichlet conditions on the magnetic field.
\code
===Number of Dirichlet sides for Hxn
0
\endcode
<li>We set the magnetic permeability in each domains where the magnetic field is approximated.
\code
===Permeability in the conductive part (1:nb_dom_H)
1.d0 2.d0
\endcode
<li>We set the conductivity in each domains where the magnetic field is approximated.
\code
===Conductivity in the conductive part (1:nb_dom_H)
1.d0 1.d0
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
3
\endcode
<li>We set the number of boundaries with Dirichlet conditions on the scalar potential and give their respective labels.
\code
===How many boundary pieces for Dirichlet BCs on phi?
1
===List of boundary pieces for Dirichlet BCs on phi
4
\endcode
<li>We set the number of interface between the magnetic field and the scalar potential and give their respective labels.
\code
===Number of interfaces between H and phi
1
===List of interfaces between H and phi
3
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
 These values of reference are in the last lines of the file <tt>debug_data_test_6</tt> in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
============================================
sphere_0.05_form.FEM, P1P2
===Reference results
1.09002132008650578E-002 L2 error on Hn
1.44935294860469559E-002 L2 error on Curl(Hn) 
9.47173125059953824E-002 L2 norm of Div(mu Hn)
3.73975917654559026E-004 H1 error on phin
\endcode


To conclude this test, we show the profile of the approximated  magnetic field \f$\textbf{H}\f$ and
 scalar potential \f$\phi\f$ at the final time.
 These figures are done in the plane \f$y=0\f$ which is 
 the union of the half plane \f$\theta=0\f$ and \f$\theta=\pi\f$.
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_06_mag_tfin.png "Magnetic field magnitude in the plane plane y=0."
    </td>
    <td align="center">
    @image html  fig_test_06_phi_tfin.png  "Scalar potential in the plane plane y=0."
    </td>
</tr>
</table>
 */
