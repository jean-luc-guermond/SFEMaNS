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
 * @page doc_debug_test_03 Test 3: Maxwell with vacuum P2-P2

<h3>Introduction</h3>
In this example, we check the correctness of SFEMaNS for a magnetic problem involving a conducting and an insulating domain. We consider Dirichlet boundary conditions.  We use P2 finite elements for the magnetic field and the scalar potential.

We solve the Maxwell equations:
\f{align}{
\begin{cases}
\partial_t (\mu^c \mathbf{H}) + \nabla \times \left(\frac{1}{\Rm \sigma}\nabla \times \mathbf{H}  \right)  =
 \nabla\times (\bu^\text{given} \times \mu^c \mathbf{H}) 
+ \nabla \times \left(\frac{1}{\Rm \sigma} \nabla\times \mathbf{j} \right) & \text{in } \Omega_1, \\
\text{div} (\mu^c \mathbf {H}) = 0  & \text{in } \Omega_1 ,\\
-\partial_t \DIV (\mu^v \GRAD( \phi)) = 0 & \text{in } \Omega_2 ,\\
\bH \times  \bn^c + \nabla \phi \times \bn^v = 0 & \text{on } \Sigma, \\
\mu^c \bH \cdot  \bn^c + \mu ^v \nabla \phi \cdot \bn^v = 0 & \text{on } \Sigma, \\
\bH \times \bn = \bH_{\text{bdy}} \times \bn & \text{on } \Gamma_1,\\
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
 the boundary data \f$\bH_{\text{bdy}}\f$ and \f$\phi_{\text{bdy}}\f$,
the initial datas \f$\bH_0\f$ and \f$\phi _0\f$. The parameter \f$\Rm\f$ is the magnetic Reynolds number.
 The parameter \f$\mu^c\f$, respectively \f$\mu^v\f$, is the magnetic permeability of the conducting region \f$\Omega_1\f$, respectively of the vacuum region \f$\Omega_2\f$.
 The parameter \f$\sigma\f$ is the conductivity in the conducting region \f$\Omega_1\f$.


<h3>Manufactured solutions</h3>
We approximate the following analytical solutions:
@f{align*}
H_r(r,\theta,z,t) &= \frac{\alpha z}{\mu^c} \left( \cos(\theta) + \frac{r}{4}\cos(2\theta)+ \frac{r^2}{9}\cos(3\theta) \right)\cos(t)
+ \frac{\beta z}{\mu^c} \left(\sin(\theta) + \frac{r}{4}\sin(2\theta)+ \frac{r^2}{9}\sin(3\theta)  \right)\cos(t),
\\ H_{\theta}(r,\theta,z,t) &= \frac{\beta z}{\mu^c}  \left( \cos(\theta) + \frac{r}{4}\cos(2\theta)+\frac{r^2}{9}\cos(3\theta) \right)\cos(t)
-\frac{\alpha z}{\mu^c} \left(\sin(\theta) + \frac{r}{4}\sin(2\theta)+\frac{r^2}{9}\sin(3\theta) \right)\cos(t),
\\ H_z(r,\theta,z,t) &= \frac{\alpha}{\mu^c} \left(r\cos(\theta) +\frac{r^2}{8}\cos(2\theta) + \frac{r^3}{27}\cos(3\theta)  \right) \cos(t)
+ \frac{\beta}{\mu^c} \left(r\sin(\theta) +\frac{r^2}{8}\sin(2\theta) + \frac{r^3}{27}\sin(3\theta)  \right) \cos(t),
\\ \phi(r,\theta,z,t) &= \frac{\alpha z}{\mu^v} \left( r \cos(\theta) + \frac{r^2}{8} \cos(2\theta) +
 \frac{r^3}{27}\cos(3\theta) \right) \cos(t) 
+ \frac{\beta z}{\mu^v} \left( r \sin(\theta) + \frac{r^2}{8} \sin(2\theta) +
 \frac{r^3}{27}\sin(3\theta) \right) \cos(t),
@f}
where \f$\alpha\f$ and \f$\beta\f$ are parameters. We also set the given velocity field as follows:
@f{align*}
u^\text{given}_r(r,\theta,z,t) &= \alpha z \left(\cos(\theta) + \frac{r}{4}\cos(2\theta) + \frac{r^2}{9}\cos(3\theta) \right)
+ \beta z \left( \sin(\theta) + \frac{r}{4}\sin(2\theta) + \frac{r^2}{9}\sin(3\theta) \right), 
\\ u^\text{given}_{\theta}(r,\theta,z,t) &= \beta z \left( \cos(\theta) + \frac{r}{4}\cos(2\theta) + \frac{r^2}{9}\cos(3\theta)\right)
-\alpha z \left( \sin(\theta) + \frac{r}{4}\sin(2\theta) + \frac{r^2}{9}\sin(3\theta) \right), 
\\ u^\text{given}_z(r,\theta,z,t) &= \alpha \left( r \cos(\theta) + \frac{r^2}{8}\cos(2\theta) + \frac{r^3}{27}\cos(3\theta)\right)
+\beta \left( r \sin(\theta) + \frac{r^2}{8}\sin(2\theta) + \frac{r^3}{27}\sin(3\theta) \right).
@f}
The source term \f$ \mathbf{j} \f$ and the boundary data \f$\bH_{\text{bdy}}, \phi_{\text{bdy}}\f$ are computed accordingly.



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
 equations are set in the file <tt>condlim_test_3.f90</tt>.
Here is a description of the subroutines and functions of interest.
<ol>
<li>First we define the parameters \f$\alpha\f$ and \f$\beta\f$ at the begining of the module
 so that every subroutines has acces to these real numbers.
\code
  REAL (KIND=8), PRIVATE  :: alpha=1.d0, beta=1.d0
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
<li>We check that the magnetic permeability is constant in the conducting domain.
\code
    IF (MAXVAL(mu_H_field) /= MINVAL(mu_H_field)) THEN
       CALL error_petsc(' BUG in condlim, mu not constant')
    END IF
\endcode
<li>We define the magnetic permeability and the radial and vertical coordinates r, z.
\code
    muH=mu_H_field(1) 
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>If the Fourier mode m is equal to 0, the magnetic field is set to zero.
\code
    IF (m==0) THEN
       vv = 0
       RETURN
    END IF
\endcode
<li>For the Fourier modes \f$m\in\{1,2,3\}\f$, we define the magnetic field depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) as follows:
\code
    IF (TYPE == 1) THEN
       vv = alpha*z*(r**(m-1))*m
    ELSEIF (TYPE == 2) THEN
       vv = beta *z*(r**(m-1))*m
    ELSEIF (TYPE ==3) THEN
       vv = beta *z*(r**(m-1))*m
    ELSEIF (TYPE == 4)  THEN
       vv =-alpha*z*(r**(m-1))*m
    ELSEIF (TYPE == 5) THEN
       vv = alpha*(r**m)
    ELSEIF (TYPE == 6) THEN
       vv = beta *(r**m)
    ENDIF
    vv = (vv/muH)*COS(t)/m**3
    RETURN
\endcode
where \f$t\f$ is the time.
</ol>
<li>The function <tt>Phiexact </tt> contains the analytical scalar potential.
 It is used to initialize the scalar potential and to impose Dirichlet boundary
 conditions on the scalar potential.
<ol>
<li>First we define the radial and vertical coordinates r, z.
\code
   r = rr(1,:)
   z = rr(2,:)
\endcode
<li>We set the scalar potential to zero.
\code
   vv(:) = 0.d0
\endcode
<li>If the Fourier mode m is equal to 0 we already defined the scalar potential correcty.
\code
   IF (m==0) RETURN
\endcode
<li>For the Fourier modes \f$m\in\{1,2,3\}\f$, we define the scalar depending of its TYPE (1 for cosine and 2 for sine) as follows:
\code
   IF  (TYPE == 1) THEN
      vv = alpha*z*(r**m)
   ELSEIF (TYPE == 2) THEN
      vv = beta *z*(r**m)
   ENDIF
   vv = (vv/mu_phi)*COS(t)/m**3
   RETURN
\endcode
where \f$t\f$ is the time.
</ol>
<li>The function <tt>Vexact</tt> is used to define the given velocity field \f$\bu^\text{given}\f$.
<ol>
<li>If the Fourier mode m is equal to 0,the velocity field is set to zero.
\code
    IF (m==0) THEN
       vv = 0
       RETURN
    END IF
\endcode
<li>We define the radial and vertical coordinates r, z.
\code
    r => H_mesh%rr(1,:)
    z => H_mesh%rr(2,:)
\endcode
<li>For the Fourier modes \f$m\in\{1,2,3\}\f$, we define the velocity field depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine).
\code
    vv(:,1) = alpha*z*(r**(m-1))*m
    vv(:,2) = beta *z*(r**(m-1))*m
    vv(:,3) = beta *z*(r**(m-1))*m
    vv(:,4) =-alpha*z*(r**(m-1))*m
    vv(:,5) = alpha*(r**m)
    vv(:,6) = beta *(r**m)
    VV = vv/m**3
    RETURN
\endcode
</ol>
<li>The function <tt>Jexact_gauss</tt> is used to define the source term \f$\textbf{j}\f$.
 Since \f$\ROT \textbf{H}=0\f$ and \f$ \bu^\text{given} \times \textbf{H}=0 \f$, this source term has to satisfy
 the relation \f$ \ROT(\frac{1}{\sigma\Rm}\textbf{j})=\partial_t (\mu \textbf{H}) \f$.
It is done by using the function Eexact_gauss as follows:
\code
   vv = -sigma* Eexact_gauss(TYPE, rr, m, mu_phi, sigma, mu_H, t)
   RETURN
\endcode
We note that we do not multiply by the magnetic Reynolds number \f$\Rm\f$ because the sigma that is given to the function <tt>Jexact_gauss</tt> has already been multiplied by \f$\Rm\f$.
<li>The function <tt>Eexact_gauss</tt> is used to define a vector, let's say \f$\textbf{D}\f$, such that
\f$ \ROT \textbf{D}=\partial_t (\mu^c \textbf{H})\f$.
<ol>
<li>We define the radial and vertical coordinates r, z.
\code
    r = rr(1)
    z = rr(2)
\endcode
<li>We set Eexact_gauss to zero.
\code
    vv = 0.d0
\endcode
<li>If the Fourier mode m is equal to 0, Eexact_gauss is already defined correcty.
\code
    IF (m == 0) RETURN
\endcode
<li>For the Fourier modes \f$m\in\{1,2,3\}\f$, we define it depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine).
\code
    IF  (TYPE == 1) THEN
       vv = 0.d0
    ELSEIF (TYPE == 2) THEN
       vv = 0.d0
    ELSEIF (TYPE ==3) THEN
       vv = alpha*(-1.d0/(m+2)*r**(m+1))
    ELSEIF (TYPE == 4)  THEN
       vv = beta *(-1.d0/(m+2)*r**(m+1))
    ELSEIF (TYPE == 5) THEN
       vv =  beta*z*(r**m)
    ELSEIF (TYPE == 6) THEN
       vv =-alpha*z*(r**m)
    ENDIF
    vv = -vv*SIN(t)/m**3
    RETURN
\endcode
As \f$\mu^c\f$ is set to one in the <tt>data</tt> file, we do not multiply by its value.
</ol>
</ol>
All the other subroutines present in the file <tt>condlim_test_3.f90</tt> are not used in this test.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file.







<h3>Setting in the data file</h3>
We describe the data file of this test.
 It is called <tt>debug_data_test_3</tt> and can be found in the following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
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
<li>We select specific Fourier modes to solve by setting:
\code
===Select Fourier modes? (true/false)
.t.
\endcode
<li>We give the list of the Fourier modes to solve.
\code
===List of Fourier modes (if select_mode=.TRUE.)
1 2 3
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
1.d-2, 100
\endcode
<li>We approximate the Maxwell equations with the variable \f$\textbf{B}=\mu\textbf{H}\f$.
\code
===Solve Maxwell with H (true) or B (false)?
.f.
\endcode
<li>We set the number of domains and their label, see the files associated to the generation of the mesh,
 where the code approximates the magnetic field.
\code
===Number of subdomains in magnetic field (H) mesh
1
===List of subdomains for magnetic field (H) mesh
1
\endcode
<li>We set the number of interface in H_mesh give their respective labels.
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
2
===List of Dirichlet sides for Hxn
2 4
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
1
===List of subdomains for magnetic potential (phi) mesh
2
\endcode
<li>We set the number of boundaries with Dirichlet conditions on the scalar potential and give their respective labels.
\code
===How many boundary pieces for Dirichlet BCs on phi?
3
===List of boundary pieces for Dirichlet BCs on phi
2 3 4
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
 These values of reference are in the last lines of the file <tt>debug_data_test_3</tt> in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
============================================
Mesh_10_form.FEM, P2P2
===Reference results
1.673900216118972E-005 L2 error on Hn     
4.334997301723562E-005 L2 error on Curl(Hn)  
2.445564215398369E-004 L2 norm of Div(mu Hn)
1.303414733032171E-005 H1 error on phin
\endcode


To conclude this test, we show the profile of the approximated  magnetic field \f$\textbf{H}\f$ and
 scalar potential \f$\phi\f$ at the final time.
 These figures are done in the plane \f$y=0\f$ which is 
 the union of the half plane \f$\theta=0\f$ and \f$\theta=\pi\f$.
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_03_mag_tfin.png "Magnetic field magnitude in the plane plane y=0."
    </td>
    <td align="center">
    @image html  fig_test_03_phi_tfin.png  "Scalar potential in the plane plane y=0."
    </td>
</tr>
</table>
 */
