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
 * @page doc_debug_test_12 Test 12: Maxwell P1 Neumann bdy

<h3>Introduction</h3>
In this example, we check the correctness of SFEMaNS for a magnetic problem involving a conducting domain. We consider Dirichlet and Neumann boundary conditions. We use P1 finite elements for the magnetic field.

We solve the Maxwell equations:
\f{align}{
\begin{cases}
\partial_t (\mu^c \mathbf{H}) + \nabla \times \left(\frac{1}{\Rm \sigma}\nabla \times \mathbf{H}  \right)  =
 \nabla\times (\bu^\text{given} \times \mu^c \mathbf{H}) 
+ \nabla \times \left(\frac{1}{\Rm \sigma} \mathbf{j}^s \right) & \text{in } \Omega_1, \\
\text{div} (\mu^c \mathbf {H}) = 0  & \text{in } \Omega_1 ,\\
\bH \times \bn = \bH_{\text{bdy}} \times \bn & \text{on } \Gamma_1,\\
 \left( \frac{1}{\Rm \sigma} \left( \ROT (\mathbf{H}) - \mathbf{j}^s  \right) - \bu^\text{given} \times \mu \mathbf{H}
 \right)  \times \bn  = \textbf{a} \times \bn & \text{on } \Gamma_2,\\
\bH_{|t=0}= \bH_0, 
\end{cases}
\f}
in the domain \f$\Omega_1= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,0.5] \times [0,2\pi) \times [0,1]\} \f$.
We also set \f$\Gamma_1= \partial \Omega_1\cap  \{ z=1\} \f$ 
and \f$\Gamma_2= \partial \Omega_1 \cap  \{ r=1\} \cap \{z=0\} \f$.
The data are the source term \f$\mathbf{j}^s\f$, the given velocity \f$\bu^\text{given}\f$,
 the boundary datas \f$\bH_{\text{bdy}}\f$, \f$\textbf{a}\f$, and
the initial data \f$\bH_0\f$. The parameter \f$\Rm\f$ is the magnetic Reynolds number.
 The parameter \f$\mu^c\f$ is the magnetic permeability of the conducting region \f$\Omega_1\f$. 
 The parameter \f$\sigma\f$ is the conductivity in the conducting region \f$\Omega_1\f$.

<h3>Manufactured solutions</h3>
We approximate the following analytical solution:
@f{align*}
H_r(r,\theta,z,t) &= \be_z,
\\ H_{\theta}(r,\theta,z,t) &= r - r^2 \sin(\theta),
\\ H_z(r,\theta,z,t) &= r z \cos(\theta).
@f}
We also set the given velocity field as follows:
@f{align*}
u^\text{given}_r(r,\theta,z,t) &= 0, 
\\ u^\text{given}_{\theta}(r,\theta,z,t) &= 0, 
\\ u^\text{given}_z(r,\theta,z,t) &= 0.
@f}
The source term \f$ \mathbf{j}^s \f$ and the boundary datas \f$\bH_{\text{bdy}}, 
\textbf{a}\f$ are computed accordingly: \f$\bE=\be_z\f$ and \f$j^s_r=-z \sin(\theta)\f$, 
\f$j^s_\theta=-z \cos(\theta)\f$, \f$b^s_z=-3 r\sin(\theta)\f$.

<h3>Generation of the mesh</h3>
The finite element mesh used for this test is named <tt>Mesh_10_form.FEM</tt> and 
 has a mesh size of \f$0.1\f$ for the P1 approximation. You can generate this mesh with the files in the following directory:
 ($SFEMaNS_MESH_GEN_DIR)/EXAMPLES/EXAMPLES_MANUFACTURED_SOLUTIONS/Mesh_10_form.
The following images show the mesh for P1 finite elements.

<table width="70%" align="center">
<tr>
    <td align="center">
    @image html  fig_Mesh_10_form.png  "Finite element mesh (P1)."
    </td>
</tr>
</table>







<h3>Information on the file <tt>condlim.f90</tt></h3>
The initial conditions, boundary conditions and the forcing term \f$\textbf{j}^s\f$ in the Maxwell
 equations are set in the file <tt>condlim_test_12.f90</tt>.
Here is a description of the subroutines and functions of interest.
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
          Hn(:,k,i) = Hexact(H_mesh, k, H_mesh%rr, list_mode(i), mu_H_field, time)
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
<li>For the Fourier modes m∈{0,1}, we define the magnetic field depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) as follows: 
\code
    IF (m==1) THEN
       IF (TYPE==4) THEN
          vv = -rr(1,:)**2
       ELSE IF (TYPE==5) THEN
          vv = rr(1,:)*rr(2,:)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (m==0) THEN
       IF (TYPE==3) THEN
          vv = rr(1,:)
       ELSE IF (TYPE==5) THEN
          vv = 1.d0
       ELSE
          vv = 0.d0
       END IF
\endcode
<li> For the other Fourier modes, the magnetic field is set to zero.
\code
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
<li>The function <tt>Jexact_gauss</tt> is used to define the source term \f$\textbf{j}^s\f$.
 Since \f$\partial_t (\mu \textbf{H})=0\f$  and \f$ \bu^\text{given} \times \textbf{H}=0 \f$, we define this source term so it satisfies the relation \f$ \ROT\textbf{j}^s=\ROT(\ROT\textbf{H} )\f$. As this term only depends of the Fourier mode \f$m=1\f$, we define it as follows:
\code
   IF (m==1) THEN
      IF (TYPE==2) THEN
         vv = -rr(2)
      ELSE IF (TYPE==3) THEN
         vv = -rr(2)
      ELSE IF (TYPE==6) THEN
         vv = -3*rr(1)
      ELSE
         vv = 0.d0
      END IF
   ELSE
      vv = 0.d0
   END IF
   RETURN
\endcode
<li>The function <tt>Eexact_gauss</tt> is used to define the boundary data \f$\textbf{a}\f$. Since  \f$ \bu^\text{given} \times \textbf{H}=0 \f$, we define this boundary data so it satisfies the relation \f$\textbf{a}= \frac{1}{\Rm \sigma} \left( \ROT (\mathbf{H}) - \mathbf{j}^s  \right)\f$.
\code
    IF (m/=0) THEN
       vv = 0
    ELSE
       IF (TYPE==5) THEN
          vv = 2.d0
       ELSE
          vv = 0.d0
       END IF
    END IF
    RETURN
\endcode
</ol>
All the other subroutines present in the file <tt>condlim_test_12.f90</tt> are not used in this test.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file.







<h3>Setting in the data file</h3>
We describe the data file of this test.
 It is called <tt>debug_data_test_12</tt> and can be found in the following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
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
<li>We solve the problem for \f$2\f$ Fourier modes.
\code
===Number of Fourier modes
2
\endcode
<li>We use \f$2\f$ processors in Fourier space.
\code
===Number of processors in Fourier space
2
\endcode
It means that each processors is solving the problem for \f$2/2=1\f$ Fourier modes.
<li>We select specific Fourier modes to solve.
\code
===Select Fourier modes? (true/false)
.t.
\endcode
<li>We give the list of the Fourier modes to solve.
\code
===List of Fourier modes (if select_mode=.TRUE.)
0 1
\endcode
We note that setting select Fourier modes to false would give the same result 
 as we select the first two Fourier modes.
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
<li>We use a time step of \f$1\f$ and solve the problem over \f$1\f$ time iterations.
\code
===Time step and number of time iterations
1.d0, 1
\endcode
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
1
===List of Dirichlet sides for Hxn
4
\endcode
We note that we did not select the boundary of index 3 that correspond to the face \f$r=0.5\f$. 
 As consequence, Neumann boundary condition are apply to this face.
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
that can be found in the following: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. 

To check the well behavior of the code, we compute three quantities:
<ol>
<li>The L2 norm of the error on the magnetic field \f$\textbf{H}\f$.
<li>The L2 norm of the error on the curl of the magnetic field \f$\textbf{H}\f$.
<li>The L2 norm of the divergence of the magnetic field \f$\textbf{B}=\mu\bH\f$.
</ol>

 They are compared to reference values to attest of the correctness of the code.  
 These values of reference are in the last lines of the file <tt>debug_data_test_12</tt> in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
============================================
Mesh_10_form.FEM, P1
===Reference results
1.3128820890830654E-003 L2 error on  Hn     
1.0628526039389356E-002 L2 error on Curl(Hn) 
3.9114181055566619E-002 L2 norm of error on Div(mu Hn)
1.                      Dummy ref
\endcode


To conclude this test, we show the profile of the approximated  magnetic field \f$\textbf{H}\f$ at the final time.
 This figure is done in the plane \f$y=0\f$ which is 
 the union of the half plane \f$\theta=0\f$ and \f$\theta=\pi\f$.
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_12_mag_tfin.png "Magnetic field magnitude in the plane plane y=0."
    </td>
</tr>
</table>
 */
