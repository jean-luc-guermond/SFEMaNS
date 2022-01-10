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
 * @page doc_debug_test_23 Test 23 Maxwell Equations with variable magnetic permeability in (r, theta, z) with jumps



<h3>Introduction</h3>
In this example, we check the correctness of SFEMaNS for stationary magnetic problem involving a conducting. The magnetic permeability is variable in \f$(r,,\theta,z)\f$ and presents a jump in \f$r=1\f$. The magnetic permeability is given. We consider Dirichlet boundary conditions. We use P2 finite elements. This set up is approximated over one time iterations with a large time step.


We solve the Maxwell equations with the magnetic field \f$\bB=\mu \bH\f$ as dependent variable:
\f{align}{
\begin{cases}
\partial_t (\mathbf{B}) + \nabla \times \left(\frac{1}{\Rm \sigma}\nabla \times( \frac{1}{\mu^c}\mathbf{B} ) \right)  =
 \nabla\times (\bu^\text{given} \times  \mathbf{B}) 
+ \nabla \times \left(\frac{1}{\Rm \sigma} \nabla\times \mathbf{j} \right) & \text{in } \Omega, \\
\text{div} (\mathbf{B}) = 0  & \text{in } \Omega ,\\
\bH_1 \times  \bn_1 +  \bH_2 \times  \bn_2 = 0 & \text{on } \Sigma_\mu,\\
\mu^c_1\bH _1 \cdot  \bn_1 +  \mu^c_2 \bH _2 \cdot  \bn_2 = 0 & \text{on } \Sigma_\mu,\\
{\bH \times \bn}_{|\Gamma} = {\bH_{\text{bdy}} \times \bn}_{|\Gamma}, &\\
\bH_{|t=0}= \bH_0, \\
\end{cases}
\f}
in the domain \f$\Omega= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,2] \times [0,2\pi) \times [0.25,1]\} \f$.
This domain is the union of two conduction domain \f$\Omega_1=\{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,1] \times [0,2\pi) \times [0.25,1]\}\f$ and \f$ \Omega_2= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [1,2] \times [0,2\pi) \times [0.25,1]\}\f$. We denote by \f$\Sigma_\mu\f$ the interface of \f$\Omega_1\f$ and \f$\Omega_2\f$.
We also set \f$\Gamma= \partial \Omega\f$.
 We define outward normals \f$\bn_1, \bn_2\f$ to the surface \f$\Sigma_\mu\f$ and the outward normals \f$\bn\f$ to the surface \f$\Gamma\f$.
The data are the source term \f$\mathbf{j}\f$, the given velocity \f$\bu^\text{given}\f$,
 the boundary data \f$\bH_{\text{bdy}}\f$ and the initial data \f$\bH_0\f$.
 The parameter \f$\Rm\f$ is the magnetic Reynolds number.
 The parameter \f$\mu^c\f$ is the magnetic permeability of the conducting region \f$\Omega\f$.
 The parameter \f$\sigma\f$ is the conductivity in the conducting region \f$\Omega\f$.

<h3>Manufactured solutions</h3>

We define the magnetic permeability in the conducting region as follows:
@f{align*}
\mu^c(r,z)=
\begin{cases}
\frac{1}{ 1 + f_{23}(r,z) \cos(m_{23}\theta) } & \text{ in } \Omega_1, \\
\frac{1+\lambda_\mu/z}{ 1 + f_{23}(r,z) \cos(m_{23}\theta)}   &\text{ in } \Omega_2, 
\end{cases}
@f}
with \f$lambda_\mu= 1\f$, \f$m_{23}=3\f$ and \f$f_{23}\f$ defined as follows:
@f{align*}
f_{23}(r,z)=\frac{1}{0.00016}\frac{\text{ratio}_\mu-1}{\text{ratio}_\mu+1} \left(r(1-r)(r-2)(z-0.25)(z-1)\right)^3,
@f}
where \f$\text{ratio}_\mu=50\f$ represents the fraction of the minimum over the maximum of the magnetic permeability in the conducting region.

We approximate the following analytical solutions:
@f{align*}
H_r(r,\theta,z,t) &= 
\begin{cases}
 r ( 1 + f_{23}(r,z) \cos(m_{23}\theta))  & \text{ in } \Omega_1, \\
 \frac{r}{1+\lambda_\mu/z} ( 1 + f_{23}(r,z) \cos(m_{23}\theta)) &\text{ in } \Omega_2, 
\end{cases}
\\ H_{\theta}(r,\theta,z,t) &=0,
\\ H_z(r,\theta,z,t) &=
\begin{cases}
 -2 z ( 1 + f_{23}(r,z) \cos(m_{23}\theta))   & \text{ in } \Omega_1, \\
0  &\text{ in } \Omega_2. 
\end{cases}
@f}
and remind that \f$\bB=\mu\bH\f$.

 We also set the given velocity field to zero.
@f{align*}
u^\text{given}_r(r,\theta,z,t) &= 0, 
\\ u^\text{given}_{\theta}(r,\theta,z,t) &= 0, 
\\ u^\text{given}_z(r,\theta,z,t) &= 0.
@f}
The source term \f$ \mathbf{j}\f$ and the boundary data \f$\bH_{\text{bdy}}\f$ are computed accordingly.


<h3>Generation of the mesh</h3>
The finite element mesh used for this test is named <tt>mesh_18_0.025.FEM</tt>. 
 The mesh size for the P1 approximation is \f$0.025\f$ in \f$\Omega\f$. You can generate this mesh with the files in the following directory:
 ($SFEMaNS_MESH_GEN_DIR)/EXAMPLES/EXAMPLES_MANUFACTURED_SOLUTIONS/mesh_18_0.05.FEM (after editing the topology files).
The following images show the mesh for P1 finite elements for the conducting and vacuum regions.

<table width="70%" align="center">
<tr>
    <td align="center">
    @image html  fig_mesh_18_0.025.png  "Finite element mesh (P1) of the conducting region."
    </td>
</tr>
</table>



<h3>Information on the file <tt>condlim.f90</tt></h3>
The initial conditions, boundary conditions and the forcing term \f$\textbf{j}\f$ in the Maxwell
 equations are set in the file <tt>condlim_test_23.f90</tt>. This condlim also uses functions of the file <tt>test_23.f90</tt> to define the magnetic permeability and its gradient. Here is a description of the subroutines and functions of interest of the file <tt>condlim_test_23.f90</tt>.
<ol>
<li>The subroutine <tt>init_maxwell</tt> initializes the magnetic field at the time\f$-dt\f$ and \f$0\f$ with \f$dt\f$ being the time step. It is set to zero as follows:
\code
    Hn1 = 0.d0
    Hn = 0.d0
    phin1 = 0.d0
    phin = 0.d0
    time=0.d0
    RETURN
\endcode
<li>The function <tt>Hexact</tt> contains the analytical magnetic field.
 It is used to impose Dirichlet boundary conditions on \f$\textbf{H}\times\textbf{n}\f$ with \f$\textbf{n}\f$ the outter normal vector.
<ol>
<li>We define a tabular id so that each nodes n, id(n) is equal to 1 or 2 depending if the node is in \f$\Omega_1\f$ or \f$\Omega_2\f$.
<ol>
<li>If rr contains the radial and vertical coordinates of each np nodes, we define id as follows.
\code
    IF (SIZE(rr,2)== H_mesh%np) THEN
       DO mm = 1, H_mesh%me !id is used to determine on which side of the interface we are
          id(H_mesh%jj(:,mm)) = H_mesh%i_d(mm) !used for initialization
       END DO
\endcode
This case would be done if we initialize the magnetic field with Hexact.
<li>If rr does not contain the information on np nodes, it means Hexact is called to set boundary conditions. In that case, Hexact is called one time per nodes. So we define id, a tabular of dimension 1, as follows:
\code
    ELSE
       IF (rr(1,1)<1) THEN !used for boundary condition
          id = 1
       ELSE
          id = 2
       END IF
    END IF
\endcode
</ol>
<li>For the Fourier mode \f$m=0\f$ and \f$m=m_{23}\f$, we define the magnetic field depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) and the subdomain considered follows:
\code
    DO n = 1, SIZE(rr,2)
       IF (m == 0) THEN
          IF (TYPE==1) THEN 
             IF(id(n)==1) THEN
                vv(n)  =  rr(1,n)
             ELSE
                vv(n)  =  rr(1,n)/(1.d0 + lambda_mu_T23/rr(2,n) )
             END IF

          ELSE IF (TYPE==5) THEN
             vv(n)  =  -2.d0*rr(2,n)
          ELSE
             vv(n) = 0.d0
          END IF
       ELSE IF (m  ==  mode_mu_T23) THEN
          IF (TYPE==1) THEN 
             IF(id(n)==1) THEN
                vv(n)   = rr(1,n)*s_test_T23(rr(1,n),rr(2,n))   
             ELSE
                vv(n)   = rr(1,n)*s_test_T23(rr(1,n),rr(2,n))/(1.d0 + lambda_mu_T23/rr(2,n))
             END IF
          ELSE IF (TYPE==5) THEN
             vv(n)  = -2.d0*rr(2,n)*s_test_T23(rr(1,n),rr(2,n))   
          ELSE
             vv(n) = 0.d0
          END IF
       ELSE
          vv = 0.d0
       END IF
    END DO
    RETURN
\endcode
where s_test_T23, respectibely lambda_mu_T23, represents the function \f$f_{23}\f$, respectively \f$lambda_\mu\f$. They are defined in the file <tt>test_23.f90</tt>.
</ol>
<li>The function <tt>Vexact</tt> is used to define the given velocity field \f$\bu^\text{given}\f$. It is set to zero.
\code
   vv = 0.d0
\endcode
<li>The function <tt>Jexact_gauss</tt> is used to define the source term \f$\textbf{j}\f$.
 We remind that \f$ \bu^\text{given}=0 \f$ and \f$\partial_t \textbf{H}=0\f$ so it is defined such that \f$\textbf{j}=\ROT \bH\f$.
<ol>
<li>First we define the cylindrical coordinates (r,z) of the node considered.
\code
   r = rr(1)
   z = rr(2)
\endcode
<li>For the Fourier mode \f$m=0\f$, the azimuthal component of \f$\textbf{j}\f$ in the domain \f$\Omega_2\f$ is defined as follows:
\code
   !DCQ: mesh_id gives us now info about what side (domain) we are.
   IF ((m==0) .AND. (mesh_id==2) .AND. (TYPE==3)) THEN   
      !J_theta
      beta=lambda_mu_T23/((1.d0 + lambda_mu_T23/z)**2 * z**2)
      vv=r*beta
\endcode
<li>For the Fourier mode \f$m=m_{23}\f$, the source term is defined of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) and the subdomain considered (mesh_id=1 in \f$\Omega_1\f$ and 2 in \f$\Omega_2\f$).
\code
   ELSE IF (m==mode_mu_T23) THEN   
      !J_r
      IF (TYPE==2) THEN
         vv  = 2.d0*m*z*b_factor_T23*r**2*( (r-1.d0)*(r-2.d0)*(z-0.25)*(z-1.d0) )**3
         !J_theta
      ELSE IF (TYPE==3) THEN 
         alpha  =6*b_factor_T23*z*( (z-0.25)*(z-1) )**3*( r* (r-1.d0)*(r-2.d0) )**2 &
              *( r*(r-1)+r*(r-2)+(r-1)*(r-2))

         IF (mesh_id==1) THEN
            vv= 3*b_factor_T23*r**4*( (r-1)*(r-2) )**3 *( (z-0.25)*(z-1) )**2&
                 *( 2*z -1 - 0.25) + alpha
         ELSE          
            beta  =lambda_mu_T23/((1.d0 + lambda_mu_T23/z)**2 * z**2)

            vv= b_factor_T23*r**4*(  (r-1.d0)*(r-2.d0) )**3 *( (z-0.25)*(z-1) )**2 &
                 *( 3*(2*z- 1 - 0.25 )/(1.d0+lambda_mu_T23/z) + beta*(z-1)*(z-0.25) ) &
                 +  alpha
         END IF
         !J_z
      ELSE   IF (TYPE==6) THEN 

         IF (mesh_id==1)  THEN
            vv =   m*s_test_T23(r,z)
         ELSE
            vv =   m*s_test_T23(r,z)/(1.d0+lambda_mu_T23/z)
         ENDIF
      ELSE
         vv = 0.d0
      END IF
   ELSE
      vv = 0.d0
   END IF
   RETURN
\endcode
</ol>
<li>The function <tt>mu_bar_in_fourier_space</tt> defines a stabilization term, denoted mu_bar, on the nodes or gauss points of the finite element mesh. It needs to be smaller than mu everywhere in the domain. It is done by using the function <tt>mu_bar_in_fourier_space_anal_T23</tt> of the file <tt>test_23.f90</tt>.
\code
    IF( PRESENT(pts) .AND. PRESENT(pts_ids) ) THEN
       vv=mu_bar_in_fourier_space_anal_T23(H_mesh,nb,ne,pts,pts_ids) 
    ELSE
       vv=mu_bar_in_fourier_space_anal_T23(H_mesh,nb,ne) 
    END IF
    RETURN
\endcode
<li>The function <tt>grad_mu_bar_in_fourier_space</tt> defines the gradient of mu_bar. It is done by using the function <tt>grad_mu_bar_in_fourier_space_anal_T2</tt> of the file <tt>test_22.f90</tt>.
\code
    vv=grad_mu_bar_in_fourier_space_anal_T23(pt,pt_id) 
    RETURN
\endcode
When the magnetic permeability depends of the time or the azimuthal direction, this function is used to define the gradient of a stabilization term mu_bar.
<li>The function <tt>mu_in_real_space</tt> defines the magnetic permeability. It is done by using the function <tt>mu_in_real_space_anal_T23</tt> of the file <tt>test_23.f90</tt>.
\code
    vv = mu_in_real_space_anal_T23(H_mesh,angles,nb_angles,nb,ne)
    RETURN
\endcode
</ol>
All the other subroutines present in the file <tt>condlim_test_23.f90</tt> are not used in this test.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file. 
 The following describes the functions of interest of the file <tt>test_23.f90</tt>.
<ol>
<li>First we define the real numbers \f$ratio_\mu\f$, \f$b\f$, \f$\lambda_\mu\f$ and the interger \f$m_{23}\f$.
\code
  REAL (kind=8),PARAMETER, PUBLIC :: ratio_mu_T23 = 50.0d0 ! the variation of mu1
  REAL (kind=8),           PUBLIC :: b_factor_T23 = (1.d0/0.00016) * (ratio_mu_T23 - 1.d0)/(ratio_mu_T23 +1.d0)
  REAL (kind=8),           PUBLIC :: lambda_mu_T23 = 1.d0
  INTEGER,                 PUBLIC :: mode_mu_T23 = 3
\endcode
<li>The function <tt>s_test_T23</tt> defines the function \f$f_{23}\f$.
\code
    vv = b_factor_T23*(  r*(r-1.d0)*(r-2.d0)*(z-0.25)*(z-1.d0)  )**3
    RETURN
\endcode
<li>The function <tt>Ds_test_T23</tt> defines the gradient of the function \f$f_{23}\f$.
<ol>
<li>We compute the r-derivative as follows:
\code
    vv(1) = b_factor_T23*((z-0.25)*(z-1.d0))**3 * &
         (  3*( r*(r-1.d0)*(r-2.d0) )**2*( r*(  (r-1)+(r-2) ) +  (r-1)*(r-2)    ) )
\endcode
<li>We compute the z-derivative as follows:
\code
    vv (2) = b_factor_T23*(  r*(r-1.d0)*(r-2.d0))**3 * &
         ( 3* (z-0.25)*(z-1.d0)  )**2 *( (z-1.d0) + (z-0.25) ) 
    RETURN
\endcode
</ol>
<li>The function <tt>mu_bar_in_fourier_space_anal_T23</tt> defines stabilization term mu_bar such that it is smaller than mu everywhere in the domain. 
<ol>
<li>First we define the radial and vertical cylindrical coordinates of each nodes. We also define a tabular local_ids that is equal to 1 if the node is in \f$\Omega_1\f$ and 2 if it's in \f$\Omega_2\f$.
\code
    IF( PRESENT(pts) .AND. PRESENT(pts_ids) ) THEN !Computing mu at  pts
       r=pts(1,nb:ne)
       z=pts(2,nb:ne)
       local_ids=pts_ids
    ELSE
       r=H_mesh%rr(1,nb:ne) !Computing mu  at nodes
       z=H_mesh%rr(2,nb:ne)
       DO m = 1, H_mesh%me
          global_ids(H_mesh%jj(:,m)) = H_mesh%i_d(m)
       END DO
       local_ids=global_ids(nb:ne)
    END IF
\endcode
<li>We define the stabilization term mu_bar.
\code
    DO n  = 1, ne - nb + 1
       s   = s_test_T23(r(n),z(n))
       IF(local_ids(n)==1) THEN
          vv(n) = 1.d0/(1.d0 + ABS(s))  !mu1_bar, DCQ: If you change mu1_bar, you have to change
          !Jexact_gauss() as well
       ELSE
          vv(n) = 1.d0/(  (1.d0 + ABS(s)))*(1.d0 + lambda_mu_T23/z(n))   !mu2_bar
       END IF
    END DO
    RETURN
\endcode
</ol>
<li>The function <tt>grad_mu_bar_in_fourier_space_anal_T23</tt> defines the gradient of mu_bar on gauss points.
<ol>
<li>We define the radial and vertical cylindrical coordinates (r,z).
\code
    r=pt(1)
    z=pt(2)
\endcode
<li>We compute \f$f_{23}(r,z)\f$ and define a function sign depending of its value.
\code
    s=s_test_T23(r,z)
    IF (s .GE. 0.d0 ) THEN
       sign =1.0
    ELSE
       sign =-1.0
    END IF
\endcode
<li>We compute the gradient of \f$f_{23}\f$ in (r,z).
\code
    tmp=Ds_test_T23(r,z)!derivative
\endcode
<li>We define the gradient of mu_bar.
\code    
    IF(pt_id(1)==1) THEN !grad_mu_1
       vv(1)= 1.d0
       vv(2)= 0.d0  
    ELSE                 !grad_mu_2
       vv(1)=1.d0 +  ( (3*r+2)*(2*lambda_mu_T23)  -  ( (2*lambda_mu_T23*(1+r)))*(3)  )  /( z*(3*r+2) )**2
       vv(2)=  (2*lambda_mu_T23*(1+r))/(3*r+2)*(-2/z**3) 
    END IF
    RETURN
\endcode
</ol>
<li>The function <tt>mu_in_real_space_anal_T23</tt> defines the magnetic permeability (depending of the node in the meridian plan and its angle).
<ol>
<li> We define a tabular id that is equal to 1 if the node is in \f$\Omega_1\f$ and 2 if it's in \f$\Omega_2\f$.
\code
    DO m = 1, H_mesh%me
       id(H_mesh%jj(:,m)) = H_mesh%i_d(m)  !DCQ: Speed Efficient but requires more memory
    END DO
\endcode
<li>We define the magnetic permeability.
\code
    DO n = nb, ne
       n_loc = n - nb + 1
       tmp   = s_test_T23(H_mesh%rr(1,n),H_mesh%rr(2,n))
       DO ang = 1, nb_angles

          IF (id(n)==1) THEN
             vv(ang,n_loc)  = 1.d0/(1.d0 + tmp*COS(mode_mu_T23*angles(ang)) )!mu_1
          ELSE
             vv(ang,n_loc)  = 1.d0/( 1.d0 + tmp*COS(mode_mu_T23*angles(ang)) ) &
                  *(  1.d0 + lambda_mu_T23/(H_mesh%rr(2,n))  )  !mu_2
          ENDIF
       END DO
    END DO
\endcode
</ol>
</ol>




<h3>Setting in the data file</h3>
We describe the data file of this test.
 It is called <tt>debug_data_test_23</tt> and can be found in the following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<ol>
<li>We use a formatted mesh by setting:
\code
===Is mesh file formatted (true/false)?
.t.
\endcode
<li>The path and the name of the mesh are specified with the two following lines:
\code
===Directory and name of mesh file
'.' 'mesh_18_0.025.FEM'
\endcode
where '.' refers to the directory where the data file is, meaning ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<li>We use two processors in the meridian section. It means the finite element mesh is subdivised in two.
\code
===Number of processors in meridian section
2
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
As a consequence, the code approximates the problem on the first four Fourier modes.
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
<li>We use a time step of \f$1\f$ and solve the problem over \f$10\f$ time iterations.
\code  
===Time step and number of time iterations
1d0, 10
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
2
===List of subdomains for magnetic field (H) mesh
1 2
\endcode
<li>We set the number of interface in H_mesh and give their respective labels.
\code
===Number of interfaces in H mesh
1
===List of interfaces in H mesh
5
\endcode
Such interfaces represent interfaces with discontinuities in magnetic permeability or interfaces between the magnetic field mesh and the temperature or the velocity field meshes.
<li>We set the number of boundaries with Dirichlet conditions on the magnetic field and give their respective labels.
\code
===Number of Dirichlet sides for Hxn
3
===List of Dirichlet sides for Hxn
2 4 3
\endcode
<li>The magnetic permeability is defined with a function of the file <tt>condlim_test_23.f90</tt>.
\code
===Is permeability defined analytically (true/false)?
.t.
\endcode
<li>We construct a stablization term, denoted mu_bar, on the gauss points by a finite element interpolation of its value on the nodes. This stabilization term is defined in the file <tt>condlim_test_23.f90</tt>.
\code
===Use FEM Interpolation for magnetic permeability (true/false)?
.t.
\endcode
<li>The magnetic permeability is variable in theta.
\code
===Is permeability variable in theta (true/false)?
.t.
\endcode
<li>We set the conductivity in each domains where the magnetic field is approximated.
\code
===Conductivity in the conductive part (1:nb_dom_H)
1.d0 1.d0
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
<li>The L2 norm of the magnetic filed \f$\bH\f$.
</ol>

 They are compared to reference values to attest of the correctness of the code.  
 These values of reference are in the last lines of the file <tt>debug_data_test_23</tt> in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
============================================
mesh_18_0.025.FEM, dt=1d0, it_max=10
===Reference results
5.117459706179962E-003   !L2 norm of div(mu (H-Hexact))
0.232311657766681        !L2 norm of curl(H-Hexact)
9.199612567417106E-003   !L2 norm of H-Hexact
4.61225833710813         !L2 norm of H
\endcode


To conclude this test, we show the profile of the approximated  magnetic field \f$\textbf{H}\f$ at the final time.
 These figures are done in the plane \f$y=0\f$ which is 
 the union of the half plane \f$\theta=0\f$ and \f$\theta=\pi\f$.
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_23_mag_tfin.png "Magnetic field magnitude in the plane plane y=0."
    </td>
</tr>
</table>

 */
