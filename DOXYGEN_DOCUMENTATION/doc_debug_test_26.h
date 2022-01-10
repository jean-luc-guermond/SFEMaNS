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
 * @page doc_debug_test_26 Test 26: Maxwell Equations with variable magnetic permeability in z and smooth jump.




<h3>Introduction</h3>
In this example, we check the correctness of SFEMaNS for stationary magnetic problem involving a conducting. The magnetic permeability is variable in z and presents a smooth jump. The magnetic permeability is given. We consider Dirichlet boundary conditions. We use P2 finite elements. This set up is approximated over one time iteration with a large time step.


We solve the Maxwell equations with the magnetic field \f$\bB=\mu \bH\f$ as dependent variable:
\f{align}{
\begin{cases}
\partial_t (\mathbf{B}) + \nabla \times \left(\frac{1}{\Rm \sigma}\nabla \times( \frac{1}{\mu^c}\mathbf{B} ) \right)  =
 \nabla\times (\bu^\text{given} \times  \mathbf{B}) 
+ \nabla \times \left(\frac{1}{\Rm \sigma} \nabla\times \mathbf{j} \right) & \text{in } \Omega, \\
\text{div} (\mathbf{B}) = 0  & \text{in } \Omega ,\\
{\bH \times \bn}_{|\Gamma} = {\bH_{\text{bdy}} \times \bn}_{|\Gamma}, &\\
\bH_{|t=0}= \bH_0, \\
\end{cases}
\f}
in the domain \f$\Omega= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,0.5] \times [0,2\pi) \times [-0.2,0.2]\} \f$. We set \f$\Gamma= \partial \Omega\f$.  We also define the outward normal \f$\bn\f$ to the surface \f$\Gamma\f$.
The data are the source term \f$\mathbf{j}\f$, the given velocity \f$\bu^\text{given}\f$,
 the boundary data \f$\bH_{\text{bdy}}\f$ and the initial data \f$\bH_0\f$.
 The parameter \f$\Rm\f$ is the magnetic Reynolds number.
 The parameter \f$\mu^c\f$ is the magnetic permeability of the conducting region \f$\Omega\f$.
 The parameter \f$\sigma\f$ is the conductivity in the conducting region \f$\Omega\f$.

<h3>Manufactured solutions</h3>

We define the magnetic permeability in the conducting region as follows:
@f{align*}
\mu^c(z)=
\begin{cases}
  \mu_{26}& \text{ for } z\leq z_1, \\
1+ (\mu_{26}-1)F_z(z)  &\text{ for } z_1 \leq z \leq z_0, \\
 1 &\text{ for } z_0 \leq z,
\end{cases}
@f}
where \f$\mu_{26}=5\f$, \f$z_0=0\f$, \f$z_1=-0.032\f$ and \f$F_z\f$ is the polynomial function of order 2 such that \f$F_z(z_0)=0\f$, \f$F_z(z_1)=1\f$ and \f$\partial_z F_z(\frac{z_0+z_1}{2})=0\f$.

We approximate the following analytical solutions:
@f{align*}
H_r(r,\theta,z,t) &= \frac{\alpha_{26}m_{26} r^{m_{26}-1}}{\mu^c(z)}\cos(m_{26}\theta),
\\ H_{\theta}(r,\theta,z,t) &=-\frac{\alpha_{26}m_{26} r^{m_{26}-1}}{\mu^c(z)}\sin(m_{26}\theta),
\\ H_z(r,\theta,z,t) &= 0,
@f}
with \f$\alpha_{26}=1\f$ and \f$m_{26}=4\f$. We also remind that \f$\bB=\mu\bH\f$.

 We also set the given velocity field to zero.
@f{align*}
u^\text{given}_r(r,\theta,z,t) &= 0, 
\\ u^\text{given}_{\theta}(r,\theta,z,t) &= 0, 
\\ u^\text{given}_z(r,\theta,z,t) &= 0.
@f}
The source term \f$ \mathbf{j}\f$ and the boundary data \f$\bH_{\text{bdy}}\f$ are computed accordingly.


<h3>Generation of the mesh</h3>
The finite element mesh used for this test is named <tt>mesh_26_0.03.FEM</tt>. 
 The mesh size for the P1 approximation is \f$0.03\f$ in \f$\Omega\f$. You can generate this mesh with the files in the following directory:
 ($SFEMaNS_MESH_GEN_DIR)/EXAMPLES/EXAMPLES_MANUFACTURED_SOLUTIONS/mesh_26_0.03.FEM.
The following images show the mesh for P1 finite elements for the conducting and vacuum regions.

<table width="70%" align="center">
<tr>
    <td align="center">
    @image html  fig_mesh_26_0.03.png  "Finite element mesh (P1) of the conducting region."
    </td>
</tr>
</table>




<h3>Information on the file <tt>condlim.f90</tt></h3>
The initial conditions, boundary conditions and the forcing term \f$\textbf{j}\f$ in the Maxwell
 equations are set in the file <tt>condlim_test_26.f90</tt>. This condlim also uses functions of the file <tt>test_26.f90</tt> to define the magnetic permeability and its gradient. Here is a description of the subroutines and functions of interest of the file <tt>condlim_test_26.f90</tt>.
<ol>
<li>The subroutine <tt>init_maxwell</tt> initializes the magnetic field at the time\f$-dt\f$ and \f$0\f$ with \f$dt\f$ being the time step. It is set to zero.
\code
    Hn1 = 0.d0
    Hn = 0.d0
    phin1 = 0.d0
    phin = 0.d0
    time =0.d0
    RETURN
\endcode
<li>The function <tt>Hexact</tt> contains the analytical magnetic field.
 It is used to impose Dirichlet boundary conditions on \f$\textbf{H}\times\textbf{n}\f$ with \f$\textbf{n}\f$ the outter normal vector.
<ol>
<li>If the Fourier mode m is not equal to \f$m_{26}\f$, we set the magnetic field to zero.
\code
    vv  = 0.d0
    IF (m .NE. test_mode_T26 ) THEN
       RETURN
    ENDIF
\endcode
<li>We define the radial and vertical coordinates \f$(r, z)\f$.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>For the Fourier mode \f$m=m_{26}\f$, we define the magnetic field depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) as follows:
\code
    DO n = 1, SIZE(rr,2)       
       tmp= alpha_T26*test_mode_T26*r(n)**(test_mode_T26-1)/mu_func_T26(r(n),z(n))
       IF (type==1) THEN
          vv(n) =  tmp
       ELSE IF (type==4) THEN
          vv(n) = -tmp
       ENDIF
    END DO
    RETURN
\endcode
where alpha_T26, test_mode_T26 are parameters defined in the file <tt>test_26.f90</tt>. The function mu_func_T26 is also defined in the file <tt>test_26.f90</tt>, it is equal to the magnetic permeability.
</ol>
<li>The function <tt>Vexact</tt> is used to define the given velocity field \f$\bu^\text{given}\f$. It is set to zero.
\code
   vv = 0.d0
   RETURN
\endcode
<li>The function <tt>Jexact_gauss</tt> is used to define the source term \f$\textbf{j}\f$.
 We remind that \f$ \bu^\text{given}=0 \f$ and \f$\partial_t \textbf{H}=0\f$ so it is defined such that \f$\textbf{j}=\ROT \bH\f$.
<ol>
<li>If the Fourier mode m is not equal to \f$m_{26}\f$, we set the \f$\textbf{j}\f$ to zero.
\code
   vv = 0.d0
   IF (m .NE. test_mode_T26) THEN
      RETURN
   ENDIF
\endcode
<li>We define the radial and vertical coordinates \f$(r,z)\f$.
\code
   r = rr(1)
   z = rr(2)
\endcode
<li>We define the z-derivative of the magnetic permeability (only depends of z) with the function <tt>Dmu_func_T26</tt> defined in the file <tt>test_26.f90</tt>.
\code
   Dmu=Dmu_func_T26(r,z)
\endcode
<li>For the Fourier mode \f$m=m_{26}\f$, we define the source term depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) as follows:
\code
   tmp=alpha_T26*test_mode_T26*r**(test_mode_T26-1)*(-1.d0/mu_func_T26(r,z)**2)
   IF (type==2) THEN
      vv =     tmp*Dmu(2)
   ELSE IF (type==3) THEN
      vv  =    tmp*Dmu(2)
   ELSE IF (type==6) THEN
      vv  =  - tmp*Dmu(1)
   ENDIF
   RETURN
\endcode
</ol>
<li>The function <tt>mu_bar_in_fourier_space</tt> defines the magnetic permeability that depends of the radial and vertical coordinates. It is done by using the function <tt>mu_bar_in_fourier_space_anal_T26</tt> of the file <tt>test_26.f90</tt>.
\code
    IF( PRESENT(pts) .AND. PRESENT(pts_ids) ) THEN
       vv=mu_bar_in_fourier_space_anal_T26(H_mesh,nb,ne,pts,pts_ids) 
    ELSE
       vv=mu_bar_in_fourier_space_anal_T26(H_mesh,nb,ne,pts) 
    END IF
    RETURN
\endcode
When the magnetic permeability depends of the time or the azimuthal direction (set "===Is permeability variable in theta (true/false)?" to true in the data file), this function is used to define a stabilization term mu_bar. It needs to be smaller than mu everywhere in the domain. The magnetic permeability is then defined with the function <tt>mu_in_real_space</tt>.
<li>The function <tt>grad_mu_bar_in_fourier_space</tt> defines the gradient of the magnetic permeability. It is done by using the function <tt>grad_mu_bar_in_fourier_space_anal_T26</tt> of the file <tt>test_26.f90</tt>.
\code
    vv=grad_mu_bar_in_fourier_space_anal_T26(pt,pt_id)
    RETURN
\endcode
When the magnetic permeability depends of the time or the azimuthal direction, this function is used to define the gradient of a stabilization term mu_bar.
</ol>
All the other subroutines present in the file <tt>condlim_test_26.f90</tt> are not used in this test.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file. 
 The following describes the functions of interest of the file <tt>test_26.f90</tt>.
<ol>
<li>First we define the real numbers \f$\mu_{26}\f$, \f$\alpha_{26}\f$ and the integer \f$m_{26}\f$. 
\code
  REAL(KIND=8),  PARAMETER:: mu_disk_T26 = 5.d0
  REAL(KIND=8),  PARAMETER:: alpha_T26 = 1.d0
  INTEGER,  PARAMETER     :: test_mode_T26 = 4;
\endcode
<li>We also define the vertical coordinates \f$z_0\f$ and \f$z_1\f$ (denoted zm_T26 and z1_T26).
\code
  REAL(KIND=8),  PARAMETER:: wjump_T26 = 0.032d0*(1.0d0)
  REAL(KIND=8),  PARAMETER:: zm_T26 = 0.d0, z1_T26 = zm_T26-wjump_T26
\endcode
We remind that the magnetic permeability is only variable when \f$z_1\leq z \leq z_0\f$. 
<li>The function <tt>smooth_jump_down_T26</tt> depends of the variable x and is defined such that it is equal to 1 for x= x0 and 0 for x=x1. Moreover its derivative in \f$x=\frac{x0+x1}{2}\f$ is zero.
\code
    a0 = x1**2*(3*x0-x1)/(x0-x1)**3; 
    a1 = -6.0*x0*x1/(x0-x1)**3; 
    a2 = (3.0*(x0+x1))/(x0-x1)**3;
    a3 = -2.0/(x0-x1)**3;

    vv = a0+a1*x+a2*x*x + a3*x*x*x
    RETURN
\endcode
<li>The function <tt> Dsmooth_jump_down_T26 </tt> is the derivative of the function <tt>smooth_jump_down_T26</tt>.
\code
    a0 = x1**2*(3*x0-x1)/(x0-x1)**3; 
    a1 = -6.0*x0*x1/(x0-x1)**3; 
    a2 = (3.0*(x0+x1))/(x0-x1)**3;
    a3 = -2.0/(x0-x1)**3;

    vv = a1+2.d0*a2*x + 3.d0*a3*x*x 
    RETURN
\endcode
<li>The function <tt>smooth_jump_up_T26</tt> depends of the variable x and is defined such that it is equal to 1 for x= x1 and 0 for x=x0. Moreover its derivative in \f$x=\frac{x0+x1}{2}\f$ is zero.
\code
    vv = 1.d0 - smooth_jump_down_T26( x,x0,x1 );
    RETURN
\endcode
<li>The function <tt> Dsmooth_jump_up_T26 </tt> is the derivative of the function <tt>smooth_jump_up_T26</tt>.
\code
    vv =  - Dsmooth_jump_down_T26( x,x0,x1 );
    RETURN
\endcode
<li>The function <tt>mu_func_T26</tt> is used to defined the magnetic permeability.
\code
    Fz=smooth_jump_up_T26(z,zm_T26,z1_T26)
    IF ( z.GE.zm_T26 ) THEN       
       vv  = 1.d0
    ELSE  IF ( z.LE. z1_T26 ) THEN
       vv  = mu_disk_T26
    ELSE   
       vv = Fz*(mu_disk_T26-1.d0) + 1.d0                    
    END IF
    RETURN
\endcode
One can note that the result of this function is exactly the magnetic permeability considerer in this test.
<li>The function <tt>Dmu_func_T26</tt> is the gradient of the function <tt>mu_func_T26</tt>.
\code
    DFz=Dsmooth_jump_up_T26(z,zm_T26,z1_T26)
    IF ( z.GE.zm_T26 ) THEN       
       vv  = 0.d0
    ELSE  IF ( z.LE. z1_T26 ) THEN
       vv  = 0.d0
    ELSE   
       vv(1)  = 0
       vv(2)  = DFz*(mu_disk_T26-1.d0)
    END IF
    RETURN
\endcode
We note that the r-derivative is zero because mu_func_T26 only depends of the z coordinate.
<li>The function <tt>mu_bar_in_fourier_space_anal_T26</tt> defines the magnetic permeability.
<ol>
<li>First we define the radial and vertical cylindrical coordinates of each nodes considered.
\code
    IF( PRESENT(pts) .AND. PRESENT(pts_ids) ) THEN !Computing mu at  pts
       r=pts(1,nb:ne)
       z=pts(2,nb:ne)
    ELSE
       r=H_mesh%rr(1,nb:ne) !Computing mu  at nodes
       z=H_mesh%rr(2,nb:ne)
    END IF
\endcode
<li>We define the magnetic permeability on the nodes considered.
\code
    DO n = 1, ne - nb + 1
       vv(n) = mu_func_T26(r(n),z(n))
    END DO
    RETURN
\endcode
</ol>
<li>The function <tt>grad_mu_bar_in_fourier_space_anal_T18</tt> defines the gradient of the magnetic permeability on gauss points.
<ol>
<li>We define the radial and vertical cylindrical coordinates.
\code
    r=pt(1)
    z=pt(2)
\endcode
<li>We define the gradient of the magnetic permeability.
\code    
    tmp=Dmu_func_T26(r,z)
    vv(1)=tmp(1)
    vv(2)=tmp(2)
\endcode
</ol>
</ol>





<h3>Setting in the data file</h3>
We describe the data file of this test.
 It is called <tt>debug_data_test_18</tt> and can be found in the following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<ol>
<li>We use a formatted mesh by setting:
\code
===Is mesh file formatted (true/false)?
.t.
\endcode
<li>The path and the name of the mesh are specified with the two following lines:
\code
===Directory and name of mesh file
'.' 'mesh_26_0.03.FEM'
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
<li>We select specific Fourier modes to solve.
\code
===Select Fourier modes? (true/false)
.t.
\endcode
<li>We give the list of the Fourier modes to solve.
\code
===List of Fourier modes (if select_mode=.TRUE.)
4
\endcode
It means the problem is only solved for the Fourier mode \f$m=4\f$.
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
1.d40 1
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
<li>We set the number of boundaries with Dirichlet conditions on the magnetic field and give their respective labels.
\code
===Number of Dirichlet sides for Hxn
3
===List of Dirichlet sides for Hxn 
4 2 3
\endcode
<li>The magnetic permeability is defined with a function of the file <tt>condlim_test_18.f90</tt>.
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
 These values of reference are in the last lines of the file <tt>debug_data_test_26</tt> in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
============================================
mesh_26_0.03.FEM, dt=1.0d40, it_max=1
===Reference results
6.232118631975291E-003     !L2 norm of div(mu (H-Hexact))
1.402063664918604E-003     !L2 norm of curl(H-Hexact)
5.555512428427833E-005     !L2 norm of H-Hexact
0.102653049435891          !L2 norm of H
\endcode


To conclude this test, we show the profile of the approximated  magnetic field \f$\textbf{H}\f$ at the final time.
 These figures are done in the plane \f$y=0\f$ which is 
 the union of the half plane \f$\theta=0\f$ and \f$\theta=\pi\f$.
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_26_mag_tfin.png "Magnetic field magnitude in the plane plane y=0."
    </td>
</tr>
</table>

*/
