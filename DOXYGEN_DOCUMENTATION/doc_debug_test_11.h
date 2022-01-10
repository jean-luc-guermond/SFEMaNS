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
 * @page doc_debug_test_11 Test 11: MHD with Temperature and precession

<h3>Introduction</h3>
In this example, we check the correctness of SFEMaNS for a thermo-magnetohydrodynamic problem with a rotating force around the vertical axis. This test uses Dirichlet boundary conditions. We note this test does not involve manufactured solutions and consists of checking four quantities, like the \f$\bL^2\f$ norm of the velocity, are the same than the values of reference.

We solve the temperature equation:
@f{align*}
\partial_t T+ \bu \cdot \GRAD T - \kappa \LAP T &= 0, \\
T_{|\Gamma} &= T_\text{bdy} , \\
T_{|t=0} &= T_0,
@f}
the Navier-Stokes equations:
@f{align*}
\partial_t\bu+\left(\ROT\bu + 2 \epsilon \textbf{e}_z \right)\CROSS\bu - \frac{1}{\Re}\LAP \bu +\GRAD p     
&= \alpha T (r\textbf{e}_r+z\textbf{e}_z) + (\ROT \textbf{H}) \times (\mu^c \textbf{H}),
\\ \DIV \bu &= 0, \\
\bu_{|\Gamma} &= \bu_{\text{bdy}} , \\
\bu_{|t=0} &= \bu_0, \\
p_{|t=0} &= p_0,
@f}
and the Maxwell equations:
@f{align*}
\partial_t (\mu^c \mathbf{H}) + \nabla \times \left(\frac{1}{\Rm \sigma}\nabla \times \mathbf{H}  \right) 
& = \nabla\times (\bu \times \mu^c \mathbf{H}), \\
\text{div} (\mu^c \mathbf {H}) &= 0   ,\\
{\bH \times \bn}_{|\Gamma} &= {\bH_{\text{bdy}} \times \bn}_{|\Gamma},\\
\bH_{|t=0}&= \bH_0.
@f}
Theses equations are solved in the domain \f$\Omega= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,1] \times [0,2\pi) \times [0,1] \; | \; R_i \leq \sqrt{r^2+z^2} \leq R_o\} \f$ with \f$R_i=7/13\f$, \f$R_o=20/13\f$ and  \f$\Gamma= \partial \Omega\f$.
 We denote by \f$\textbf{e}_r\f$, respectively \f$\textbf{e}_z\f$, the unit vector in the radial direction, respectively vertical direction. We solve these equation in the mantle frame, meaning we set homogeneous Dirichlet boundary conditions of the velocity field. It explains the presence of the term \f$ 2 \epsilon \textbf{e}_z \CROSS\bu  \f$ in the Navier-Stokes equations where \f$\epsilon\f$ is the angular velocity.
The data are the boundary datas \f$T_\text{bdy}, \bu_{\text{bdy}}, \bH_{\text{bdy}}\f$,
 the initial data \f$T_0, \bu_0, p_0, \bH_0\f$.
 The parameters are the thermal diffusivityr \f$\kappa\f$, the kinetic Reynolds number \f$\Re\f$, the thermal gravity number \f$\alpha\f$, the magnetic Reynolds number \f$\Rm\f$, the magnetic permeability \f$\mu^c\f$ and the conductivity  \f$\sigma\f$ of the fluid.




<h3>Manufactured solutions</h3>
As mentionned earlier this test does not involve manufactured solutions. As a consequence, we do not consider specific source term and only initialize the variables to approximate.

 The temperature is iniatialized as follows:
@f{align*}
T(r,\theta,z, t=0) & = \frac{R_i R_o}{\sqrt{r^2 + z^2}} - R_i + \frac{21}{\sqrt{12920 \pi}}
(1- 3x^2 +3x^4 -x^6)\sin(
\varphi)^4 \sin(4\theta) ,
@f}
where \f$ x = 2 \sqrt{r^2+z^2} - R_i - R_o \f$ and \f$\varphi=atan2(r,z)\f$ is the polar angle of the spherical coordinate system.

The initial velocity field and pressure are set to zero.
@f{align*}
 u_r(r,\theta,z,t=0) &= 0, 
\\ u_{\theta}(r,\theta,z,t=0) &= 0, 
\\ u_z(r,\theta,z,t=0) &= 0, 
\\ p(r,\theta,z,t=0) &= 0.
@f}

The magnetic field is iniatialized as follows:
@f{align*}
\\ H_r(r,\theta,z,t=0) &= 2 \sqrt{10^{-4}} \left( B_1(r,\theta,z) \sin(\varphi) + B_2(r,\theta,z) \cos(\varphi) \right),
\\ H_{\theta}(r,\theta,z,t=0) &= 2 \sqrt{10^{-4}} \sin(2\varphi)\frac{15}{8\sqrt{2}} 
\sin\left(\pi(\sqrt{r^2+z^2} - R_i)\right),
\\ H_z(r,\theta,z,t=0) &= 2 \sqrt{10^{-4}}\left( B_1(r,\theta,z) \cos(\varphi) - B_2(r,\theta,z) \sin(\varphi) \right),
@f}
where we define \f$\varphi\f$ as previously while \f$ B_1\f$ and \f$B_2\f$ are defined as follows:
@f{align*}
B_1(r,\theta,z) &= \cos(\varphi)\frac{5}{8\sqrt{2}} \frac{\left( -48 R_i R_o + (4 R_o + R_i(4+3 R_o)) 6 \sqrt{r^2+z^2}
- 4(4+3(R_i+R_o))(r^2+z^2) + 9 (r^2+z^2)^{3/2}\right)}{\sqrt{r^2+z^2}},
 \\ B_2(r,\theta,z) &= \sin(\varphi)\frac{-15}{4\sqrt{2}}
\frac{(\sqrt{r^2+z^2}-R_i)(\sqrt{r^2+z^2}-R_o)(3\sqrt{r^2+z^2}-4)}{\sqrt{r^2+z^2}}.
@f}
The boundary datas \f$T_\text{bdy}, \bu_{\text{bdy}}, \bH_{\text{bdy}}\f$ are computed accordingly.


<h3>Generation of the mesh</h3>

The finite element mesh used for this test is named <tt>Mesh_BENCH1_20.FEM</tt>. 
 The mesh size is \f$0.05\f$ for the P1 approximation. You can generate this mesh with the files in the following directory:
 ($SFEMaNS_MESH_GEN_DIR)/EXAMPLES/EXAMPLES_MANUFACTURED_SOLUTIONS/Mesh_BENCH1_20.
The following image show the mesh for P1 finite elements.
<table width="70%" align="center">
<tr>
    <td align="center">
    @image html  fig_Mesh_BENCH1_20.png  "Finite element mesh (P1)."
    </td>
</tr>
</table>


<h3>Information on the file <tt>condlim.f90</tt></h3>
The initial conditions, boundary conditions and the forcing terms are set in the file <tt>condlim_test_11.f90</tt>.
Here is a description of the subroutines and functions of interest.
<ol>
<li>First we define the numbers \f${test11}_{Ri}\f$, \f${test11}_{Ro}\f$, \f${test11}_{Rossby}\f$ and \f$\pi\f$ at the begining of the module so that every subroutines has access to these real numbers.
\code
  REAL(KIND=8), PRIVATE :: test11_Ri=7/13d0, test11_Ro=20/13d0, test11_Rossby=1.d-4
  REAL(KIND=8) :: pi=ACOS(-1.d0)
\endcode
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
<li>The subroutine <tt>init_temperature</tt> initializes the temperature at the time\f$-dt\f$ and \f$0\f$ with \f$dt\f$ being the time step. It is done by using the function temperature_exact as follows:
\code
    time = 0.d0
    DO i= 1, SIZE(list_mode)
       mode = list_mode(i) 
       DO j = 1, 2 
          tempn_m1(:,j,i) = temperature_exact(j, mesh%rr, mode, time-dt)
          tempn   (:,j,i) = temperature_exact(j, mesh%rr, mode, time)
       ENDDO
    ENDDO
\endcode
<li>The function <tt>vv_exact</tt> contains the analytical velocity field.
 It is used to initialize the velocity field and to impose Dirichlet boundary
 conditions on the velocity field. It is set to zero.
\code
    vv = 0.d0
    RETURN
\endcode
<li>The function <tt>pp_exact</tt> contains the analytical pressure.
 It is used to initialize the pressure and is set to zero.
\code
    vv=0.d0
    RETURN
\endcode
<li>The function <tt>temperature_exact</tt> contains the analytical temperature. It is used to initialize the temperature and to impose Dirichlet boundary condition on the temperature.
<ol>
<li>We construct the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li> For the Fourier modes \f$m\in \{0,4\}\f$ we define the temperature depending of its TYPE (1 for cosine and 2 for sine) as follows:
\code
    IF (m==0 .OR. m==4) THEN 
       DO n = 1, SIZE(rr,2)
          rho(n)=SQRT(rr(1,n)**2+rr(2,n)**2)
          theta=ATAN2(rr(1,n),rr(2,n))
          x=2*rho(n) - test11_Ri - test11_Ro
          A(n)=(21/SQRT(17920*Pi))*(1-3*x**2+3*x**4-x**6)*SIN(theta)**4 
       END DO
       IF (m==0 .AND. TYPE==1) THEN
          vv=( test11_Ri*test11_Ro/rho)- test11_Ri
       ELSE IF (m==4 .AND. TYPE==1) THEN
          vv= A
       ELSE 
          vv = 0.d0
       END IF
\endcode
<li> For the other Fourier modes, we set the temperature to zero.
\code
    ELSE
       vv = 0.d0
    END IF
    RETURN
\endcode
</ol>
<li>The function <tt>source_in_temperature</tt> computes the source term, denoted \f$f_T\f$ in previous tests, of the temperature equation. As it is not used in this test, we set it to zero.
\code
    vv = 0.d0
    RETURN
\endcode
<li>The function <tt>source_in_NS_momentum</tt> computes the source term \f$\alpha T(r\textbf{e}_r+z\textbf{e}_z) \f$ of the Navier-Stokes equations.
<ol>
<li>We construct the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>We compute the term  \f$\alpha T r\textbf{e}_r\f$ (component radial so TYPE 1 and 2).
\code    
    IF (TYPE==1) THEN
       vv = inputs%gravity_coefficient*opt_tempn(:,1,i)*r
    ELSE IF (TYPE==2) THEN
       vv = inputs%gravity_coefficient*opt_tempn(:,2,i)*r
\endcode
<li>We compute the term  \f$\alpha T z\textbf{e}_z\f$ (component vertical so TYPE 5 and 6).
\code
    ELSE IF (TYPE==5) THEN
       vv = inputs%gravity_coefficient*opt_tempn(:,1,i)*z
    ELSE IF (TYPE==6) THEN
       vv = inputs%gravity_coefficient*opt_tempn(:,2,i)*z
\endcode
<li>The azimuthal component of this source term is set to zero.
\code
    ELSE
       vv = 0.d0
    END IF
    RETURN
\endcode
</ol>
<li>The subroutine <tt>init_maxwell</tt> initializes the magnetic field
at the time\f$-dt\f$ and \f$0\f$ with \f$dt\f$ being the time step. We remind the initial magnetic field only depends of the Fourier mode \f$m=0\f$. Firslty we define it for spherical coordinates (Brho, Btheta, Bphi) so we can later define for cylindrical coordinates depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine).
\code
    Hn  = 0.d0
    DO i = 1, SIZE(list_mode)
       IF (list_mode(i) /= 0) CYCLE
       normalization = 2*SQRT(test11_Rossby)
       DO n=1,SIZE(H_mesh%rr,2)
          rho   = SQRT(H_mesh%rr(1,n)**2+H_mesh%rr(2,n)**2)
          theta(n) = ATAN2(H_mesh%rr(1,n),H_mesh%rr(2,n))
          Brho(n)  = COS(theta(n))*5.d0/(SQRT(2d0)*8d0)*(-48d0*test11_Ri*test11_Ro &
               + (4*test11_Ro+test11_Ri*(4+3*test11_Ro))*6*rho &
               -4*(4+3*(test11_Ri+test11_Ro))*rho**2+9*rho**3)/rho
          Btheta(n)= SIN(theta(n))*(-15d0/(SQRT(2d0)*4d0))*((rho-test11_Ri)*(rho-test11_Ro)*(3*rho-4))/rho
          Bphi(n)  = SIN(2*theta(n))*(15d0/(SQRT(2d0)*8d0))*SIN(Pi*(rho-test11_Ri)) 
       END DO
       Hn(:,1,i) = normalization*(Brho*SIN(theta) + Btheta*COS(theta))
       Hn(:,3,i) = normalization*Bphi
       Hn(:,5,i) = normalization*(Brho*COS(theta) - Btheta*SIN(theta))
    END DO
    Hn1 = Hn
    time = 0.d0
    phin =0.d0  
    phin1=0.d0
    RETURN
\endcode
The scalar potential \f$\phi\f$ is initialized but we note it is not used in this example. The variable Brho and Btheta plays the role of the variable \f$B_1\f$ and \f$B_2\f$ involves in the above definition of the magnetic field for a time \f$t=0\f$. One can check that the radial and vertical component of this initial magnetic field is zero on the border of the domain (if you set rho to test11_Ri and test11_Ro you will get a Brho and Btheta equal to zero).
<li>The function <tt>Hexact</tt> contains the analytical magnetic field.
 It is used to impose Dirichlet boundary
 conditions on \f$\textbf{H}\times\textbf{n}\f$ with \f$\textbf{n}\f$ the outter normal vector.
\code
    vv=0.d0
    RETURN
\endcode
<li>The function <tt>Jexact_gauss</tt> is used to define the source term \f$\textbf{j}\f$ that is not used in this test. So we set it to zero.
\code
   vv=0.d0
   RETURN
\endcode
</ol>
All the other subroutines present in the file <tt>condlim_test_11.f90</tt> are not used in this test.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file.




<h3>Setting in the data file</h3>
We describe the data file of this test.
It is called <tt>debug_data_test_11</tt> and can be found in the following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<ol>
<li>We use a formatted mesh by setting:
\code
===Is mesh file formatted (true/false)?
.t.
\endcode
<li>The path and the name of the mesh are specified with the two following lines:
\code
===Directory and name of mesh file
'.' 'Mesh_BENCH1_20.FEM'
\endcode
where '.' refers to the directory where the data file is, meaning ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<li>We use two processors in the meridian section. It means the finite element mesh is subdivised in two. To do so, we write:
\code
===Number of processors in meridian section
2
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
0 4 8 
\endcode
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
<li>We use a time step of \f$0.02\f$ and solve the problem over \f$20\f$ time iterations.
\code
===Time step and number of time iterations
2d-2, 20
\endcode
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
2
===List of boundary pieces for full Dirichlet BCs on velocity
2 4
\endcode
<li>We set the kinetic Reynolds number \f$\Re\f$.
\code
===Reynolds number
1000d0
\endcode
<li>We want to take into account a rotating force of axis \f$\textbf{e}_z\f$ and angular velocity \f$\epsilon\f$ in the wall frame. So we add the term \f$ 2\epsilon \textbf{e}_z \times \bu \f$ in the left hand side of the Navier-Stokes equations. Such features is already programmed in SFEMaNS via the following lines in your data file.
<ol>
<li>We set the option precession to true.
\code
===Is there a precession term (true/false)?
.t.
\endcode
It adds the term \f$ 2\epsilon \left(\sin(\alpha\pi)\textbf{e}_x+\cos(\alpha\pi)\textbf{e}_z\right) \times \bu \f$ in the left hand side of the Navier-Stokes equations with \f$\textbf{e}_x\f$ the normal vector associated to the x cartesian coordinate.
<li>We set the precession rate \f$\epsilon\f$ and the precesion angle \f$\alpha\f$.
\code
===Precession rate
1.d0
===Precession angle over pi
0.d0
\endcode
</ol>
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
<li>We solve the temperature equation (in the same domain than the Navier-Stokes equations).
\code
===Is there a temperature field?
.t.
\endcode
<li>We set the number of domains and their label, see the files associated to the generation of the mesh,
 where the code approximated the temperature equation.
\code
===Number of subdomains in temperature mesh
1
===List of subdomains for temperature mesh
1
\endcode
<li>We set the thermal diffusivity \f$\kappa\f$.
\code
===Diffusivity coefficient for temperature
1.d-3
\endcode
<li>We set the thermal gravity number \f$\alpha\f$.
\code
===Non-dimensional gravity coefficient
0.065 
\endcode
<li>We set the number of boundaries with Dirichlet conditions on the temperature and give their respective labels.
\code
===How many boundary pieces for Dirichlet BCs on temperature?
2
===List of boundary pieces for Dirichlet BCs on temperature
2 4
\endcode
<li>We give information on how to solve the matrix associated to the time marching of the temperature.
<ol>
<li>
\code
===Maximum number of iterations for temperature solver
100
\endcode
<li>
\code
===Relative tolerance for temperature solver
1.d-6
===Absolute tolerance for temperature solver
1.d-10
\endcode
<li>
\code
===Solver type for temperature (FGMRES, CG, ...)
GMRES
===Preconditionner type for temperature solver (HYPRE, JACOBI, MUMPS...)
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
5000d0  
\endcode
<li>We set stabilization coefficient for the divergence of the magnetic field and the penalty of the Dirichlet and interface terms.
\code
===Stabilization coefficient (divergence)
1.d0
===Stabilization coefficient for Dirichlet H and/or interface H/H
100.d0
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
</ol>

<h3> Outputs and value of reference </h3>

The outputs of this test are computed with the file <tt>post_processing_debug.f90</tt> 
that can be found in the following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. 

To check the well behavior of the code, we compute four quantities:
<ol>
<li>The H1 norm of the velocity.
<li>The L2 norm of the magnetic field.
<li>The L2 norm of the pressure.
<li>The L2 norm of the temperature
</ol>
These quantities are computed at the final time \f$t=0.4\f$.
 They are compared to reference values to attest of the correctness of the code.  
 These values of reference are in the last lines of the file <tt>debug_data_test_11</tt> in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
============================================
Mesh_BENCH1_20.FEM, P2
===Reference results
0.14529939453854082      H1 norm of velocity    
0.16031055031353644      L2 norm of magnetic field  
1.47953318917485640E-002 L2 norm of pressure
1.1061039638796786       L2 norm of temperature
\endcode


To conclude this test, we show the profile of the approximated pressure, velocity magnitude,
 temperature and magnetic field magnitude at the final time.
 These figures are done in the plane \f$y=0\f$ which
 is the union of the half plane \f$\theta=0\f$ and \f$\theta=\pi\f$.
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_11_pre_tfin.png "Pressure in the plane plane y=0."
    </td>
    <td align="center">
    @image html  fig_test_11_vel_tfin.png  "Velocity magnitude in the plane plane y=0."
    </td>
</tr>
</table>
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_11_temp_tfin.png "Temperature in the plane plane y=0."
    </td>
    <td align="center">
     @image html  fig_test_11_mag_tfin.png "Magnetic field magnitude in the plane plane y=0."
    </td>
</tr>
</table>
 */
