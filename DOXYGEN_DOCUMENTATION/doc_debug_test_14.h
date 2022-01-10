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
 * @page doc_debug_test_14 Test 14: ARPACK for Maxwell and eigen values problem

<h3>Introduction</h3>
In this example, we check the correctness of SFEMaNS for an eigenvalue problem of a magnetic set up. The set up involves a conducting domain only. We consider Dirichlet boundary conditions. We use P2 finite elements for the magnetic field.

We approximate the first five eigenvalues (with the largest real part) of the Maxwell equations:
\f{align}{
\begin{cases}
\partial_t (\mu^c \mathbf{H}) + \nabla \times \left(\frac{1}{\Rm \sigma}\nabla \times \mathbf{H}  \right)  =
 \nabla\times (\bu^\text{given} \times \mu^c \mathbf{H})  & \text{in } \Omega_1, \\
\text{div} (\mu^c \mathbf {H}) = 0  & \text{in } \Omega_1 ,\\
\bH \times \bn = \bH_{\text{bdy}} \times \bn & \text{on } \Gamma_1,\\
\bH_{|t=0}= \bH_0,
\end{cases}
\f}
in the conducting domain \f$\Omega_1= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in [0,1] \times [0,2\pi) \times [-1,1]\} \f$ with \f$\Gamma_1= \partial \Omega_1\f$.
The data are the given velocity \f$\bu^\text{given}\f$,
 the boundary data \f$\bH_{\text{bdy}}\f$,
the initial datas \f$\bH_0\f$. The parameter \f$\Rm\f$ is the magnetic Reynolds number.
 The parameter \f$\mu^c\f$ is the magnetic permeability of the conducting region \f$\Omega_1\f$.
 The parameter \f$\sigma\f$ is the conductivity in the conducting region \f$\Omega_1\f$.

<h3>Manufactured solutions</h3>
We consider the following solutions:
@f{align*}
H_r(r,\theta,z,t) &= A \frac{\sin(n_0z)}{r}(R_i-r)^2 (R_o-r)^2\cos(\theta),
\\ H_{\theta}(r,\theta,z,t) &= -2 A \sin(n_0 z)(r-R_i)(r-R_o)(2r-R_o-R_i) \sin(\theta),
\\ H_z(r,\theta,z,t) &= 0,
@f}
where \f$n_0=1\f$, \f$R_i=1\f$, \f$R_o=2\f$ and \f$A=10^{-3}\f$.

We also set the given velocity field as follows:
@f{align*}
u^\text{given}_r(r,\theta,z,t) &= - \text{amp}_{MND} \frac{\pi}{h} r (1-r^2)(1+2r) \cos(\frac{2\pi z}{h}), 
\\ u^\text{given}_{\theta}(r,\theta,z,t) &= \text{amp}_{MND}\frac{4\pi}{h} r (1-r) \sin^{-1}(z)   ,
\\ u^\text{given}_z(r,\theta,z,t) &= \text{amp}_{MND} (1-r)(1+r+0.5 r^2)\sin(\frac{2\pi z}{h}),
@f}
where \f$\text{amp}_{MND}=1\f$ and \f$h=2\f$.



<h3>Generation of the mesh</h3>
The finite element mesh used for this test is named <tt>VKS_MND_form_10.FEM</tt> and 
 has a mesh size of \f$0.1\f$ for the P1 approximation. You can generate this mesh with the files in the following directory:
 ($SFEMaNS_MESH_GEN_DIR)/EXAMPLES/EXAMPLES_MANUFACTURED_SOLUTIONS/VKS_MND_form_10.
The following images show the mesh for P1 finite elements.

<table width="70%" align="center">
<tr>
    <td align="center">
    @image html  fig_mesh_VKS_MND_form_10.png  "Finite element mesh (P1) of the conducting region."
    </td>
</tr>
</table>





<h3>Information on the file <tt>condlim.f90</tt></h3>
The initial conditions, boundary conditions and the forcing term \f$\textbf{j}\f$ in the Maxwell
 equations are set in the file <tt>condlim_test_14.f90</tt>.
Here is a description of the subroutines and functions of interest.
<ol>
<li>We define the real number \f$\pi\f$.
\code
  REAL(KIND=8) :: pi=ACOS(-1.d0)
\endcode
<li>The subroutine <tt>init_maxwell</tt> initializes the magnetic field at the time\f$-dt\f$ and \f$0\f$ with \f$dt\f$ being the time step. It is done  as follows:
<ol>
<li>We set the magnetic fields to zero.
\code
    time  = 0.d0
    Hn    = 0.d0
    Hn1   = 0.d0
    phin  = 0.d0
    phin1 = 0.d0
\endcode
We note the scalar potential, not used in this test, is also initialized. 
<li>We define the parameters \f$n_0, R_i, R_0\f$ and \f$A\f$.
\code
    n0 = 1
    Ri = 1.d0
    Ro = 2.d0
    A = 1d-3
\endcode
<li>For the Fourier mode \f$m=1\f$, we define the magnetic fields depending of their TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine).
\code
    IF (H_mesh%me /= 0) THEN
       DO i = 1, SIZE(list_mode)
          IF (list_mode(i) == 1) THEN
             Hn1(:,1,i) = A*(SIN(n0*H_mesh%rr(2,:))/H_mesh%rr(1,:))*(Ri-H_mesh%rr(1,:))**2 &
                  * (Ro-H_mesh%rr(1,:))**2
             Hn(:,1,i)  = A*(SIN(n0*H_mesh%rr(2,:))/H_mesh%rr(1,:))*(Ri-H_mesh%rr(1,:))**2 &
                  * (Ro-H_mesh%rr(1,:))**2
             Hn1(:,4,i) = -2*A*SIN(n0*H_mesh%rr(2,:))*(H_mesh%rr(1,:)-Ri)&
                  *(H_mesh%rr(1,:)-Ro)*(2.d0*H_mesh%rr(1,:)-Ro-Ri)
             Hn(:,4,i)  = -2*A*SIN(n0*H_mesh%rr(2,:))*(H_mesh%rr(1,:)-Ri)&
                  *(H_mesh%rr(1,:)-Ro)*(2.d0*H_mesh%rr(1,:)-Ro-Ri)
          END IF
       END DO
    END IF
    RETURN
\endcode
We note that Hn and Hn1 are equal (the problem is time independent).
</ol>
<li>The function <tt>Hexact</tt> is used to impose Dirichlet boundary
 conditions on \f$\textbf{H}\times\textbf{n}\f$ with \f$\textbf{n}\f$ the outter normal vector. It is set to zero as follows:
\code
    vv = 0.d0
    RETURN
\endcode
<li>The function <tt>Vexact</tt> is used to define the given velocity field \f$\bu^\text{given}\f$.
<ol>
<li>We set the parameters \f$\epsilon, h, \zeta, \text{amp}_{MND}\f$ and \f$\text{amp}_{fl}\f$ used in the definition of the velocity field.
\code
    REAL(KIND=8) :: eps=1.d-5, height=2.d0, zeta=30.d0, amp_MND=1.d0, amp_fl=0.d0
\endcode
<li>We set the velocity field to zero.
\code    
    vv = 0.d0
\endcode
<li>For the Fourier mode \f$m=0\f$, we define the velocity field depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine).
\code
    IF (m==0) THEN
       virgin = .TRUE.
       DO mjl = 1, H_mesh%me
          !IF (H_mesh%i_d(mjl)/= 4) CYCLE
          !We are in the sodium
          DO njl = 1, H_mesh%gauss%n_w
             i =  H_mesh%jj(njl,mjl)
             IF (.NOT.virgin(i)) CYCLE
             virgin(i) = .FALSE.

             rr = H_mesh%rr(1,i)
             zz = H_mesh%rr(2,i)

             vv(i,1) = amp_MND*(-(PI/height)*rr*(1.d0-rr)**2*(1.d0+2.d0*rr)*COS(2.d0*PI*zz/height))
             vv(i,3) = amp_MND*(4.d0*height/PI)*rr*(1.d0-rr)*ASIN(zz)
             IF (zz .LE. (-eps)) THEN ! On est en bas
                vv(i,3) = vv(i,3)+ amp_fl*rr*SIN(PI*rr)*(1.-2.*zeta*(zz+1.)**2)*EXP(-zeta*(zz+1.)**2)
             ELSE IF (zz .GE. eps) THEN ! On est en haut
                vv(i,3) = vv(i,3)+ amp_fl*rr*SIN(PI*rr)*(1.-2.*zeta*(zz-1.)**2)*EXP(-zeta*(zz-1.)**2)
             ENDIF
             vv(i,5) = amp_MND*(1.-rr)*(1.+rr-5.*rr**2)*SIN(2.*PI*zz/height)
!!$              vel_loc(i) = sqrt(vv(i,1)**2 + vv(i,3)**2 + vv(i,5)**2)
          END DO !njl
       END DO
    END IF
\endcode
We note that the part multiphy by amp_fl is equal to zero. As a result, the parameters \f$\epsilon\f$, \f$\zeta\f$ and \f$\text{amp}_{nl}\f$ could be removed from this function.
</ol>
<li>The function <tt>Jexact_gauss</tt> is used to define the source term \f$\textbf{j}\f$ that is not used in this test. It is set to zero.
\code
   vv = 0.d0
   RETURN
\endcode
</ol>
All the other subroutines present in the file <tt>condlim_test_14.f90</tt> are not used in this test.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file.







<h3>Setting in the data file</h3>
We describe the data file of this test.
 It is called <tt>debug_data_test_14</tt> and can be found in the following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<ol>
<li>We use a formatted mesh by setting:
\code
===Is mesh file formatted (true/false)?
.t.
\endcode
<li>The path and the name of the mesh are specified with the two following lines:
\code
===Directory and name of mesh file
'.' 'VKS_MND_form_10.FEM'
\endcode
where '.' refers to the directory where the data file is, meaning ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<li>We use one processor in the meridian section. It means the finite element mesh is not subdivised.
\code
===Number of processors in meridian section
1
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
<li>We do not select specific Fourier modes to solve.
\code
===Select Fourier modes? (true/false)
.f.
\endcode
<li>We give the list of the Fourier modes to solve.
\code
===Problem type: (nst, mxw, mhd, fhd)
'mxw'
===Restart on velocity (true/false)
.f.
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
<li>We use a time step of \f$0.01\f$ and solve the problem over one time iteration.
\code
===Time step and number of time iterations
1.d-2, 1
\endcode
<li>To solve an eigenvalue problem with arpack we need to set the following parameters.
<ol>
<li>We use arpack.
\code
===Do we use Arpack?
.t.
\endcode
<li>We set the number of eigenvalues to compute.
\code
===Number of eigenvalues to compute
5
\endcode
<li>We set the maximum number of iterations Arpack can do.
\code
===Maximum number of Arpack iterations
3000
\endcode
<li>We set the tolerance.
\code
===Tolerance for Arpack
1.d-3
\endcode
<li>We give information on which eigenvalues we are looking for.
\code
===Which eigenvalues ('LM', 'SM', 'SR', 'LR' 'LI', 'SI')
'LR'
\endcode
We note that "LM", "SM", "SR", "LR", "LI" and "SI" respectively means largest magnitude, smallest magnitude, smallest real part, largest real part, largest imaginary part and smallest imaginary part.
<li>We generate vtu files with the following option.
\code
===Create 2D vtu files for Arpack? (true/false)
.t.
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
3
===List of Dirichlet sides for Hxn
2 4 3
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
80.d0
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
<li>The real part of the first eigenvalue of the mode 0.
<li>The relative divergence of the first eigenvalue of the mode 0.
<li>The real part of the first eigenvalue of the mode 1.
<li>The relative divergence of the first eigenvalue of the mode 1.
</ol>

 They are compared to reference values to attest of the correctness of the code.  
 These values of reference are in the last lines of the file <tt>debug_data_test_14</tt> in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
============================================
Quick test with VKS_MND_form_10.FEM, dt=1.d-2, tol=1.d-3, 5 eigenvalues
===Reference results
-4.97410008173271648E-002  Real part eigenvalue 1, mode 0
7.64643293777484551E-003   |div(Bn)|_L2/|Bn|_H1, eigenvalue 1, mode 0
4.71269276486770920E-003   Real part eigenvalue 1, mode 1
1.38537382005007436E-002   |div(Bn)|_L2/|Bn|_H1, eigenvalue 1, mode 1
\endcode

 */
