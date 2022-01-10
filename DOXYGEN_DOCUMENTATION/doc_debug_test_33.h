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
 * @page doc_debug_test_33 Test 33: Ferrohydrodynamic

<h3>Introduction</h3>
In this example, we check the correctness of SFEMaNS for a ferrohydrodynamic problem. Maxwell, Navier-Stokes and temperature equations are solved. The domain presents a solid and a fluid subdomain. A force known as the Kelvin force is added to the source term of Navier-Stokes. It represents the action of the magnetic field on the magnetic fluid. Boussinesq force is also present. Maxwell equations are reduced to the magnetostatic equations due to the quasi-static approximation for the magnetic field. The fluid and solid regions do not have the same thermal properties. The volumetric heat capacity and the thermal conductivity are different. This case is usefull to study the cooling by convection and magnetoconvection of a solid by a ferrofluid.

The domain of computation is \f$\Omega= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in (0,1) \times [0,2\pi) \times (0,1)\} \f$. We note \f$\Gamma= \partial \Omega\f$. It is composed of a solid and a fluid subdomain:
\f{align*}{
\overline{\Omega} = \overline{\Omega_s} \cup \overline{\Omega_f}.\\
\f}
The subdomains are defined followingly: \f$\Omega_s = \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in (0,1/2) \times [0,2\pi) \times (0,1)\} \f$ and \f$\Omega_f= \{ (r,\theta,z) \in {R}^3 : (r,\theta,z) \in (1/2,1) \times [0,2\pi) \times (0,1)\} \f$. We note \f$\Gamma_f = \partial \Omega_f\f$.

We solve the temperature equations:
\f{align*}{
\begin{cases}
c\partial_t T+ c\tilde{\bu} \cdot \GRAD T - \DIV (\lambda \GRAD T) &= f_T, \\
T_{|\Gamma} &= T_\text{bdy} , \\
T_{|t=0} &= T_0,
\end{cases}
\f}
in the domain \f$\Omega\f$. The extended velocity is defined by \f$\tilde{\bu} = \bu\f$ in \f$\Omega_f\f$ and \f$0\f$  in \f$\Omega_s\f$. The volumetric heat capacity \f$c\f$ and the thermal conductivity \f$\lambda\f$ are piecewise constant functions of space: \f$c = c_f\f$ in \f$\Omega_f\f$ and \f$c_s\f$ in \f$\Omega_s\f$, with \f$c_f = 2 c_s\f$; \f$\lambda = \lambda_f\f$ in \f$\Omega_f\f$ and \f$\lambda_s\f$ in \f$\Omega_s\f$, with \f$\lambda_s = 10 \lambda_f\f$. We recall that these parameters are dimensionless. 

We solve the Navier-Stokes equations:
\f{align*}{
\begin{cases}
\partial_t\bu+\left(\ROT\bu\right)\CROSS\bu - \frac{1}{\Re}\LAP \bu +\GRAD p     
&= \alpha T \textbf{e}_z + \chi(T) \GRAD \left( \frac{H^2}{2} \right) + \bef,
\\ \DIV \bu &= 0, \\
\bu_{|\Gamma_f} &= \bu_{\text{bdy}} , \\
\bu_{|t=0} &= \bu_0, \\
p_{|t=0} &= p_0,
\end{cases}
\f}
in the domain \f$\Omega_f\f$. We denote by \f$\textbf{e}_z\f$ the unit vector in the vertical direction. The function of the temperature \f$\chi\f$ in the Kelvin force is defined by
\f{align*}{
\chi(T) = T^2
\f}

We solve the magnetostatic equations:
\f{align*}{
\begin{cases}
\ROT \bH &= \bJ,\\ 
\DIV (\mu \bH) &= 0, \\
\bH \CROSS \bn &= \bH_{\text{bdy}} \CROSS \bn, \text{ on } \Gamma,\\
\bH_{|t=0} &= \bH_0, \\
\end{cases}
\f}
in the domain \f$\Omega\f$. The permeability \f$\mu\f$ is constant in the whole domain. Solving magnetostatic equations is possible by using a very small electrical conductivity in the whole domain (see data file description). The other terms of the induction equation are thus neglectable in that case.

The data are the source terms \f$f_T\f$, \f$\bef\f$ and \f$\bJ\f$, the boundary data \f$T_\text{bdy}\f$, \f$\bu_{\text{bdy}}\f$ and \f$\bH_{\text{bdy}}\f$, the initial data \f$T_0\f$, \f$\bu_0\f$, \f$p_0\f$ and \f$\bH_0\f$. The parameters are the volumetric heat capacities \f$c_s\f$ and \f$c_f\f$, the thermal conductivities \f$\lambda_s\f$ and \f$\lambda_f\f$, the kinetic Reynolds number \f$\Re\f$, the thermal gravity number \f$\alpha\f$ and the permeability \f$\mu\f$. 

<h3>Manufactured solutions</h3>

We approximate the following analytical solutions:
@f{align*}{
T(r,\theta,z,t) & = r^2(r-r_0)^2\sin(z)(1+\cos(\theta))\cos(t),
\\ u_r(r,\theta,z,t) &= -(r-r_0)^2\cos(z)(1+\cos(\theta))\cos(t),
\\ u_{\theta}(r,\theta,z,t) &= (r-r_0)^2\cos(z)(1+\cos(\theta))\cos(t),
\\ u_z(r,\theta,z,t) &=  \frac{r-r_0}{r} \sin(z) ((3r-r_0)(1+\cos(\theta))+(r-r_0)\sin(\theta))\cos(t),
\\ p(r,\theta,z,t) &= 0,
\\ H_r(r,\theta,z,t) &= - r^3 \exp(z) (1 + \sin(\theta) \cos(t),
\\ H_\theta(r,\theta,z,t) &= r^3 \exp(z) (1 + \sin(\theta)) \cos(t),
\\ H_z(r,\theta,z,t) &= r^2 \exp(z) (4 - \cos(\theta) + 4 \sin(\theta)) \cos(t),
@f}
where \f$r_0 = 1/2\f$ is the limit between solid and fluid regions. The velocity is the curl of a vector field, it is thus divergence free. Same thing for the magnetic field.

The source terms \f$f_T, \bef\f$, \f$\bJ\f$ and the boundary data \f$T_\text{bdy}, \bu_{\text{bdy}}\f$, \f$\bH_{\text{bdy}}\f$ are computed accordingly.

<h3>Generation of the mesh</h3>

The finite element mesh used for this test is named <tt>SOLID_FLUID_10.FEM</tt> and 
 has a mesh size of \f$0.1\f$ for the P1 approximation. 
You can generate this mesh with the files in the following directory:
 ($SFEMaNS_MESH_GEN_DIR)/EXAMPLES/EXAMPLES_MANUFACTURED_SOLUTIONS/SOLID_FLUID_10.
The following image shows the mesh for P1 finite elements.
<table width="100%" align="center">
<tr>
    <td align="center">
    @image html  fig_SOLID_FLUID_10.png  "Finite element mesh (P1)."
    </td>
</tr>
</table>

<h3>Information on the file <tt>condlim.f90</tt></h3>

The initial conditions, boundary conditions and the forcing terms are set in the file <tt>condlim_test_33.f90</tt>.
Here is a description of the subroutines and functions of interest.
<ol>
<li>The subroutine <tt>init_velocity_pressure</tt> initializes the velocity field
 and the pressure at the times \f$-dt\f$ and \f$0\f$ with \f$dt\f$ being the time step.
 This is done by using the functions vv_exact and pp_exact as follows:
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
<li>The subroutine <tt>init_temperature</tt> initializes the temperature at the times \f$-dt\f$ and \f$0\f$ with \f$dt\f$ the time step. 
This is done by using the function temperature_exact as follows:
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
 conditions on the velocity field.
<ol>
<li>The limit between fluid and solid region is defined:
\code
    REAL(KIND=8)                                      :: r0 = 0.5d0
\endcode
<li>We construct the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>We define the velocity field depending of its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) and of its mode m as follows:
\code
    IF (TYPE==1) THEN
       IF ((m==0).OR.(m==1)) THEN
          vv = -(r-r0)**2*cos(t)*cos(z)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE==3) THEN
       IF ((m==0).OR.(m==1)) THEN
          vv = (r-r0)**2*cos(t)*cos(z)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE==5) THEN
       IF ((m==0).OR.(m==1)) THEN
          vv = (r-r0)*cos(t)*sin(z)/r * (3*r-r0)
       ELSE
          vv = 0.d0
      END IF
    ELSE IF (TYPE==6) THEN
       IF (m==1) THEN
          vv = (r-r0)*cos(t)*sin(z)/r * (r-r0)
       ELSE
          vv = 0.d0
       END IF
    ELSE
       vv = 0.d0
    END IF
\endcode
where \f$t\f$ is the time. It is important to specify the null types or modes to avoid nonsense results.
</ol>
<li>The function <tt>pp_exact</tt> contains the analytical pressure. It is used to initialize the pressure. The pressure is equal to 0 (for every type or mode).
\code
    vv = 0.d0
\endcode
<li>The function <tt>temperature_exact</tt> contains the analytical temperature. It is used to initialize the temperature and to impose Dirichlet boundary condition on the temperature.
<ol>
<li>The limit between fluid and solid region is defined:
\code
    REAL(KIND=8)                                      :: r0 = 0.5d0
\endcode
<li>We construct the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>We define the temperature depending on its TYPE (1 and 2 for cosine and sine) and on its mode as follows:
\code
    IF ((TYPE==1).AND.((m==0).OR.(m==1))) THEN
       vv = r**2*(r-r0)**2*sin(z)*cos(t)
    ELSE
       vv = 0.d0
    END IF
\endcode
</ol>
<li>The function <tt>Hexact</tt> contains the analytical magnetic field. It is used to initialize the magnetic field and enforce the boundary conditions.
<ol>
<li>We construct the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>We define the magnetic field depending on its TYPE (1 and 2 for the component radial cosine and sine, 3 and 4 for the component azimuthal cosine and sine, 5 and 6 for the component vertical cosine and sine) and on its mode as follows:
\code
IF (TYPE == 1) THEN
       IF (m == 0) THEN
          vv = - exp(z) * r**3 * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 2) THEN
       IF (m == 1) THEN
          vv = - exp(z) * r**3 * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 3) THEN
       IF (m == 0) THEN
          vv = exp(z) * r**3 * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 4) THEN
       IF (m == 1) THEN
          vv = exp(z) * r**3 * cos(t)
       ELSE 
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 5) THEN
       IF (m == 0) THEN
          vv = 4 * exp(z) * r**2 * cos(t)
       ELSE IF (m == 1) THEN
          vv = - exp(z) * r**2 * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE
       IF (m == 1) THEN
          vv = 4 * exp(z) * r**2 * cos(t)
       ELSE
          vv = 0.d0
       END IF
    END IF
\endcode
</ol>
<li>The function <tt>source_in_temperature</tt> computes the source term \f$f_T\f$ of the temperature equations.
<ol>
<li>Arrays c and lambda for the value of the volumetric heat capacity and thermal conductivity at each node must be declared.
\code
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv, r, z, c, lambda
\endcode
<li>The limit between fluid and solid region is defined:
\code
    REAL(KIND=8)                                      :: r0 = 0.5d0
\endcode
<li>We construct the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>The c array is filled based on the data file. The solid volumetric heat capacity is used for the region \f$r \le r_0\f$ and the fluid volumetric heat capacity is used in the fluid region \f$r > r_0\f$. 
\code
DO i=1,SIZE(rr,2)
       IF (r(i).LE.r0) THEN
          c(i) = inputs%vol_heat_capacity(1)
       ELSE
          c(i) = inputs%vol_heat_capacity(2)
       END IF
    END DO
\endcode
<li>The lambda array is filled based on the data file. The solid conductivity is used for the region \f$r \le r_0\f$ and the fluid conductivity is used in the fluid region \f$r > r_0\f$. The temperature_diffusivity input contains the thermal conductivity when a volumetric heat capacity is used.
\code
DO i=1,SIZE(rr,2)
       IF (r(i).LE.r0) THEN
          lambda(i) = inputs%temperature_diffusivity(1)
       ELSE
          lambda(i) = inputs%temperature_diffusivity(2)
       END IF
    END DO
\endcode
<li>The source term \f$ f_T = c\partial_t T+ c\tilde{\bu} \cdot \GRAD T - \DIV(\lambda\GRAD T) \f$ is defined in two parts. Firstly, we define the part \f$ c\partial_t T - \DIV (\lambda\GRAD T) \f$:
\code
    IF (TYPE==1) THEN
       IF (m==0) THEN
          vv = - (-(r**4 + 18*r*r0 - 2*r**3*r0 - 4*r0**2 + r**2*(-16.d0 + r0**2))*lambda*Cos(t) &
               + c*r**2*(r - r0)**2*Sin(t))*Sin(z)
       ELSE IF (m==1) THEN
          vv = -(-(r**4 + 16*r*r0 - 2*r**3*r0 - 3*r0**2 + r**2*(-15.d0 + r0**2))*lambda*Cos(t) &
               + c*r**2*(r - r0)**2*Sin(t))*Sin(z)
       ELSE
          vv = 0.d0
       END IF
    ELSE
       vv = 0.d0
    END IF
\endcode
Secondly, we add the part \f$ c\tilde{\bu} \cdot \GRAD T \f$, which is different from 0 only in the fluid:
\code
    IF (TYPE==1) THEN
       IF (m==0) THEN
          DO i=1,size(rr,2)
             IF (r(i)>r0) THEN
                vv(i) = vv(i) + c(i)*(-3*r(i)*(r(i) - r0)**4*Cos(t)**2*Sin(2*z(i)))/4.d0
             END IF
          END DO
       ELSE IF (m==1) THEN
          DO i=1,size(rr,2)
             IF (r(i)>r0) THEN
                vv(i) = vv(i) - c(i)*r(i)*(r(i) - r0)**4*Cos(t)**2*Sin(2*z(i))
             END IF
          END DO
       ELSE IF (m==2) THEN
          DO i=1,size(rr,2)
             IF (r(i)>r0) THEN
                vv(i) = vv(i) - c(i)/lambda(i)*(r(i)*(r(i) - r0)**4*Cos(t)**2*Sin(2*z(i)))/4.d0
             END IF
          END DO
       END IF
    END IF
\endcode
In the second part, the array is filled cell by cell because we have to test if the node is in the fluid region.
</ol>
<li>The function <tt>source_in_NS_momentum</tt> computes the source term \f$\alpha T \textbf{e}_z+\bef\f$ of the Navier-Stokes equations.
<ol>
<li>The coefficient \f$\alpha\f$ of the Boussinesq force is declared and the limit between fluid and solid region is defined:
\code
    REAL(KIND=8)                                      :: alpha, r0 = 0.5d0
\endcode
<li>The coefficient \f$\alpha\f$ is defined based on the data file.
\code
    alpha = inputs%gravity_coefficient
\endcode
<li>We construct the radial and vertical coordinates r, z.
\code
    r = rr(1,:)
    z = rr(2,:)
\endcode
<li>We construct the first part of the source term containing the Boussinesq force \f$\alpha T \be_z\f$ and the piece of \f$\bef\f$ that cancels it.
\code
    IF (TYPE==5) THEN 
       vv = alpha*(opt_tempn(:,1,i) - temperature_exact(1,rr,mode,time))
    ELSE IF (TYPE==6) THEN
       vv = alpha*(opt_tempn(:,2,i) - temperature_exact(2,rr,mode,time))
    ELSE 
       vv = 0.d0
    END IF
\endcode
The array opt_temp is an input argument and produces the Boussinesq force. The temperature_exact function is defined in the same file and leads to the cancelling term of \f$\bef\f$.
<li>The part of the source term \f$ \partial_t\bu+\left(\ROT\bu\right)\CROSS\bu - \frac{1}{\Re}\LAP \bu +\GRAD p\f$ is then added. It depends on the TYPE (1-6) and the mode (0-2).
\code
IF (TYPE==1) THEN
       IF (mode==0) THEN
          vv = vv -(2*r*(r**4 - 2*r**3*r0 + r0**2 + r**2*(-3 + r0**2))*Cos(time)*Cos(z) + &
               (r - r0)*Re*Cos(time)**2*(14*r**3 - 5*r**2*r0 - 5*r*r0**2 + 2*r0**3 + &
               (9*r**5 - 21*r**4*r0 + 5*r*r0**2 - 2*r0**3 + r**3*(-14 + 15*r0**2) + r**2*(5*r0 - 3*r0**3))*Cos(2*z)) - &
               2*r**3*(r - r0)**2*Re*Cos(z)*Sin(time))/(2.*r**3*Re)
       ELSE IF (mode==1) THEN
          vv = vv +(-(r*(r**4 - 2*r*r0 - 2*r**3*r0 + 2*r0**2 + r**2*(-2 + r0**2))*Cos(time)*Cos(z)) - &
               (3*r**2 - 4*r*r0 + r0**2)*Re*Cos(time)**2*(3*r**2 - r0**2 + (2*r**4 - 4*r**3*r0 + r0**2 + &
               r**2*(-3 + 2*r0**2))*Cos(2*z)) + r**3*(r - r0)**2*Re*Cos(z)*Sin(time))/ (r**3*Re)
       ELSE IF (mode==2) THEN
          vv = vv -((r - r0)*Cos(time)**2*(4*r**2 - r*r0 - r0**2 + (3*r**4 - 7*r**3*r0 + r0**2 + &
               r**2*(-4 + 5*r0**2) + r*(r0 - r0**3))*Cos(2*z)))/(2.*r**2)
       END IF
    ELSE IF (TYPE==2) THEN
       IF (mode==1) THEN
          vv = vv + ((r - r0)**2*Cos(time)*(-2*r*Cos(z) + Re*Cos(time)*(r**4 - r*r0 - 2*r**3*r0 + &
               r0**2 + r**2*(-3 + r0**2) + (3*r**2 + r*r0 - r0**2)*Cos(2*z))))/(r**3*Re)
       ELSE IF (mode==2) THEN
          vv = vv + ((r - r0)**2*Cos(time)**2*(r**4 - r*r0 - 2*r**3*r0 + r0**2 + r**2*(-3 + r0**2) + &
               (3*r**2 + r*r0 - r0**2)*Cos(2*z)))/(2.*r**3)
       END IF
    ELSE IF (TYPE==3) THEN
       IF (mode==0) THEN
          vv = vv + (-3*r*(r - r0)**3*(3*r - r0)*Re*Cos(time)**2 + 2*(r**4 - 2*r**3*r0 + r0**2 + &
               r**2*(-3 + r0**2))*Cos(time)*Cos(z) - 2*r**2*(r - r0)**2*Re*Cos(z)*Sin(time))/(2.*r**2*Re)
       ELSE IF (mode==1) THEN
          vv = vv -(-2*r*(r**4 - 2*r*r0 - 2*r**3*r0 + 2*r0**2 + r**2*(-2 + r0**2))*Cos(time)*Cos(z) + &
               (r - r0)**3*(3*r - r0)*Re*Cos(time)**2*(1 + 4*r**2 - Cos(2*z)) + &
               2*r**3*(r - r0)**2*Re*Cos(z)*Sin(time))/(2.*r**3*Re)
       ELSE IF (mode==2) THEN
          vv = vv + ((r - r0)**3*(3*r - r0)*Cos(time)**2*(-1 - r**2 + Cos(2*z)))/(2.*r**3)
       END IF
    ELSE IF (TYPE==4) THEN
       IF (mode==1) THEN
          vv = vv + ((r - r0)**2*Cos(time)*(-4*r*Cos(z) + Re*Cos(time)*((-3*r + r0)**2 + &
               (2*r**4 + 6*r*r0 - 4*r**3*r0 - r0**2 + r**2*(-9 + 2*r0**2))*Cos(2*z))))/(2.*r**3*Re)
       ELSE IF (mode==2) THEN
          vv = vv + ((r - r0)**2*Cos(time)**2*(4*r - 2*r0 + (r*(-4 + (r - r0)**2) + 2*r0)*Cos(2*z)))/(2.*r**2)
       END IF
    ELSE IF (TYPE==5) THEN
       IF (mode==0) THEN
          vv = vv + (((3*r**4 - 4*r**3*r0 - r0**2 + r**2*(-3 + r0**2))*Cos(time) + &
               r*(r - r0)**2*(3*r**4 - r*r0 - 6*r**3*r0 + 2*r0**2 + r**2*(-4 + 3*r0**2))*Re*Cos(time)**2*Cos(z) - &
               r**2*(3*r**2 - 4*r*r0 + r0**2)*Re*Sin(time))*Sin(z))/(r**3*Re)
       ELSE IF (mode==1) THEN
          vv = vv + ((3*r**3 - 4*r0 - 4*r**2*r0 + r*r0**2)*Cos(time)*Sin(z) - r*(3*r**2 - 4*r*r0 + r0**2)*Re*Sin(time)*Sin(z) + &
               ((r - r0)**2*(4*r**4 - 2*r*r0 - 8*r**3*r0 + 3*r0**2 + r**2*(-5 + 4*r0**2))*Re*Cos(time)**2*Sin(2*z))/2.)/(r**2*Re)
       ELSE IF (mode==2) THEN
          vv = vv + ((r - r0)**2*(r**4 - r*r0 - 2*r**3*r0 + r0**2 + r**2*(-1 + r0**2))*Cos(time)**2*Sin(2*z))/(2.*r**2)
       END IF
    ELSE IF (TYPE==6) THEN
       IF (mode==1) THEN
          vv = vv -((((-r**3 + 2*r0 + 2*r**2*r0 - r*r0**2)*Cos(time) + 4*r*(r - r0)**3*Re*Cos(time)**2*Cos(z) + &
               r*(r - r0)**2*Re*Sin(time))*Sin(z))/(r**2*Re))
       ELSE IF (mode==2) THEN
          vv = vv -(((r - r0)**3*Cos(time)**2*Sin(2*z))/r)
       END IF
    END IF
\endcode
<li>The last part of the source term  \f$- \chi(T) \GRAD \left( \frac{H^2}{2} \right)\f$ is then added. It depends on the TYPE (1-6) and the mode (0-4).
\code
IF (TYPE == 1) THEN
       IF (mode == 0) THEN
          vv = vv - 0.25 * exp(2*z) * r**7 * (215.d0 + 51*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 1) THEN
          vv = vv - 5 * exp(2*z) * r**7 * (11.d0 + 3*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 2) THEN
          vv = vv + 14 * exp(2*z) * r**7 * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 3) THEN
          vv = vv + exp(2*z) * r**7 * (19.d0 + 3*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 4) THEN
          vv = vv + 0.75 * exp(2*z) * r**7 * (5.d0 + r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 2) THEN
       IF (mode == 1) THEN
          vv = vv - 3 * exp(2*z) * r**7 * (24.d0 + 5*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 2) THEN
          vv = vv - 4 * exp(2*z) * r**7 * (13.d0 + 3*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 3) THEN
          vv = vv - exp(2*z) * r**7 * (8.d0 + 3*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 4) THEN
          vv = vv + 2 * exp(2*z) * r**7 * (r - r0)**4 * cos(time)**4 * sin(z)**2 
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 3) THEN
       IF (mode == 0) THEN
          vv = vv - exp(2*z) * r**7 * (15.d0 + 2*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 1) THEN
          vv = vv - 0.5 * exp(2*z) * r**7 * (48.d0 + 7*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 2) THEN
          vv = vv - 2 * exp(2*z) * r**7 * (5.d0 + r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 3) THEN
          vv = vv - 0.5 * exp(2*z) * r**9 * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 4) THEN
          vv = vv + exp(2*z) * r**7 * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 4) THEN
       IF (mode == 1) THEN
          vv = vv - 0.5 * exp(2*z) * r**7 * (25.d0 + 2*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 2) THEN
          vv = vv - 0.25 * exp(2*z) * r**7 * (61.d0 + 6*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2  
       ELSE IF (mode == 3) THEN
          vv = vv - 0.5 * exp(2*z) * r**7 * (17.d0 + 2*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2  
       ELSE IF (mode == 4) THEN
          vv = vv - 0.125 * exp(2*z) * r**7 * (15.d0 + 2*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 5) THEN
       IF (mode == 0) THEN
          vv = vv - 0.125 * exp(2*z) * r**8 * (215.d0 + 34*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 1) THEN
          vv = vv - 2.5 * exp(2*z) * r**8 * (11.d0 + 2*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 2) THEN
          vv = vv + 7 * exp(2*z) * r**8 * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 3) THEN
          vv = vv + 0.5 * exp(2*z) * r**8 * (19.d0 + 2*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 4) THEN
          vv = vv + 0.125 * exp(2*z) * r**8 * (15.d0 + 2*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 6) THEN
       IF (mode == 1) THEN
          vv = vv - exp(2*z) * r**8 * (36.d0 + 5*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 2) THEN
          vv = vv - 2 * exp(2*z) * r**8 * (13.d0 + 2*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 3) THEN
          vv = vv - exp(2*z) * r**8 * (4.d0 + r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 4) THEN
          vv = vv + exp(2*z) * r**8 * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE
          vv = 0.d0
       END IF
    END IF
\endcode
</ol>
<li>The function <tt>Jexact_gauss</tt> computes the source term \f$\bJ\f$ of the magnetostatic equations.
<ol>
<li>We construct the radial and vertical coordinates r, z.
\code
    r = rr(1)
    z = rr(2)
\endcode
<li>The current density is defined by type and mode.
\code
IF (TYPE == 1) THEN
       IF (m == 0) THEN
          vv = - exp(z) * r**3 * cos(t)
       ELSE IF (m == 1) THEN
          vv = 4 * exp(z) * r * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 2) THEN
       IF (m == 1) THEN
          vv = - exp(z) * r * (-1.d0 + r**2) * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 3) THEN
       IF (m == 0) THEN
          vv = - exp(z) * r * (8.d0 + r**2) * cos(t)
       ELSE IF (m == 1) THEN
          vv = 2 * exp(z) * r * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 4) THEN
       IF (m == 1) THEN
          vv = - exp(z) * r * (8.d0 + r**2) * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 5) THEN
       IF (m == 0) THEN
          vv = 4 * exp(z) *r**2 * cos(t)
       ELSE IF (m == 1) THEN
          vv = exp(z) * r**2 * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE 
       IF (m == 1) THEN
          vv = 4 * exp(z) * r**2 * cos(t)
       ELSE
          vv = 0.d0
       END IF
    END IF
\endcode
</ol>
<li>The function <tt>chi_coeff_law</tt> computes the coefficient \f$\chi(T)\f$ of the Kelvin force. It takes a real \f$T\f$ and gives the real \f$\chi(T)\f$. In this case \f$\chi(T) = T^2\f$ so we implement:
\code
vv = temp**2
\endcode
</ol>
All the other subroutines present in the file <tt>condlim_test_33.f90</tt> are not used in this test.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file.


<h3>Setting in the data file</h3>

We describe the data file of this test.
 It is called <tt>debug_data_test_33</tt> and can be found in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<ol>
<li>We use a formatted mesh by setting:
\code
===Is mesh file formatted (true/false)?
.t.
\endcode
<li>The path and the name of the mesh are specified with the two following lines:
\code
===Directory and name of mesh file
'.' 'SOLID_FLUID_10.FEM'
\endcode
where '.' refers to the directory where the data file is, meaning ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
<li>We use one processor in the meridian section. It means the finite element mesh is not subdivised.
\code
===Number of processors in meridian section
1
\endcode
<li>We solved the problem for \f$5\f$ Fourier modes.
\code
===Number of Fourier modes
5
\endcode
The Fourier modes are not detailed so the first three modes \f$0,1,2,3,4\f$ are solved.
<li>We use \f$5\f$ processors in Fourier space.
\code
===Number of processors in Fourier space
5
\endcode
It means that each processors is solving the problem for \f$5/5=1\f$ Fourier modes.
<li>We approximate the ferrohydrodynamic equations by setting:
\code
===Problem type: (nst, mxw, mhd, fhd)
'fhd'
\endcode
<li>We do not restart the computations from previous results.
\code
===Restart on velocity (true/false)
.f.
===Restart on magnetic field (true/false)
.f.
\endcode
It means the computation starts from the time \f$t=0\f$.
<li>We use a time step of \f$0.01\f$ and solve the problem over \f$100\f$ time iterations.
\code
===Time step and number of time iterations
1.d-2 100
\endcode
<li>We set the number of domains and their label, see the files associated to the generation of the mesh,
 where the code approximates the Navier-Stokes equations.
\code
===Number of subdomains in Navier-Stokes mesh
1
===List of subdomains for Navier-Stokes mesh
2
\endcode
<li>We set the number of boundaries with Dirichlet conditions on the velocity field and give their respective labels.
\code
===How many boundary pieces for full Dirichlet BCs on velocity?
4
===List of boundary pieces for full Dirichlet BCs on velocity
2 3 4 5
\endcode
<li>We set the kinetic Reynolds number \f$\Re\f$.
\code
===Reynolds number
1.d0
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
<li>We solve the temperature equation.
\code
===Is there a temperature field?
.t.
\endcode
<li>We set the coefficient \f$\alpha\f$ of Boussinesq force.
\code
===Non-dimensional gravity coefficient
1.d0
\endcode
<li>We set the number of domains and their label, see the files associated to the generation of the mesh,
 where the code approximates the temperature equation.
\code
===Number of subdomains in temperature mesh
2
===List of subdomains for temperature mesh
1 2
\endcode
<li>We set the volumetric heat capacity \f$c\f$ and the thermal conductivity \f$\lambda\f$.
\code
===Volumetric heat capacity (1:nb_dom_temp)
1.d0 2.d0
===Thermal conductivity (1:nb_dom_temp)
10.d0 1.d0
\endcode
<li>We set the number of boundaries with Dirichlet conditions on the velocity and give their respective labels.
\code
===How many boundary pieces for Dirichlet BCs on temperature?
3
===List of boundary pieces for Dirichlet BCs on temperature
2 4 5
\endcode
<li>Setting the interfaces between regions where only the temperature is solved and regions where velocity and temperature are solved is not required for 'fhd' problems.
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
2
===List of subdomains for magnetic field (H) mesh
1 2
\endcode
<li>We set the number of interfaces in H_mesh and give their respective labels.
\code
===Number of interfaces in H mesh
1
===List of interfaces in H mesh
3
\endcode
There is no permeability jump but it is necessary to indicate that there is an interface in H_mesh in order to impose boundary conditions on the velocity on this interface.
<li>We set the number of boundaries with Dirichlet conditions on the magnetic field and give their respective labels.
\code
===Number of Dirichlet sides for Hxn
3
===List of Dirichlet sides for Hxn
2 4 5
\endcode
<li>We set the magnetic permeability in each domains where the magnetic field is approximated.
\code
===Permeability in the conductive part (1:nb_dom_H)
1.d0 1.d0
\endcode
<li>We set the conductivity in each domains where the magnetic field is approximated.
\code
===Conductivity in the conductive part (1:nb_dom_H)
1.d-20 1.d-20
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
<li>We set the stabilization coefficients for the divergence of the magnetic field and the penalty of the Dirichlet and interface terms.
\code
===Stabilization coefficient (divergence)
1.d0
===Stabilization coefficient for Dirichlet H and/or interface H/H
1.d0
\endcode
We note that these coefficients are usually set to \f$1\f$.
<li>We set the number of domains and their label, see the files associated to the generation of the mesh,
 where the code approximates the potential scalar: none, only an H formulation is used.
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
that can be found in: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. 

To check the well behavior of the code, we compute four quantities:
<ol>
<li>The L2-norm of error on u divided by the L2-norm of u exact.
<li>The L2-norm of error on p.
<li>The L2-norm of error on T divided by L2-norm of T exact.
<li>The L2-norm of error on H divided by L2-norm of H exact.
</ol>
These quantities are computed at the final time \f$t=1\f$.
 They are compared to reference values to attest of the correct behavior of the code.  
 These values of reference are in the last lines of the file <tt>debug_data_test_33</tt> in ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
============================================
(SOLID_FLUID_10.FEM)
===Reference results
1.410510605417536E-005  L2-norm of error on u / L2-norm of u exact
3.825370980418075E-004  L2-norm of error on p
7.607510571699082E-005  L2-norm of error on T / L2-norm of T exact
8.605845065931233E-005  L2-norm of error on H / L2-norm of H exact
\endcode

To conclude this test, we display the profile of the approximated pressure, velocity magnitude, temperature and magnetic field magnitude at the final time.
 These figures are done in the plane \f$y=0\f$ which
 is the union of the half plane \f$\theta=0\f$ and \f$\theta=\pi\f$.
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_33_pre_tfin.png "Pressure in the plane y=0."
    </td>
</tr>
</table>
<table width="100%" align="center">
<tr>
    <td align="center">
    @image html  fig_test_33_vel_tfin.png  "Velocity magnitude in the plane y=0."
    </td>
</tr>
</table>
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_33_temp_tfin.png "Temperature in the plane y=0."
    </td>
</tr>
</table>
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_33_mag_tfin.png "Magnetic field magnitude in the plane y=0."
    </td>
</tr>
</table>
 */
