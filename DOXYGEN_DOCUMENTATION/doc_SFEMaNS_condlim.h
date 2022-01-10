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
 * @page doc_SFEMaNS_condlim Fortran file condlim.f90

In this section we describe the subroutines and the functions of the file <tt>condlim.f90</tt>. They allow to impose initial conditions, boundaries conditions and the source terms of the equations considered. We remind that the code SFEMaNS can approximate the Navier-Stokes equation, the temperature equation and the Maxwell equations for various kind of problems (variable density, solid obstacle in the fluid domain, variable magnetic permeability). A template of the file <tt>condlim.f90</tt> is available in the following directory: <tt>($SFEMaNS_DIR)/TEMPLATE</tt>.

<h1>Settings for the Navier-Stokes equations</h1>

<h2>General case</h2>
The code SFEMaNS can approximate the solutions of the Navier-Stokes equations in an axisymmetric domain \f$\Omega_{\bu}\f$. For a problem with constant density, constant viscosity and without solid obstacle, they are written as follows:
@f{align*}
\partial_t\bu+\left(\ROT\bu\right)\CROSS\bu - \frac{1}{\Re}\LAP \bu +\GRAD p     &=\bef,
\\ \DIV \bu &= 0, \\
\bu_{|\Gamma_1} &= \bu_{\text{bdy}} , \\
{\bu \cdot \textbf{n}}_{|\Gamma_2} &= 0 , \\
\bu_{|t=0} &= \bu_0, \\
p_{|t=0} &= p_0.
@f}
where \f$\textbf{n}\f$ is the outward normals of \f$\partial\Omega_{\bu}\f$. We also define \f$\Gamma_1\f$ and \f$\Gamma_2\f$ such that \f$\Gamma_1 \cup \Gamma_2=\partial \Omega_{\bu}\f$ and \f$ \Gamma_1 \cap \Gamma_2 = \emptyset \f$. The parameter \f$\Re\f$ is the kinetic Reynolds number.

Remark: for magnetohydrodynamic problem, the Lorentz force is added in the right hand-side of the Navier-Stokes equations. It is equal to \f$(\ROT\bH)\times \bB\f$.

To approximate such problems with the code SFEMaNS, the following subroutines/functions need to be edited. All required information to edit them are provided in the links below.
<table width="90%" align="center">
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_init_velocity_pressure  initializes the velocity field, the pressure and the pressure increment.
    </td>
</tr>
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_source_in_NS defines the source term \f$\textbf{f}\f$.
    </td>
</tr>
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_vv_exact defines the Dirichlet boundary condition \f$\bu_{\text{bdy}}\f$. It can also be used by the subroutine <tt>init_velocity_pressure</tt> to initialize the velocity field.
    </td>
</tr>
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_pp_exact can be used by the subroutine <tt>init_velocity_pressure</tt> to initialize the pressure. We note that Dirichlet boundary conditions on the pressure are not implemented.
    </td>
</tr>
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_extension_vel defines the velocity field in a solid domain where the Navier-Stokes equations are not approximated. It is used when the domains of the temperature or the magnetic field are larger than the domain where the Navier-Stokes equations are approximated.
    </td>
</tr>
</table>



<h2>Problem with variable density and/or dynamical viscosity</h2>
The code SFEMaNS can approximate the solutions of the Navier-Stokes equations with variable density and viscosity in an axisymmetric domain \f$\Omega_{\bu}\f$. We assume that the immiscible fluids are stratified. For each interface betweem two fluids, we introduce a level set function, denoted by \f$\varphi\f$. It is  equal to 0 in one fluid and 1 in the other. The interface is represented by \f$\varphi^{-1}(\{1/2\})\f$. The level set is the solution of the following equation in \f$\Omega_{\bu}\f$:
@f{align*}
\partial_t \varphi + \bu \cdot \GRAD \varphi = f_\varphi.
@f}
For a two phase field, we recontruct the density and dynamical viscosity as follows:
@f{align*}
\rho=\rho_1+(\rho_2-\rho_1) F(\varphi),
\\ \eta=\eta_1+(\eta_2-\rho_1) F(\varphi),
@f}
with \f$\eta_i\f$ and \f$\rho_i\f$ are the dynamical viscosities and densities of the two fluids. The function \f$F(\phi)\f$ is either the identity (called linear reconstruction) or a piece-wise polynomial function (called reg reconstruction).

We remind that the Navier-Stokes equations are approximated with
 the momentum \f$\bm=\rho \bu\f$ as dependent variable. We refer
 to the section \ref doc_intro_SFEMaNS_possibilities_nst_4
 for more information on the mathematical formulation of such problems.

In addition to the previous functions described above, the following functions of the file <tt>condlim.f90</tt> need to be edited. All required information to edit them are provided in the links below. 

<table width="90%" align="center">
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_init_level_set initializes the level set(s).
    </td>
</tr>
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_source_in_level_set defines the source term \f$f_\varphi \f$.
    </td>
</tr>
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_level_set_exact can be used by the subroutine <tt>init_level_set</tt> to initialize the level set(s) (the boundary conditions on the level set follow from the ones on the velocity).
    </td>
</tr>
</table>


<h2>Problem with non axisymmetric solid obstacle</h2>
The code SFEMaNS can consider solid non axisymmetric obstacle \f$\Omega_\text{obs}\f$
 in the domain \f$\Omega_{\bu}\f$ where the Navier-Stokes equations are approximated.
 The velocity in \f$\Omega_\text{obs}\f$ can be non zero and the domain
 \f$\Omega_\text{obs}\f$ can be time dependent. Such feature are made possible thanks
 to a penalty method. We refer to the section
 \ref doc_intro_SFEMaNS_possibilities_nst_3
 for more information on the mathematical formulation of such problems.

To consider such problems, the following functions needs to be set properly.

<table width="90%" align="center">
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_penal_in_real_space defines the penalty function \f$\chi\f$ that is equal to 1 in the fluid domain \f$\Omega_{\bu} \setminus \Omega_\text{obs} \f$ and 0 in the domain \f$\Omega_\text{obs}\f$. We note that the penalty function can be continuous (have it go from the values 0 to 1 in a few mesh cells around the interface between the solid and the fluid domain).
    </td>
</tr>
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_imposed_velocity_by_penalty defines the velocity in the domain \f$\Omega_\text{obs}\f$.
    </td>
</tr>
</table>

<h2>Problem with ferrofluid</h2>
The code SFEMaNS can consider ferrofluid. For such problems, the action of the magnetic field on the velocity field is reported by the Kelvin force. This force is equal to \f$g(T) \GRAD \bH^2 \f$ where \f$\bH\f$ is the magnetic field, \f$T\f$ is the temperature and \f$g(T)\f$ is a coefficient that depends of the temperature. We note that the definition of this coefficient can depend of the problem considered.

To consider such problems, the following function needs to be edited.
<table width="90%" align="center">
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_kelvin_force_coeff defines the coefficient \f$g(T)\f$ depending of the temperature.
    </td>
</tr>
</table>



<h1>Settings for the temperature equation</h1>

The code SFEMaNS can approximate the solution of the temperature equation in an axisymmetric domain \f$\Omega_{T}\f$. These equations are written as follows:
@f{align*}
\partial_t T+ \bu \cdot \GRAD T - \kappa \LAP T &= f_T, \\
T_{|\Gamma_1} &= T_\text{bdy} , \\
\partial_{\textbf{n}}T_{|\Gamma_2} &= 0, \\
T_{|t=0} &= T_0,
@f}
where \f$\textbf{n}\f$ is the outward normals of \f$\partial\Omega_T\f$. We also define \f$\Gamma_1\f$ and \f$\Gamma_2\f$ such that \f$\Gamma_1 \cup \Gamma_2=\partial \Omega_T\f$ and \f$ \Gamma_1 \cap \Gamma_2 = \emptyset \f$. The parameter \f$\kappa\f$ is the thermal diffusivity. We note that for thermohydrodynamic problems, we assume that \f$\Omega_{\bu} \subset \Omega_T\f$.

To approximate such problems with the code SFEMaNS, the following subroutines/functions need to be edited. All required information to edit them are provided in the links below.
<table width="90%" align="center">
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_init_temperature initializes the temperature.
    </td>
</tr>
<tr valign=top>
    <td align="left">
     @subpage doc_SFEMaNS_condlim_source_in_temperature defines the source term \f$f_T\f$.
    </td>
</tr>
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_temperature_exact defines the Dirichlet boundary condition \f$T_{\text{bdy}}\f$. It can be used by the subroutine <tt>init_temperature</tt> to initialize the temperature.
    </td>
</tr>
</table>



<h1>Settings for the Maxwell equations</h1>


<h2>General case</h2>
The code SFEMaNS can approximate the solutions of the  Maxwell equations in an axisymmetric domain \f$\Omega_\text{mxw}\f$. The domain is the union of an axisymmetric conducting domain \f$\Omega_{\bH}\f$ and an axisymmetric insulating domain \f$\Omega_\phi\f$. The Maxwell equations are written as follows:
\f{align}{
\begin{cases}
\partial_t (\mu^c \mathbf{H}) + \nabla \times \left(\frac{1}{\Rm \sigma}\nabla \times \mathbf{H}  \right)  =
 \nabla\times (\bu^\text{given} \times \mu^c \mathbf{H}) 
 + \nabla \times \left(\frac{1}{\Rm \sigma} \nabla\times \mathbf{j} \right) & \text{in } \Omega_{\text{H}}, \\
\text{div} (\mu^c \mathbf {H}) = 0  & \text{in } \Omega_{\text{H}} ,\\
 - \partial_t \DIV (\mu^v \GRAD( \phi)) = 0 & \text{in } \Omega_\phi ,\\ 
\bH \times  \bn^c + \nabla \phi \times \bn^v = 0 & \text{on } \Sigma, \\
\mu^c \bH \cdot  \bn^c + \mu ^v \nabla \phi \cdot \bn^v = 0 & \text{on } \Sigma, \\
\bH \times \bn = \bH_{\text{bdy}} \times \bn & \text{on } \Gamma_{\bH}^1,\\
 \left( \frac{1}{\Rm \sigma} \left( \ROT (\mathbf{H}) - \mathbf{j}  \right) - \bu \times \mu \mathbf{H}
 \right)  \times \bn   = \textbf{a} \times \bn &  \text{on } \Gamma_{\bH}^2,\\
\phi = \phi_{\text{bdy}}  & \text{on } \Gamma_\phi,\\
\bH_{|t=0}= \bH_0, \\
\phi_{|t=0}= \phi _0,
\end{cases}
\f}
where \f$\Sigma= \Omega_{\bH} \cap \Omega_\phi \f$, \f$\Gamma_{\bH}^1 \cup \Gamma_{\bH}^2= \partial \Omega_{\bH}\setminus \Sigma \f$, \f$\Gamma_{\bH}^1 \cap \Gamma_{\bH}^2=\emptyset\f$ and \f$\Gamma_\phi= \partial \Omega_\phi\setminus \Sigma \f$. We denote by \f$\textbf{n}\f$ the outward normal vector to \f$\Gamma_{\bH}^1 \cup \Gamma_{\bH}^2\f$. The function \f$\textbf{j}\f$ is a source term to define. The parameter \f$\Rm\f$ is the magnetic Reynolds number. We also denoted by \f$\sigma\f$ the conductivity in \f$\Omega_{\bH}\f$ and by \f$\mu^c\f$ the magnetic permeability in in \f$\Omega_{\bH}\f$. They both are piecewise constant. The magnetic permeability in the insulating domain is denoted \f$\mu^v\f$. The velocity \f$\f$ comes from the Navier-Stokes equations or functions of the file <tt>condlim.f90</tt>. We remind that the scalar potential satisfies the relation \f$\bH=\GRAD \phi\f$ in \f$\Omega_\phi\f$. We note that for magnetohydrodynamic problems, we assume that \f$\Omega_{\bu} \subset \Omega_T \subset \Omega_\text{mxw}\f$.

To approximate such problems with the code SFEMaNS, the following subroutines/functions need to be edited. All required information to edit them are provided in the links below.
<table width="90%" align="center">
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_init_maxwell initializes the magnetic field \f$\bH\f$ and the scalar potential \f$\phi\f$.
    </td>
</tr>
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_Hexact defines the Dirichlet boundary condition \f$\bH_{\text{bdy}}\f$. It can be used by the subroutine <tt>init_maxwell</tt> to initialize the magnetic field.
    </td>
</tr>
<tr valign=top>
    <td align="left">
  @subpage doc_SFEMaNS_condlim_Phiexact defines the Dirichlet boundary condition \f$\phi_{\text{bdy}}\f$. It can be used by the subroutine <tt>init_maxwell</tt> to initialize the scalar potential.
    </td>
</tr>
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_Vexact defines the velocity field \f$\bu_\text{given}\f$ that is a time independent function. It is only used when the solutions of the Navier-stokes equations are not approximated.
    </td>
</tr>
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_Jexact_gauss defines the source term \f$\textbf{j}\f$.
    </td>
</tr>
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_Eexact_gauss defines the boundary data \f$\textbf{a}\f$.
    </td>
</tr>
</table>


<h2>Quasi-static approximation</h2>
The code SFEMaNS allows to approximate the Maxwell equation in a quasi-static regime
 with the magnetic field \f$\bB=\mu\bH\f$ as dependent variable. The total magnetic
 field \f$\bB\f$ is written \f$\bB_b+\bb\f$ where \f$\bB_b\f$ represent a basic state
 that is time independent. The Maxwell equations are approximated with \f$\bb\f$
 where the time derivative and the term
 \f$\nabla\times (\bu^\text{given} \times \mathbf{b})\f$ are disregarded.
 The term \f$\nabla\times (\bu^\text{given} \times \mathbf{B}_b)\f$ is still
 computed and that the source term \f$\textbf{j}\f$ has to report the action of the term
 \f$\ROT( \frac{1}{\sigma\Rm}\ROT \bB_b)\f$. The equations can be written as follows.
\f{align*}{
 \ROT \left( \frac{1}{\sigma\Rm} \ROT \textbf{b} \right)
  = \nabla\times (\bu^\text{given} \times \mathbf{B}_b)
  + \ROT\left( \frac{1}{\sigma\Rm}\textbf{j}\right).
\f}

Remark: the Lorentz force in the Navier-Stokes equations is taken equal to \f$(\ROT\bH_b)\times \bb + (\ROT \bh)\times \bB_b\f$ with \f$\bh =\mu \bb\f$ and \f$\bH_b=\mu \bB_b\f$.

To approximate a quasi-static regime of the Maxwell equation you need to set properly the following function:
<table width="90%" align="center">
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_H_B_quasi_static defines the base states \f$\bB_b\f$ and \f$\bH_b\f$.
    </td>
</tr>
</table>
 All required information to edit them are provided in the above link.



<h2>Problem with a given magnetic permeability variable in \f$(r,\theta,z,t)\f$</h2>
The code SFEMaNS can take into account magnetic permeability that are variable in
 \f$(r,\theta,z,t)\f$ with \f$t\f$ being the time. It is done by approximating
 the Maxwell equations with the magnetic field \f$\bB\f$ as dependent variable.
 Moreover, the diffusion term \f$\ROT (\frac{1}{\sigma\Rm}\ROT(\frac{\bB}{\mu}))\f$
 is made explicit. The numerical scheme is stabilized with a term that involves a
 function \f$\overline{\mu}(r,z)\f$ that has to be smaller than \f$\mu(r,\theta,z,t)\f$
 for all \f$(r,\theta,z,t)\f$. 
 We refer to the section \ref doc_intro_SFEMaNS_possibilities_mxw_2
 for more information on the mathematical formulation of such problems. 

Remark:
<ol>
<li>If the magnetic permeability does not depend of the time and the azimuthal coordinate \f$\theta\f$, the term \f$\ROT( \ROT (\frac{1}{\mu}\bB))\f$ can be treated implicitly and no stabilization term is added.
<li>If the magnetic permeability 
is a piece wise constant function with re (one value on each axisymmetric domains where the 
</ol>

To approximate such problems with the code SFEMaNS, the following subroutines/functions need to be edited. All required information to edit them are provided in the links below. 
<table width="90%" align="center">
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_mu_bar_in_fourier_space. If the magnetic permeability does not depend of \f$(\theta,t)\f$, this function defines the magnetic permeability \f$\mu\f$. If \f$\mu\f$ depends of the azimuthal direction \f$\theta\f$ and/or the time \f$t\f$, this function defines a function \f$\overline{\mu}(r,z)\f$. It satisfies the relation \f$\overline{\mu}(r,z)\leq \mu(r,\theta,z,t)\f$ for all \f$(r,\theta,z,t)\f$ considered.
    </td>
</tr>
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_grad_mu_bar_in_fourier_space defines the gradient in the radial and vertical directions of the function mu_bar_in_fourier_space described above.
    </td>
</tr>
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_mu_in_real_space defines the magnetic permeability in \f$\Omega_\text{mxw}\f$ when its depends of the time and/or the azimuthal direction.
    </td>
</tr>
</table>

<h2>Multiphase problem with variable electrical conductivity</h2>
The code SFEMaNS can approximate the Navier-Stokes equations for stratified
 immiscible fluids. The electrical conductivity of these fluids may be different.
 As a consequence, the electrical conductivity is variable in \f$(r,\theta,z,t)\f$
 with t being the time. To approximate the Maxwell equations for such problems,
 the diffusion term \f$\ROT (\frac{1}{\sigma\Rm}\ROT(\frac{\bB}{\mu}))\f$ is
 made explicit.  The numerical scheme is stabilized with a term that involves
 a function \f$\overline{\sigma}(r,z)\f$ that has to be smaller than
 \f$\sigma(r,\theta,z,t)\f$ for all \f$(r,\theta,z,t)\f$. 
 We refer to the section \ref doc_intro_SFEMaNS_possibilities_nst_4
 for more information on the mathematical formulation of such problems.

Remark:
<ol>
<li>The magnetic permeability \f$\mu\f$ needs to be constant.
<li>The Maxwell equations can either be approximated with \f$\bB\f$ or \f$\bH\f$ as  dependent variable.
</ol>

To approximate such problems with the code SFEMaNS, the following function needs to be edited.
<table width="90%" align="center">
<tr valign=top>
    <td align="left">
    @subpage doc_SFEMaNS_condlim_sigma_bar_in_fourier_space. It defines the stabilization term \f$\overline{\sigma}\f$ such that \f$\overline{\sigma}(r,z)\leq \sigma\f$ for all \f$(r,\theta,z,t)\f$ considered.
    </td>
</tr>
</table>
 */
