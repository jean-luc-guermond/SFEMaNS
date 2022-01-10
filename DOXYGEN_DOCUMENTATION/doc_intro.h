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
 * @page doc_intro SFEMaNS presentation

This section starts with a presentation of the equations that
 are implemented in SFEMaNS.


@section doc_intro_equations Equations considered by SFEMaNS

The following equations are implemented in SFEMaNS.
<ol>
<li>The incompressible Navier-Stokes equations. In a domain \f$\Omega\f$, these
 equations are written as follows:
@f{align*}
\partial_t\bu+\left(\ROT\bu\right)\CROSS\bu 
- \frac{1}{\Re}\LAP \bu +\GRAD p     &=\bef,
\\ \DIV \bu &= 0,
@f}
with \f$\bu\f$ the velocity field, \f$p\f$ the pressure, \f$\Re\f$
 the kinetic Reynolds number and \f$\bef\f$ a source term.

Remark: One can also consider the abovee equations with variable density and 
viscosity as described in @subpage doc_intro_SFEMaNS_possibilities_nst_4

<li>The heat equation. In a domain \f$\Omega\f$, these equations
 are written as follows:
@f{align*}
 C \partial_t T+ \DIV(T \bu) - \DIV (\lambda \GRAD T) &= f_T,
@f}
with \f$T\f$ the temperature, \f$\bu\f$ the velocity field,
 \f$C\f$ the volumetric heat capacity, \f$\lambda\f$ the thermal
 conducitivty and \f$f_T\f$ a source term.

<li>The Maxwell equations. In a conducting domain \f$\Omega_c\f$,
 these equations are written as follows:
@f{align*}
\partial_t (\mu^c \bH^c) + \nabla \times \left(\frac{1}{\Rm \sigma}
\nabla \times \bH^c  \right)  =
 \nabla\times (\bu \times \mu^c \bH^c) 
+ \nabla \times \left(\frac{1}{\Rm \sigma}\mathbf{j}^s \right), \\
\text{div} (\mu^c \bH^c) = 0 ,
@f}
with \f$\bH^c\f$ the magnetic field, \f$\bu\f$ the velocity field,
 \f$\textbf{j}^s\f$ a source term, \f$\mu^c\f$ the magnetic permeability,
 \f$\sigma\f$ the electrical conductivity and \f$\Rm\f$ the magnetic
 Reynolds number. If the magnetic permeability is discontinuous across
 a surface denoted \f$\Sigma_\mu\f$, the following equations have to
 be satisfied on \f$\Sigma_\mu\f$:
@f{align*}
\bH^c_1 \times  \bn_1 +  \bH^c_2 \times  \bn_2 = 0,\\
\mu^c_1\bH^c_1 \cdot  \bn_1 +  \mu^c_2 \bH^c_2 \cdot  \bn_2 = 0 ,\\
@f}
where \f$\bn_1, \bn_2\f$ are outward normal to the surface
 \f$\Sigma_\mu\f$. \f$\bn_1\f$ points from $\f$\Omega_1\f$ to $\f$\Omega_2\f$ and
\f$\bn_1=-\bn_2\f$.  

<li> If the conductivity in a subdomain of the magnetic domain is zero,
we adopt a potential representation for the magnetic field. Let \f$\Omega_v\f$ be 
the insulating domain in question (assumed to be simply connected).  \f$\Omega_v\f$ is
henceforth referred to as vacuum. The Maxwell equations are then written as follows in  \f$\Omega_v\f$:
@f{align*}
-\mu^v \partial_t \LAP \phi = 0 ,
@f}
where \f$\phi\f$ the scalar potential such that \f$\bH=\GRAD \phi\f$ in the vacuum.
 The following continuity conditions across the interface
 \f$\Sigma=\Omega_c \cap \Omega_v\f$ have to be satisfied:
@f{align*}
\bH^c \times  \bn^c + \nabla \phi \times \bn^v = 0 , \\
\mu^c \bH^c \cdot  \bn^c + \mu ^v \nabla \phi \cdot \bn^v = 0 ,
@f}
with \f$\bn^c\f$ and \f$\bn^v\f$ the outward normals to
 the surface \f$\Sigma\f$.
</ol>

Remark: the above equations are supplemented by initial
 and boundaries conditions.

@section doc_intro_frame_work Key structural assumptions


@subsection doc_intro_frame_work_domain_axi Domain geometry and axisymmetric hypothesis
The code SFEMaNS uses cylindrical coordinates \f$(r,\theta,z)\f$ and assumes
a priori that the computational domain is axi-symmetric. The approximation method 
consists of using a Fourier
 decomposition in the azimuthal direction and Lagrange finite elements in the meridian
 section.

@subsection doc_intro_frame_work_domain_decomp Domain decomposition and simply connected insulating sub-domain hypothesis

The computational domain \f$\Omega\f$ is divided into the following three sub-domains:
<ol>
<li>A conducting fluid domain, denoted by \f$\Omega_{c,f}\f$, where the conductivity,
 permeability, viscosity and density of the fluid are constant and positive.
<li>A conducting solid domain, denoted by \f$\Omega_{c,s}\f$, where the velocity 
 of the solid is imposed. This sub-domain is assumed to be a finite union of disjoint 
 solid axisymmetric domains \f$\Omega_{c,s}^i\f$ with positive constant 
 conductivity \f$\sigma_i\f$ 
 and permeability \f$\mu_i\f$. We denote by \f$I\f$ the set that
 contains the integers \f$i\f$.
<li>An insulating domain, called vacuum and denoted by \f$\Omega_v\f$,
 where the electrical conductivity \f$\sigma\f$ is zero and the relative 
 magnetic permeability \f$\mu^v\f$ is 1 by default (but can be set to an other value by the user).
</ol>

The insulating sub-domain \f$\Omega_v\f$ is assumed to be simply
 connected so the magnetic field are written \f$\bH=\GRAD\phi\f$.
 The scalar potential \f$\phi\f$ can be proved to be the solution of
 the following equation in  \f$\Omega_v\f$:
@f{align*}
-\mu^v \partial_t \LAP \phi = 0.
@f}

Remarks:
<ol>
<li>The Navier-stokes equations are approximated in the fluid domain
 \f$\Omega_{c,f}\f$.
<li>The heat equation is approximated in a domain
 \f$\Omega_T=\Omega_{c,f} \cup (\underset{j\in J}{\cup}\Omega_{c,s}^j)\f$
 with \f$J\f$ a set included in \f$I\f$.
<li>The magnetic field \f$\bH^c\f$ is approximated in
 \f$\Omega_c=\Omega_{c,f} \cup (\underset{i\in I}{\cup}\Omega_{c,s}^i)\f$.
<li>The scalar potential \f$\phi\f$ is approximated in \f$\Omega_v\f$.
</ol>

@subsection doc_intro_num_app SFEMaNS's possibilities

The following set ups can be considered by the code SFEMaNS:
<ol>
<li>Hydrodynamics (NST). The Navier-Stokes equations are approximated.
 Thermal effect can also be considered.
<li>Magnetism (MXW). The Maxwell equations are approximated with a given velocity field.
<li>Magnetohydrodynamics (MHD). The Navier-Stokes and the Maxwell equations
 are approximated. The Lorentz force \f$\textbf{f}_\text{L}=
 (\ROT \bH) \times (\mu\bH)\f$ is added to the
 source term \f$\textbf{f}\f$ of the Navier-Stokes equations.
 Thermal effect can also be considered when solving the temperature equation.
<li>Ferrohydrodynamics (FHD). All of the above equations are approximated.
 The Kelvin force \f$\textbf{f}_\text{fhd}=g(T) \GRAD(\frac{\bH^2}{2})\f$,
 with \f$g(T)\f$ a user-defined scalar function, is added to the
 source term \f$\textbf{f}\f$ of the Navier-Stokes equations.
 The action of the magnetic field on the temperature can be taken into account.
 The change of temperature can also be taken into account in the Maxwell equations.
</ol>

The following extensions are also available in SFEMaNS but require a good knowledge of the code to be used properly:
<ol>
<li>Stabilization method, called entropy viscosity, for problems with large kinetic Reynolds numbers.
<li>Non axisymmetric geometry.
<li>Multiphase flow problem with variable density, dynamical viscosity and electrical conductivity.
<li>Magnetic permeability depending of the time and the azimuthal direction. We note the
 variation in \f$\theta\f$ is smooth while jumps in the \f$(r,z)\f$ direction can be considered.
<li>Quasi-static approximation of the MHD equations.
</ol>

The approximation methods of the above setting are described in the section
 \ref doc_num_app. The use is referred to this
<a href='doc_SFEMaNS_condlim.html'><code>section</code> </a>
for more details on the quasi-static approximation of the MHD equations.


*/

@section doc_intro_frame_work_num_app Numerical approximation


@subsection doc_intro_fram_work_num_app_Fourier_FEM Fourier/Finite element representation
The SFEMaNS code uses a hybrid Fourier/Finite element formulation.
 The Fourier decomposition allows to approximate the problem’s 
 solutions for each Fourier mode independently, modulo nonlinear
 terms that are made explicit. The variables are then approximated 
 on a meridian section of the domain with a finite element method.

The numerical approximation of a function \f$f\f$ is written in the following generic form:
@f{align*}
f(r,\theta,z,t)=f_h^{0,\cos}(r,z,t) + 
 \sum_{m=1}^M f_h^{m,\cos} \cos(m\theta) + f_h^{m,\sin} \sin(m\theta),

@f}
with \f$(r,\theta,z)\f$ the cylindrical coordinates, \f$t\f$ the time and \f$M\f$ 
 the number of Fourier modes considered. The unknown \f$f_h^{m,\cos}\f$ and
 \f$f_h^{m,\sin}\f$ can be approximated independtly modulo the computation 
 of nonlinear terms. Introducing the functions \f$\cos_m = \cos(m\theta)\f$,
 \f$\sin_m = \sin(m\theta)\f$ and a basis functions \f$(\phi_j)_{j \in J}\f$
 of the finite element space of the meridian section results in
 \f$(\phi_j \cos_m)_{j\in J, m \in [|0,M|]} \cup (\phi_j \sin_m)_{j\in J, m \in [|1,M|]}\f$
 being a basis of the space of approximation.

@subsection doc_intro_fram_work_num_app_space Space of approximation
 We set the number of Fourier mode to \f$M\f$. We define the meridian sections
 \f$\Omega_{c,f}^{2D}\f$, \f$\Omega_{T}^{2D}\f$,
 \f$\Omega_{c}^{2D}\f$, \f$\Omega_v^{2D}\f$ and \f$\Omega^{2D}\f$ of
 the domains \f$\Omega_{c,f}\f$, \f$\Omega_{T}\f$, \f$\Omega_{c}\f$,
 \f$\Omega_v\f$ and \f$\Omega\f$.
 We also consider \f$\left\{ \mathcal{T}_h \right\}_{h > 0}\f$
 a family of meshes of the meridian plane \f$\Omega^{2D}\f$
 composed of disjoint triangular cells \f$K\f$ of diameters at most \f$h\f$.
 For a given \f$h>0\f$, we assume that the mesh \f$\mathcal{T}_h\f$
 can be restricted to each sub domain of interests. These sub-meshes are denoted 
 \f$\mathcal{T}_h^{c,f}\f$, \f$\mathcal{T}_h^{T}\f$, \f$\mathcal{T}_h^{c}\f$
 and \f$\mathcal{T}_h^{v}\f$. The approximation of the solutions of the
 Navier-Stokes, heat and Maxwell equations either involve \f$\mathbb{P}_1\f$ or
 \f$\mathbb{P}_2\f$ Lagrange finite elements. The following defines
 the space of approximation used for each dependent variable.

<ol>
<li>The velocity field \f$\bu\f$ and the pressure \f$p\f$ are 
 respectively approximated in the following spaces: 
\f{align*}{
\textbf{V}_{h}^\bu := \left\{   \textbf{v} =\sum\limits_{k=-M}^M \textbf{v}_h^k (r,z) 
e^{ik \theta} ; \textbf{v}_h^k \in  \textbf{V}_{h}^{\bu,2D} \text{, }
 \overline{\textbf{v}_h^k}=\textbf{v}_h^{-k} \text{, } -M \leq k \leq M \right\} ,\\
S_{h}^p := \left\{    q_h= \sum\limits_{k=-M}^M q_h^k (r,z) e^{i k \theta} ;
 q_h^k \in S_{h}^{p,2D} \text{, } \overline{q_h^k}=q_h^{-k} \text{, }
 -M \leq k \leq M  \right\},
\f}
where we introduce the following finite element space: 
\f{align*}{
\textbf{V}_{h}^{\bu,2D}  : = \left\{ \textbf{v}_h \in C^0(\overline{\Omega_{c,f}^{2D}}) ;
 \textbf{v}_h|_K \in \mathbb{P}_2^6 \text{ } \forall K \in \mathcal{T}_h^{c,f}
     \right\} ,\\
S_{h}^{p,2D} : = \left\{ q_h \in   C^0(\overline{\Omega_{c,f}^{2D}}) ;
 q_h|_K \in \mathbb{P}_1^2 \text{ }  \forall K \in \mathcal{T}_h^{c,f}  \right\} .
\f}
We also introduce the subspace \f$\textbf{V}_{h,0}^{2D}\f$
 of \f$\textbf{V}_{h}^\bu\f$ composed of vector field that are
 zero on the boundaries of \f$\Omega_{c,f}\f$.
<li>The temperature \f$T\f$ is approximated in the following space:
\f{align*}{
S_h^T := \left\{    q_h= \sum\limits_{k=-M}^M q_h^k (r,z) e^{i k \theta} ;
 q_h^k \in S_{h}^{T,2D} \text{, } \overline{q_h^k}=q_h^{-k} \text{, }
 -M \leq k \leq M  \right\},
\f}
where we introduce the following finite element space: 
\f{align*}{
S_{h}^{T,2D} : = \left\{ q_h \in   C^0(\overline{\Omega_{T}^{2D}}) ;
 q_h|_K \in \mathbb{P}_2^2 \text{ }  \forall K \in \mathcal{T}_h^{T}  \right\} .
\f}
<li>The magnetic field \f$\bH^c\f$ and the scalar potentiel \f$\phi\f$
 are respectively approximated in the following spaces:
\f{align*}{
\textbf{V}_h^{\bH^c} :=
\left\{
\textbf{b}= 
\sum\limits_{k=-M}^M \textbf{b}_h^k (r,z) e^{ik \theta} ;
 \textbf{b}_h^k \in  \textbf{V}_h^{\bH^c,2D}, \; -M \leq k \leq M \right\},
\\
S_h^{\phi} :=
\left\{
\varphi= 
\sum\limits_{k=-M}^M \varphi_h^k (r,z) e^{ik \theta} ;
 \varphi_h^k \in  S_{h}^{\phi,2D}, \; -M \leq k \leq M
\right\},
\f}
where we introduce the following finite element spaces: 
\f{align*}{
\textbf{V}_{h}^{\bH^c, 2D} : = \left\{ \textbf{b}_h \in 
C^0(\overline{\Omega_{c}^{2D}});
 \textbf{b}_h|_K \in \mathbb{P}_{l_\bH}^6 \text{ }
 \forall K \in \mathcal{T}_h^{c}     \right\} ,\\
S_{h}^{\phi, 2D}  : = \left\{ \varphi_h \in   C^0(\overline{\Omega_{v}^{2D}}) ;
 \varphi_h|_K \in \mathbb{P}_{l_\phi}^2 \text{ }  \forall K \in \mathcal{T}_h^v ,  \right\} 
\f}
with \f$l_\bH\in\{1,2\}\f$ and \f$l_\phi\in\{1,2\}\f$ such that \f$l_\phi \geq l_\bH\f$.
</ol>

@subsection doc_intro_fram_work_num_app_time_marching Time marching

To present the time marching, we introduce a time step \f$\tau\f$ and denote
 by \f$f^n\f$ the approximation of a function at the time \f$t_n=n \tau\f$.
 When approximating the Navier-Stokes, heat and Maxwell equations, the time marching
 can be summarized by the four following steps:
<ol>
<li>Initialization of the temperature \f$T^0, T^1\f$, velocity fields \f$\bu^0, \bu^1\f$ , the dynamical pressure \f$p^0, p^1\f$ , the magnetic field \f$\bH^0, \bH^1\f$ and the scalar potential \f$\phi^0, \phi^1\f$.
<li>Approximation of \f$T^{n+1}\f$ after the nonlinear terms are computed with extroplation involving \f$T^n, T^{n-1}, \bu^n\f$ and \f$\bu^{n-1}\f$.
<li>Approximation of \f$\bu^{n+1}\f$ and \f$p^{n+1}\f$ after the nonlinear terms are computed with extroplation involving \f$\bu^n, \bu^{n-1}, \bH^n\f$ and \f$\bH^{n-1}\f$.
<li>Approximation of \f$\bH^{n+1}\f$ and \f$\phi^{n+1}\f$ after the nonlinear terms are computed with extroplation involving \f$\bu^{n+1}, \bH^n\f$ and \f$\bH^{n-1}\f$.
</ol>


@section doc_intro_SFEMaNS_weak_form_extensions Weak formulation and extensions

This section introduces the weak formulations implemented in SFEMaNS and
 additional features/extensions of the code. The notations introduced
 previously, such as the domain of approximation for each equations or
 the time step \f$\tau\f$, are unchanged.

@subsection doc_intro_SFEMaNS_possibilities_nst Hydrodynamic setting


@subsubsection doc_intro_SFEMaNS_possibilities_nst_1 Approximation of the Navier-Stokes equations

The approximation of the Navier-Stokes equations is based on a
 rotational form of the prediction-correction projection method
 detailed in <a href='http://www.ams.org/journals/mcom/2004-73-248/S0025-5718-03-01621-1/'>
<code>Guermond and Shen (2004)</code></a>. As the code SFEMaNS
 approximates the predicted velocity, a penalty method of the
 divergence of the velocity field is also implemented.

 The method proceeds as follows:
<ol>
<li>Initialization of the velocity field, the pressure
 and the pressure increments.
<li>For \f$n\geq0\f$ let \f$\bu^{n+1}\f$, that
 matches the Dirichlet boundary conditions of the
 problem, be the solutions of the following formulation for all
  \f$\textbf{v}\f$ in \f$\textbf{V}_{h,0}^\bu\f$:
\f{equation}{
\label{eq:SFEMaNS_weak_from_NS_1}
\int_{\Omega_{c,f}} \frac{3}{2 \tau} \textbf{u}^{n+1} \cdot \textbf{v}
 + \frac{1}{\Re} \GRAD \textbf{u}^{n+1} :  \GRAD \textbf{v} 
 + \text{c}_\text{div} h^{-1} \DIV\bu^{n+1} \DIV \bv= 
 - \int_{\Omega_{c,f}} ( \frac{-4 \textbf{u}^n
 + \textbf{u}^{n-1}}{2 \tau}
  + \GRAD ( p^n +\frac{4\psi^n -  \psi^{n-1}}{3} ) ) \cdot \textbf{v}  \\
  + \int_{\Omega_{c,f}} ( \textbf{f}^{n+1} - (\ROT \textbf{u}^{*,n+1} )
 \times \textbf{u}^{*,n+1} ) \cdot \textbf{v} ,
\f}
where \f$h\f$ is the local mesh size, \f$\text{c}_\text{div}\f$
 is a penalty coefficent,
 \f$\textbf{u}^{*,n+1}=2\textbf{u}^n-\textbf{u}^{n-1}\f$ and 
\f$\GRAD \textbf{u} : \GRAD \textbf{v} = \sum_{i,j} \partial_i u_j \partial_i v_j\f$.
We remind that  \f$\Re\f$ is the kinetic Reynolds number
 and \f$\textbf{f}\f$ a source term.
<li>Computation of \f$\psi^{n+1}\f$ and \f$\delta^{n+1}\f$ 
 solutions in \f$S_h^p\f$ of:
\f{equation}{
\label{eq:SFEMaNS_weak_from_NS_2}
\int_{\Omega_{c,f}} \GRAD \psi^{n+1} \cdot \GRAD q 
= \frac{3}{2 \tau} \int_{\Omega_{c,f}} \textbf{u}^{n+1} \cdot \GRAD q,
\f}
\f{equation}{
\label{eq:SFEMaNS_weak_from_NS_3}
\int_{\Omega_{c,f}} q \delta^{n+1} = \int_{\Omega_{c,f}} q \DIV \textbf{u} ^{n+1},
\f}
 for all \f$q\f$ in \f$S_h^p\f$. 
<li>The pressure is updated as follows:
\f{equation}{
\label{eq:SFEMaNS_weak_from_NS_4}
 p^{n+1} = p^n + \psi^{n+1} - \frac{1}{\Re} \delta^{n+1} .
\f}
</ol>

@subsubsection doc_intro_SFEMaNS_possibilities_nst_2 Entropy viscosity for under resolved computation

Hydrodynamic problems with large kinetic Reynolds number 
 introduce extremely complex flows. Approximating all of 
 the dynamics's scales of such problems is not always possible
 with present computational ressources. To address this problem,
 a nonlinear stabilization method called entropy viscosity is 
 implemented in SFEMaNS. This method has been introduced by
 <a href='http://link.springer.com/article/10.1007%2Fs10915-010-9445-3'>
 <code>Guermond et al. (2011)</code></a>. It consists of introducing an artifical
 viscosity, denoted \f$\nu_{E}\f$, that is taken proportional
 to the default of equilibrium of an energy equation.

This implementation of this method in SFEMaNS can be summarized
 in the three following steps:
<ol>
<li>Define the residual of the Navier-Stokes at
 the time \f$t_n\f$ as follows:
\f{equation}{
\textbf{Res}_\text{NS}^n:=
\frac{\bu^n-\bu^{n-2}}{ 2 \tau}
-\frac{1}{\Re} \LAP \bu^{n-1}
+ \ROT (\bu^{*,n-1}) \times \bu^{*,n-1}
  + \GRAD p^{n-1} -\textbf{f}^{n-1} ,
\f}
<li>Compute the entropy viscosity on each mesh cell K as follows:
\f{equation}{
\label{eq:SFEMaNS_NS_entropy_viscosity}
\nu_{E|K}^{n}:=\min\left(c_\text{max} h \|\bu^{n-1}\|_{\bL^\infty(K)},
 c_\text{e} h^2 \frac{\|\textbf{Res}_\text{NS}^n \cdot
 \bu^{n-1}\|_{\bL^\infty(K)}}{\|\bu^{n-1}\|_{\bL^\infty(K)}^2}\right),
\f}
with \f$h\f$ the local mesh size of the cell K, 
 \f$c_\text{max}=\frac{1}{2}\f$ for \f${\mathbb P}_1\f$ 
finite elements and \f$c_{\max}=\frac{1}{8}\f$ for \f${\mathbb P}_2\f$
 finite elements. The coefficient \f$c_\text{e}\f$ is a tunable
 constant in the interval \f$(0,1)\f$. It is generally set to one.
<li>When approximating \f$\bu^{n+1}\f$, the term
 \f$-\DIV(\nu_{E}^n \GRAD \bu^n)\f$
 is added in the left handside of the Navier-Stokes equations.
</ol>

Thus defined, the entropy viscosity is expected to be smaller
 than the consistency error in the smooth regions. In regions
 with large gradients, the entropy viscosity switches to the first
 order viscosity \f$\nu_{\max|K}^n:=c_\text{max} h_K \|\bu^{n-1}\|_{\bL^\infty(K)}\f$.
 Note that \f$\nu_\max^n\f$ corresponds to the artifical viscosity
 induces by first order up-wind scheme in the finite difference
 and finite volume litterature.

Remark: To facilitate the explicit treatment of the entropy viscosity, 
the following term can be added in the left handside of the Navier-Stokes
equations:
\f{equation}{
\label{eq:SFEMaNS_NS_LES_c1}
- \DIV( c_1 h \GRAD (\bu^{n+1}-\bu^{*,n+1})).
\f}
with \f$h\f$ the local mesh size and \f$c_1\f$ is a tunable constant.
 The coefficient \f$c_1\f$ should be of the same order of
 \f$c_\text{max} \|\bu\|_{\bL^\infty(\Omega_{c,f})}\f$.

@subsubsection doc_intro_SFEMaNS_possibilities_nst_3 Extension to non axisymmetric geometry

A penalty method of  <a href='http://www.sciencedirect.com/science/article/pii/S0168927407000815'>
<code>Pasquetti et al. (2008)</code></a>) is implemented so 
 the code SFEMaNS can report of the presence of non axisymmetric
 solid domain in \f$\Omega_{c,f}\f$. Such solid domains can either
 be driving the fluid or represents an obtacle to the fluid motion
 when their velocity is zero. The domain \f$\Omega_{c,f}\f$, where
 the Navier-Stokes equations are approximated, is splitted into a 
 fluid domain \f$\Omega_\text{fluid}\f$ and a solid
 domain \f$\Omega_\text{obs}\f$. These sub domains can be non
 axisymmetric and time dependent. The penalty method introduces
 a penalty function \f$\chi\f$. It is used to force the velocity
 field approximated by the Navier-Stokes equations to
 match the given velocity field of the solid in \f$\Omega_\text{obs}\f$.
 This penalty function is defined as follows:
\f{equation}{
\label{eq:SFEMaNS_NS_penal_1}
\chi = 
\left\{
\begin{array}{c}
 1 \text{ in } \Omega_\text{fluid}, \\
 0 \text{ in } \Omega_\text{obs}.
\end{array}
\right.
\f}

The velocity field is updated as follows:
\f{equation}{
\label{eq:SFEMaNS_NS_penal_2}
\frac{3\bu^{n+1}}{2\tau} 
- \frac{1}{\Re} \LAP \bu^{n+1} 
=   
-\GRAD p^{n} 
+  \chi^{n+1} \left(\frac{4\bu^n - \bu^{n-1}}{2\tau} 
- \GRAD( \frac{4\psi^n-\psi^{n-1}}{3}) \right) 
\\
+ \chi^{n+1} \left( 
- ( \ROT \bu^{*,n+1} ) \times\bu^{*,n+1} 
 + \textbf{f}^{n+1} \right)
 + (1 - \chi^{n+1}) \frac{3\bu^{n+1}_\text{obst}}{2\tau},
\f}
 with \f$\bu_\text{obs}\f$ the given velocity of the solid obstacle. 
 As the subdomains \f$\Omega_\text{fluid}\f$ and \f$\Omega_\text{obs}\f$
 can be time dependent so is the penalty function \f$\chi\f$.
 Note that the original scheme is recovered where \f$\chi=1\f$.

Remark: the correction and update of the pressure is not modified.

@subsubsection doc_intro_SFEMaNS_possibilities_nst_4 Extension to multiphase flow problem

The code SFEMaNS can approximate two phase flow problems.
 The governing equations are written as follows:
\f{equation}{
\label{eq:SFEMaNS_NS_multiphase_1}
\partial_t \rho + \DIV( \textbf{m}) = 0,
\f}
\f{equation}{
\label{eq:SFEMaNS_NS_multiphase_2}
\partial_t(\textbf{m}) 
+ \DIV(\textbf{m}{\otimes}\bu)
- \frac{2}{\Re} \DIV(\eta \varepsilon(\bu)) 
= -\GRAD p +  \textbf{f},
\f}
\f{equation}{
\label{eq:SFEMaNS_NS_multiphase_3}
\DIV \bu = 0,
\f}
where \f$\rho\f$ is the density, \f$\eta\f$ the dynamical viscosity,
 \f$\bm=\rho\bu\f$ the momentum, \f$\textbf{f}\f$ a forcing term and 
\f$\varepsilon(\bu)=\GRAD^s \bu = \frac12 (\GRAD \bu +(\GRAD \bu)^\sf{T})\f$.
 The densities, respestively dynamical viscosities, of the two fluids are denoted
 \f$\rho_0\f$ and \f$\rho_1\f$, respectively \f$\eta_0\f$ and \f$\eta_1\f$.


The approximation method is based on the following ideas.
<ol>
<li>Use of a level set method to follow the interface evolution.
 The method consists of approximating \f$\varphi\f$ that takes
 value in \f$[0,1]\f$ solution of:
\f{align*}{
\partial_t \varphi + \bu \cdot \GRAD \varphi=0.
\f}
The level set is equal to 0 in a fluid and 1 in the other fluid. 
 The interface between the fluid is represented by \f$\varphi^{-1}({1/2})\f$.
<li>Use the momentum as dependent variable for the Navier-Stokes equations. 
The mass matrix becomes time independent and can be treated with pseudo-spectral method.
<li>Rewritte the diffusive term \f$- \frac{2}{\Re} \DIV(\eta \varepsilon(\bu)) \f$ as follows:
\f{align*}{
- \frac{2}{\Re} \DIV(\eta \varepsilon(\bu)) = 
- \frac{2}{\Re} \DIV(\overline{\nu} \varepsilon(\bm))
- \left( \frac{2}{\Re} \DIV(\eta \varepsilon(\bu))
- \frac{2}{\Re} \DIV(\overline{\nu} \varepsilon(\bm)) \right)
\f}
with \f$\overline{\nu}\f$ a constant satisfying \f$\overline{\nu}\geq \frac{\eta}{\rho}\f$.
 The first term is made implicit while the second is treated explicitly.
 The resulting stiffness matrix is time independent and does not involve nonlinearity.
<li>The level set and Navier-Stokes equations are stabilized with the same entropy viscosity.
 For each mesh cell \f$K\f$ and each time iteration \f$n\f$,
 the entropy viscosity \f$\nu_E\f$ is defined as follows:
\f{align*}{
\nu_{E|K}^{n}:=\min\left(c_\text{max} h \|\bu^{n-1}\|_{\bL^\infty(K)},
 c_\text{e} h^2 \frac{\ \max( 
 |\textbf{Res}_\text{NS}^n \cdot  \bu^{n-1}\|_{\bL^\infty(K)}, 
|\text{Res}_\rho^n \|\bu{n-1}\|^2|
}{\|\bu^{n-1}\|_{\bL^\infty(K)}\|\bm^{n-1}\|_{\bL^\infty(K)}}
\right),
\f}
where 
\f{align*}{
\textbf{Res}_\text{NS}^n=
\frac{\bm^n-\bm^{n-2}}{ 2 \tau}
-\frac{1}{\Re} \DIV (\eta^{n-1}\epsilon(\bu^{n-1}))
+ \DIV(\bm^n{\otimes}\bu^n) + \GRAD p^{n-1} -\textbf{f}^{n-1} ,
\f}
and
\f{align*}{
\text{Res}_\rho^n= \frac{\rho^n-\rho^{n-2}}{ 2 \tau}
+ \DIV (\bm^{m-1}).
\f}
To facilitate the explicit treatment of the entropy viscosity, 
the term \f$- \DIV( c_1 h \GRAD (\bu^{n+1}-\bu^{n}))\f$, respectively
 \f$-\DIV( c_1 h \GRAD (\varphi^{n+1}-\varphi^n))\f$, can be added
 in the left handside of the Navier-Stokes, respectively of level set equation.
<li>A compression term that allows the level set to not get flatten over time
 iteration is added. It consists of adding the following term in the right
 handside of the level set equation:
\f{align*}{
 - \DIV \left(c_\text{comp}\nu_E h^{-1} \varphi(1-\varphi)\frac{\GRAD\varphi}{\|\varphi\|}\right).
\f}
The coefficient \f$c_\text{comp}\f$ a tunable constant in \f$[0,1]\f$.
  We generally set \f$c_\text{comp}=1\f$.
</ol>

After initialization, the first time order algorithm (BDF1) proceeds as follows:
<ol>
<li>Compute \f$\varphi^{n+1}\f$ solution of
\f{align}{
\frac{\varphi^{n+1}-\varphi^n}{\tau}  =  - \bu^n \cdot \GRAD \varphi^n
 + \DIV \left(
 \nu_E^n\GRAD \varphi^n
 - c_\text{comp} \nu_E^n h^{-1}  \varphi^n(1-\varphi^n)\frac{\GRAD\varphi^n}{\|\varphi^n\|} 
   \right).
\f}
<li>Reconstruct \f$\rho^{n+1}\f$ and \f$\eta^{n+1}\f$ as follows:
\f{align*}{
\rho^{n+1} =  \rho_0 + (\rho_1 - \rho_0) F(\varphi^{n+1}), \qquad
\eta  =  \eta_0 + (\eta_1 - \eta_0) F(\varphi^{n+1}),
\f}
where \f$F\f$ is either equal to the identity,
 \f$F(\varphi)=\varphi\f$, or a piecewise ponylomial function defined by:
\f{align*}{
F(\varphi) = 
\begin{cases}
0 & \text{if $\varphi - 0.5\le -c_{\text{reg}}$}, \\
\frac12 \left(1+\frac{(\varphi-0.5)((\varphi-0.5)^2 - 3 c_{\text{reg}}^2)}{-2c^3_{\text{reg}}}\right)
& \text{if $|\varphi - 0.5| \le c_{\text{reg}}$}, \\
1 & \text{if $c_{\text{reg}} \le \varphi - 0.5$}.
\end{cases} 
\f}
The tunable coefficient \f$c_\text{reg}\f$ lives in \f$[0,0.5]\f$. We generally set \f$c_\text{reg}=0.5\f$.
<li>Compute \f$\bm^{n+1}\f$ solution of:
\f{align}{
\frac{\bm^{n+1}-\bm^n}{\tau} - \frac{2\overline{\nu}}{\Re}\DIV(\epsilon(\bm^{n+1})-\epsilon(\bm^n))
= \frac{2}{\Re}\DIV( \eta^n\epsilon(\bu^n))
- \DIV(\bm^n\times\bu^n)
- \GRAD(p^n+\phi^n)
+\textbf{f}^{n+1}.
\f}
<li>Update the pressure as follows:
<ol>
<li>Solve \f$\phi^{n+1}\f$ solution of
\f{align*}{
- \LAP \phi^{n+1} = \frac{\rho_\text{min}}{\tau} \DIV \bu^{n+1},
\f}
with \f$\rho_\text{min}=\min(\rho_0,\rho_1)\f$ and \f$\bu^{n+1}=\frac{1}{\rho^{n+1}}\bm^{n+1}\f$.
<li>Set \f$p^{n+1}=p^n + \phi^{n+1} - \frac{\eta_\text{min}}{\Re} \DIV \bu^{n+1}\f$
 with \f$\eta_\text{min}=\min(\eta_0,\eta_1)\f$.
</ol>
</ol>


Remarks:
<ol>
<li>This method can be used to approximate problems with
 a stratification or an inclusion of \f$n\geq 3\f$ fluids.
 One level set is approximated per interface between two
 fluids. The fluids properties are reconstructed with
 recursive convex combinations. 
<li>MHD multiphase problems with variable electrical conductivity
 between the fluids can also be considered. The electrical
 conductivity in the fluid is reconstructed with the level set
 the same way the density and the dynamical viscosity are.
 The magnetic field \f$\bH^{n+1}\f$ is updated as follows:
\f{align}{
\frac{3\bH^{n+1}-4\bH^n+\bH^{n-1}}{2\tau} 
 + \ROT \left( \frac{1}{\overline{\sigma}\Rm}
 \ROT ( \bH^{n+1}-\bH^{*,n+1}) \right)
=
 - \ROT\left( \frac{1}{\sigma\Rm} \ROT \bH^{*,n+1} \right)
 + \ROT (\bu^{n+1}\times \mu^c \bH^{*,n+1})
 + \ROT \left( \frac{1}{\sigma\Rm} \textbf{j}^{n+1} \right)
\f}
with \f$\bH^{*,n+1}=2\bH^{n+1}-\bH^n\f$ and \f$\overline{\sigma}\f$ a
 function depending of the radial and vertical
 coordinates \f$(r,z)\f$ such that
 \f$\overline{\sigma}(r,z)\leq \sigma(r,\theta,z,t)\f$ for
 all \f$(r,\theta,z,t)\f$ considered.
</ol>



@subsection doc_intro_SFEMaNS_possibilities_temp Heat equation's weak formulation

The heat equations is approximated as follows.
<ol>
<li>Initialization of the temperature.
<li>For all \f$n\geq0\f$ let \f$T^{n+1}\f$, that matches the
 Dirichlet boundary conditions of the problem, be the solution
 of the following formulation for all \f$v\in S_h^T\f$:
\f{equation}{
\label{eq:SFEMaNS_weak_form_temp}
\int_{\Omega_T} \frac{3 C }{2 \tau}T^{n+1} v
 + \lambda \GRAD T^{n+1} \cdot \GRAD v 
 =  - \int_{\Omega_T} \left( \frac{4 T^n -T^{n-1}}{2 \tau} 
    - \DIV (T^{*,n+1} \bu^{*,n+1}) + f_T^{n+1}\right) v, 
\f}
where \f$T^{*,n+1}=2 T^n - T^{n-1}\f$. We remind that \f$C\f$ is
 the volumetric heat capacity, \f$\lambda\f$ the thermal conductivty
 and \f$f_T\f$ a source term.
</ol>



@subsection doc_intro_SFEMaNS_possibilities_mxw Magnetic setting

The code SFEMaNS uses \f$\bH^1\f$ conforming Lagrange finite element to approximate
 the magnetic field. As a consequence, the zero divergence condition on the
 magnetic field cannot be enforced by standard penalty technique for
 non-smooth and non-convex domains. 
 To overcome this obstacle, a method inspired of
 <a href='http://www.ams.org/journals/mcom/2011-80-276/S0025-5718-2011-02464-6/'>
<code>Bonito and Guermond (2011)</code></a>
 has been implemented. This method consists of introducting a
 magnetic pressure denoted \f$ p_\text{m}\f$ and proceeds as follows.
<ol>
<li>Add the term \f$-\mu^c \GRAD p_\text{m}^c\f$ in the right handside
 of the magnetic field \f$\bH^c\f$ equation where \f$p_\text{m}^c\f$
 if the solution in \f$\Omega^c\f$ of:
\f{align*}{
- \DIV( h_\text{loc}^{2(1-\alpha)} \GRAD p_m^{c,n+1} ) &=
 - \DIV( \mu^c \bH^{c,n+1}) ,  
 \\
p_m^{c,n+1}|_{\partial \Omega_c} &= 0, 
\f}
where \f$h\f$ is the local mesh size and \f$\alpha\f$ a
 constant parameter in \f$[0.6,0.8]\f$.
<li>Add the term \f$ -\DIV(\mu^v \GRAD p_\text{m}^v)\f$ in the right handside
 of the scalar potential \f$\phi\f$ equation where  \f$p_\text{m}^v\f$
 is the solution in \f$\Omega^v\f$ of:
\f{align*}{
\LAP p_m^{v,n+1} = \LAP \phi^{n+1}, \\
\GRAD p_m^{v,n+1} \cdot \textbf{n} |_{\partial \Omega_v} = 0. 
\f}
</ol>

We note that the magnetic pressure can be eliminated from the equation
 of the scalar potential \f$\phi\f$. We refer to
 <a href='http://www.sciencedirect.com/science/article/pii/S0021999111002749'>
 <code>Guermond et al. (2011)</code></a> for more details.
 The approximation space used
 to approximate \f$ p_\text{m}^c\f$ is the following:
\f{align*}{
S_h^{p_\text{m}^c} :=
\left\{
\varphi= 
\sum\limits_{k=-M}^M \varphi_h^k (r,z) e^{ik \theta} ;
 \varphi_h^k \in  S_{h}^{p_\text{m}^c,2D}, \; -M \leq k \leq M
\right\},
\f}
where we introduce the following finite element space: 
\f{align*}{
S_{h}^{p_\text{m}^c, 2D}  : = \left\{ \varphi_h \in   C^0(\overline{\Omega_{c}^{2D}}) ;
 \varphi_h|_K \in \mathbb{P}_1^2 \text{ }  \forall K \in \mathcal{T}_h^c ,  \right\}.
\f}

In addition, an interior penalty method is used to enforce the continuity conditions
 across the interfaces \f$\Sigma_\mu\f$ and \f$\Sigma\f$. We refer to
 <a href='http://www.sciencedirect.com/science/article/pii/S0021999106002944'>
 <code>Guermond et al. (2007)</code></a> for more details.

@subsubsection doc_intro_SFEMaNS_possibilities_mxw_1 Approximation of the Maxwell equations with H

The Maxwell equations are approximated as follows:
<ol>
<li>Initialization of the magnetic field \f$\bH^c\f$, the scalar potential \f$\phi\f$ and the magnetic pressure \f$p_\text{m}^c\f$.
<li>For all \f$n\geq 1\f$, computation of \f$(\bH^{c,n+1},\phi^{n+1},p_\text{m}^{c,n+1})\f$
 solutions of the following formulation for all \f$b\in \bV_h^{\bH^c} \f$,
 \f$\varphi\in S_h^{\phi}\f$
 and \f$q\in S_h^{p_\text{m}^c} \f$.
\f{align*}{
& \int_{\Omega_c}\mu^c \frac{D\bH^{c,n+1}}{\Delta t}\SCAL \bb
  +\int_{\Omega_c} \frac{1}{\sigma R_m} \ROT  \bH ^{c,n+1}\cdot \ROT \bb 
  +\int_{\Omega_v} \muv\frac{\GRAD D\phi^{n+1}}{\Delta t}\SCAL \GRAD\varphi 
  +\int_{\Omega_v} \muv\GRAD\phi^{n+1}\SCAL \GRAD\varphi -
  \int_{\partial\Omega_v} \muv\varphi \bn\SCAL \GRAD \phi^{n+1}\\
& + \beta_1\left(\int_{\Omega_c} \mu^c\GRAD p_\text{m}^{c,n+1}\SCAL\bb
  - \int_{\Omega_c} \mu^c\bH^{c,n+1}\SCAL\GRAD q +
  \int_{\Omega_c} h^{2(1-\alpha)}\GRAD p_\text{m}^{c,n+1}\SCAL \GRAD q  
  + \int_{\Omega_c}
  h^{2\alpha}\DIV (\mu^c \bH^{c,n+1} )\DIV (\mu^c  \bb)\right)\\
& +\int_{\Sigma_{\mu}} \left \{  \frac{1}{\sigma R_m} \ROT {\bH ^{c,n+1}}  \right \}
  \cdot  \left ( { \bb_1}\times \bn_1^c + { \bb_2}\times \bn_2^c\right )\\
& +\beta_3 \int_{\Sigma_{\mu}} h^{-1}  \left( { \bH_1^{c,n+1}}\times \bn_1^c
  + {\bH_2^{c,n+1}}\times \bn_2^c\right ) \SCAL   \left ( { \bb_1}\times \bn_1^c + { \bb_2}\times \bn_2^c\right )\\
& +\beta_1 \int_{\Sigma_{\mu}} h^{-1} \left({ \mu^c_1\bH_1^{c,n+1}}\cdot \bn_1^c 
  + {\mu^c_2 \bH_2^{c,n+1}}\cdot \bn_2^c\right ) \SCAL   \left ( {\mu^c_1}{ \bb_1}\cdot \bn_1^c
  + {\mu^c_2}{ \bb_2}\cdot \bn_2^c\right )\\
& +\int_{\Sigma} \frac{1}{\sigma R_m} \ROT {\bH ^{c,n+1}} \cdot
  \left( { \bb }\times  \bn^c +  \nabla \varphi ^{n+1}\times \bn^v\right)
  + \beta_2  \int_\Sigma h^{-1} \left( {\bH^{c,n+1}}\CROSS \bn_1^c
  + {\GRAD \phi^{n+1}}\CROSS \bn_2^c\right )  \SCAL (\bb\CROSS \bnc +
  \GRAD\varphi\CROSS \bnv)\\
& + \beta_1 \int_\Sigma h^{-1} \left( { \mu^c\bH ^{c,n+1}}\cdot \bn_1^c
  + {\GRAD \phi^{n+1}}\cdot \bn_2^c\right )  \SCAL ({\mu^c}\bb\cdot \bnc +
  \GRAD\varphi \cdot \bnv)\\
& + \int _{\Gamma_c}  \frac{1}{\sigma R_m} \ROT   \bH ^{c,n+1} \cdot ( \bb  \CROSS \bnc)
  + \beta _3\left(
  \int_{\Gamma_c} h^{-1} \left( { \bH^{c,n+1}}\CROSS \bn^c \right )  \SCAL (\bb\CROSS \bnc)
  \right )\\
& =\\
& \int_{\Omega_c} \left( \frac{1}{\sigma R_m}\bj^s + \bu^{n+1} \times  \mu^c \bH^{*,n+1} \right )
  \cdot \ROT \bb  
 + \int _{\Sigma_{\mu}} \left \{  \frac{1}{\sigma R_m}\bj^s  +
  \bu^{n+1} \times \mu^c \bH^{*,n+1}  \right \} \cdot
  \left( { \bb_1}\times \bn_1^c + { \bb_2}\times \bn_2^c\right )\\
& +\int_{\Sigma}\left (  \frac{1}{\sigma R_m} \bj^s + \bu^{n+1} \times \mu^c \bH^{*,n+1}
  \right)\cdot \left ( { \bb }\times  \bn^c +  \nabla \varphi \times \bn^v\right)
  +\int_{\Gamma_c}(\ba \times \bn) \cdot \left ({\bb} \times \bn \right) + \int_{\Gamma_v}
  (\ba \times \bn) \cdot (\nabla \varphi \times \bn)\\
& + \int_{\Gamma_c} \left ( \frac{1}{\sigma R_m}\bj^s   + \bu^{n+1} \times 
  \mu^c \bH^{*,n+1} \right )\cdot ( \bb  \CROSS \bnc)
  +\beta_3 \int_{\Gamma_c} h^{-1} 
  \left( {\bH}_\text{bdy}^{c,n+1}\CROSS \bn^c \right)  \SCAL (\bb\CROSS \bnc) ,
\f}
where we set \f$D\bH^{c,n+1}=\dfrac{3\bH^{c,n+1}-4\bH^{c,n}+\bH^{c,n-1}}{2}\f$, 
\f$D\phi^{c,n+1}=\dfrac{3\phi^{c,n+1}-4\phi^{c,n}+\phi^{c,n-1}}{2}\f$,
 \f$\bH^{*,n+1}=2\bH^{c,n}-\bH^{c,n-1}\f$. We use the operator \f$\{.\}\f$ defined by
 \f$\{f\}=\frac{f_1+f_2}{2}\f$ on the interface \f$\Sigma_\mu\f$.
 The constants \f$\beta_1, \beta_2\f$ and \f$\beta_3\f$ are penalty coefficients.
 They are normalized by \f$(\sigma\Rm)^{-1}\f$ so their value can be set to one
 in the data file. The function \f$\bH_\text{bdy}^{c}\f$ is a user function
 used to impose Dirichlet boundary conditions on the surface
 \f$\Gamma_c=\partial\Omega_c \setminus \Sigma\f$.

</ol>





@subsubsection doc_intro_SFEMaNS_possibilities_mxw_2 Extension to magnetic permeability variable in time and azimuthal direction

The use of a Fourier decomposition in the azimuthal direction leads us to use
 the magnetic field \f$\bB^c=\mu\bH^c\f$ as dependent variable of the Maxwell equations
 in the conducting domain. The mass matrix becomes time independent and can be computed with pseudo-spectral methods.
 To get a time independent stiffness matrix that does not involve nonlinearity, the diffusive term
 \f$\ROT \left(\frac{1}{\sigma\Rm} \ROT\frac{\bB^c}{\mu} \right)\f$ is rewritten as follows:
\f{align*}{
   \ROT \left( \frac{1}{\sigma\Rm} \ROT \frac{\bB^c}{\mu} \right) = 
   \ROT \left( \frac{1}{\sigma\Rm \overline{\mu}} \ROT\frac{\bB^c}{\mu} \right)
 + \ROT \left( \frac{1}{\sigma\Rm} \ROT ((\frac{1}{\mu}-\frac{1}{\overline{\mu}})\bB^c) \right)
\f}
with \f$\overline{\mu}\f$ a function depending of the radial and vertical
 coordinates \f$(r,z)\f$ such that \f$\overline{\mu}(r,z)\leq \mu(r,\theta,z,t)\f$ for
 all \f$(r,\theta,z,t)\f$ considered. The first term is then made implicit while
 the term involving \f$\frac{1}{\mu}\f$ is treated explicitly.


Under the previous notations and assuming,
\f{align*}{
 \bar{\mu} &= \mu \quad \text{on} \quad \Sigma, \\
\bar{\mu} &= \mu \quad \text{on}  \quad \Sigma_{\mu},\\
\bar{\mu} &= \mu \quad \text{on}  \quad \Gamma_c,
\f}
the Maxwell equations are approximated as follows.
<ol>
<li>Initialization of the magnetic field \f$\bB^c\f$, the scalar potential \f$\phi\f$ and the magnetic pressure \f$p_\text{m}^c\f$.
<li>For all \f$n\geq 1\f$, computation of \f$(\bB^{c,n+1},\phi^{n+1},p_\text{m}^{c,n+1})\f$
 solutions of the following formulation for all \f$b\in \bV_h^{\bH^c} \f$,
 \f$\varphi\in S_h^{\phi}\f$
 and \f$q\in S_h^{p_\text{m}^c} \f$.
\f{align*}{
& \int_{\Omega_c}\frac{D\bB^{c,n+1}}{\Delta t}\SCAL \bb
  + \int _{\Omega_c} \frac{1}{\sigma R_m} \ROT  \frac{\bB ^{c,n+1}}{\overline{\mu^c}}\cdot \ROT \bb 
  + \int_{\Omega_v} \muv\frac{\GRAD D\phi^{n+1}}{\Delta t}\SCAL \GRAD\varphi 
  + \int_{\Omega_v} \muv\GRAD\phi^{n+1}\SCAL \GRAD\varphi 
  - \int_{\partial\Omega_v} \muv\varphi \bn\SCAL \GRAD \phi^{n+1}\\
& + \beta_1\left(\int_{\Omega_c} \overline{\mu^c}\GRAD p_\text{m}^{c,n+1}\SCAL\bb
  - \int_{\Omega_c} \bB^{c,n+1}\SCAL\GRAD q + \int_{\Omega_c} h^{2(1-\alpha)}
  \GRAD p_\text{m}^{c,n+1}\SCAL \GRAD q  
  + \int_{\Omega_c} h^{2\alpha}
  \overline{\mu^c} \DIV \bB^{c,n+1} \DIV   \bb \right)\\
& +\int _{\Sigma_{\mu}} \left\{ \frac{1}{\sigma R_m} 
  \ROT \frac{\bB ^{c,n+1}}{\overline{\mu^c}}   \right \}
  \cdot  \left ( { \bb_1}\times \bn_1^c + { \bb_2}\times \bn_2^c\right )\\
& +\beta_3 \int_{\Sigma_{\mu}} h^{-1} \left(
  \frac{\bB_1^{c,n+1}}{\overline{\mu^c}_1}\times \bn_1^c + \frac{\bB_2^{c,n+1}}{\overline{\mu^c_2}}\times \bn_2^c
  \right) \SCAL   \left ( { \bb_1}\times \bn_1^c + { \bb_2}\times \bn_2^c\right )\\
& +\beta_1 \int_{\Sigma_{\mu}} h^{-1}
  \left( {\bB_1^{c,n+1}}\cdot \bn_1^c + {\bB_2^{c,n+1}}\cdot \bn_2^c\right)
  \SCAL   \left( \overline{\mu^c_1}{ \bb_1}\cdot \bn_1^c + \overline{\mu^c_2}{ \bb_2}\cdot \bn_2^c\right )\\
& +\int _{\Sigma} \frac{1}{\sigma R_m} \ROT \frac{\bB ^{c,n+1}}{\overline{\mu^c}} \cdot
  \left( {\bb }\times  \bn^c +  \nabla \varphi \times \bn^v\right)
  + \beta_2 \int_{\Sigma}  h^{-1}
  \left( \frac{\bB^{c,n+1}}{\overline{\mu^c}}\CROSS \bn_1^c + {\GRAD \phi ^{n+1}}\CROSS \bn_2^c\right)
  \SCAL (\bb\CROSS \bnc +  \GRAD\varphi\CROSS \bnv)\\
& + \beta_1 \int_{\Sigma} h^{-1}
  \left( {\bB ^{c,n+1}}\cdot \bn_1^c + {\GRAD \phi ^{n+1}} \cdot \bn_2^c\right)
  \SCAL \left(\overline{{\mu^c}}\bb\cdot \bnc +
  \GRAD\varphi \cdot \bnv \right )\\
& + \int_{\Gamma_c}  \frac{1}{\sigma R_m} \ROT \frac{\bB ^{c,n+1}}{\overline{\mu^c}}
  \cdot ( \bb  \CROSS \bnc) + \beta_3\left( \int_{\Gamma_c}  h^{-1}
  \left( \frac{ \bB_\text{bdy}^{c,n+1}}{\overline{\mu^c}}\CROSS \bn^c \right )  \SCAL (\bb\CROSS \bnc)
  \right)\\
& =\\
& \int_{\Omega_c} \frac{1}{\sigma R_m} \ROT  (\langle \overline{\mu^c},{\mu^c}
  \rangle {\bB^{*,n+1}})\cdot \ROT \bb \\
& \int_{\Omega_c} \left( \frac{1}{\sigma R_m}\bj^s + \bu^{n+1} \times 
  \bB^{*,n+1} \right )\cdot \ROT \bb
  + \int_{\Sigma_{\mu}} \left \{  \frac{1}{\sigma R_m}\bj^s  
  + \bu^{n+1} \times \bB^{*,n+1}  \right \} \cdot  
  \left( { \bb_1}\times \bn_1^c + { \bb_2}\times \bn_2^c\right )\\
& +\int_{\Sigma}\left(  \frac{1}{\sigma R_m} \bj^s
  + \bu^{n+1} \times \bB^{*,n+1}  \right )
  \cdot \left ( { \bb }\times  \bn^c +  \nabla \varphi \times \bn^v\right)
 +\int_{\Gamma_c}(\ba \times \bn) \cdot \left ({\bb} \times \bn \right)
  + \int_{\Gamma_v}(\ba \times \bn) \cdot (\nabla \varphi \times \bn)\\
& + \int_{\Gamma_c} \left ( \frac{1}{\sigma R_m}\bj^s  + \bu^{n+1}
  \times  \bB^{*,n+1} \right )\cdot ( \bb  \CROSS \bnc)
  +\beta_3 \int_{\Gamma_c} h^{-1} \left({\bB}_\text{bdy}^{c,n+1}\CROSS \bn^c \right )  \SCAL (\bb\CROSS \bnc),
\f}
where we set \f$\bB^{*,n+1}=2\bB^n-\bB^{n-1}\f$ and  
 \f$\langle \overline{\mu^c},{\mu^c}\rangle=\frac{1}{\overline{\mu^c}}- \frac{1}{\mu^c}\f$.
 The function \f$\bB_\text{bdy}^{c}\f$ is a user function  used to impose Dirichlet
  boundary conditions on the surface \f$\Gamma_c=\partial\Omega_c \setminus \Sigma\f$.



 */
