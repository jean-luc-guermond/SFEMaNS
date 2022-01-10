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
 * @page MaxwellScheme_wH Maxwell Numerical Scheme using H

Complete scheme using \f$\bH\f$,


\f{align*}{
&\bH^c |_{t=0}  = \bH_0 ^c, \quad \phi |_{t=0}= \phi _0,\\
& \int_{\Omegac}\mu^c \frac{D\bH^{c,n+1}}{\Delta t}\SCAL \bb
+ \int_{\Omegav} \muv\frac{\GRAD D\phi^{n+1}}{\Delta t}\SCAL \GRAD\varphi 
+\int _{\Omega_c} \frac{1}{\sigma R_m} \ROT  \bH ^{c,n+1}\cdot \ROT \bb 
\\
& +\int _{\Sigma_{\mu}} \left \{  \frac{1}{\sigma R_m} \ROT {\bH ^{c,n+1}}  \right \} \cdot  \left ( { \bb_1}\times \bn_1^c + { \bb_2}\times \bn_2^c\right )\\
& +\beta_3 \sum_{F\in
    \Sigma_{\mu }} h_F^{-1}\int_{F}   \left ( { \bH_1^{c,n+1}}\times \bn_1^c + {\bH_2^{c,n+1}}\times \bn_2^c\right ) \SCAL   \left ( { \bb_1}\times \bn_1^c + { \bb_2}\times \bn_2^c\right )\\
& +\beta_1 \sum_{F\in
    \Sigma_{\mu }} h_F^{-1}\int_{F}   \left ( { \mu^c_1\bH_1^{c,n+1}}\cdot \bn_1^c + {\mu^c_2 \bH_2^{c,n+1}}\cdot \bn_2^c\right ) \SCAL   \left ( {\mu^c_1}{ \bb_1}\cdot \bn_1^c + {\mu^c_2}{ \bb_2}\cdot \bn_2^c\right )\\
&+\int _{\Sigma} \frac{1}{\sigma R_m} \ROT {\bH ^{c,n+1}} \cdot \left ( { \bb }\times  \bn^c +  \nabla \varphi ^{n+1}\times \bn^v\right)
+ \beta_2 \sum_{F\in \Sigma} h_F^{-1}\int_{F}  \left ( {\bH^{c,n+1}}\CROSS \bn_1^c + {\GRAD \phi^{n+1}}\CROSS \bn_2^c\right )  \SCAL (\bb\CROSS \bnc +
  \GRAD\varphi\CROSS \bnv)\\
&+ \beta_1 \sum_{F\in \Sigma} h_F^{-1}\int_{F}  \left ( { \mu^c\bH ^{c,n+1}}\cdot \bn_1^c + {\GRAD \phi^{n+1}}\cdot \bn_2^c\right )  \SCAL ({\mu^c}\bb\cdot \bnc +
  \GRAD\varphi \cdot \bnv)\\
& + \beta_1\left(\int_\Omegac \color{blue}{\mu^c}\GRAD p\SCAL\bb
- \int_\Omegac \color{blue}{\mu^c}\bH^{c,n+1}\SCAL\GRAD q + \sum_{K\in\calF_h^c}
\int_{K^{3D}} h_K^{2(1-\alpha)}\GRAD p\SCAL \GRAD q  
+ \sum_{K\in\calF_h^c} \int_{K^{3D}}
h_K^{2\alpha}\DIV (\color{blue}{\mu^c} \bH^{c,n+1} )\DIV (\color{blue}{\mu^c}  \bb)\right)\\
&+\int_{\Omega_v} \muv\GRAD\phi^{n+1}\SCAL \GRAD\varphi -
\int_{\partial\Omega_v} \muv\varphi \bn\SCAL \GRAD \phi^{n+1}\\
&+  
\int _{\Gamma^c _1}  \frac{1}{\sigma R_m} \ROT   \bH ^{c,n+1} \cdot ( \bb  \CROSS \bnc)
+ \beta _3\left (
\sum_{F\in \Gamma ^c _1} h_F^{-1}\int_{F}  \left ( { \bH^{c,n+1}}\CROSS \bn^c \right )  \SCAL (\bb\CROSS \bnc)
\right )\\
&=\\
& \int _{\Omega_c} \left( \frac{1}{\sigma R_m}\bj^s + \tilde {\bu} \times  \mu^c \bH^* \right )\cdot \ROT \bb  + \int _{\Sigma_{\mu}} \left \{  \frac{1}{\sigma R_m}\bj^s  + \tilde {\bu} \times \mu^c \bH^*  \right \} \cdot  \left ( { \bb_1}\times \bn_1^c + { \bb_2}\times \bn_2^c\right )\\
&+\int _{\Sigma}\left (  \frac{1}{\sigma R_m} \bj^s + \tilde {\bu} \times \mu^c \bH^*  \right )\cdot \left ( { \bb }\times  \bn^c +  \nabla \varphi \times \bn^v\right)
+\int _{\Gamma ^c _2}(\ba \times \bn) \cdot \left ({\bb} \times \bn \right) + \int_{\Gamma _v}(\ba \times \bn) \cdot (\nabla \varphi \times \bn)\\
&+ \int _{\Gamma^c _1} \left ( \frac{1}{\sigma R_m}\bj^s   + \tilde {\bu} \times  \mu^c \bH^* \right )\cdot ( \bb  \CROSS \bnc)
+\beta _3 \left ( \sum_{F\in \Gamma ^c _1} h_F^{-1}\int_{F}  \left ( {\bH}^{\text{given}}\CROSS \bn^c \right )  \SCAL (\bb\CROSS \bnc) \right).
\f}

<h3>The Code</h3>

 The implementation of the above scheme is implemented in the function <code> \link  update_maxwell_with_h::maxwell_decouple_with_h\endlink </code> which is in the file <b> <code>\link maxwell_update_time_with_H.f90\endlink</code></b>. We first describe what  <code> fortran </code> functions corresponds  to the terms of the  left hand side of the  scheme above. Then we do the same for the right hand side.

<h4>Left Hand Side</h4>
<ol>

<li> The implicit part  of 

\f{align*}{

& \int_{\Omegac}\mu^c \frac{D\bH^{c,n+1}}{\Delta t}\SCAL \bb
+ \int_{\Omegav} \muv\frac{\GRAD D\phi^{n+1}}{\Delta t}\SCAL \GRAD\varphi 
\f}

and the terms,
\f{align*}{

&
+\int _{\Omega_c} \frac{1}{\sigma R_m} \ROT  \bH ^{c,n+1}\cdot \ROT \bb \\
&+\int _{\Sigma} \frac{1}{\sigma R_m} \ROT {\bH ^{c,n+1}} \cdot \left ( { \bb }\times  \bn^c +  \nabla \varphi ^{n+1}\times \bn^v\right)
+ \beta_2 \sum_{F\in \Sigma} h_F^{-1}\int_{F}  \left ( {\bH^{c,n+1}}\CROSS \bn_1^c + {\GRAD \phi^{n+1}}\CROSS \bn_2^c\right )  \SCAL (\bb\CROSS \bnc +
  \GRAD\varphi\CROSS \bnv)\\
&+ \beta_1 \sum_{F\in \Sigma} h_F^{-1}\int_{F}  \left ( { \mu^c\bH ^{c,n+1}}\cdot \bn_1^c + {\GRAD \phi^{n+1}}\cdot \bn_2^c\right )  \SCAL ({\mu^c}\bb\cdot \bnc +
  \GRAD\varphi \cdot \bnv)\\
& + \beta_1\left(\int_\Omegac \mu^c\GRAD p\SCAL\bb
- \int_\Omegac \mu^c\bH^{c,n+1}\SCAL\GRAD q + \sum_{K\in\calF_h^c}
\int_{K^{3D}} h_K^{2(1-\alpha)}\GRAD p\SCAL \GRAD q  
+ \sum_{K\in\calF_h^c} \int_{K^{3D}}
h_K^{2\alpha}\DIV ({\mu^c} \bH^{c,n+1} )\DIV ({\mu^c}  \bb)\right)\\
&+\int_{\Omega_v} \muv\GRAD\phi^{n+1}\SCAL \GRAD\varphi -
\int_{\partial\Omega_v} \muv\varphi \bn\SCAL \GRAD \phi^{n+1}\\
\f}
are computed in the function <code>\link update_maxwell_with_h::mat_h_p_phi_maxwell\endlink</code>.

<li> The Dirichlet terms,
\f{align*}{
& 
\int _{\Gamma^c _1}  \frac{1}{\sigma R_m} \ROT   \bH ^{c,n+1} \cdot ( \bb  \CROSS \bnc)
+ \beta _3\left (
\sum_{F\in \Gamma ^c _1} h_F^{-1}\int_{F}  \left ( { \bH^{c,n+1}}\CROSS \bn^c \right )  \SCAL (\bb\CROSS \bnc)
\right ),\\
\f}


are computed in the function <code> \link update_maxwell_with_h::mat_dirichlet_maxwell\endlink</code>.

<li> The implicit terms  on the interface \f$ \Sigma_{\mu}\f$,

\f{align*}{
& \int _{\Sigma_{\mu}} \left \{  \frac{1}{\sigma R_m} \ROT {\bH ^{c,n+1}}  \right \} \cdot  \left ( { \bb_1}\times \bn_1^c + { \bb_2}\times \bn_2^c\right )\\
& +\beta_3 \sum_{F\in
    \Sigma_{\mu }} h_F^{-1}\int_{F}   \left ( { \bH_1^{c,n+1}}\times \bn_1^c + {\bH_2^{c,n+1}}\times \bn_2^c\right ) \SCAL   \left ( { \bb_1}\times \bn_1^c + { \bb_2}\times \bn_2^c\right )\\
& +\beta_1 \sum_{F\in
    \Sigma_{\mu }} h_F^{-1}\int_{F}   \left ( { \mu^c_1\bH_1^{c,n+1}}\cdot \bn_1^c + {\mu^c_2 \bH_2^{c,n+1}}\cdot \bn_2^c\right ) \SCAL   \left ( {\mu^c_1}{ \bb_1}\cdot \bn_1^c + {\mu^c_2}{ \bb_2}\cdot \bn_2^c\right )\\
\f}

are computed in the function <code>\link update_maxwell_with_h::mat_maxwell_mu\endlink</code>.
</ol>

We now describe the <code> fortran </code> functions concerning the right hand side.

<h4>Right Hand Side</h4>

<ol>

<li> The source (explicit) terms  of 
\f{align*}{
& \mu^c \frac{D\bH^{c,n+1}}{\Delta t}\SCAL \bb
+  \muv\frac{\GRAD D\phi^{n+1}}{\Delta t}\SCAL \GRAD\varphi 

\f}

are computed at the nodes in <code>\link update_maxwell_with_h::maxwell_decouple_with_h\endlink</code>. The integral of this sources are
computed  in the  function <code>\link update_maxwell_with_h::courant_int_by_parts \endlink</code>.

<li>The terms,
\f{align*}{
& \int _{\Omega_c} \left( \frac{1}{\sigma R_m}\bj^s + \tilde {\bu} \times  \mu^c \bH^* \right )\cdot \ROT \bb  &+\int _{\Sigma}\left (  \frac{1}{\sigma R_m} \bj^s + \tilde {\bu} \times \mu^c \bH^*  \right )\cdot \left ( { \bb }\times  \bn^c +  \nabla \varphi \times \bn^v\right)
\f}

are computed  in the  function <code>\link update_maxwell_with_h::courant_int_by_parts \endlink</code>.

<li>The term on the interface \f$ \Sigma_{\mu}\f$,
\f{align*}{
 \int _{\Sigma_{\mu}} \left \{  \frac{1}{\sigma R_m}\bj^s  + \tilde {\bu} \times \mu^c \bH^*  \right \} \cdot  \left ( { \bb_1}\times \bn_1^c + { \bb_2}\times \bn_2^c\right ),\\
\f}

is computed  in the  function <code>\link update_maxwell_with_h::surf_int \endlink</code>.

<li>The boundary terms,
\f{align*}{
&\int _{\Gamma ^c _2}(\ba \times \bn) \cdot \left ({\bb} \times \bn \right) + \int_{\Gamma _v}(\ba \times \bn) \cdot (\nabla \varphi \times \bn),\\
\f}

are computed  in the  function <code>\link update_maxwell_with_h::surf_int \endlink</code>.



<li>The Dirichlet boundary terms,
\f{align*}{
& \int _{\Gamma^c _1} \left ( \frac{1}{\sigma R_m}\bj^s   + \tilde {\bu} \times   \mu^c \bH^* \right )\cdot ( \bb  \CROSS \bnc)
 +\beta _3 \left ( \sum_{F\in \Gamma ^c _1} h_F^{-1}\int_{F}  \left ( {\bH}^{\text{given}}\CROSS \bn^c \right )  \SCAL (\bb\CROSS \bnc) \right),
\f}

are computed  in the  function <code>\link update_maxwell_with_h::rhs_dirichlet\endlink</code>.

 



</ol>





 */
