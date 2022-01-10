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
 * @page doc_SFEMaNS_condlim_source_in_NS The function source_in_NS_momentum


It is used to define the source term \f$\textbf{f}\f$ of the Navier-Stokes
 equations. We refer to the section
 \ref doc_intro_SFEMaNS_possibilities_nst_1 for more information
 on the formulation of the Navier-Stokes equations in SFEMaNS. 

This function defines the source term for one given Fourier mode, one given component (radial cosine, radial sine, azimuthal cosine, azimuthal sine, vertical cosine or vertical sine) on all the nodes of the finite element mesh.  We denote by vv_mesh the finite element mesh used to approximate the velocity field.

<h3>Inputs and outputs</h3>
The inputs of this function are the following:
<ol>
<li><tt>TYPE</tt> is the component of the source term that we compute (radial cosine, radial sine, etc.). It is an integer between one and six.
<li><tt>rr</tt> is a real valued tabular that contains two columns with dimensions (2,vv_mesh\%np). The tabular rr(1,:) contains the radial cylindrical coordinate of all the nodes of the finite element mesh vv_mesh. Respectively, rr(2,:) contains the vertical coordinates of these nodes.
<li><tt>mode</tt> is the Fourier mode considered. It is an integer.
<li><tt>i</tt> is the label that the processor associates to the Fourier mode considered. It is an integer.
<li><tt>time</tt> is the time at which the source term is computed. It is a real number.
<li><tt>opt_density</tt> is the density. This input is optional.
<li><tt>opt_tempn</tt> is the temperature. This input is optional.
</ol>

The output of this function is a real valued tabular vv of dimension SIZE(rr,2).

Remarks:
<ol>
<li>We note that opt_density and opt_tempn are optional because solving the Navier-Stokes equations does not imply to solve the temperature or the level set equations.
<li>The density is available when the following parameter is set to true in the data:
\code
===Is there a level set?
.t.
\endcode
<li>The temperature is available when the following parameter is set to true in the data file:
\code
===Is there a temperature field?
.t.
\endcode
</ol>

<h3>Exemple</h3>
Here is an exemple where the source term is equal to \f$-\rho \textbf{e}_z\f$ with  \f$\rho\f$ being the density and \f$\textbf{e}_z\f$ being the unit vector in the vertical direction.

\code
    IF (PRESENT(opt_density)) THEN
       IF (TYPE==5) THEN
          vv = -opt_density(:,1,i)
       ELSE IF (TYPE==6.AND.mode>0) THEN
          vv = -opt_density(:,2,i)
       ELSE
          vv = 0.d0
       END IF
    ELSE
       CALL error_petsc('problem in source_in_NS_momentum: set level set to true in data to have a density')
    END IF
    RETURN
\endcode

Remarks:
<ol>
<li>We check that the variable opt_density is present. If not, it means the data file is not correctly set: it does not use a level set and a variable density. As a consequence, the code is stopped with the function <tt>error_petsc</tt>.
<li>The source term only depends of the vertical direction (TYPE 5 and 6).
<li>When we consider the vertical cosine part of the source term, it involve the cosine part of the density. Same goes when considering sine part.
<li>The vertical sine part of the source term is set to zero when the Fourier mode considered is zero. It is due to the fact that \f$\sin(0\theta)=0\f$ for all \f$\theta\in[0,2\pi]\f$.
</ol>


We refer to the sections \ref doc_debug_test and
\ref doc_example_physical_pb for more examples.

 */
