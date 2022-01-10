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
 * @page doc_SFEMaNS_condlim_Hexact The function Hexact

It is used to define the boundary condition \f$\bH_\text{bdy}\f$. It can also be used by the function <tt>init_maxwell</tt> so initial conditions match the boundary conditions.

This function defines boundary conditions for one given Fourier mode, one given component (radial cosine, radial sine, azimuthal cosine, azimuthal sine, vertical cosine or vertical sine) on a number of nodes, denoted nb_node, of the finite element mesh.  We denote by H_mesh the finite element mesh used to approximate the magnetic field.

<h3>Inputs and outputs</h3>
The inputs of this function are the following:
<ol>
<li><tt>H_mesh</tt> is the finite element mesh used to approximate the magnetic field.
<li><tt>TYPE</tt> is the component of the source term that is computed (radial cosine, radial sine, etc.). It is an integer between one and six.
<li><tt>rr</tt> is a real valued tabular that contains two columns with dimensions (2,nb_node). The tabular rr(1,:) contains the radial cylindrical coordinate of each nodes considered. Respectively, rr(2,:) contains the vertical coordinates of these nodes. We note that the integer nb_node is generally equal to the total number of node H_mesh\%np or the number of nodes where Dirichlet boundary conditions are applied.
<li><tt>m</tt> is the Fourier mode \f$m\f$ considered. It is an integer.
<li><tt>mu_H_field</tt> is a list of real that contains the magnetic permeability of each conducting domains.
<li><tt>t</tt> is the time at which this term is computed. It is a real number.
</ol>

The output of this function is a real valued tabular vv of dimension SIZE(rr,2).


<h3>Exemple</h3>
Here is an exemple where we set \f$\bH_\text{bdy}= z\textbf{e}_r + 3r\sin(2\theta)\textbf{e}_z\f$ with \f$\textbf{e}_r\f$, respectively \f$\textbf{e}_z\f$, the unit vector in the radial direction, respectively in the vertical direction.

The corresponding code lines are written as follows.
\code
    IF (TYPE==1.AND.m==0) THEN
       vv = rr(2,:)
    ELSE IF (TYPE==6.AND.m==2) THEN
       vv = 3.d0*rr(1,:)
    ELSE
       vv = 0.d0
    END IF
    RETURN
\endcode

Remark: The magnetic field \f$\mu\textbf{H}\f$ has a zero divergence so the divergence of the boundary data \f$\bH_\text{bdy}\f$ has to be zero.

We refer to the sections \ref doc_debug_test and
\ref doc_example_physical_pb for more examples.

 */
