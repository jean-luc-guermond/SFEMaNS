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
 * @page doc_SFEMaNS_condlim_Jexact_gauss The function Jexact_gauss


It is used to define the source term \f$\textbf{j}\f$. It represents the action of an externally imposed distribution of current. It can also be used when approximated a perturbation of a basic state solution of the Maxwell equations.

This function defines the source term for one given Fourier mode, one given component (radial cosine, radial sine, azimuthal cosine, azimuthal sine, vertical cosine or vertical sine) on given node or gauss point of the finite element mesh.

<h3>Inputs and outputs</h3>
The inputs of this function are the following:
<ol>
<li><tt>TYPE</tt> is the component of the source term that is computed (radial cosine, radial sine, etc.). It is an integer between one and six.
<li><tt>rr</tt> is a list of two real numbers. The tabular rr(1) contains the radial cylindrical coordinate of the node or the gauss point considered. Respectively, rr(2) contains the vertical coordinates of this node or gauss point.
<li><tt>m</tt> is the Fourier mode \f$m\f$ considered. It is an integer.
<li><tt>mu_phi</tt> is the magnetic perbeability in the insulating region.
<li><tt>sigma</tt> is the magnetic Reynolds number multiplied by the electrical conductivity of the domain that contains the node/gauss point considered.
<li><tt>mu_H_field</tt> is the magnetic permeability of the domain that contains the node/gauss point considered.
<li><tt>t</tt> is the time at which this term is computed. It is a real number.
<li><tt>mesh_id</tt> is the label of the domain where this source term is computed. It is an interger.
<li><tt>opt_B_ext</tt> is the component of the Fourier mode <tt>m</tt> of the magnetic field extrapolated at time <tt>t</tt>. This input optional.
</ol>

The output of this function is a real number vv.


<h3>Exemple</h3>
Here is an exemple where we approximate the pertubation of the basic state \f$\bH_b=r\textbf{e}_\theta\f$ of the Maxwell equations. We denote \f$\textbf{e}_\theta\f$ the unit vector in the azimuthal direction and  \f$\textbf{e}_z\f$ the unit vector in the vertical direction. As a consequence, the source term has to satisfy the relation \f$\textbf{j}=-\ROT \bH_b\f$. So we  set \f$\textbf{j}= -2 \textbf{e}_z\f$.

The corresponding code lines are written as follows.
\code
    IF (TYPE==5.AND.m==0) THEN
       vv = -2.d0
    ELSE
       vv = 0.d0
    END IF
    RETURN
\endcode

We refer to the sections \ref doc_debug_test and
\ref doc_example_physical_pb for more examples.


 */
