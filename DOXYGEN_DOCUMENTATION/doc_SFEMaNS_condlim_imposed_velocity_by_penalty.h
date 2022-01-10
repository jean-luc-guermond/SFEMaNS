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
 * @page doc_SFEMaNS_condlim_imposed_velocity_by_penalty The function imposed_velocity_by_penalty




It is used to define the velocity in the solid obstacle denoted
 \f$\bu_\text{obs}\f$. We note this velocity does not depends of
 the cylindrical coordinate \f$\theta\f$. We refer to the section
 \ref doc_intro_SFEMaNS_possibilities_nst_3
 for more information on the formulation of the Navier-Stokes
 equations in SFEMaNS. 

This function defines a velocity for the Fourier mode zero, for all components (radial cosine, radial sine, azimuthal cosine, azimuthal sine, vertical cosine or vertical sine) on a number of nodes, denoted nb_node, of the finite element mesh. We denote by vv_mesh the finite element mesh used to approximate the velocity field.

<h3>Inputs and outputs</h3>
The inputs of this function are the following:
<ol>
<li><tt>rr</tt> is a real valued tabular that contains two columns with dimensions (2,nb_node). The tabular rr(1,:) contains the radial cylindrical coordinate of each nodes considered. Respectively, rr(2,:) contains the vertical coordinates of these nodes. We note that the integer nb_node is generally equal to the total number of node vv_mesh\%np.
<li><tt>t</tt> is the time at which this term is computed. It is a real number.
</ol>

The output of this function is a real valued tabular vv with two columns of dimension (SIZE(rr,2),6).


<h3>Exemple</h3>
Here is an exemple where we impose a solid rotation. It means that we set \f$\bu=r\textbf{e}_\theta\f$ on the boundary of the domain where \f$\textbf{e}_\theta\f$ is the unit vector in the azimuthal direction.

The corresponding code lines are written as follows.
\code
    vv(:,1)   = rr(1,:)
    vv(:,2:6) = 0.d0
    RETURN
\endcode


We refer to the sections \ref doc_debug_test (see test 24) and
\ref doc_example_physical_pb for more examples.


 */
