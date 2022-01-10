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
 * @page doc_SFEMaNS_condlim_source_in_level_set The function source_in_level_set

It is used to define the source term \f$f_\varphi\f$ of the level
 set equations. We refer to the section
\ref doc_intro_SFEMaNS_possibilities_nst_4
 for more information on the formulation of the Navier-Stokes
 equations in SFEMaNS for multiphase problem. We remind there
 is one source term per level set.

This function defines the source term for one given interface, one given Fourier mode, one given component (cosine or sine) on all the nodes of the finite element mesh. We denote by level_set_mesh the finite element mesh used to approximate the level set.

<h3>Inputs and outputs</h3>
The inputs of this function are the following:
<ol>
<li><tt>TYPE</tt> is the component of the source term that is computed (cosine or sine). It is an integer between one and two.
<li><tt>rr</tt> is a real valued tabular that contains two columns with dimensions (2,level-set_mesh\%np). The tabular rr(1,:) contains the radial cylindrical coordinate of all the nodes of the finite element mesh level_set_mesh. Respectively, rr(2,:) contains the vertical coordinates of these nodes.
<li><tt>m</tt> is the Fourier mode \f$m\f$ considered. It is an integer.
<li><tt>interface_nb</tt> is the number of the interface considered. It is an integer.
<li><tt>t</tt> is the time at which the source term is computed. It is a real number.
</ol>

The output of this function is a real valued tabular vv of dimension SIZE(rr,2) which is equal to the number of nodes in level_set_mesh.

Remark:
<ol>
<li>For physical applications, these source terms should be set to zero. Indeed the transport of the interface by the velocity is an advection equation. So it does not involve source term. 
<li>This function has been added to consider manufactured solutions that are not solution of the advection equation.
</ol>


<h3>Exemple</h3>
Here is an exemple where the source term is set to zero.
\code
    vv = 0.d0
    RETURN
\endcode

We refer to the sections \ref doc_debug_test and
\ref doc_example_physical_pb for more examples.

 */
