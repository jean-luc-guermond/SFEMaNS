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
 * @page doc_SFEMaNS_condlim_sigma_bar_in_fourier_space The function sigma_bar_in_fourier_space

It is used to define a stabilization term \f$\overline{\sigma}(r,z)\f$
 for multiphase problems where the electrical conductivity varies
 between the fluids considered.  We refer to the section
 \ref doc_intro_SFEMaNS_possibilities_nst_4
 for more information on the formulation of the Navier-Stokes
 equations in SFEMaNS. 

This function defines a scalar function depending only of the radial and vertical coordinates. It is defined on the nodes of the finite element mesh.

<h3>Inputs and outputs</h3>
The input of this function the mesh <tt>H_mesh</tt> where the magnetic field is approximated.

The output of this function is a real valued tabular vv of dimension SIZE(H_mesh%rr,2). It is equal to the number of nodes, H_mesh\%np, of the finite element mesh used to approximate the magnetic field.

Remark:
<ol>
<li>The electrical conducitivites of the fluids are set in the data as follows.
\code
===Conductivity of fluid 0, fluid 1, ...
1.d0 2.d0
\endcode
<li>These electrical conductivities are stored in the variable inputs\%sigma_fluid. It is a real valued tabular of dimension the number of fluids considered.
</ol>

<h3>Exemple</h3>
Here is an exemple where we set the stabilization term \f$\overline{\sigma}\f$ to half of the  minimum of the fluids electrical conductivities.

The corresponding code lines are written as follows.
\code
    vv=0.5d0*MINVAL(inputs%sigma_fluid)
    RETURN
\endcode

We refer to the sections \ref doc_debug_test and
\ref doc_example_physical_pb for more examples.

 */
