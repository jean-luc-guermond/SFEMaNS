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
 * @page doc_SFEMaNS_condlim_kelvin_force_coeff The function kelvin_force_coeff

It is used to define the coefficient \f$g(T)\f$ of the Kelvin force
 for ferrofluid problems. This coefficient only depends of the temperature.
 We refer to the section \ref doc_intro_num_app
 for more information on the formulation of the Navier-Stokes
 equations in SFEMaNS for ferrofluid problem. 

This function defines a scalar coefficient for a given temperature
 defined in the physical space.

<h3>Inputs and outputs</h3>
The inputs of this function is a real denoted temp. It represents the value of the temperature at a given point of the 3D Navier-Stokes domain (radial, azimuthal and vertical cylindrical coordinates are fixed).

The output of this function is a real vv.

Remark: 
<ol>
<li>The input temperature is given in the physical space for a coordinate \f$(r,\theta,z)\f$. It means you can use product, division, etc.
<li>The output vv is also defined in the physical space and not in the Fourier space. Meaning it is defined for the coordinates \f$(r,\theta,z)\f$ where the given temperature is defined.
</ol>

<h3>Exemple</h3>
Here is an exemple where the following coeffcient for the Kelvin force is considered: \f$g(T)=1+ |T| + 3 T^2\f$.

The corresponding code lines are written as follows.
\code
    vv = 1.d0 + ABS(temp) + 3.d0*temp**2
    RETURN
\endcode

We refer to the sections \ref doc_debug_test (see test 33) and
\ref doc_example_physical_pb for more examples.


 */
