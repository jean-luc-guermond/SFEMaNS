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
 * @page doc_SFEMaNS_condlim_grad_mu_bar_in_fourier_space The function grad_mu_bar_in_fourier_space


It is used to define the gradient in the radial and vertical direction
 of the function <tt>mu_bar_in_fourier_space</tt>. We remind that
 <tt>mu_bar_in_fourier_space</tt> is  either a magnetic permeability
 \f$\mu(r,z)\f$ or a stabilization term \f$\overline{\mu}(r,z)\f$
 when the magnetic permeability depends of the time or the azimuthal
 direction \f$\theta\f$.  We refer to the section
 \ref doc_intro_SFEMaNS_possibilities_mxw_2
 for more information on the formulation of the Maxwell equations in SFEMaNS. 

This function defines a gradient in the radial and vertical direction
 for one given gauss point of the finite element mesh.

<h3>Inputs and outputs</h3>
The inputs of this function are the following:
<ol>
<li><tt>pts</tt> is real valued tabular of dimension two. pts(1) is the radial cylindrical coordinate of the gauss point considered. pts(2) is the vertical cylindrical coordinate of the gauss point considered.
<li><tt>pts_id</tt> is the label of the domain that contains the gauss point considered. 
</ol>

The output of this function is a real valued tabular vv of dimension two.

<h3>Exemple</h3>
Here is an exemple where we consider the magnetic permeability \f$\mu(r,z)=1 + r z^2\f$.

The corresponding code lines are written as follows.
\code
    vv(1)=pts(2)**2
    vv(1)=2.d0*pts(1)*pts(2)
    RETURN
\endcode


We refer to the sections \ref doc_debug_test (see tests 17, 18, 22, 23
 or 29) and \ref doc_example_physical_pb for more examples.


 */
