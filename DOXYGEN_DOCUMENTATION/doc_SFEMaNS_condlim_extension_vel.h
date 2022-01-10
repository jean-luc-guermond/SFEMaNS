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
 * @page doc_SFEMaNS_condlim_extension_vel  The function extension_velocity



It is used to extend the velocity field approximated with the Navier-Stokes equations to a velocity field on a larger domain. It is only used when the temperature or the magnetic field equations are approximated in addition of the Navier-Stokes equations.

This function defines a velocity field for one given Fourier mode, one given component (radial cosine, radial sine, azimuthal cosine, azimuthal sine, vertical cosine or vertical sine) on all the nodes of the finite element mesh considered.

<h3>Inputs and outputs</h3>
The inputs of this function are the following:
<ol>
<li><tt>TYPE</tt> is the component of the source term that is computed (radial cosine, radial sine, etc.) It is an integer between one and six.
<li><tt>H_mesh</tt> is the mesh where the velocity field is defined. It has H_mesh\%np nodes and the radial and vertical cylindrical coordinates are in the tabular H_mesh\%rr.
<li><tt>mode</tt> is the Fourier mode considered. It is an integer.
<li><tt>t</tt> is the time at which this term is computed. It is a real number.
</ol>

The output of this function is a real valued tabular vv of dimension H_mesh\%np.

Remark: The extension_velocity defines a velocity on the whole temperature or magnetic field domain. However only its value outside the Navier-Stokes domain are used.


<h3>Exemple</h3>
Here is an exemple where we extend the velocity field of Navier-Stokes domain to a region where the velocity field is \f$\bu=(r-z)\cos(t)\textbf{e}_\theta\f$. We by \f$\textbf{e}_\theta\f$ the unit vector in the azimuthal direction. 

The corresponding code lines are written as follows.
\code
    IF (TYPE==3.AND.m==0)
       vv = (H_mesh%rr(1,:)-H_mesh%rr(2,:))COS(t)
    ELSE
       vv = 0.d0
    END IF
    RETURN
\endcode

Remark:
<ol>
<li>The divergence of the velocity has to be zero.
<li>The velocity on the interface between the Navier-Stokes domains and the othert domains needs to match. It means that the outputs of the function <tt>vv_exact</tt> and <tt>extension_velocity</tt> have to match on the interfaces.
</ol>


We refer to the sections \ref doc_debug_test and
\ref doc_example_physical_pb for more examples.


 */
