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
 * @page doc_SFEMaNS_condlim_H_B_quasi_static The function H_B_quasi_static

It is used to define the basic states \f$\bB_b\f$ and \f$\bH_b\f$. 

This function defines a time independent magnetic field for one given Fourier mode on all the nodes of the finite element mesh.  We denote by H_mesh the finite element mesh used to approximate the magnetic field.

<h3>Inputs and outputs</h3>
The inputs of this function are the following:
<ol>
<li><tt>char_h_b</tt> is a character of one letters. It is equal to B or H and allows to know which magnetic field is computed.
<li><tt>rr</tt> is a real valued tabular that contains two columns with dimensions (2,H_mesh\%np). The tabular rr(1,:) contains the radial cylindrical coordinate of all the nodes of the finite element mesh H_mesh. Respectively, rr(2,:) contains the vertical coordinates of these nodes.
<li><tt>m</tt> is the Fourier mode \f$m\f$ considered. It is an integer.
</ol>

The output of this function is a real valued tabular vv of dimension (SIZE(rr,2),6).

<h3>Exemple</h3>
Here is an exemple where we set \f$\bH_\text{b}= r\textbf{e}_\theta\f$ with \f$\textbf{e}_\theta\f$ the unit vector in the azimuthal direction. Moreover, we assume that the Maxwell equations are approximated on one domain. The magnetic permeability of the domain is avaible via the data inputs\%mu_H(1).

The corresponding code lines are written as follows.
\code
    IF (char_h_b=='B') THEN
       IF (TYPE==3.AND.m==0) THEN
          vv = inputs\%mu_H(1)*rr(1,:)
       ELSE
          vv = 0.d0
       END IF
    ELSE
       IF (TYPE==3.AND.m==0) THEN
          vv = rr(1,:)
       ELSE
          vv = 0.d0
       END IF
    END IF
    RETURN
\endcode


We refer to the section \ref doc_example_physical_pb for more examples.

 */
