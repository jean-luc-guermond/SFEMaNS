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
 * @page doc_SFEMaNS_condlim_mu_bar_in_fourier_space The function mu_bar_in_fourier_space

It is used to define either a magnetic permeability \f$\mu(r,z)\f$
 or a stabilization term \f$\overline{\mu}(r,z)\f$ when the magnetic
 permeability depends of the time or the azimuthal direction \f$\theta\f$.
  We refer to the section 
 \ref  doc_intro_SFEMaNS_possibilities_mxw_2
 for more information on the formulation of the Maxwell equations in SFEMaNS. 

This function defines a scalar function depending only of the radial and vertical coordinates. The function is either defined on a number of nodes or one gauss points of the finite element mesh.

<h3>Inputs and outputs</h3>
The inputs of this function are the following:
<ol>
<li><tt>H_mesh</tt> is the mesh where the magnetic field is approximated.
<li><tt>nb</tt> is the label of the first node or gauss point considered. It is an integer.
<li><tt>ne</tt> is the label of the last node or gauss point considered. It is an integer.
<li><tt>pts</tt> is a real valued tabular that contains two columns with dimensions (2,ne-nb+1). The tabular pts(1,:) contains the radial cylindrical coordinate of each nodes considered. Respectively, pts(2,:) contains the vertical coordinates of these nodes.
<li><tt>pts_id</tt> is the label of the domain that contains the gauss points considered. This input is optional, it is only used when doing the computation on gauss points.
</ol>

The output of this function is a real valued tabular vv of dimension ne-nb+1.

Remark: The input <tt>pts_id</tt> is only used when the magnetic permeability presents a discontinuity. Such discontinuities occur on an interface between two domains where the magnetic field is approximated. The presence of a jump is set in the data file with the following lines:
\code
===Number of interfaces in H mesh
??
===List of interfaces in H mesh
??
\endcode
The jump is treated with a penalty method.

<h3>Exemple</h3>
Here is an exemple where we consider the magnetic permeability \f$\mu(r,z)=1 + r z^2\f$.

The corresponding code lines are written as follows.
\code
    IF (PRESENT(pts)) THEN
       DO  n = 1, SIZE(rr,2)
          vv(n)=1.d0+pts(1,n)*pts(2,n)**2
       END DO
    ELSE
       DO  n = 1, SIZE(rr,2)
          r=H_mesh%rr(1,nb+n-1)
          z=H_mesh%rr(2,nb+n-1)
          vv(n)=1.d0 + r*z**2
       END DO
    END IF
    RETURN
\endcode


We refer to the sections \ref doc_debug_test (see tests 17, 18, 22, 23 or 29)
 and \ref doc_example_physical_pb for more examples.


 */
