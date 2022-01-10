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
 * @page doc_SFEMaNS_condlim_mu_in_real_space The function mu_in_real_space


It is used to define a magnetic permeability depending of the time and all
 of the space direction \f$(r,\theta, z)\f$. If the permeability magnetic
 does not depend of the time and the azimuthal direction, this function
 is not used. We refer to the section
 \ref doc_intro_SFEMaNS_possibilities_mxw_2 
 for more information on the formulation of the Maxwell equations in SFEMaNS. 

This function defines a scalar function in the physical space for a number
 of angles, denoted nb_angles, on a number of nodes of the finite element mesh.

<h3>Inputs and outputs</h3>
The inputs of this function are the following:
<ol>
<li><tt>H_mesh</tt> is the mesh where the magnetic field is approximated.
<li><tt>angles</tt> is the list of the angles where the penalty function is computed. These reals numbers are in the interval \f$[0,2\pi[\f$.
<li><tt>nb_angles</tt> is the number of angles considered. It is an interger.
<li><tt>nb</tt> is the label of the first node considered. It is an integer.
<li><tt>ne</tt> is the label of the last node considered. It is an integer.
<li><tt>time</tt> is the time at which this term is computed. It is a real number.
</ol>

The output of this function is a real valued tabular vv with two columns of dimension (nb_angles,ne-nb+1).


<h3>Exemple</h3>
Here is an exemple where the following magnetic permeability is considered: \f$\mu(r,\theta,z,t) = 1 + r^2 + |z|(1+\cos(t-\theta))\f$.

The corresponding code lines are written as follows.
\code
    DO i = 1, nb_angles
       DO n = 1, ne-nb+1
          r=H_mesh%rr(1,nb+n)
          theta=angles(i)
          z=H_mesh%rr(2,nb+n)
          vv(i,n)=1.d0 + r**2 + ABS(z)*(1.d0+COS(time-theta))
       END DO
    END DO
    RETURN
\endcode

Remark: If the problem involves a magnetic perbeability with discontinuities,
 it is necessary to use a continuous magnetic permeabilty with a sharp
 gradient inthe regions where the discontinuities occur.


We refer to the sections \ref doc_debug_test (see test 22, 23 and 29) and
\ref doc_example_physical_pb for more examples.


 */
