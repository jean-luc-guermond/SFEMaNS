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
 * @page doc_SFEMaNS_condlim_pp_exact The function pp_exact

It can be used by the function <tt>init_velcoity_pressure</tt> to initialize the pressure.

This function defines a pressure for one given Fourier mode, one given component (cosine or sine) on a number of nodes, denoted nb_node, of the finite element mesh. We denote by pp_mesh the finite element mesh used to approximate the velocity field.

<h3>Inputs and outputs</h3>
The inputs of this function are the following:
<ol>
<li><tt>TYPE</tt> is the component of the source term that is computed (cosine or sine). It is an integer between one and two.
<li><tt>rr</tt> is a real valued tabular that contains two columns with dimensions (2,nb_node). The tabular rr(1,:) contains the radial cylindrical coordinate of each nodes considered. Respectively, rr(2,:) contains the vertical coordinates of these nodes. We note that the integer nb_node is generally equal to the total number of node pp_mesh\%np.
<li><tt>m</tt> is the Fourier mode considered. It is an integer.
<li><tt>t</tt> is the time at which thei initialization is computed. It is a real number.
</ol>

The output of this function is a real valued tabular vv of dimension SIZE(rr,2).


<h3>Exemple</h3>
Here is an exemple where we use this function to initalize a pressure to \f$p(r,\theta,z,t=0)=r^2 + 3z\sin(\theta)\f$ 

The corresponding code lines are written as follows.
\code
    IF (TYPE==1.AND.m==0) THEN
       vv = rr(1,:)**2
    ELSE IF (TYPE==2.AND.m==1) THEN
       vv = 3.d0*rr(2,:)
    ELSE
       vv = 0.d0
    END IF
    RETURN
\endcode


We refer to the sections \ref doc_debug_test and
\ref doc_example_physical_pb for more examples.



 */
