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
 * @page doc_SFEMaNS_condlim_level_set_exact  The function level_set_exact

It can be used by the function <tt>init_level_set</tt> to initialize the level set functions.

This functions defines a level set for one given interface, one given Fourier mode,  one given component (cosine or sine) on a number of nodes, denoted nb_node, of the finite element mesh. We denote by level_set_mesh the finite element mesh used to approximate the level set.


<h3>Inputs and outputs</h3>
The inputs of this function are the following:
<ol>
<li><tt>interface_nb</tt> is the label of the interface considered. It is an integer between one and the number of fluids minus one.
<li><tt>TYPE</tt> is the component of the level set is computed (cosine or sine). It is an integer between one and two.
<li><tt>rr</tt> is a real valued tabular that contains two columns with dimensions (2,nb_node). The tabular rr(1,:) contains the radial cylindrical coordinate of each nodes considered. Respectively, rr(2,:) contains the vertical coordinates of these nodes. We note that the integer nb_node is generally equal to the total number of node level_set_mesh\%np.
<li><tt>m</tt> is the Fourier mode \f$m\f$ considered. It is an integer.
<li><tt>t</tt> is the time at which this term is computed. It is a real number.
</ol>

The output of this function is a real valued tabular vv of dimension SIZE(rr,2).

<h3>Exemple</h3>
Here is an exemple where we consider a problem with one level set (so two fluids). The initial interface is \f$z=1\f$. We define a level set \f$\varphi\f$ that depends of the vertical coordinates as follows:
\f{align}
\varphi(z)= \frac{1}{2}(1+\tanh(\frac{z-1}{0.05})).
 \f}
One can note that \f$\varphi\f$ is close to zero for z smaller than 0.9, close to one for z larger than 1.1 and equal to 0.5 for z equal to 1.

The corresponding code lines are written as follows.
\code
    IF (interface_nb==1) THEN
       IF (TYPE==1.AND.m==0) THEN
          vv = 0.5d0*(1.d0+TANH((z-1.d0)/0.05d0))
       ELSE
          vv = 0.d0
       END IF
    ELSE
       CALL error_petsc('problem in level_set_exact: there is more than one level set')
    END IF
    RETURN
\endcode


We refer to the sections \ref doc_debug_test and
\ref doc_example_physical_pb for more examples.


 */
