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
 * @page doc_SFEMaNS_condlim_penal_in_real_space The function penal_in_real_space

It is used to define the penalty function \f$\chi\f$. We remind this function
 is a scalar function equal to 0 in the solid obstacle and 1 in the fluid
 region. We refer to the section 
 \ref doc_intro_SFEMaNS_possibilities_nst_3
 for more information on the formulation of the Navier-Stokes
 equations in SFEMaNS. 

This function defines a scalar function in the physical space for a number of angles, denoted nb_angles, on a number of nodes of the finite element mesh. We denote by mesh the finite element mesh used to approximate the velocity field.

<h3>Inputs and outputs</h3>
The inputs of this function are the following:
<ol>
<li><tt>mesh</tt> is the mesh where the velocity field is approximated.
<li><tt>angles</tt> is the list of the angles where the penalty function is computed. These reals numbers are in the interval \f$[0,2\pi[\f$.
<li><tt>nb_angles</tt> is the number of angles considered. It is an interger.
<li><tt>nb</tt> is the label of the first node considered. It is an integer.
<li><tt>ne</tt> is the label of the last node considered. It is an integer.
<li><tt>rr_gauss</tt> is a real valued tabular that contains two columns with dimensions (2,ne-nb+1). The tabular rr(1,:) contains the radial cylindrical coordinate of each nodes considered. Respectively, rr(2,:) contains the vertical coordinates of these nodes. We note that we usually use ne=1 and nb=mesh\%np.
<li><tt>time</tt> is the time at which this term is computed. It is a real number.
</ol>

The output of this function is a real valued tabular vv with two columns of dimension (nb_angles,ne-nb+1).


<h3>Exemple</h3>
Here is an exemple where the solid obsctable is the domain \f$\Omega_\text{obs} = \{(r,\theta,z)\in [0,1]\times[0,\pi]\times[-0.1,0.1] \}\f$.

The corresponding code lines are written as follows.
\code
    DO i = 1, nb_angles
       DO n = 1, SIZE(rr,2)
          r=rr(1,n)
          theta=angles(i)
          z=rr(2,n)
          IF ((r.LE.1).AND.(theta.LE.ACOS(-1.d0)).AND.(ABS(z).LE.0.1d0)) THEN
             vv(i,n)=0.d0
          ELSE
             vv(i,n)=1.d0
          END IF
       END DO
    END DO
    RETURN
\endcode
We note that it is better to use a smooth penalty function and
 not a discontinuous one as above.

We refer to the sections \ref doc_debug_test  (see test 24) and
\ref doc_example_physical_pb for more examples.


 */
