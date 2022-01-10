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
 * @page doc_SFEMaNS_condlim_source_in_temperature The function source_in_temperature

It is used to define the source term \f$f_T\f$ of the temperature equation. We refer to the section \ref doc_intro_SFEMaNS_possibilities_temp for more information on the formulation of the temperature equation in SFEMaNS. 

This function defines the source term for one given Fourier mode, one given component (cosine or sine) on all the gauss points of the finite element mesh. We denote by temp_mesh the finite element mesh used to approximate the temperature.

<h3>Inputs and outputs</h3>
The inputs of this function are the following:
<ol>
<li><tt>TYPE</tt> is the component of the source term that is computed (cosine or sine). It is an integer between one and two.
<li><tt>rr</tt> is a real valued tabular that contains two columns with dimensions (2,temp_mesh\%gauss\%l_G*temp_mesh\%me). We remind that me is the number of cells of the mesh and that l_G is the number of gauss points per cell. The tabular rr(1,:) contains the radial cylindrical coordinate of each gauss points of temp_mesh. Respectively, rr(2,:) contains the vertical coordinates of these gauss points.
<li><tt>m</tt> is the Fourier mode \f$m\f$ considered. It is an integer.
<li><tt>t</tt> is the time at which the source term is computed. It is a real number.
</ol>


The output of this function is a real valued tabular vv of dimension SIZE(rr,2) which is equal to the number of gauss points in temp_mesh.



<h3>Exemple</h3>
Here is an exemple where the source term \f$f_T\f$ is set to \f$ z + 2r\cos(\theta) + rz \sin(2\theta)\f$.
\code
    IF (TYPE==1.AND.m==0) THEN
       vv = rr(2,:)
    ELSE IF (TYPE==1.AND.m==1) THEN
       vv = 2.d0*rr(1,:)
    ELSE IF  (TYPE==2.AND.m==2) THEN
       vv = r(1,:)*rr(2,:)
    ELSE
       vv = 0.d0
    END IF
    RETURN
\endcode

We refer to the sections \ref doc_debug_test and
\ref doc_example_physical_pb for more examples.

 */
