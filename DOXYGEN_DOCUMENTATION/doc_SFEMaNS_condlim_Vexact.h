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
 * @page doc_SFEMaNS_condlim_Vexact The function Vexact



It is used to define the velocity field \f$\bu_\text{given}\f$ when the Navier-Stokes equations are not approximated. Moreover, in the data file the following line has to be set to false:
\code
===Restart on velocity (true/false)
\endcode
If it set to true, the velocity is read from a suite file of the Navier-Stokes equations and not this function.


This function defines a time independent velocity field for one given Fourier mode on all of the nodes of the finite element mesh.

<h3>Inputs and outputs</h3>
The inputs of this function are the following:
<ol>
<li><tt>m</tt> is the Fourier mode \f$m\f$ considered. It is an integer.
<li><tt>H_mesh</tt> is the finite element mesh used to approximate the magnetic field.

</ol>

The output of this function is a real valued tabular vv of dimension (H_mesh\%np,6).


<h3>Exemple</h3>
Here is an exemple where we set \f$\bu_\text{given}=r\textbf{e}_\theta\f$  where \f$\textbf{e}_\theta\f$ is the unit vector in the azimuthal direction.


The corresponding code lines are written as follows.
\code
    IF (TYPE==3.AND.m==0) THEN
       vv = rr(1,:)
    ELSE
       vv = 0.d0
    END IF
    RETURN
\endcode

Remark: The divergence of the velocity field has to be zero.

We refer to the sections \ref doc_debug_test and
\ref doc_example_physical_pb for more examples.

 */

