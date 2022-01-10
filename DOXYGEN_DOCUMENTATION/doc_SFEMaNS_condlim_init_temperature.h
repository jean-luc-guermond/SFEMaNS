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
 * @page doc_SFEMaNS_condlim_init_temperature The subroutine init_temperature


It is used to initialize the temperature. 

<h3>Inputs and outputs</h3>
The inputs of this function are the following:
<ol>
<li><tt>mesh</tt> is the finite element mesh used to approximate the temperature.
<li><tt>dt</tt> is the time step.
<li><tt>list_mode</tt> is a list of integers which contains the Fourier modes approximated.
</ol>
As the mesh can be subdivised, we note that the mesh depends of the processor considered when doing parallel computing. Same goes for the list of Fourier mode <tt>list_mode</tt>.

The outputs of this function are the following:
<ol>
<li><tt>time</tt> is the time when the computations starts. 
<li><tt>tempn_m1</tt> is the temperature at the time <tt>time-dt</tt>.
<li><tt>tempn</tt> is the temperature at the time <tt>time</tt>.
</ol>
We note that the temperature is a real valued scalar so its format is real value a tabular of three columns of dimension (mesh\%np,2,SIZE(list_mode)). We remind that mesh\%np is the number of nodes of the finite element mesh.

<h3>Exemple</h3>
Here is an exemple where we use the function temperature_exact (so initial conditions satisfy the boundary conditions).
\code
    time = 0.d0
    DO i= 1, SIZE(list_mode)
       mode = list_mode(i) 
       DO j = 1, 2 
          tempn_m1(:,j,i) = temperature_exact(j, mesh%rr, mode, time-dt)
          tempn   (:,j,i) = temperature_exact(j, mesh%rr, mode, time)
       ENDDO
    ENDDO
    RETURN
\endcode
We note that the integer mode, i and j need to be declared. It is done by adding the  following line in the declaration of the function <tt>init_temperature</tt>:
\code
    INTEGER                                    :: mode, i, j 
\endcode


We refer to the sections \ref doc_debug_test and
\ref doc_example_physical_pb for more examples.



 */
