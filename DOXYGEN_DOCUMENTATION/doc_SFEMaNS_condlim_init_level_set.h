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
 * @page doc_SFEMaNS_condlim_init_level_set  The subroutine init_level_set

It is used to initialize the level set involved for multiphase problem with variable density or dynamical viscosity. 

<h3>Inputs and outputs</h3>
The inputs of this function are the following:
<ol>
<li><tt>pp_mesh</tt> is the finite element mesh used to approximate the level set.
<li><tt>dt</tt> is the time step.
<li><tt>list_mode</tt> is a list of integers which contains the Fourier modes approximated.
</ol>
As the mesh can be subdivised, we note that  pp_mesh depends of the processor considered when doing parallel computing. Same goes for the list of Fourier mode <tt>list_mode</tt>.

The outputs of this function are the following:
<ol>
<li><tt>time</tt> is the time when the computations starts. 
<li><tt>level_set_m1</tt> is the level set at the time <tt>time-dt</tt>.
<li><tt>level_set</tt> is the level set at the time <tt>time</tt>.
</ol>
We note that the level set is a list of real valued scalar functions (one scalar per interface between fluids). So its format is a real valued tabular of four columns of dimension (inputs\%nb_fluid -1,pp_mesh\%np,2,SIZE(list_mode)). We remind that pp_mesh\%np is the number of nodes of the finite element mesh. The integer inputs\%nb_fluid is the number of immiscible fluids (only used for multiphase problem).

<h3>Exemple</h3>
Here is an exemple where we use the function level_set_exact to initialize the level set.
\code
    time = 0.d0
    DO i= 1, SIZE(list_mode)
       mode = list_mode(i) 
       DO j = 1, 2
          !===level_set
          DO n = 1, inputs%nb_fluid -1
             level_set_m1(n,:,j,i) = level_set_exact(n,j,vv_mesh%rr,mode,time-dt)  
             level_set   (n,:,j,i) = level_set_exact(n,j,vv_mesh%rr,mode,time)
          END DO
       END DO
    END DO
    RETURN
\endcode
We note that the integer mode, i, j and n need to be declared. It is done by adding the  following line in the declaration of the function <tt>init_temperature</tt>:
\code
    INTEGER                                    :: mode, i, j 
\endcode

We refer to the sections \ref doc_debug_test and
\ref doc_example_physical_pb for more examples.

 */
