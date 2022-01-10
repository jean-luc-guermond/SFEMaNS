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
 * @page doc_SFEMaNS_condlim_init_velocity_pressure The subroutine init_velocity_pressure

It is used to initialize the velocity field, the pressure and the pressure increment. 

<h3>Inputs and outputs</h3>
The inputs of this function are the following:
<ol>
<li><tt>mesh_f</tt> is the finite element mesh used to approximate the velocity field.
<li><tt>mesh_c</tt> is the finite element mesh used to approximate the pressure.
<li><tt>dt</tt> is the time step.
<li><tt>list_mode</tt> is a list of integers which contains the Fourier modes approximated.
</ol>
As the mesh can be subdivised, we note that  mesh_f and mesh_c depend of the processor considered when doing parallel computing. Same goes for the list of Fourier mode <tt>list_mode</tt>.

The outputs of this function are the following:
<ol>
<li><tt>time</tt> is the time when the computations starts. 
<li><tt>un_m1</tt> is the velocity field at the time <tt>time-dt</tt>.
<li><tt>un</tt> is the velocity field at the time <tt>time</tt>.
<li><tt>pn_m1</tt> is the pressure at the time <tt>time-dt</tt>.
<li><tt>pn</tt> is the pressure at the time <tt>time</tt>.
<li><tt>phin_m1</tt> is the increment pressure at the time <tt>time-dt</tt>.
<li><tt>phin</tt> is the increment pressure at the time <tt>time</tt>.
</ol>

Remarks:
<ol>
<li>For a velocity field with a zero divergence, the pressure increments satisfy the relations: phin=pn-pn_m1 and phin_m1=pn_m1-pn_m2.
<li>The velocity field is a vector so its format is a real valued tabular of three columns with dimension (mesh_f\%np,6,SIZE(list_mode)). We remind that mesh_f\%np is the number of nodes of the finite element mesh <tt>mesh_f</tt>.
<li>The pressure and increment pressure are scalar so their format is a real valued tabular of three columns of dimension (mesh_c\%np,2,SIZE(list_mode)). We remind that mesh_c\%np is the number of nodes of the finite element mesh <tt>mesh_p</tt>.
</ol>

<h3>Exemple</h3>
Here is an exemple where we use the function vv_exact and pp_exact to initialize the velocity field and the pressire. Thus, the initial conditions satisfy the boundary conditions.
\code
   time = 0.d0
    DO i= 1, SIZE(list_mode)
       mode = list_mode(i) 
       DO j = 1, 6 
          !===velocity
          un_m1(:,j,i) = vv_exact(j,mesh_f%rr,mode,time-dt)  
          un   (:,j,i) = vv_exact(j,mesh_f%rr,mode,time)
       END DO
       DO j = 1, 2
          !===pressure
          pn_m2(:)       = pp_exact(j,mesh_c%rr,mode,time-2*dt)
          pn_m1  (:,j,i) = pp_exact(j,mesh_c%rr,mode,time-dt)
          pn     (:,j,i) = pp_exact(j,mesh_c%rr,mode,time)
          phin_m1(:,j,i) = pn_m1(:,j,i) - pn_m2(:)
          phin   (:,j,i) = Pn   (:,j,i) - pn_m1(:,j,i)
       ENDDO
    ENDDO
    RETURN
\endcode
We note that the integers mode, i and j have to be declared. Same for the real valued tabular pn_m2 whose dimension is mesh_c\%np. It is done by adding the two following lines in the declaration of the function <tt>init_velocity_pressure</tt>:
\code
    INTEGER                                    :: mode, i, j 
    REAL(KIND=8), DIMENSION(mesh_c%np)         :: pn_m2
\endcode

We refer to the sections \ref doc_debug_test and
\ref doc_example_physical_pb for more examples.


 */
