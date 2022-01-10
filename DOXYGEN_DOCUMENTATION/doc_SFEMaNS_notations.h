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
 * @page doc_SFEMaNS_notations Notations in SFEMaNS





@section doc_notation_gen General notations




@subsection doc_notation_gen_vect_scal Declaration of scalars and vectors


Fourier modes and nodes associated to one processor.

np = number of nodes per processors.
list_mode = list of the Fourier modes that the processors is dealing with.


Declaration of scalar and vector.

dimension scalar = (np, 2, SIZE(list_mode))
dimension vector = (np, 6, SIZE(list_mode))



@subsection doc_notation_gen_comm Communication between processors

Communication between processors.

comm_one_d_ns
comm_one_d_temp
comm_one_d




@section doc_notation_derived_data_mesh The derived data type mesh

<tt>def_type_mesh.f90</tt>

me
np
jj
rr
hloc

mes
nps
jjs

Gauss:
n_w
l_G
hloc_gauss
ww
dw
rj








 */
