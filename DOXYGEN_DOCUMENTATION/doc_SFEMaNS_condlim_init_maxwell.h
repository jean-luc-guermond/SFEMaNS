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
 * @page doc_SFEMaNS_condlim_init_maxwell The subroutine init_maxwell

It is used to initialize the magnetic field \f$\bH\f$ and the scalar potential \f$\phi\f$.

<h3>Inputs and outputs</h3>
The inputs of this function are the following:
<ol>
<li><tt>H_mesh</tt> is the finite element mesh used to approximate the magnetic field.
<li><tt>phi_mesh</tt> is the finite element mesh used to approximate the scalar potential.
<li><tt>dt</tt> is the time step.
<li><tt>mu_H_field</tt> is a list of real. It contains the magnetic permeability of each conducting domains.
<li><tt>mu_phi</tt> is the magnetic permeability in the insulating domain where \f$\phi\f$ is approximated.
<li><tt>list_mode</tt> is a list of integers which contains the Fourier modes approximated.
</ol>
As the mesh can be subdivised, we note that  H_mesh and phi_mesh depend of the processor considered when doing parallel computing. Same goes for the list of Fourier mode <tt>list_mode</tt>.

The outputs of this function are the following:
<ol>
<li><tt>time</tt> is the time when the computations starts. 
<li><tt>Hn_m1</tt> is the magnetic at the time <tt>time-dt</tt>.
<li><tt>Hn</tt> is the magnetic field at the time <tt>time</tt>.
<li><tt>phin_m1</tt> is the scalar potential at the time <tt>time-dt</tt>.
<li><tt>phin</tt> is the scalar potential at the time <tt>time</tt>.
</ol>

Remarks:
<ol>
<li>The magnetic field is a vector so its format is a real valued tabular of three columns with dimension (H_mesh\%np,6,SIZE(list_mode)). We remind that H_mesh\%np is the number of nodes of the finite element mesh <tt>H_mesh</tt>.
<li>The scalar potential is a scalar so its format is a real valued tabular of three columns of dimension (phi_mesh\%np,2,SIZE(list_mode)). We remind that phi_mesh\%np is the number of nodes of the finite element mesh <tt>phi_mesh</tt>.
</ol>

<h3>Exemple</h3>
Here is an exemple where we use the function Hexact and Phiexact (so initial conditions satisfy the boundary conditions).
\code
    time = -dt
    DO k=1,6
       DO i=1, SIZE(list_mode)
          Hn1(:,k,i) = Hexact(H_mesh,k, H_mesh%rr, list_mode(i), mu_H_field, time)
          IF (inputs%nb_dom_phi>0) THEN
             IF (k<3) THEN
                phin1(:,k,i) = Phiexact(k, phi_mesh%rr, list_mode(i) , mu_phi, time)
             ENDIF
          END  IF
       ENDDO
    ENDDO

    time = time + dt
    DO k=1,6
       DO i=1, SIZE(list_mode)
          Hn(:,k,i) = Hexact(H_mesh,k, H_mesh%rr, list_mode(i), mu_H_field, time)
          IF (inputs%nb_dom_phi>0) THEN
             IF (k<3) THEN
                phin(:,k,i) = Phiexact(k, phi_mesh%rr, list_mode(i), mu_phi, time)
             ENDIF
          END  IF
       ENDDO
    ENDDO
    RETURN
\endcode
We note that the integers i and k have to be declared. It is done by adding the two following lines in the declaration of the function <tt>init_maxwell</tt>:
\code
    INTEGER                                    :: i, k
\endcode


We refer to the sections \ref doc_debug_test and
\ref doc_example_physical_pb for more examples.

 */
