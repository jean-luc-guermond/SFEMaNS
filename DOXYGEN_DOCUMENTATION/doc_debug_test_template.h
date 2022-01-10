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
 * @page doc_debug_test_xx Test xx: 

<h3>Introduction</h3>




<h3>Manufactured solutions</h3>




<h3>Generation of the mesh</h3>


<table width="70%" align="center">
<tr>
    <td align="center">
    @image html  fig_Mesh_10_form.png  "Finite element mesh (P1)."
    </td>
</tr>
</table>







<h3>Information on the file <tt>condlim.f90</tt></h3>




All the other subroutines present in the file <tt>condlim_test_xx.f90</tt> are not used in this test.
We refer to the section \ref doc_SFEMaNS_condlim for a description of all the subroutines of the condlim file.







<h3>Setting in the data file</h3>
We describe the data file of this test.
 It is called <tt>debug_data_test_xx</tt> and can be found in the directory ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.




<h3> Outputs and value of reference </h3>

The outputs of this test are computed with the file <tt>post_processing_debug.f90</tt> 
that can be found in: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. 

To check the well behavior of the code, we compute four quantities:
<ol>
<li> 
<li> 
<li> 
<li> 
</ol>

 They are compared to reference values to attest of the correct behavior of the code.  
 These values of reference are in the last lines of the file <tt>debug_data_test_xx</tt> in ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC. They are equal to:
\code
values of reference
\endcode


To conclude this test, we display the profile of ??? at the final time.
 The two figures are done in the plane \f$y=0\f$ which is 
 the union of the half plane \f$\theta=0\f$ and \f$\theta=\pi\f$.
<table width="100%" align="center">
<tr>
    <td align="center">
     @image html  fig_test_01_pre_tfin.png "Pressure in the plane plane y=0."
    </td>
    <td align="center">
    @image html  fig_test_01_vel_tfin.png  "Velocity magnitude in the plane plane y=0."
    </td>
</tr>
</table>
 */
