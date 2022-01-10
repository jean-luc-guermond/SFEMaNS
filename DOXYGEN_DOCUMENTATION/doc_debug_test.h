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
 * @page doc_debug_test Examples with manufactured solutions

In this section we describe numerous examples that mainly involve
 manufactured solutions on hydrodynamic, magnetic and MHD set ups.
These tests have been implemented in the code SFEMaNS so users can
 check that the installation of  SFEMaNS and the required softwares
 (PETSC, ARPACK, etc.) have been done correctly.
 These tests can be run with the shell <tt>debug_SFEMaNS_template</tt>.
 We refer to the section \ref doc_install_sfemans_check_install for
 more information on how to check the correct installation of the
 code SFEMaNS.

We note that these tests don't read the information in 
 <tt>condlim.f90</tt>, <tt>read_user_data.f90</tt> and <tt>main.f90</tt>:
<ol>
<li>The initial/boundary conditions and forcing term are read in a condlim file specific for each test.
 They can be found in the following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC/. They are named 
 <tt>condlim_test_nb_test</tt> with nb_test the number of the test.
<li>Each test uses its own data that are in the following directory: 
 ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC/. They are named <tt>debug_data_nb_test</tt> with nb_test the number of the test.
<li>As these tests compute similar outputs, the postprocessing is done via the same fortran file
 post_processing_debug.f90. It can be found in the following directory: ($SFEMaNS_DIR)/MHD_DATA_TEST_CONV_PETSC.
</ol>

The description of each of these tests follows this structure:
<ol>
<li> Introduction: describe the purpose of the test and the equations to solve.
<li> Manufactured solutions: introduce the solutions considered.
<li> Generation of the mesh: information on the mesh and its topology.
<li> Information on the file <tt>condlim</tt>: describe the implementation of initial/boundary conditions and source terms.
<li> Setting in the data file: description of the <tt>data</tt> file.
<li> Outputs and value of reference: describe outputs considered and the values to recover with SFEMaNS.
</ol>

Here is the list of the tests. Their respective title provide information on 
the equations involved.
<table width="90%" align="center">
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_01
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_02
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_03
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_04
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_05
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_06
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_07
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_08
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_09
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_10
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_11
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_12
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_13
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_14
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_15
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_16
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_17
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_18
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_19
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_20
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_21
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_22
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_23
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_24
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_25
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_26
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_27
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_28
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_29
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_30
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_31
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_32
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_33
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_34
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_35
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_36
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_37
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_38
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_39
    </td>
</tr>
<tr valign=top>
    <td align="left">
@subpage doc_debug_test_40
    </td>
</tr>
</table>

 */
