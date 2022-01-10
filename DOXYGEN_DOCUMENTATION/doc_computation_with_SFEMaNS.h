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
 * @page doc_computation_with_SFEMaNS Computations with SFEMaNS

To generate an executable of the code SFEMaNS and use it to do computations, you need the following files:
<ol>
<li><tt>make.inc</tt> that defines environment variables (path to PETSc, FFTW, ARPACK) and compilation options.
<li><tt>my_make</tt> that allows to generate a file called <tt>makefile</tt>. It contains the tree of the code (which module need which modules, etc.).
<li><tt>main.f90</tt> that is used to set the desired outputs with the subroutine <tt>my_post_processing</tt>.
<li><tt>condlim.f90</tt> that is used to set the initial conditions, boundary conditions and source term of the problem considered.
<li><tt>read_user_data.f90</tt> that is used to add new variables in the code. For instance, such variable can be used in <tt>condlim.f90</tt> to change the amplitude of a forcing term. 
<li><tt>data</tt> that gives various informations on the problem considered.
</ol>

Templates of each of the above files are available in the following directory: <tt>$SFEMaNS_DIR/TEMPLATE</tt>.

 We refer to this <a href='doc_installation.html#doc_install_sfemans'><code>section</code> </a> for more details on how to edit the file <tt>make.inc</tt>.

This section describes the fortran files and the data file needed to do computations with SFEMaNS. It is splitted into the following subsections.

<table width="90%" align="center">
<tr valign=top>
    <td width="30%" align="left">
    @subpage doc_SFEMaNS_notations
    </td>
    <td align="left">
    Description of the notations used in the code SFEMaNS.
    </td>
</tr>
<tr valign=top>
    <td width="30%" align="left">
    @subpage doc_SFEMaNS_data
    </td>
    <td align="left">
    Description of the <tt>data</tt> file.
    </td>
</tr>
<tr valign=top>
    <td width="30%" align="left">
    @subpage doc_SFEMaNS_read_user_data
    </td>
    <td align="left">
    Description of the file <tt>read_user_data.f90</tt>. It generates a derived data type called, denoted user, so one can add new variables in the code.
    </td>
</tr>
<tr valign=top>
    <td width="30%" align="left">
    @subpage doc_SFEMaNS_condlim
    </td>
    <td align="left">
    Description of the file <tt>condlim.f90</tt>. It is used to set the initial conditions, boundary conditions and the forcing terms of the equations considered.
    </td>
</tr>
<tr valign=top>
    <td width="30%" align="left">
    @subpage doc_SFEMaNS_main
    </td>
    <td align="left">
    Description of the file <tt>main.f90</tt>. We focus on the subroutine <tt>my_post_processing</tt> which computes the desired outputs (divergence, energy, 2D and 3D plots, etc.).
    </td>
</tr>
</table>


 */

<tr valign=top>
    <td width="30%" align="left">
    @subpage doc_SFEMaNS_tools_post_processing
    </td>
    <td align="left">
    Information on existing post-processing tools (display evolution of energy mode by mode, display energy spectrum at final time, etc.).
    </td>
</tr>
</table>


 */
