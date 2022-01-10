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
* @mainpage Documentation


The objective of this documentation is to help users and developers to use all the features of SFEMaNS and add new ones. The documentation is divided in three sections.

<h2>Required information to compute with SFEMaNS</h2>
<table width="90%" align="center">
<tr valign=top>
    <td width="30%" align="left">
    @subpage doc_intro
    </td>
    <td align="left">
    General information on SFEMaNS (equations solved, computational domain, and key features of the code).
    </td>
</tr>
<tr valign=top>
    <td width="30%" align="left">
    @subpage doc_num_app
    </td>
    <td align="left">
    Information on the approximation techniques.
    Description of the algorithms and weak formulations that are implemented.
    </td>
</tr>
<tr valign=top>
    <td width="30%" align="left">
    @subpage doc_installation
    </td>
    <td align="left">
    Download and check installation of SFEMaNS.
    </td>
</tr>
<tr valign=top>
    <td width="30%" align="left">
    @subpage doc_mesh_generator
    </td>
    <td align="left">
    Creation of grid with a mesh generator.
    </td>
</tr>
<tr valign=top>
    <td width="30%" align="left">
    @subpage doc_computation_with_SFEMaNS
    </td>
    <td align="left">
    Description of the files that the user must provide to actually compute something (condlim.f90, main.f90, read_user_data.f90 and data). Templates are provided. Information on available post-processing tools.
    </td>
</tr>
<tr valign=top>
    <td width="30%" align="left">
    @subpage doc_mesh_interpol
    </td>
    <td align="left">
    Description of the interpolation of restart files on a mesh M1 to restart files on a mesh M2.
    </td>
</tr>
</table>


<h2>Examples of computations in various settings</h2>
<table width="90%" align="center" >
<tr valign=top>
    <td width="30%" align="left">
    @subpage doc_debug_test
    </td>
    <td align="left">
    Examples showing how to use SFEMaNS with manufactured solutions. These examples are used to check the correctness of the installation of SFEMaNS and the various softwares that are invoked by SFEMaNS (PETSc, FFTW3, ARPACK).
    </td>
</tr>
<tr valign=top>
    <td width="30%" align="left">
    @subpage doc_example_physical_pb
    </td>
    <td align="left">
    Examples showing how to use SFEMaNS in various physical settings (hydrodynamic, magnetism, MHD). These examples been published in referred journals.
    </td>
</tr>
</table>



<h2>Additional information to add terms in equations or implement new equation</h2>
<table width="90%" align="center" >
<tr valign=top>
    <td width="30%" align="left">
    in construction
    </td>
    <td align="left">
    Theoretical reminder on weak formulations and Gaussian quadratures.
  Description of the fortran files that assemble matrices (left-handside) and vectors (right-handside) induced by the weak formulation of the MHD equations.
    </td>
</tr>
<tr valign=top>
    <td width="30%" align="left">
    in construction
    </td>
    <td align="left">
    Computation of nonlinear terms with SFEMaNS via Fast Fourier Transform (FFTW3 package).
    </td>
</tr>
</table>
*/

<tr valign=top>
    <td width="30%" align="left">
    in construction
    </td>
    <td align="left">
    This section follows the structure/tree of the code and proposes a description of the fortran files
 involved in the initialization process and in the approximation of an equation
 (temperature, mass, Navier-Stokes and Maxwell).
    </td>
</tr>
</table>
*/
