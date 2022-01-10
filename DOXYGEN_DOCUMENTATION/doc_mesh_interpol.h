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
 * @page doc_mesh_interpol Mesh Interpolation and restart files

The code SFEMaNS allows to interpolate restart files (from the Navier-Stokes equations, the temperature equation and the Maxwell equations) on a different finite element mesh. It is done with an executable called <tt>i.exe</tt> and a shell script called <tt>mesh_interpol</tt>. To interpolate restart files, the following files are required:
<ol>
<li><tt>make.inc</tt> that defines environment variable (path to PETSc, FFTW, ARPACK) and compilation options. 
<li><tt>my_make_mesh_interpol</tt> that allows to generate a file called makefile.
<li><tt>main_mesh_interpol.f90</tt>, <tt>mesh_interpolation.f90</tt> and <tt>tools_interpol.f90</tt>. They are the fortran files that generate the executable i.exe.
<li><tt>mesh</tt><tt>_interpolation</tt> is a shell script. It allows to turn restart files on a finite element mesh M1 into restart files on a finite element mesh M2.
<li><tt>data_interpol</tt> that gives various informations on the interpolation to do.
<li>the suite files to interpol on a different finite element mesh.
<li>the mesh_part file associated to the suite files. It contains information on how the meridian section has been partitionned for parallel computing.
</ol>

Templates of each of the above files (modulo the suite and mesh_part files) are available in the following directory: $SFEMaNS_DIR/MESH_INTERPOLATION. A template of the file <tt>make.inc</tt> is available in the following directory: $SFEMaNS_DIR/TEMPLATE. We note that only the files <tt>mesh_interpol</tt> and <tt>data_interpol</tt> require modifications.

The following describes how to generate the executable <tt>i.exe</tt> and modify the files <tt>mesh_interpol</tt> and <tt>data_interpol</tt>.


@section doc_mesh_inter_exec_gene Generation of the executable i.exe

First create a copy of the directory $SFEMaNS_DIR/MESH_INTERPOLATION in the directory where you plan to do the interpolation. It can be done as follows:
\code
cp -rf  $SFEMaNS_DIR/MESH_INTERPOLATION MESH  MY_MESH_INTERPOL
\endcode

Then generate the executable <tt>i.exe</tt> as follows.
<ol>
<li>Go into the directory MY_MESH_INTERPOL.
\code
cd MY_MESH_INTERPOL
\endcode
<li>Copy the <tt>make.inc</tt> template in your directory as follows:
\code
cp $SFEMaNS_DIR/TEMPLATE/make.inc_template make.inc
\endcode
<li>Edit the file <tt>make.inc</tt>. We refer to the section 
\ref doc_install_sfemans_check_install for more details.
 We note that you can use the make.inc used to generate
 the executable a.exe.
<li>Generate a makefile with the command:
\code
./my_make_mesh_interpol
\endcode
<li>Compile the executable i.exe with the command:
\code
make i.exe -j 4
\endcode
The option -j 4 means 4 processors are used to do the compilation. It can be removed or used with a different number of processors.
</ol>

Remark: the shell <tt>debug_SFEMaNS_template</tt>, used to check the correct intallation of the code, generates an executable i.exe. Indeed by following the instructions of this <a href='doc_installation.html#doc_install_sfemans_check_install'><code>section</code></a>, an executable i.exe is available in the following directory: MY_APPLICATIONS_SFEMaNS_v5.3/CHECK_INSTALLATION/INTERPOLATION_MESH. It can be used here to avoid generating a new executable i.exe.



@section doc_mesh_inter_restart_file Interpollation of restart files


@subsection doc_mesh_inter_mesh_inter The shell mesh_interpolation

Edit the shell <tt>mesh</tt><tt>_interpolation</tt> as follows.
<ol>
<li>Uncomment the following line.
\code
#export MY_MPIRUN=mpirun
\endcode
<li>Define the variable MY_MPIRUN depending of your work environment (mpirun, mpiexec, etc.).
</ol>

@subsection doc_mesh_inter_data The file data_interpol

Edit the file <tt>data_interpol</tt> as follows.
<ol>
<li>Give information on the initial restart and mesh_part files.
<ol>
<li>Set the number of processors used in meridian section.
\code
===Number of processors in meridian section (Input)
1
\endcode
<li>Set the initial problem type approximated.
\code
===Problem type (Input): (nst, mxw, mhd, fhd)
'nst'
\endcode
<li>Set if the temperature equation was approximated. It informs the code of the presence of restart files for the temperature field.
\code
===Is there an input temperature field?
.t.
\endcode
<li>The path and the name of the initial mesh are specified with the two following lines:
\code
===Directory and name of input mesh file
'.' 'SOLID_FLUID_10.FEM'
\endcode
<li>Set the format of the initial mesh.
\code
===Is input mesh file formatted (true/false)?
.t.
\endcode
</ol>
<li>Give information on the interpolated restart and mesh_part files.
<ol>
<li>Set the desired number of processors in meridian section.
\code
===Number of processors in meridian section (Output)
5
\endcode
<li>The interpolated restart files are going to be used on which problem type.
\code
===Problem type (Output): (nst, mxw, mhd, fhd)
'nst'
\endcode
<li>Set if restart files of the temperatured field are generated.
\code
===Is there an output temperature field?
.t.
\endcode
<li>The path and the name of the new mesh are specified with the two following lines:
\code
===Directory and name of output mesh file
'.' 'SOLID_FLUID_20.FEM'
\endcode
<li>Set the format of the new mesh. 
\code
===Is output mesh file formatted (true/false)?
.t.
\endcode
</ol>
<li>Give information of the interpolation process.
<ol>
<li>Set if the initial datas are interpolated on a new mesh.
\code
===Should data be interpolated on new mesh? (True/False)
.t.
\endcode
<li>Set the number of time steps where the files are interpolated.
\code
===How many files (times steps) should be converted?
1
\endcode
<li>Set the index of the first time step.
\code
===What is the starting index? (suite*_I*index.mesh*), index is an integer in [1,999]
1
\endcode
<li>Set the number of Fourier modes present in the initial restart files.
\code
===Number of Fourier modes
3
\endcode
If the initial restart file are defined on a specific list of Fourier modes, add the following lines.
<ol>
<li>Set if a specific list of Fourier is used.
\code
===Select Fourier modes? (true/false)
\endcode
<li>Give the list of Fourier modes considered.
\code
===List of Fourier modes (if select_mode=.TRUE.)
\endcode
</ol>
</ol>
<li>Set periodic condition (if any).
<ol>
<li>Set the number of pair of boundaries that are periodic.
\code
===How many pieces of periodic boundary?
1
\endcode
<li>Give the label of the boundaries and the vector that lead to the first boundary to the second one.
\code
===Indices of periodic boundaries and corresponding vectors
4 2 .0d0 1.d0
\endcode
We note that we need as much as lines as the number of pairs of boundaries with periodic condition.
</ol>
<li>Give information on the Navier-Stokes mesh (if any).
<ol>
<li>Set the number of domains and their label where the code approximates the Navier-Stokes equations.
\code
===Number of subdomains in Navier-Stokes mesh
1
===List of subdomains for Navier-Stokes mesh
2
\endcode
<li>Set if there is a level (for multiphase problem).
\code
===Is there a level set?
.f.
\endcode
<li>If a level set is present, set the number of fluids considered.
\code
===How many fluids?
2
\endcode
</ol>
<li>Give information on the temperature mesh (if any).
<ol>
<li>Set the number of domains and their label where the code approximates the temperature equation.
\code
===Number of subdomains in temperature mesh
2
===List of subdomains for temperature mesh
1 2
\endcode
<li>Set the number of interface between the temperature and the velocity field domains and give their respective labels.
\code
===Number of interfaces between velocity and temperature only domains (for nst applications)
1
===List of interfaces between velocity and temperature only domains (for nst applications)
3
\endcode
It is only used when there is no magnetic field domain.
</ol>
<li>Give information on the magnetic field mesh (if the Maxwell equations are involved).
<ol>
<li>Set the number of domains and their label where the code approximates the Maxwell equations with the magnetic field.
\code
===Number of subdomains in magnetic field (H) mesh
0
===List of subdomains for magnetic field (H) mesh
1
\endcode
<li>Set the number of interface in H_mesh and give their respective labels.
\code
===Number of interfaces in H mesh
0
===List of interfaces in H mesh
0
\endcode
Such interfaces represent interfaces with discontinuities in magnetic permeability or interfaces between the magnetic field mesh and the temperature or the velocity field meshes.
<li>Set the type of finite element used to approximate the magnetic field.
\code
===Type of finite element for magnetic field
2
\endcode
</ol>
<li>Give information on the scalar potential mesh (if the Maxwell equations are involved).
<ol>
<li>Set the number of domains and their label where the code approximates the Maxwell equations with the scalar potential.
\code
===Number of subdomains in magnetic potential (phi) mesh
0
===List of subdomains for magnetic potential (phi) mesh
2
\endcode
<li>Set the number of interface between the magnetic field and the scalar potential and give their respective labels.
\code
===Number of interfaces between H and phi
1
===List of interfaces between H and phi
5
\endcode
<li>Set the type of finite element used to approximate the scalar potential.
\code
===Type of finite element for scalar potential
2
\endcode
</ol>
<li>You can generate files to check the correctness of the interpolation.
\code
===Check construction with plotmtv? (True/False)
.f.
\endcode
These files have the extension ".plt" and can be open with plotmtv.
</ol>



@subsection doc_mesh_inter_process Generation of interpolated restart files

To interpolate your restart files on the new mesh, you need to copy the following files in your directory:
<ol>
<li>The restart files.
<li>The file mesh_part associated to the restart files. This file contains information on how the initial mesh has been partitionned for parallel computing. It depends of the number of processors used per meridian section.
</ol>

Then type the following command:
\code
./mesh_interpolation
\endcode
It generates restart files on the new finite element mesh and a mesh_part file associated to these interpolated restart files. We note that to use these files to restart a computation on this new mesh, the number of the iteration of the restart files (_I001, _I002, etc.) has to be removed.

Remark: the initial restart files are erased during the interpolation process. So always keep a copy of them in an other directory.
 */
