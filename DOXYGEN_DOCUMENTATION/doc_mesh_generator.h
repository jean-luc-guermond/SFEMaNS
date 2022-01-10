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
 * @page doc_mesh_generator Mesh Generator


@section Mesh_Gene_Download Download

First download the mesh generator of SFEMaNS by typing the command:
\code
wget  http://math.tamu.edu/~guermond/DOWNLOADS/MESH_GENERATOR_SFEMaNS.tar.gz
\endcode

Extract this file in the directory where you want to install the mesh generator:
\code
tar -xvf MESH_GENERATOR_SFEMaNS.tar.gz
\endcode
It creates a directory called MESH_GENERATOR_SFEMaNS.
 To compile the executables that allow to generate P1-P2 finite element meshes, type the following commands:
\code
cd MESH_GENERATOR_SFEMaNS
make clean
make all
\endcode
It generates three executables called maill.exe, colle_p1p1.exe and symmetrize.exe. 
 Note that you may have to modify the file <tt> make.inc </tt> if you don't have access to mpi.f90.
 For pratical use, we advise the user to create aliases by typing:
\code
alias SFEMANS_MESH_GEN_DIR="path to directory of the mesh generator"
alias maill.exe = "($SFEMANS_MESH_GEN_DIR)/Code/maill.exe"
alias colle_p1p2.exe = "($SFEMANS_MESH_GEN_DIR)/Collage/colle_p1p2.exe"
alias symmetrize.exe = "($SFEMANS_MESH_GEN_DIR)/Collage/symmetrize.exe"
\endcode
Such aliases should be written in your .bashrc (or equivalent) to avoid redefining them everytime you log in.


@section Mesh_Gene_Visualization Visualization

When generating a mesh, as detailled below, many files are created so one can check if the mesh is generated correctly.
These files, that have the extension ".plt",  can be visualized with the executable plotmtv with the command:
\code
plotmtv file_to_visualize
\endcode
This executable and its documentation are available in  ($SFEMANS_MESH_GEN_DIR)/PLOTMTV.
 For practical use, you can create an alias by typing:
\code
alias plotmtv="($SFEMANS_MESH_GEN_DIR)/PLOTMTV/plotmtv"
\endcode


Plotmtv offers various mode of visualization (zooming, rotation, 2D and 3D rendering, etc.).
You can also create postscript file by clicking on <tt> print to file</tt>. 
Such a postscript file is always named  <tt>dataplot.ps</tt> so don't forget to modify its name after creation.
 Note that you can also edit files with the extension ".plt".


@section Mesh_Gene_create_mesh Generate a mesh

We now describe the procedure to be followed to create mesh of a
two-dimensional domain of your choice. To do that you need three files named: 
<table width="90%" align="center" >
<tr valign=top>
    <td width="20%" align="left">
     <tt>my_project.in </tt> to set the boundary of the mesh.
    </td>
</tr>
<tr valign=top>
    <td width="20%" align="left">
     <tt>topology.name_of_mesh </tt> to set the topology of the mesh at the boundary.
    </td>
</tr>
<tr valign=top>
    <td width="20%" align="left">
     <tt>backgrid.name_of_mesh </tt> to set specific topology at points inside the domain.
    </td>
</tr>
</table>
These files allow to create a grid with P1 elements with the following command:
\code
 maill.exe < my_project.in
\endcode


In order to create a grid with P1-P2 elements, you need a file named <tt>data_colle</tt> and type the following command:
\code
colle_p1p2.exe
\endcode

Eventually the P1-P2 mesh can be symmetrized with respect to a given line. It is done by editing the file <tt>data_symmetrize</tt> and typing the command:
\code
symmetrize.exe
\endcode

The following sections describe the above procedure on a template example. Before going on, you need to go to your work directory and type the following command:
\code
cp ($SFEMANS_MESH_GEN_DIR)/EXAMPLES/TEMPLATE/* .
\endcode
It copies the following files <tt>my_example.in</tt>, <tt>topology.example</tt>, <tt>backgrid.example</tt>, <tt>data_colle</tt> and <tt>data_symmetrize</tt> in your directory. The following explains how to edit these files and how to use the above executables correctly.


@subsection Mesh_Gene_P1_mesh Generation of a mesh P1

We describe how to modify the files   <tt>my_example.in</tt>, <tt>topology.example</tt> and <tt>backgrid.example </tt> that
 allows to generate a mesh P1. We note that such a mesh can't be used by the code SFEMaNS which requires a P1-P2 mesh. 
 Such a P1-P2 mesh is generated later using this P1 mesh.

@subsubsection  Mesh_Gene_my_project_file  The my_project.in file

The first step consists of defining the boundary of your domain.
Draw your domain on a piece of paper. Divide the boundary into
elementary paths. Starting from \f$1\f$, number the paths
composing the outer boundary going clockwise. Do the same
for the path extremities. The path extremities are called vertices. 
If any is present, continue to number the
inner boundaries and corresponding inner vertices
by going anti-clockwise. The rule is that the domain must always be
on your right-hand side as you go along boundaries.

The information on paths and vertices is stored in
<tt>my_project.in</tt>. For instance copy <tt>my_example.in</tt> to <tt>my_project.in</tt>
as follows
\code
cp  my_example.in my_project.in
\endcode

<ol>
<li> The first line of  <tt>my_project.in</tt>
is composed of the name of the project and an integer called integer_domain_index that, for the time being, you
set to \f$1\f$. It will be used later when assigning domain where each equations need to be solved in the finite element code.

name_of_mesh  integer_domain_index

<li> The second line contains an integer that is
equal to the number of paths that composes the boundary of the domain,
say \f$n_p\f$.

<li> The following \f$2 n_p\f$ lines describe the paths.
Four types of paths are possible:
<ol>
<li> Line segments:  line.
<li> Circle pieces:  circle.
<li> Ellipse pieces:  ellipse.
<li> Paths tabulated by the user in a data file.
</ol>

<li>  If you want to specify a line segment, proceed as follows by typing:

 line  integer_bdy_index

\f$x_1\f$ \f$y_1\f$ \f$x_2\f$ \f$y_2\f$

Here "integer_bdy_index" is an integer that you assign to the path. It will
be used latter when assigning boundary conditions in the finite element code.
For the time being you can always set this number to \f$1\f$ or to whatever value you want.

\f$x_1\f$ is the \f$x\f$-coordinate of the first vertex of the path.
\f$y_1\f$ is the \f$y\f$-coordinate of the first vertex of the path.
\f$x_2\f$ and \f$y_2\f$ are the \f$x\f$- and \f$y\f$-coordinates of the
second vertex of the path.

<li>  If you want to specify a portion of circle, proceed as follows by typing:

 circle  integer_bdy_index

\f$x_c\f$ \f$y_c\f$ \f$r\f$ \f$\theta_1\f$ \f$\theta_2\f$

where \f$x_c\f$ and \f$y_c\f$ are the \f$x\f$- and \f$y\f$-coordinates of 
the center of the circle; \f$r\f$ is the radius of the circle;
\f$\theta_1\f$ is the angle, measured is degrees, of the first vertex and \f$\theta_2\f$
is the angle of the second vertex. The origin of the angles is the horizontal 
half line starting at the circle center and pointing to the right.

<li>  If you want to specify a piece of ellipse, proceed as follows by typing:

 ellipse  integer_bdy_index

\f$x_c\f$ \f$y_c\f$ \f$a_x\f$  \f$a_y\f$ \f$\theta_1\f$ \f$\theta_2\f$

where \f$x_c\f$ and \f$y_c\f$ are the \f$x\f$- and \f$y\f$-coordinates of 
the center of the ellipse; \f$a_x\f$ is the first semi-axis and
\f$a_y\f$ is the second semi-axis;
\f$\theta_1\f$ is the angle, measured is degrees, of the first vertex and \f$\theta_2\f$
is the angle of the second vertex. The origin of the angles is the horizontal 
half line starting at the center of the ellipse and pointing to the right.

<li> If you want to specify a path tabulated in a file, proceed as follows by typing:

 data  integer_bdy_index

 file_containing_my_data

 where file_containing_my_data is a file written with the following format

\f[
\begin{matrix}
\text{number_of_points} &  &\\ 
x_1  &                           &  y_1 \\
x_2  &                           &  y_2 \\
\vdots &                         & \vdots \\
x_{\text{number_of_points}} & &   y_{\text{number_of_points}}
\end{matrix}
\f]
The first line contains an integer that is equal to the number of points that are 
tabulated. The other lines contain the coordinates of the points composing the path.
The first point is the first vertex of the path
and the last point is the second vertex. The points must be ordered so that as one 
progresses along the path the domain is on the right-hand side.
</ol>


@subsubsection Mesh_Gene_topology_file The topology.name_of_mesh file

The second step of the work consists of defining the mesh size at the boundary of the domain.
This is done by specifying the corresponding data in the file  <tt>topology.name_of_mesh</tt>.
The first part of the name, i.e.  topology is fixed, whereas
the second part, i.e.  name_of_mesh in the present case, is specified
by the user in the first line in the  <tt>my_project.in file</tt>.

To prepare your own file
copy   <tt>topology.example</tt> in <tt>topology.name_of_mesh</tt> as follows
\code
cp  topology.example topology.name_of_mesh
\endcode
and edit topology.name_of_mesh as follows.

<ol>
<li> The first line of the file is a comment that explains 
the meaning of the two integers that are on the second line.

<li> The second line contains two integers.
The first one is the number of paths composing the boundary
and the second one is the number of vertices.
These two numbers must be equal to \f$n_p\f$.
<li> The third line is a blank line.
<li> Then there are packets of lines
limited by the words  BEGIN and  END and these packets are separated by one
blank line. There are as many packets as paths, i.e. \f$n_p\f$ in all.
The  BEGIN and END lines are comments.
<ol>
<li> The first line after  BEGIN is a comment recalling 
the meaning of the integers in the second line.
<li> The second line contains four integers.
The first integer is the number of the paths. 
This number must be identical to the order of 
appearance of the path in  my_project.in.
The second integer says at how many curvilinear abscissas 
one will give the value of the mesh size; let \f$n_h\f$ be this integer.
The third integer is the number of the first vertex of the path and the 
fourth integer is the number of the second vertex of the path.
<li> The third line is a comment to remind the meaning of the following lines.
<li> Then there are \f$n_h\f$ lines before reaching the END line.
Each line is composed of one integer and two real numbers. Always set the first integer to \f$1\f$. The second number is the curvilinear abscissa
where one wants to enforce the value of the mesh size. The curvilinear abscissa \f$s\f$ is comprised
between \f$0\f$ and \f$1\f$. One must always start with \f$s=0\f$ and finish with \f$s=1\f$.
The third number is the value of the mesh size that one wants to set.
</ol>
</ol>

@subsubsection Mesh_Gene_backgrid_file The backgrid.name_of_mesh file

It happens that one needs to to refine the grid
inside the domain to increase the accuracy of the computation.
This can be done by using the file <tt>backgrid.name_of_mesh</tt>. 
The first part of the name, i.e.  backgrid is fixed, whereas
the second part, i.e.  name_of_mesh in the present case, is specified
by the user in the first line of the <tt>my_project.in</tt> file.
This files must exists even if the user does not want to refine the grid inside the domain.

Copy the file  <tt>backgrid.example</tt> in  <tt>backgrid.name_of_mesh</tt> as follows:
\code
cp  backgrid.example backgrid.name_of_mesh
\endcode
Then edit  <tt>backgrid.name_of_mesh</tt> for your own purposes.
<ol>
<li> The first line of the file is a comment recalling the meaning of the integer in second line.
<li> The second line contains an integer that is equal to the number of points
inside the domain where the user wants to enforce the mesh size, say \f$n_r\f$.
The number \f$n_r\f$ must be set to zero if the user does not want to refine the mesh.
<li> The third line has to be the following:

         X        Y      H 

It is a comment recalling the structure of the following lines.
<li> The \f$n_r\f$ lines thereafter contain three real numbers each, say
\f[
x \qquad y \qquad h
\f]
\f$x\f$ and \f$y\f$ are the coordinate of the point where one wants to refine the mesh.
\f$h\f$ is the mesh size that one wants to impose.
</ol>

We note that if \f$n_r=0\f$, the mesh generator only reads the first three lines of the file <tt> backgrid.name_of_mesh</tt>.

@subsubsection Mesh_Gene_maill The executable maill.exe

Once you have modified the three files described above, you can type the following command:
\code
 maill.exe < my_project.in
\endcode
It will create a mesh with P1 finite element called <tt>FEM.name_of_mesh</tt>.
 Moreover a file called <tt>gridplot.name_of_mesh</tt> is also created.
 It can be visualized with plotmtv to check the correct construction of the mesh.


@subsection Mesh_Gene_create_P2_mesh Generation of a P1-P2 mesh

To create a P1-P2 finite element mesh that the code SFEMaNS can used,
 you need to use the file named <tt>data_colle</tt> and the executable <tt>colle_p1p1.exe</tt>. 
After this step, you can symmetrize your P1-P2 mesh with respect to a given line with a file called <tt>data_symmetrize</tt>
 and the executable <tt>symmetrize.exe</tt>. We draw the attention of the reader that the symmetrization step is optional.

@subsubsection Mesh_Gene_data_colle_file The data_colle file 

This file allows to define the number of subdomains, meaning P1 meshes, that you are you going to glue together so you can
generate a global P1-P2 mesh that will be used by the code SFEMaNS.
<ol>
<li> The first line is a logical (.t. to generate formatted mesh and .f. to generate unformatted mesh)
<li> The second line is the number of P1-mesh to glue.
<li> The third line contains the path and the name of the first P1 mesh.
<li> If you want to create P1-P2 mesh that is the union of two or more P1 mesh, you need to add the three following lines for each P1 mesh to add:
<ol>
<li> The first line is an integer \f$\text{nb}_\text{inter}\f$ that represents the number of interfaces to keep between the two meshes.
<li> The second line contains \f$\text{nb}_\text{inter}\f$ integer that are the index of the interfaces to keep
 (defined in the files <tt>*.in</tt>).
<li> The third line contains the path and the name of the P1 mesh to add.
</ol>
</ol>

Here is an example of a file <tt> data_colle</tt> that generates a formatted P1-P2 mesh by gluing
  two P1-meshes. We assume that theses meshes are in your working directory and are called
  <tt>FEM.name_of_mesh1</tt> and <tt>FEM.name_of_mesh2</tt>.
 Moreover we decide to keep one interface of index 3.
\code
.t.    ! .t.= formatted  .f.=unformatted
2      ! nb of subdomain
. FEM.name_of_mesh1
1     ! nb of interfaces to keep
3    ! index of interfaces
. FEM.name_of_mesh2
\endcode

@subsubsection Mesh_Gene_colle_p1p12 The executable colle_p1p2.exe
Once you modified your data_colle, you can type the following command:
\code
colle_p1p2.exe
\endcode
To the question " Attention to cavities: number of cavities = ", type the number of cavities present in your mesh.
A cavity is a hole in the domain; it is not part of your mesh.

The executable generates a P1-P2 mesh called COLLE.FEM that can be used by the code SFEMaNS.
 Moreover it also generates files that can be visualized with plotmtv such as colle.plt (final P1 mesh), colle_dom1.plt (final P2 mesh).
 We note the file sides.plt, that represent the final P1 mesh, also displays the boundary and the interface you conserved.
 When visualizing sides.plt with plotmtv, you can click on "3D plot" to see the index of all of the boundaries and interfaces 
 (previously defined in the file <tt>my_project.in</tt>).


@subsubsection Mesh_Gene_symmetrize Optional symmetrization step
The last feature of SFEMaNS's mesh generator is to symmetrize a P1-P2 mesh with respect to a given line.
First you need to modify the file <tt>data_symmetrize</tt> that you previously copied in your work directory.
<ol>
<li> The first line contains a logical (.t. if the unsymmetrized mesh is formatted and .f. if it is unformatted).
<li> The second line contains the directory of the mesh to symmetrize and the name of the mesh to symmetrize.
<li> The third line contains the three coefficient a, b and c where ax+by+c=0 is the equation of the symmetry line.
<li> The fourth (and last) line contains a logical (.t. to keep the interface between sub-meshes, .f. else).
</ol>

Once you modified the file <tt>data_symmetrize</tt>, you can type the following command:
\code
symmetrize.exe
\endcode
It will generate a symmetrized P1-P2 mesh called <tt>symm_COLLE.FEM</tt> and "*.plt" files 
(such as sides_p2.plt, i_d_p2.plt, etc.) that can be visualized with plotmtv.



@subsection Mesh_Gene_Example Example 

The goal of this example is to use all the executables described above to generate a P1-P2 symmetrized mesh.
We plan to construct a mesh that correspond to a cylinder with radius R=1 and height H=2 that is
 included in a sphere of center (0,0) and of radius 11.
 The symmetrization will be done with respect to the line y=0.
 The coordinates (x,y) correspond to the cylindrical coordinates (r,z) of the cylinder, meaning r for the radius and the z for the height in the direction of the revolution axis of the cylinder. 
 All the files required to do such a mesh can be found in the directory:

<tt> ($SFEMaNS_MESH_GEN_DIR)/EXAMPLES/EXAMPLE_DOC </tt> .


@subsubsection Mesh_Gene_Example_step1 Generation of the P1 meshes

First we set the boundaries of the cylinder in the file <tt>rect_in.in </tt> as follows:
\code
rect_in 1
4
line 1
0. 0. 0. 1.
line 2
0. 1. 1. 1.
line 3
1. 1. 1. 0.
line 4
1. 0. 0. 0.
\endcode
Since we plan to symmetrize the mesh with respect to the line \f$y=0\f$, this files represent a meridian plane of the cylinder's upper half.

Then we choose a mesh size of \f$0.05\f$ in P1 by creating the file  <tt>topology.rect_in</tt> as follows:
\code
NV      NE
 4       4

BEGIN
       E       N      BV      EV
       1       2       1       2
   IDATA   SDATA   HDATA
       1     0.0    0.05
       1     1.0    0.05
END

BEGIN
       E       N      BV      EV
       2       2       2       3
   IDATA   SDATA   HDATA
       1     0.0   0.05
       1     1.0   0.05
END

BEGIN
       E       N      BV      EV
       3       2       3       4
   IDATA   SDATA   HDATA
       1     0.0    0.05
       1     1.0    0.05
END

BEGIN
       E       N      BV      EV
       4       2       4       1
   IDATA   SDATA   HDATA
       1     0.0   0.05
       1     1.0   0.05
END
\endcode
Moreover we impose a mesh size of \f$ 0.025\f$ inside the domain at the point \f$(x,y)=(0.5,0.5)\f$
 by setting the <tt>backgrid.rect_in</tt> as follows:
\code
        N
        1
        X        Y      H
       0.5     0.5      0.025
\endcode
We can now generate the P1 mesh with the command:
\code
maill.exe < rect_in.in
\endcode
It creates a file called <tt> FEM.rect_in </tt> that contains the information of the upper cylinder P1 mesh. 
 In addition, a file called gridplot.rect_in is generated and can be visualized with plotmtv.

We now repeat the process for the spherical domain. While doing so, we need to ensure that the mesh size (topology) of 
 the sphere and the cylinder coincide on their interfaces and also that these interfaces have the same index.
 So we create a the file <tt>circle_ext.in</tt> as follows:
\code
circle_ext 2
5
circle 5
0 0  11 90.d0  0.d0
line 4
11 0. 1. 0
line 3
1.0 0. 1.0 1.0
line 2
1.0 1.0 0.0 1.0
line 1
0.0 1.0 0 11.0
\endcode
where we keep the index 2 and 3 for the interfaces between the cylinder and the sphere. Note that the lines 1 and 4 have the same
 index that in <tt>rect.in</tt> but they are different boundaries.

 Regarding the topology of the sphere, we use a mesh size of \f$1\f$ on the outter boundary, meaning on \f$ x^2+y^2=11^2\f$.
 Moreover, we enforce a mesh size of \f$0.05\f$ for \f$ x\leq 2 \f$ and \f$ y \leq 2\f$,
 which correspond to \f$ 10 \%\f$ of the arc lenght. 
 The resulting file <tt>topology.circle_ext</tt> is written as follows:
\code
NV      NE
 5       5

BEGIN
       E       N      BV      EV
       1       2       1       2
   IDATA   SDATA   HDATA
       1     0.0    1.0
       1     1.0    1.0
END

BEGIN
       E       N      BV      EV
       2       3       2       3
   IDATA   SDATA   HDATA
       1    0.00   1.0
       1    0.90   0.05
       1    1.00   0.05
END

BEGIN
       E       N      BV      EV
       3       2       3       4
   IDATA   SDATA   HDATA
       1     0.0    0.05
       1     1.0    0.05
END

BEGIN
       E       N      BV      EV
       4       2       4       5
   IDATA   SDATA   HDATA
       1    0.00   0.05
       1    1.00   0.05
END

BEGIN
       E       N      BV      EV
       5       3       5       1
   IDATA   SDATA   HDATA
       1    0.00   0.05
       1    0.10   0.05
       1    1.00   1.00
END
\endcode

We do not impose specific refinement inside the domain so we define the file <tt>backgrid.circle_ext</tt> as follows:
\code
        N
        0
        X        Y      H
\endcode

We create the P1 mesh by executing the following commands:
\code
maill.exe < circle_ext.in
\endcode
It generates a P1 mesh called <tt>FEM.circle_ext </tt> and a plotmtv file called <tt>backgrid.circle_ext</tt>.

@subsubsection Mesh_Gene_Example_step2 Generation of the P1-P2 global mesh

To create a global P1-P2 mesh of the domain, we define the following <tt>data_colle</tt>:
\code
.t.    ! .t.= formatted  .f.=unformatted
2      ! nb of subdomain
. FEM.rect_in
2      ! nb of interfaces to keep
2 3    ! index of interfaces
. FEM.circle_ext
\endcode
where we specified that we want to keep the interfaces of index 2 and 3 between the cylinder and the circle
 (so one can later impose boundary or interface conditions).

 We can now execute the following commands to generate the global P1-P2 mesh:
\code
colle_p1p2.exe
\endcode
We answer 0 to the question: "Attention to cavities: number of cavities = ".
It generates a P1-P2 mesh called <tt>COLLE.FEM </tt> that can be used by the code SFEMaNS.
The correct generation of the mesh can be checked by visualizing various outputs (such as side_p1.plt, sides.plt, etc.).
with plotmtv.


Eventually we symmetrize the mesh with respect to the line \f$ y=0 \f$ with the following <tt> data_symmetrize</tt>:
\code
.t. !    .t.= formatted  .f.=unformatted
'.' 'COLLE.FEM'     ! directory of the mesh to symmetrize   and  name of the mesh to symmetrize
0.d0 1.d0 .0d0      ! a,b,c,  where ax+by+c=0 is the equation of symmetry line
.t.                 ! keep sides set to true
\endcode
and by typing the command:
\code
symmetrize.exe
\endcode
It generates a symmetrized P1-P2 mesh of <tt>COLLE.FEM</tt> called <tt>symm_COLLE.FEM </tt>. 
 The symmetrized mesh can be visualized with the files <tt>sides_p1.plt</tt>,  <tt>i_d_p1.plt</tt>, etc.


It concludes this example. We note other examples of mesh generation are present in <tt> ($SFEMaNS_MESH_GEN_DIR)/EXAMPLES</tt>.
 We refer to the sections \ref doc_debug_test and
 \ref doc_example_physical_pb for more details.


 */
