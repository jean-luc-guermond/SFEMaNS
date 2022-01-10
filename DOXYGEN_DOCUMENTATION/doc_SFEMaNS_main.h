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
 * @page doc_SFEMaNS_main Fortran file main.f90

This section describes the fortran file <tt>main.f90</tt>. It focuses on the subroutine <tt>my_post_processing</tt> that computes the desired outputs. We note that the other subroutines of the file <tt>main.f90</tt> should not be modified and are not described here. A template of the file <tt>main.f90</tt> is available in the following directory: <tt>($SFEMaNS_DIR)/TEMPLATE</tt>.


Remark: If ARPACK is not installed, the following lines of the file <tt>main.f90</tt> have to be commented.
<ol>
<li>
\code
    USE arpack_mhd
\endcode
<li>
\code
     CALL solver_arpack_mhd(comm_one_d,H_mesh,phi_mesh,&
          inputs%dt,list_mode,mu_H_field)
\endcode
</ol>
Apart from these lines, we insist on the fact that only the subroutine <tt>my_post_processing</tt> should be modified.

This section is splitted into four subsections. First, the structure of the subroutine <tt>my_post_processing</tt> is described. Then information on functions that compute classic quantities (like L2 norm, divergence, etc.) are given. The third subsection describes how to generate 2D and 3D visualization files for Paraview. Eventually an example is provided.


@section doc_SFEMaNS_main_my_post_processing The subroutine my_post_processing

The subroutine <tt>my_post_processing</tt> is called after each time iterations. We denote by <tt>tn</tt> the time after the time iteration n. This subroutine has access to the variables defined at the begining of the file <tt>main.f90</tt>.  Here is a list of the variables that are meaningfull to compute outputs:
<ol>
<li><tt>pp_mesh</tt> is the finite element mesh used to approximate the pressure and the level set.
<li><tt>vv_mesh</tt> is the finite element mesh used to approximate the velocity field.
<li><tt>pn</tt> is the pressure at time tn.
<li><tt>un</tt> is the velocity field at time tn.
<li><tt>level_set</tt> is the level set at time tn (for multiphase computation).
<li><tt>density</tt> is the density at time tn (for multiphase computation).
<li><tt>H_mesh</tt> is the finite element mesh used to approximate the magnetic field.
<li><tt>phi_mesh</tt> is the finite element mesh used to approximate the scalar potential.
<li><tt>Hn</tt> and <tt>Bn</tt> are the magnetic fields at time tn. We remind that \f$\mu\textbf{H}=\bB\f$ with \f$\mu\f$ the magnetic permeability.
<li><tt>sigma_field</tt> is the list of the electrical conductivity of each domains where the magnetic field is approximated.
<li><tt>mu_h_field</tt> is the list of the magnetic permeability of each domains where the magnetic field is approximated.
<li><tt>temp_mesh</tt> is the finite element mesh used to approximate the temperature.
<li><tt>temperature</tt> is the temperature at time tn.
<li><tt>temperature_diffusivity_field</tt> is the list of the temperature diffusivity of each domains where the temperature is approximated.
<li><tt>m_max_c</tt> is the number of Fourier modes approximated by each processors.
<li><tt>list_mode</tt> is the list of the Fourier mode approximated by the processor. We note that m_max_c=SIZE(list_mode).
<li><tt>time</tt> is the time tn.
<li><tt>comm_one_d_ns</tt> is the communicator on the Navier-Stokes equations domain.
<li><tt>comm_one_d_temp</tt> is the communicator on the temperature equation domain.
<li><tt>comm_one_d</tt> is the communicator on the whole domain. It is equal to the Maxwell equations domain if they are approximated.
</ol>


The structure of the subroutine <tt>my_post_processing</tt> is the following:
<ol>
<li>First, the ranks of each processors are defined. It allows to have one processor writing outputs that depends of a specific Fourier mode or region of the finite element mesh. It can be done as follows.
\code
    !===Check ranks
    IF (vv_mesh%me /=0) THEN
       CALL MPI_Comm_rank(comm_one_d_ns(1), rank_ns_S, ierr)
       CALL MPI_Comm_rank(comm_one_d_ns(2), rank_ns_F, ierr)
    ELSE
       rank_ns_s = -1
       rank_ns_f = -1
    END IF
    CALL MPI_Comm_rank(comm_one_d(1), rank_S, ierr)
    CALL MPI_Comm_rank(comm_one_d(2), rank_F, ierr)
\endcode
Remarks:
<ol>
<li>The global rank of each processor, defined earlier in the <tt>main.f90</tt>, is denoted <tt>rank</tt>. It is usefull when writing an output that does not depend of a specific Fourier mode or a specific region (like the total kinetic energy).
<li>Comunicators have a dimension equal to two. The first dimension refers to a communication between subdomains of the finite element mesh for a fixed list of Fourier modes. The second dimension refers to a communication between different list of Fourier modes on a fixed subdomain of the finite element mesh.
</ol>

<li>The numerical stability of the computation is checked after each time iteration. It is done as follows.
\code
    !===Check divergence of fields
    IF (inputs%check_numerical_stability) THEN
       IF (inputs%type_pb=='nst' .OR. inputs%type_pb=='mhd' .OR. inputs%type_pb=='fhd') THEN
          norm = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un)
       ELSE 
          norm = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn)
       END IF
       IF (norm>1.d2) THEN
          CALL error_petsc('From my_post_processing: numerical unstability')
       END IF
    END IF
\endcode
To enable this feature, the following lines are required  in the data file:
\code
===Check numerical stability (true/false)
.t.
\endcode
It is set to false by default.
<li>Verbose and user's outputs are computed every inputs\%freq_en time iterations.
 The parameter inputs\%freq_en is set in the data as follows.
\code
===Frequency to write restart file
10
\endcode
This part of the subroutine is framed by the following lines code:
\code
    !===Put your postprocessing stuff here
    IF (MOD(it,inputs%freq_en) == 0) THEN
       (verbose and user ouputs)    
    END IF ! end freq_en
\endcode
<li>Generation of 2D and 3D visualization files of the approximated variables.
 This files can be visualized with the software Paraview. They are generated
 every inputs\%freq_plot time iterations. The parameter inputs\%freq_plot is
 set in the data as follows.
\code
===Frequency to create plots
10
\endcode
This part of the subroutine is framed by the following lines code:
\code
    IF (MOD(it,inputs%freq_plot) == 0) THEN
       (3D and 2D plot generation)
    END IF ! end freq_plot
\endcode
</ol>

Remark: The first two parts of the subroutine <tt>my_post_processing</tt> don't
 need to be modified. Only the section verbose/users ouputs and the section for
 the visualization files requires modifications depending of the outputs desired.

The verbose ouputs are the average computational time per time iterations,
 the CFL and the divergence of the velocity and magnetic fields. 
 Here are the required information to have these value computed 
 every inputs\%freq_en time  iterations.
<ol>
<li>In the verbose/user ouputs section, add the following lines:
\code
       !===Verbose
       CALL write_verbose(rank)
\endcode
<li>Add the following lines in the data file:
\code
===Verbose timing? (true/false)
.t.
===Verbose divergence? (true/false)
.t.
===Verbose CFL? (true/false)
.t.
\endcode
Note that computing the divergence does not imply to compute the CFL.
 These verbose outputs are independent of each other.
 They are all set to false by default.
</ol>

@section doc_SFEMaNS_main_tools_user_output Tools to compute user outputs

The subroutine <tt>my_post_processing</tt> uses the module <tt>tn_axi</tt>.
 It gives acces to four functions that can compute various quantity like L2
 norm, H1 norm, scalar product of vectors. Here is a description of these
 functions where we denote by norm a real number.
<ol>
<li>The function <tt>dot_product_SF</tt> computes the  scalar product of
 two vectors. It is called as follows:
\code
norm=dot_product_SF(comm_one_d, H_mesh, list_mode, Hn, Bn) 
\endcode
where comm_one_d and H_mesh are respectively the communicator and the
 finite element mesh where the vectors Hn and Bn are defined.
<li>The function <tt>norm_SF</tt> can compute the L2, H1, sH1 norm of a
 scalar or a vector. It can also computes the L2 norm of the divergence
 and the curl of a vector. It is called as follows:
\code
          norm = norm_SF(comm_one_d_ns, 'L2', vv_mesh, list_mode, un)
\endcode
where comm_one_d_ns and vv_mesh are respectively the communicator and
 the finite element mesh where un is defined. The characters 'L2' can
 be switch to H1, sH1, div and curl.
<li>The function norm_S can compute the L2, H1, sH1 norm of one Fourier
 component of a scalar or a vector. It can also computes the L2 norm
 of the divergence and the curl of one Fourier component of a vector.
 It is called as follows:
\code
          norm = norm_S(comm_one_d, 'L2', vv_mesh, list_mode(i:i), un(:,:,i:i))
\endcode
where comm_one_d_ns and vv_mesh are respectively the communicator and
 the finite element mesh where un is defined. The  characters 'L2' can
 be switch to H1, sH1, div and curl. The integer i is an integer between
 1 and SIZE(list_mode). It represents the label of the Fourier component considered.
<li>The function norm_S_L1_zero_mode computes the L1 norm of the Fourier
 mode zero of a scalar or a vector. It is called as follows:
\code
norm= norm_S_L1_zero_mode(comm_one_d_ns, pp_mesh, list_mode, pn)
\endcode
where comm_one_d_ns and pp_mesh are respectively the communicator and
 the finite element mesh where pn is defined.
</ol>

In addition, the file <tt>main.f90</tt> contains two subroutines that
 compute a drag force for problem with solid obstacle and the level set
 conservation for multiphase problem. The following gives information
 on these two subroutines.

The subroutine <tt>FORCES_AND_MOMENTS</tt> computes the drag force of
 a fluid driven by a vertical velocity field in the presence of a solid obstacle.
<ol>
<li>The presence of the solid obstacle is reported with a penalty function
 \f$\chi\f$. It is equal to one in the fluid and zero in the solid.
 The obstacle velocity is denoted \f$\bu_\text{obs}\f$. These functions
 are both defined in the file <tt>condlim.f90</tt>. They are only used
 if the following parameters are set to true in the data file:
\code
===Use penalty in NS domain (true/false)?
.t.
===Use nonzero velocity in solids (true/false)?
.t.
\endcode
These parameters are set to false by default. We refer to this
 <a href='doc_SFEMaNS_condlim.html'><code>section</code></a> for more information.
<li>The drag force can be shown to be equal to
 \f$-\frac{2}{\pi}\int \frac{(1-\chi)(\bu-\bu_\text{obs})\cdot \textbf{e}_z}{dt} \f$
 with \f$dt\f$ being the time step and \f$\textbf{e}_z\f$ being the unit vector
 in the vertical direction.
<li>The inputs are the time, the finite element mesh for the velocity field,
 the communicator for the Navier-Stokes domain, the list of Fourier mode and
 the velocity field.
<li>The output is written in the file <tt>fort.12</tt>.
<li>This subroutine is called as follows:
\code
          IF (inputs%if_compute_momentum_pseudo_force) THEN
             CALL FORCES_AND_MOMENTS(time,vv_mesh,comm_one_d_ns,list_mode,un)
          END IF
\endcode
The parameter inputs\%if_compute_momentum_pseudo_force is set to true in the data as follows:
\code
===Compute z momentum (true/false)?
.t.
\endcode
The default value is false.
</ol>

The subroutine <tt>compute_level_set_conservation</tt> computes the relative
 error on the level set conservation.
<ol>
<li>For each level set, it computes the term 
 \f$\frac{|\int_\Omega (\varphi_{t=tn}-\varphi_{t=0})|}{\int_\Omega \varphi_{t=0}}\f$.
<li> The inputs are the time, the finite element of the level set,
 the communicator for the Navier-Stokes domain, the list of Fourier mode and the level set.
<li>The output is written in the file <tt>fort.97</tt>.
<li>This subroutine is called as follows:
\code
          IF (inputs%if_level_set) THEN
             CALL compute_level_set_conservation(time,pp_mesh,comm_one_d_ns,list_mode,level_set)
          END IF
\endcode
</ol>


@section doc_SFEMaNS_main_tools_visualization Tools for visualization with Paraview

The code SFEMaNS can generate 3D and 2D visualization files of real
 valued scalar function and real valued vector function.  These files
 are generated with the subroutines <tt>vtu_3d</tt> and <tt>make_vtu_file_2D</tt>.
 We note that the 2D plots are generated for each Fourier component of the scalar/vector function. 

These subroutines generate files with the extension ".vtu" and ".pvd".
 The "vtu" files are generated every inputs\%freq_plot time iterations.
 We note that one "vtu" file is generated per processor in meridian section.
 They contains the informations of the variable to visualize at that time
 of the computation. One file with the extension ".pvd" is also generated.
 It contains an overview of all the "vtu" file generated. When opening the
 "pvd" file with the software Paraview, it allows to have access of the
 visualization of the variable on the whole domain at different times 
 without loading multiple files.

The file <tt>main.f90</tt> provided in the TEMPLATE directory of SFEMaNS
 generates 3D visualization files for the following variable:
<ol>
<li>the velocity field,
<li>the pressure,
<li>the level set(s),
<li>the density,
<li>the temperature,
<li>the magnetic field,
<li>the scalar potential.
</ol>
It can also generate 2D files of the above variables when uncommenting the following lines:
\code
!!$    !===VTU 2d======================================================================
!!$    CHARACTER(LEN=3)   :: st_mode
!!$    CHARACTER(LEN=200) :: header
!!$    CHARACTER(LEN=3)   :: name_of_field
\endcode
and the section between the following lines:
\code
       !===Generation of 2D plots for each Fourier mode (uncomment if wanted)
       (generation 2D files)
       !===End Generation of 2D plots for each Fourier mode (uncomment if wanted)
\endcode

The subroutines vtu_3d generates 3D visualization files.
<ol>
<li>It is called as follows:
\code
   CALL vtu_3d(var, 'var_mesh', 'File_name', 'Var', what, opt_it=it_plot)
\endcode
<li>The inputs are the following:
<ol>
<li><tt>var</tt> is the variable to plot (like un, pn, Hn, etc.).
<li><tt>var_mesh</tt> is the name of the mesh where var is defined
 (vv_mesh, pp_mesh, temp_mesh, H_mesh or phi_mesh).
<li><tt>File_name</tt> is the name of the pvd file generated.
 The name of the vtu files also start with File_name. However they
 also present information on the subdomain they represent (_S001, _S002, etc.).
<li><tt>'Var'</tt> is a character of three letters. It is used in
 Paraview when you want to display the value of the variable <tt>var</tt>.
<li><tt>what</tt> is equal to 'new' or 'old'. If it is 'new',
 the subroutine generates a pvd file else it only updates the pvd file.
<li><tt>opt_it</tt> is an optional integer. It allows to create vtu
 time every inputs\%freq_plot time iterations by completing the name
 of the vtu files with "_I001", "_I002", etc.
</ol>
<li>The pvd output is named "File_name.pvd". The vtu ouputs have
 the following format "File_name_S001_I001.vtu".
</ol>

The subroutine make_vtu_2D generates 2D visualization files.
<ol>
<li>It is called as follows:
\code
      CALL make_vtu_file_2D(comm_one_d_var(1), var_mesh, header, var(:,:,i), name_of_field, what, opt_it=it_plot)
\endcode
<li>The inputs are the following:
<ol>
<li><tt>comm_one_d_var(1)</tt> is the comunicator associated to the variable var. Only the first dimension is given. It represents communication between subsections of the finite element mesh (the Fourier component is fixed).
<li><tt>var_mesh</tt> is the name of the mesh where var is defined (vv_mesh, pp_mesh, temp_mesh, H_mesh or phi_mesh).
<li><tt>header</tt> is the name of the pvd file. The name of the vtu files also start with header. However they also present information on the subdomain they represent (_S001, _S002, etc.).
<li><tt>var(:,:,i)</tt> is Fourier component of the variable to plot. We note that the integer <tt>i</tt> is only the label of the Fourier component, the Fourier mode is equal to list_mode(i).
<li><tt>name_of_field</tt> is a character of three letters. It is used in Paraview when you want to display the value of the variable <tt>var</tt>.
<li><tt>what</tt> is equal to 'new' or 'old'. If it is new, the subroutine generates a pvd file else it only updates the pvd file.
<li><tt>opt_it</tt> is an optional integer. It allows to create vtu time every inputs\%freq_plot time iterations by completing the name of the vtu files with "_I001", "_I002", etc.
</ol>
<li>The pvd output is named "header.pvd". The vtu ouputs have the following format "header_S001_I001.vtu".
</ol>


@section doc_SFEMaNS_main_example Example

We give a description of the subroutine <tt>my_post_processing</tt> of the template file <tt>main.f90</tt>. This file can be found in the following directory: ($SFEMaNS_DIR)/TEMPLATE. We refer to the section  <a href='doc_example_physical_pb.html'><code>Examples on physical problems</code></a> for more examples.

<ol>
<li>Modules are added at the begining of the subroutine <tt>my_post_processing</tt>.
\code
    USE sub_plot
    USE chaine_caractere
    USE tn_axi
    USE boundary
    USE sft_parallele
    USE verbose
\endcode
<li>The declaration of local variables is done as follows.
\code
    IMPLICIT NONE
    INTEGER,                             INTENT(IN) :: it
    REAL(KIND=8)                                    :: err, norm
    INTEGER                                         :: i, it_plot
    CHARACTER(LEN=3)                                :: what
    INTEGER                                         :: rank_S, rank_F 
    INTEGER                                         :: rank_ns_S, rank_ns_F
    REAL(KIND=8), DIMENSION(pp_mesh%np, 2, SIZE(list_mode)) :: level_1
!!$    !===VTU 2d======================================================================
!!$    CHARACTER(LEN=3)   :: st_mode
!!$    CHARACTER(LEN=200) :: header
!!$    CHARACTER(LEN=3)   :: name_of_field
\endcode
<li>The ranks of each processors are defined.
\code
    !===Check ranks
    IF (vv_mesh%me /=0) THEN
       CALL MPI_Comm_rank(comm_one_d_ns(1), rank_ns_S, ierr)
       CALL MPI_Comm_rank(comm_one_d_ns(2), rank_ns_F, ierr)
    ELSE
       rank_ns_s = -1
       rank_ns_f = -1
    END IF
    CALL MPI_Comm_rank(comm_one_d(1), rank_S, ierr)
    CALL MPI_Comm_rank(comm_one_d(2), rank_F, ierr)
\endcode
<li>The numerical stability of the computation is checked after each time iteration.
\code
    !===Check divergence of fields
    IF (inputs%check_numerical_stability) THEN
       IF (inputs%type_pb=='nst' .OR. inputs%type_pb=='mhd' .OR. inputs%type_pb=='fhd') THEN
          norm = norm_SF(comm_one_d_ns, 'L2', vv_mesh, list_mode, un)
       ELSE 
          norm = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn)
       END IF
       IF (norm>1.d2) THEN
          CALL error_petsc('From my_post_processing: numerical unstability')
       END IF
    END IF
\endcode
<li>The verbose and user's outputs are computed every inputs\%freq_en time iterations. 
\code
    !===Put your postprocessing stuff here
    IF (MOD(it,inputs%freq_en) == 0) THEN
\endcode
<ol>
<li>The verbose outputs are computed as follows.
\code
       !===Verbose
       CALL write_verbose(rank)
\endcode
<li>The user ouputs that involve the Navier-Stokes equations and the temperature equation are computed if these equations are approximated. It is done by checking the type of problem approximated.
\code
       IF (inputs%type_pb=='nst' .OR. inputs%type_pb=='mhd' .OR. inputs%type_pb=='fhd') THEN
\endcode
<ol>
<li>If the parameter inputs\%if_compute_momentum_pseudo_force is set to true is the data, the drag force is computed. It is written in the file fort.12.
\code
          IF (inputs%if_compute_momentum_pseudo_force) THEN
             !===Compute the term -2/pi*integral((1-chi)*(u-u_solid)/dt)
             !===chi is the penalty function, u_solid the velocity of the solid, dt the time step
             !==Output written in the file fort.12
             CALL FORCES_AND_MOMENTS(time,vv_mesh,comm_one_d_ns,list_mode,un)
          END IF
\endcode
<li>The divergence of the velocity field is computed as follows. It is written is the file fort.31.
\code
          err = norm_SF(comm_one_d, 'div', vv_mesh, list_mode, un)
          norm = norm_SF(comm_one_d, 'H1', vv_mesh, list_mode, un)
          IF (rank == 0) THEN
             !===Divergence of velocity field
             WRITE(31,*) time, err/norm
          END IF
\endcode
<li>The L2 norm of each Fourier component m of the velocity field is written is the file fort.xxx with xxx=100+m. It uses the function norm_S as follows.
\code
          DO i=1,SIZE(list_mode)
             norm = norm_S(comm_one_d, 'L2', vv_mesh, list_mode(i:i), un(:,:,i:i))
             IF (rank_ns_S == 0) THEN
                !===L2 norm of Fourier mode list_mode(i) of velocity field un
                WRITE(100+list_mode(i),*) time, norm 
             END IF
          END DO
\endcode
<li>The L2 norm of the velocity is computed as follows. It is written in the file fort.98.
\code
          err = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un)
          norm = norm_SF(comm_one_d, 'sH1', vv_mesh, list_mode, un)
          IF (rank == 0) THEN
             WRITE(98,*) 'norm L2 of velocity', time, err
             WRITE(*,*) 'norm L2 of velocity', time, err
             WRITE(*,*) 'semi norm H1 of velocity', time, norm
          END IF
\endcode
<li>The L2 norm of the pressure is computed as follows.
\code
          err = norm_SF(comm_one_d, 'L2', pp_mesh, list_mode, pn)
          IF (rank == 0) THEN
             WRITE(*,*) 'norm L2 of pressure', time, err
          END IF
\endcode
<li>If inputs\%if_level_set is set to true (multiphase problem), the conservation of the level set is computed as follows. It is written in the file fort.97.
\code
          IF (inputs%if_level_set) THEN
             !===Compute the term integral(level_set-level_set_t=0)/integral(level_set_t=0)
             !===Output written in file fort.97
             CALL compute_level_set_conservation(time,pp_mesh,comm_one_d_ns,list_mode,level_set)
          END IF
\endcode
<li>If inputs\%if_temperature is set to true, the L2 norm of the temperature is computed as follows.
\code
          IF (inputs%if_temperature) THEN
             err = norm_SF(comm_one_d, 'L2', temp_mesh, list_mode, temperature)
             IF (rank == 0) THEN
                WRITE(*,*) 'norm L2 of temperature', time, err
             END IF
          END IF
\endcode
</ol>
<li>End of the user outputs involving the Navier-Stokes equations and the temperature equation.
\code
       END IF ! end nst or mhd or fhd
\endcode
<li>The user outputs involving the Maxwell equations are computed if these equations are approximated. It is done by checking the type of problem approximated.
\code
       IF (inputs%type_pb/='nst') THEN
\endcode
<ol>
<li>The L2 norm of the magnetic field \f$\bH\f$ is computed as follows. It is written in the file fort.41.
\code
          err = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn)
          IF (rank == 0) THEN
             !===L2 norm of magnetic field
             WRITE(41,*) time, err
          END IF
\endcode
<li>The divergence and the L2 norm of the magnetic field \f$\bB\f$ are computed as follows. They are written in the files fort.51 and fort.52.
\code
          err = norm_SF(comm_one_d, 'div', H_mesh, list_mode, Bn)
          norm = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Bn)
          IF (rank == 0) THEN
             !===L2 norm of div(Bn)
             WRITE(51,*) time, err, err/norm
             WRITE(52,*) time, err, norm
             WRITE(*,*) 'norm L2 of magnetic field', time, norm
          END IF
\endcode
<li>The L2 norm of each Fourier component m of the magnetic field \f$\bH\f$ is written is the file fort.xxx with xxx=200+m. It uses the function norm_S as follows. 
\code
          DO i=1,SIZE(list_mode)
             norm = norm_S(comm_one_d, 'L2', H_mesh, list_mode(i:i), Hn(:,:,i:i))
             IF (rank_S == 0) THEN
                !===L2 norm of Fourier mode list_mode(i) of magnetic field Hn
                WRITE(200+list_mode(i),*) time, norm
             END IF
          END DO
\endcode
</ol>
<li>End of the user ouputs involving the Maxwell equations.
\code
       END IF ! end /=nst
\endcode
</ol>
<li>End of the verbose and user ouputs.
\code
    END IF ! end freq_en
\endcode
<li>The visualization files are generated every inputs\%freq_plot time iterations. 
\code
    IF (MOD(it,inputs%freq_plot) == 0) THEN
\endcode
<ol>
<li>The inputs what and it_plot of the subroutines <tt>vtu_3d</tt> and <tt>make_vtu_file_2D</tt>  are defined as follows.
\code
       !===Plot whatever you want here
       IF (it==inputs%freq_plot) THEN
          what = 'new'
       ELSE
          what = 'old'
       END IF
       it_plot = it/inputs%freq_plot
\endcode
<li>The generation of the 3D files is done with the subroutine vtu_3d as follows.
\code
       !===Generation of 3D plots
       IF (inputs%if_level_set) THEN
          level_1=level_set(1,:,:,:)
          CALL vtu_3d(density, 'vv_mesh', 'Density', 'density', what, opt_it=it_plot)
          CALL vtu_3d(level_1, 'pp_mesh', 'Level_1', 'level_1', what, opt_it=it_plot)
          IF (SIZE(level_set,1).GE.2) THEN
             level_1=level_set(2,:,:,:)
             CALL vtu_3d(level_1, 'pp_mesh', 'Level_2', 'level_2', what, opt_it=it_plot)
          END IF
       END IF
       IF (inputs%type_pb/='mxw' .AND. inputs%type_pb/='mxx') THEN
          CALL vtu_3d(un, 'vv_mesh', 'Velocity', 'vel', what, opt_it=it_plot)
          CALL vtu_3d(pn, 'pp_mesh', 'Pressure', 'pre', what, opt_it=it_plot)
       END IF
       IF (inputs%if_temperature) THEN
          CALL vtu_3d(temperature, 'temp_mesh', 'Temperature', 'temp', what, opt_it=it_plot)
       END IF
       IF (inputs%type_pb/='nst') THEN
          CALL vtu_3d(Hn, 'H_mesh', 'MagField', 'mag', what, opt_it=it_plot)
          IF (inputs%nb_dom_phi>0) THEN
             CALL vtu_3d(phin, 'phi_mesh', 'ScalPot', 'phi', what, opt_it=it_plot)
          END IF
       END IF
       !==End generation of 3D plots
\endcode
<li>The generation of the 2D files is commented. It uses the subroutine make_vtu_file_2D.
\code
       !===Generation of 2D plots for each Fourier mode (uncomment if wanted)
       !===Proceed as follows to make 2D plots in the Fourier space (using Hn for instance)
       !===what = 'new' if first image, what= 'old' for later images  (CHARACTER(LEN=3)   :: what)
       !===WRITE(st_mode,'(I3)') list_mode(i)                         (CHARACTER(LEN=3)   :: st_mode)
       !===header = 'Hn_'//'_mode_'//trim(adjustl(st_mode))           (CHARACTER(LEN=200) :: header)
       !===name_of_field = 'Hn' (for instance)                        (CHARACTER(LEN=3)   :: name_of_field)
       !===CALL make_vtu_file_2D(comm_one_(1), H_mesh, header, Hn(:,:,i), name_of_field, what, opt_it=1)
!!$       IF (inputs%if_level_set) THEN
!!$          !===2D plots for each mode of the first level set
!!$          DO i = 1, m_max_c
!!$             WRITE(st_mode,'(I3)') list_mode(i)
!!$             header = 'Ln_'//'mode_'//trim(adjustl(st_mode)) 
!!$             name_of_field = 'Ln'
!!$             CALL make_vtu_file_2D(comm_one_d_ns(1), pp_mesh, header, level_1(:,:,i), name_of_field, what, opt_it=it_plot)
!!$             header = 'Dn_'//'mode_'//trim(adjustl(st_mode))
!!$             name_of_field = 'Dn'
!!$             CALL make_vtu_file_2D(comm_one_d_ns(1), vv_mesh, header, density(:,:,i), name_of_field, what, opt_it=it_plot)
!!$          END DO
!!$       END IF
!!$       IF (inputs%type_pb/='mxw' .AND. inputs%type_pb/='mxx') THEN
!!$          !===2D plots for each mode of the velocity field and the pressure
!!$          DO i = 1, m_max_c
!!$             WRITE(st_mode,'(I3)') list_mode(i)
!!$             header = 'Vn_'//'mode_'//trim(adjustl(st_mode))
!!$             name_of_field = 'Vn'
!!$             CALL make_vtu_file_2D(comm_one_d_ns(1), vv_mesh, header, un(:,:,i), name_of_field, what, opt_it=it_plot)
!!$             WRITE(st_mode,'(I3)') list_mode(i)
!!$             header = 'Pn_'//'mode_'//trim(adjustl(st_mode))
!!$             name_of_field = 'Pn'
!!$             CALL make_vtu_file_2D(comm_one_d_ns(1), pp_mesh, header, pn(:,:,i), name_of_field, what, opt_it=it_plot)
!!$          END DO
!!$       END IF
!!$       IF (inputs%if_temperature) THEN
!!$          !===2D plots for each mode of the temperature
!!$          DO i = 1, m_max_c
!!$             WRITE(st_mode,'(I3)') list_mode(i)
!!$             header = 'Tn_'//'_mode_'//trim(adjustl(st_mode))
!!$             name_of_field = 'Tn'
!!$             CALL make_vtu_file_2D(comm_one_d(1), temp_mesh, header, temperature(:,:,i), name_of_field, what, opt_it=it_plot)
!!$          END DO
!!$       END IF
!!$       IF (inputs%type_pb/='nst') THEN
!!$          !===2D plots for each mode of the magnetic field and the scalar potential
!!$          DO i = 1, m_max_c
!!$             WRITE(st_mode,'(I3)') list_mode(i)
!!$             header = 'Hn_'//'_mode_'//trim(adjustl(st_mode))
!!$             name_of_field = 'Hn'
!!$             CALL make_vtu_file_2D(comm_one_d(1), H_mesh, header, Hn(:,:,i), name_of_field, what, opt_it=it_plot)
!!$             IF (inputs%nb_dom_phi>0) THEN
!!$                WRITE(st_mode,'(I3)') list_mode(i)
!!$                header = 'Phin_'//'_mode_'//trim(adjustl(st_mode))
!!$                name_of_field = 'Phin'
!!$                CALL make_vtu_file_2D(comm_one_d(1), phi_mesh, header, phin(:,:,i), name_of_field, what, opt_it=it_plot)
!!$             END IF
!!$          END DO
!!$       END IF
       !===End Generation of 2D plots for each Fourier mode (uncomment if wanted)
\endcode
</ol>
<li>End of the section that generates visualization files.
\code
    END IF ! end freq_plot
\endcode
</ol>

Remark: If more than 100 Fourier modes are used, a conflict arises between different energy files. Indeed, the L2 norm of the Fourier component i of the velocity is  written in the files fort.xxx with xxx=100+i. The same is done for the magnetic field with the files fort.yyy with yyy=200+i. This problem can be overcome by switching 200 to a larger number. One can also defines proper name for the energy outputs file.

 */
