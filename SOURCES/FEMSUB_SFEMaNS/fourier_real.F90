MODULE fourier_to_real_for_vtu
  USE my_util
  USE dyn_line
  USE def_type_mesh
  USE sub_plot
  USE input_data
#include "petsc/finclude/petsc.h"
  USE petsc
  PUBLIC :: sfemans_initialize_postprocessing
  PUBLIC :: transfer_fourier_to_real, make_vtu_file_2D
  PUBLIC :: vtu_3d
  PUBLIC :: compute_rot_h
  TYPE(petsc_csr_LA),  PUBLIC  :: vizu_grad_phi_LA
  TYPE(petsc_csr_LA),  PUBLIC  :: vizu_rot_h_LA
  TYPE(petsc_csr_LA),  PUBLIC  :: vizu_rot_u_LA
  TYPE(mesh_type), POINTER     :: mesh_grad_phi
  TYPE(mesh_type), POINTER     :: mesh_rot_h
  PRIVATE
  TYPE(mesh_type), TARGET      :: vv_mesh_3d
  TYPE(mesh_type), TARGET      :: pp_mesh_3d
  TYPE(mesh_type), TARGET      :: H_mesh_3d
  TYPE(mesh_type), TARGET      :: phi_mesh_3d
  TYPE(mesh_type), TARGET      :: temp_mesh_3d
  TYPE(mesh_type), TARGET      :: conc_mesh_3d
  TYPE(mesh_type), TARGET      :: vv_mesh_p2_iso_p4
  TYPE(mesh_type), TARGET      :: vv_local_mesh
  TYPE(mesh_type), TARGET      :: pp_local_mesh
  TYPE(mesh_type), TARGET      :: H_local_mesh
  TYPE(mesh_type), TARGET      :: phi_local_mesh
  TYPE(mesh_type), TARGET      :: temp_local_mesh
  TYPE(mesh_type), TARGET      :: conc_local_mesh
  TYPE(dyn_int_line), DIMENSION(:), POINTER :: vv_loc_to_glob
  TYPE(dyn_int_line), DIMENSION(:), POINTER :: pp_loc_to_glob
  TYPE(dyn_int_line), DIMENSION(:), POINTER :: H_loc_to_glob
  TYPE(dyn_int_line), DIMENSION(:), POINTER :: phi_loc_to_glob
  TYPE(dyn_int_line), DIMENSION(:), POINTER :: temp_loc_to_glob
  TYPE(dyn_int_line), DIMENSION(:), POINTER :: conc_loc_to_glob
  INTEGER,            DIMENSION(:), POINTER :: fourier_list_mode
  INTEGER                                   :: fourier_nb_angle
  INTEGER                                   :: fourier_width
  INTEGER                                   :: fourier_rank
  INTEGER                                   :: fourier_nb_procs
  LOGICAL,  SAVE                            :: if_first_grad_phi=.TRUE.
  LOGICAL,  SAVE                            :: if_first_rot_h=.TRUE.
  LOGICAL,  SAVE                            :: if_first_rot_u=.TRUE.
  INTEGER, DIMENSION(:), POINTER            :: vizu_list_mode
  MPI_Comm                                  :: fourier_communicator
  MPI_Comm, DIMENSION(2)                    :: vizu_grad_curl_comm
  Mat                                       :: mat_grad_phi
  Mat                                       :: mat_rot_h
  Mat                                       :: mat_rot_u
  Vec                                       :: vx_phi, vb_phi, vx_phi_ghost
  KSP                                       :: ksp_grad_phi
  KSP                                       :: ksp_rot_h
  KSP                                       :: ksp_rot_u
CONTAINS

  SUBROUTINE sfemans_initialize_postprocessing(comm_one_d, vv_mesh_in, pp_mesh, H_mesh, phi_mesh, temp_mesh, &
       conc_mesh,list_mode_in, opt_nb_plane)
    USE input_data
    USE sub_plot
    IMPLICIT NONE
    TYPE(mesh_type), POINTER                   :: vv_mesh
    TYPE(mesh_type), TARGET,        INTENT(IN) :: vv_mesh_in
    TYPE(mesh_type),                INTENT(IN) :: pp_mesh
    TYPE(mesh_type), TARGET                    :: H_mesh
    TYPE(mesh_type), TARGET                    :: phi_mesh
    TYPE(mesh_type),                INTENT(IN) :: temp_mesh
    TYPE(mesh_type),                INTENT(IN) :: conc_mesh
    INTEGER,          DIMENSION(:), TARGET :: list_mode_in
    INTEGER, OPTIONAL,              INTENT(IN) :: opt_nb_plane
    INTEGER                                    :: code
    PetscErrorCode                             :: ierr
    PetscMPIInt                                :: rank, nb_procs, petsc_rank
    MPI_Comm, DIMENSION(2)                     :: comm_one_d
    CALL MPI_Comm_rank(PETSC_COMM_WORLD,petsc_rank,ierr)
    CALL MPI_Comm_rank(comm_one_d(2), rank, ierr)
    CALL MPI_Comm_size(comm_one_d(2), nb_procs, ierr)
    fourier_communicator = comm_one_d(2)
    fourier_rank = rank
    fourier_nb_procs = nb_procs
    ALLOCATE(fourier_list_mode(SIZE(list_mode_in)))
    fourier_list_mode = list_mode_in
    !IF ((.NOT.inputs%if_just_processing) .AND. (inputs%freq_plot .GT. inputs%nb_iteration)) THEN
    !   !===Not preparation of vizu meshes
    !   RETURN
    !END IF
    IF (inputs%type_fe_velocity .GE. 3) THEN
       !===divide each cell of vv_mesh into 4 cells
       CALL divide_mesh_into_four_subcells(vv_mesh_in,vv_mesh_p2_iso_p4)
       vv_mesh => vv_mesh_p2_iso_p4
    ELSE
       vv_mesh => vv_mesh_in
    END IF
    CALL divide_mesh(comm_one_d(2), vv_mesh, vv_local_mesh, vv_loc_to_glob)
    CALL divide_mesh(comm_one_d(2), pp_mesh, pp_local_mesh, pp_loc_to_glob)
    CALL divide_mesh(comm_one_d(2), H_mesh,  H_local_mesh,  H_loc_to_glob)
    CALL divide_mesh(comm_one_d(2), phi_mesh, phi_local_mesh, phi_loc_to_glob)
    CALL divide_mesh(comm_one_d(2), temp_mesh, temp_local_mesh, temp_loc_to_glob)
    CALL divide_mesh(comm_one_d(2), conc_mesh, conc_local_mesh, conc_loc_to_glob)
    CALL prepare_3d_grids(vv_mesh,  vv_local_mesh,  vv_mesh_3d, petsc_rank, opt_nb_plane)
    CALL prepare_3d_grids(pp_mesh,  pp_local_mesh,  pp_mesh_3d, petsc_rank, opt_nb_plane)
    CALL prepare_3d_grids(H_mesh,   H_local_mesh,   H_mesh_3d, petsc_rank,  opt_nb_plane)
    CALL prepare_3d_grids(phi_mesh, phi_local_mesh, phi_mesh_3d, petsc_rank, opt_nb_plane)
    CALL prepare_3d_grids(temp_mesh, temp_local_mesh, temp_mesh_3d, petsc_rank, opt_nb_plane)
    CALL prepare_3d_grids(conc_mesh, conc_local_mesh, conc_mesh_3d, petsc_rank, opt_nb_plane)
    mesh_rot_h => H_mesh
    mesh_grad_phi => phi_mesh
    CALL MPI_COMM_DUP(comm_one_d(1),vizu_grad_curl_comm(1) , code)
    CALL MPI_COMM_DUP(comm_one_d(2),vizu_grad_curl_comm(2) , code)
    vizu_list_mode => list_mode_in
  END SUBROUTINE sfemans_initialize_postprocessing

  SUBROUTINE prepare_3d_grids(mesh, local_mesh, mesh_3d, rank, opt_nb_plane)
    USE dyn_line
    USE def_type_mesh
    USE input_data
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type),                INTENT(IN) :: mesh
    TYPE(mesh_type),                INTENT(IN) :: local_mesh
    TYPE(mesh_type),                INTENT(OUT):: mesh_3d
    INTEGER, OPTIONAL,              INTENT(IN) :: opt_nb_plane
    INTEGER,                        INTENT(IN) :: rank
    LOGICAL,      DIMENSION(:),     POINTER    :: virgin
    INTEGER,      DIMENSION(:),     POINTER    :: renumber
    REAL(KIND=8)                               :: pi, dtheta, theta
    INTEGER                                    :: i, k, m, n, nb_angle, count, dim, np1, &
         gap, m_max, m_max_pad

    !===Deal with empty mesh
    IF (mesh%me==0) THEN
       mesh_3d%me = 0
       mesh_3d%np = 0
       mesh_3d%gauss%n_w = 0
       ALLOCATE(mesh_3d%jj(0,0))
       ALLOCATE(mesh_3d%rr(0,0))
       RETURN
    END IF

    !===Compute gap
    IF (inputs%select_mode) THEN
       IF (SIZE(inputs%list_mode_lect)==1) THEN
          gap = inputs%list_mode_lect(1) + 1
       ELSE
          gap = inputs%list_mode_lect(2) - inputs%list_mode_lect(1)
       END IF
    ELSE
       gap = 1
    END IF
    fourier_width = (inputs%m_max/fourier_nb_procs)*gap

    !===Compute nb_angle
    m_max = gap*inputs%m_max
    IF (PRESENT(opt_nb_plane)) THEN
       IF (opt_nb_plane> 2*m_max-1) THEN
          m_max_pad = (opt_nb_plane+1)/2
       ELSE
          m_max_pad = m_max
       END IF
    ELSE
       m_max_pad = m_max
    END IF
    nb_angle = 2*m_max_pad-1
    fourier_nb_angle = nb_angle

    !===Create 3d mesh
    dim = 3
    pi = ACOS(-1.d0)
    dtheta = 2*pi/nb_angle
    IF (mesh%gauss%n_w==3) THEN
       mesh_3d%np = local_mesh%np*nb_angle
       mesh_3d%me = local_mesh%me*nb_angle
       ALLOCATE(mesh_3d%rr(dim,mesh_3d%np))
       count = 0
       DO k = 1, nb_angle
          theta = (k-1)*dtheta
          DO i = 1, local_mesh%np
             count = count + 1
             mesh_3d%rr(1,count) = local_mesh%rr(1,i)*COS(theta)
             mesh_3d%rr(2,count) = local_mesh%rr(1,i)*SIN(theta)
             mesh_3d%rr(3,count) = local_mesh%rr(2,i)
          ENDDO
       END DO
       ALLOCATE(mesh_3d%jj(6,mesh_3d%me)) !===P1 3d VTK wedge 13
       mesh_3d%gauss%n_w = 6
       count = 0
       DO k = 1, nb_angle-1
          DO m = 1, local_mesh%me
             count = count + 1
             mesh_3d%jj(1:3,count) = local_mesh%jj(1:3,m)+(k-1)*local_mesh%np
             mesh_3d%jj(4:6,count) = local_mesh%jj(1:3,m)+k*local_mesh%np
          END DO
       END DO
       DO m = 1, local_mesh%me
          count = count + 1
          mesh_3d%jj(1:3,count) = local_mesh%jj(1:3,m)+(k-1)*local_mesh%np
          mesh_3d%jj(4:6,count) = local_mesh%jj(1:3,m)
       END DO
    ELSE IF (mesh%gauss%n_w==6) THEN
       ALLOCATE(virgin(local_mesh%np), renumber(local_mesh%np))
       virgin = .TRUE.
       renumber = -1
       count = 0
       DO m = 1, local_mesh%me
          DO n = 1, 3
             i = local_mesh%jj(n,m)
             IF (.NOT.virgin(i)) CYCLE
             virgin(i) = .FALSE.
             count = count + 1 !===Count P1 nodes on local_mesh
             renumber(i) = count
          END DO
       END DO
       np1 = count
       mesh_3d%np = (local_mesh%np+np1)*nb_angle
       mesh_3d%me = local_mesh%me*nb_angle
       ALLOCATE(mesh_3d%rr(dim,mesh_3d%np))
       ALLOCATE(mesh_3d%jj(15,mesh_3d%me)) !===P2 3d VTK wedge 26
       mesh_3d%gauss%n_w = 15
       count = 0
       DO k = 1, nb_angle-1
          DO m = 1, local_mesh%me
             count = count + 1
             mesh_3d%jj(1:3,count) = local_mesh%jj(1:3,m) + (k-1)*(local_mesh%np+np1)
             mesh_3d%jj(4:6,count) = local_mesh%jj(1:3,m) + k*(local_mesh%np+np1)
             mesh_3d%jj(7,count)   = local_mesh%jj(6,m)   + (k-1)*(local_mesh%np+np1)
             mesh_3d%jj(8,count)   = local_mesh%jj(4,m)   + (k-1)*(local_mesh%np+np1)
             mesh_3d%jj(9,count)   = local_mesh%jj(5,m)   + (k-1)*(local_mesh%np+np1)
             mesh_3d%jj(10,count)  = local_mesh%jj(6,m)   + k*(local_mesh%np+np1)
             mesh_3d%jj(11,count)  = local_mesh%jj(4,m)   + k*(local_mesh%np+np1)
             mesh_3d%jj(12,count)  = local_mesh%jj(5,m)   + k*(local_mesh%np+np1)
             mesh_3d%jj(13:15,count) = renumber(local_mesh%jj(1:3,m)) &
                  + (k-1)*(local_mesh%np+np1) + local_mesh%np
          END DO
       END DO
       k = nb_angle
       DO m = 1, local_mesh%me
          count = count + 1
          mesh_3d%jj(1:3,count) = local_mesh%jj(1:3,m) + (k-1)*(local_mesh%np+np1)
          mesh_3d%jj(4:6,count) = local_mesh%jj(1:3,m)
          mesh_3d%jj(7,count)   = local_mesh%jj(6,m)   + (k-1)*(local_mesh%np+np1)
          mesh_3d%jj(8,count)   = local_mesh%jj(4,m)   + (k-1)*(local_mesh%np+np1)
          mesh_3d%jj(9,count)   = local_mesh%jj(5,m)   + (k-1)*(local_mesh%np+np1)
          mesh_3d%jj(10,count)  = local_mesh%jj(6,m)
          mesh_3d%jj(11,count)  = local_mesh%jj(4,m)
          mesh_3d%jj(12,count)  = local_mesh%jj(5,m)
          mesh_3d%jj(13:15,count) = renumber(local_mesh%jj(1:3,m)) &
               + (k-1)*(local_mesh%np+np1) + local_mesh%np
       END DO
       !===Create 3d mesh
       DO k = 1, nb_angle
          theta = (k-1)*dtheta
          DO i = 1, local_mesh%np
             n = (k-1)*(local_mesh%np+np1) + i
             mesh_3d%rr(1,n) = local_mesh%rr(1,i)*COS(theta)
             mesh_3d%rr(2,n) = local_mesh%rr(1,i)*SIN(theta)
             mesh_3d%rr(3,n) = local_mesh%rr(2,i)
             IF (virgin(i)) CYCLE !===This is a vertex
             n = (k-1)*(local_mesh%np+np1) + local_mesh%np + renumber(i)
             mesh_3d%rr(1,n) = local_mesh%rr(1,i)*COS(theta+dtheta/2)
             mesh_3d%rr(2,n) = local_mesh%rr(1,i)*SIN(theta+dtheta/2)
             mesh_3d%rr(3,n) = local_mesh%rr(2,i)
          ENDDO
       END DO
       !===Cleanup
       DEALLOCATE(virgin,renumber)
    ELSE IF (mesh%gauss%n_w==10) THEN
       IF (rank==0) THEN
          WRITE(*,*) 'VTK files not done for P3 meshes'
       END IF
    ELSE
       CALL error_petsc('Bug in prepare_3d_grids: mesh%gauss%n_w= not correct')
    END IF
  END SUBROUTINE prepare_3d_grids

  SUBROUTINE vtu_3d(comm_one_d,field, name_of_mesh, header, name_of_field, what, opt_it, opt_grad_curl, opt_2D, opt_mesh_in)
    USE dyn_line
    USE def_type_mesh
    USE vtk_viz
    USE sft_parallele
    USE input_data
    USE my_util
    USE tn_axi
    IMPLICIT NONE
    TYPE(mesh_type),  OPTIONAL                   :: opt_mesh_in
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN), TARGET :: field
    REAL(KIND=8), DIMENSION(:,:,:),   POINTER    :: field_in
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE,  TARGET     :: field_tmp
    CHARACTER(*),                     INTENT(IN) :: name_of_mesh, header, name_of_field, what
    INTEGER, OPTIONAL,                INTENT(IN) :: opt_it
    CHARACTER(*), OPTIONAL,           INTENT(IN) :: opt_grad_curl
    LOGICAL, OPTIONAL,                INTENT(IN) :: opt_2D
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE  :: field_viz
    !===VTU 2d======================================================================
    CHARACTER(LEN=3)   :: st_mode
    CHARACTER(LEN=200) :: header_2D
    CHARACTER(LEN=3)   :: name_of_field_2D
    INTEGER            :: i
    REAL(KIND=8)       :: norm
    INTEGER            :: rank_S, rank_F, nb_S
    !===Declare PETSC===============================================================
    PetscErrorCode :: ierr
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d


    IF (PRESENT(opt_grad_curl)) THEN !=== compute rot_u for P3
       IF ((opt_grad_curl == 'curl_u').AND.(PRESENT(opt_mesh_in))) THEN
          IF (opt_mesh_in%gauss%n_w==10 .AND. inputs%type_fe_velocity==3) THEN
             ALLOCATE(field_viz(SIZE(field,1), SIZE(field,2), SIZE(field,3)))
             CALL compute_rot_u(field, field_viz, opt_mesh_in) !=== only opt_mesh_in
          END IF
       END IF
    END IF

    IF (PRESENT(opt_mesh_in)) THEN
       IF (opt_mesh_in%gauss%n_w==10 .AND. inputs%type_fe_velocity==3) THEN
          ALLOCATE(field_tmp(vv_mesh_p2_iso_p4%np,SIZE(field,2),SIZE(field,3)))
          DO i = 1, SIZE(field,3)
             IF (.NOT.(PRESENT(opt_grad_curl))) THEN !=== visualization velocity
                CALL interpolate(field(:,:,i),opt_mesh_in,field_tmp(:,:,i),vv_mesh_p2_iso_p4) !===field_tmp=u
             ELSE IF (opt_grad_curl == 'curl_u') THEN !=== visualization curl
                CALL interpolate(field_viz(:,:,i),opt_mesh_in,field_tmp(:,:,i),vv_mesh_p2_iso_p4) !===field_tmp=curl_u
             ELSE
                CALL error_Petsc('Bug in vtu_3d: P3 and wrong opt_grad_curl')
             END IF
          END DO
          IF (ALLOCATED(field_viz)) DEALLOCATE(field_viz)
          field_in => field_tmp !=== if P3: links to u or rot_u on vv_mesh_p2_iso_p
       ELSE
          field_in => field !=== not P3 but (present(opt_mesh(in)))
       END IF
    ELSE
       IF (name_of_mesh == 'vv_mesh' .AND. inputs%type_fe_velocity==3) THEN
          CALL error_Petsc('Bug in vtu_3d: P3 for velocity BUT opt_mesh_in not present')
       ELSE
          field_in => field !=== not present (opt_mesh_in)
       ENDIF
    END IF

    CALL MPI_Comm_rank(comm_one_d(1), rank_S, ierr)
    CALL MPI_Comm_rank(comm_one_d(2), rank_F, ierr)
    CALL MPI_Comm_size(comm_one_d(1), nb_S, ierr)

    IF (PRESENT(opt_grad_curl)) THEN
       IF (opt_grad_curl == 'grad') THEN !===grad phi
          ALLOCATE(field_viz(SIZE(field_in,1), 3*SIZE(field_in,2), SIZE(field_in,3)))
          CALL compute_grad_phi(field_in, field_viz)
       ELSE IF ((opt_grad_curl == 'curl_h').AND.(PRESENT(opt_mesh_in))) THEN !===curl h
          ALLOCATE(field_viz(SIZE(field_in,1), SIZE(field_in,2), SIZE(field_in,3)))
          CALL compute_rot_h(field_in, field_viz)
          norm = norm_SF(comm_one_d, 'L2', opt_mesh_in, vizu_list_mode, field_viz)
          IF (rank_S == 0) THEN
             WRITE(102,*) norm
          END IF
       ELSE IF ((opt_grad_curl == 'curl_u').AND.(PRESENT(opt_mesh_in))) THEN !===curl u
          IF (opt_mesh_in%gauss%n_w==10 .AND. inputs%type_fe_velocity==3) THEN !===P3
             ALLOCATE(field_viz(SIZE(field_in,1), SIZE(field_in,2), SIZE(field_in,3)))
             field_viz = field_in  !===field_in = curl u
          ELSE
             ALLOCATE(field_viz(SIZE(field_in,1), SIZE(field_in,2), SIZE(field_in,3)))
             CALL compute_rot_u(field_in, field_viz, opt_mesh_in)  !===field_in = u
          END IF
       ELSE
          CALL error_petsc('Bug in vtu_3d: name_of_opt_grad_curl is not correct')
       END IF
    ELSE
       ALLOCATE(field_viz(SIZE(field_in, 1), SIZE(field_in,2), SIZE(field_in,3)))
       field_viz = field_in
    END IF

    IF (PRESENT(opt_it)) THEN
       CALL vtu_3d_base(field_viz, name_of_mesh, header, name_of_field, what, opt_it)
    ELSE
       CALL vtu_3d_base(field_viz, name_of_mesh, header, name_of_field, what)
    END IF

    IF (PRESENT(opt_2D)) THEN
       IF (opt_2D) THEN
          DO i = 1, SIZE(vizu_list_mode)
             WRITE(st_mode,'(I3)') vizu_list_mode(i)  !=== (CHARACTER(LEN=3) :: st_mode)
             header_2D = 'J_'//'mode_'//trim(adjustl(st_mode)) !=== (CHARACTER(LEN=200) :: header)
             name_of_field_2D = 'J' !===(for instance) (CHARACTER(LEN=3)   :: name_of_field)
             IF (PRESENT(opt_it)) THEN
                CALL make_vtu_file_2D(vizu_grad_curl_comm(1), mesh_rot_h, &
                     header_2D, field_viz(:,:,i), name_of_field_2D, what, opt_it=opt_it)
             ELSE
                CALL make_vtu_file_2D(vizu_grad_curl_comm(1), mesh_rot_h, &
                     header_2D, field_viz(:,:,i), name_of_field_2D, what)
             END IF
          END DO
       END IF
    END IF

    DEALLOCATE(field_viz)
    IF (ALLOCATED(field_tmp)) DEALLOCATE(field_tmp)
  END SUBROUTINE vtu_3d

  SUBROUTINE vtu_3d_base(field_in, name_of_mesh, header, name_of_field, what, opt_it)
    USE dyn_line
    USE def_type_mesh
    USE vtk_viz
    USE sft_parallele
    USE input_data
    USE my_util
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN) :: field_in
    CHARACTER(*),                     INTENT(IN) :: name_of_mesh, header, name_of_field, what
    INTEGER, OPTIONAL,                INTENT(IN) :: opt_it
    TYPE(mesh_type),                  POINTER    :: mesh_3d
    TYPE(mesh_type),                  POINTER    :: local_mesh
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE  :: field_out_FFT
    REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE  :: field_out
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE  :: field_out_xyz_FFT
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE  :: dist_field
    LOGICAL,      DIMENSION(:),       POINTER    :: virgin
    INTEGER,      DIMENSION(:),       POINTER    :: renumber
    INTEGER,      DIMENSION(:),       POINTER    :: list_mode
    TYPE(dyn_int_line), DIMENSION(:), POINTER    :: loc_to_glob
    INTEGER                                      :: width, mode_min
    INTEGER, DIMENSION(1)                        :: loc
    REAL(KIND=8)                                 :: pi, dtheta, theta
    INTEGER                                      :: i, k, m, n, np_b, np_e, np_loc, np_loc_max, &
         type, nb_angle, count, dim, np1, nb_procs, rank, it

    !===Initialization
    width    = fourier_width
    nb_procs = fourier_nb_procs
    rank     = fourier_rank
    nb_angle = fourier_nb_angle
    list_mode => fourier_list_mode
    ALLOCATE(loc_to_glob(nb_procs))
    IF (name_of_mesh=='vv_mesh') THEN
       DO n = 1, nb_procs
          loc_to_glob(n)%DIL => vv_loc_to_glob(n)%DIL
       END DO
       mesh_3d     => vv_mesh_3d
       local_mesh  => vv_local_mesh
    ELSE IF (name_of_mesh=='pp_mesh') THEN
       DO n = 1, nb_procs
          loc_to_glob(n)%DIL => pp_loc_to_glob(n)%DIL
       END DO
       mesh_3d     => pp_mesh_3d
       local_mesh  => pp_local_mesh
    ELSE IF (name_of_mesh=='H_mesh') THEN
       DO n = 1, nb_procs
          loc_to_glob(n)%DIL => H_loc_to_glob(n)%DIL
       END DO
       mesh_3d     => H_mesh_3d
       local_mesh  => H_local_mesh
    ELSE IF (name_of_mesh=='phi_mesh') THEN
       DO n = 1, nb_procs
          loc_to_glob(n)%DIL => phi_loc_to_glob(n)%DIL
       END DO
       mesh_3d     => phi_mesh_3d
       local_mesh  => phi_local_mesh
    ELSE IF (name_of_mesh=='temp_mesh') THEN
       DO n = 1, nb_procs
          loc_to_glob(n)%DIL => temp_loc_to_glob(n)%DIL
       END DO
       mesh_3d     => temp_mesh_3d
       local_mesh  => temp_local_mesh
       !SB 06/05/2022
    ELSE IF (name_of_mesh=='conc_mesh') THEN
       DO n = 1, nb_procs
          loc_to_glob(n)%DIL => conc_loc_to_glob(n)%DIL
       END DO
       mesh_3d     => conc_mesh_3d
       local_mesh  => conc_local_mesh
    ELSE
       CALL error_petsc('Bug in vtu_3d: name_of_mesh is not correct')
    END IF

    !===Compute np_loc_max
    np_loc_max = -1
    DO n = 1, nb_procs
       np_loc_max = MAX(np_loc_max,MAXVAL(loc_to_glob(n)%DIL))
    END DO
    ALLOCATE(dist_field(nb_procs*np_loc_max,SIZE(field_in,2),width))

    !===Navier subdomain is divided into fourier_rank parts
    !===This loop ditributes fourier modes to the sub-meshes (1 to fourier_rank)
    dist_field = 0.d0
    mode_min = rank*width
    DO i = 1, width
       IF (MINVAL(ABS(list_mode - (mode_min+i-1)))/=0) CYCLE
       loc = MINLOC(ABS(list_mode - (mode_min+i-1)))
       DO type = 1, SIZE(field_in,2)
          DO n = 1, nb_procs
             np_b = (n-1)*np_loc_max + 1
             np_loc = SIZE(loc_to_glob(n)%DIL)
             np_e = np_b + np_loc - 1
             dist_field(np_b:np_e,type,i) = field_in(loc_to_glob(n)%DIL,type,loc(1))
          END DO
       END DO
    END DO
    !===Compute field_out_FFT in real space on local_mesh
    CALL FFT_PAR_REAL(fourier_communicator, dist_field, field_out_FFT, opt_nb_plane=nb_angle)

    ALLOCATE(field_out_xyz_FFT(SIZE(field_out_FFT,1), SIZE(field_out_FFT,2), SIZE(field_out_FFT,3)))
    IF (SIZE(field_in,2)==2) THEN
       field_out_xyz_FFT=field_out_FFT
    ELSE
       pi = ACOS(-1.d0)
       dtheta = 2*pi/nb_angle
       DO k = 1, nb_angle
          theta = (k-1)*dtheta
          field_out_xyz_FFT(k,1,:) = field_out_FFT(k,1,:)*COS(theta)-field_out_FFT(k,2,:)*SIN(theta)
          field_out_xyz_FFT(k,2,:) = field_out_FFT(k,1,:)*SIN(theta)+field_out_FFT(k,2,:)*COS(theta)
          field_out_xyz_FFT(k,3,:) = field_out_FFT(k,3,:)
       END DO
    END IF

    IF (SIZE(field_in,2)==2) THEN
       dim = 1
    ELSE IF (SIZE(field_in,2)==6) THEN
       dim = 3
    ELSE
       CALL error_petsc('Bug in vtu_3d: SIZE(field_in,2) not correct')
    END IF
    IF (nb_angle /= SIZE(field_out_xyz_FFT,1)) THEN
       CALL error_petsc('Bug in vtu_3d: nb_angle /= SIZE(field_out_xyz_FFT,1')
    END IF

    !===Reconstruct field_out in real space on mesh_3d
    pi = ACOS(-1.d0)
    dtheta = 2*pi/nb_angle
    IF (local_mesh%gauss%n_w==3) THEN
       ALLOCATE(field_out(dim,mesh_3d%np))
       count = 0
       DO k = 1, nb_angle
          DO i = 1, local_mesh%np
             count = count + 1
             field_out(:,count) = field_out_xyz_FFT(k,:,i)
          ENDDO
       END DO

    ELSE IF (local_mesh%gauss%n_w==6) THEN
       ALLOCATE(virgin(local_mesh%np), renumber(local_mesh%np))
       virgin = .TRUE.
       renumber = -1
       count = 0
       DO m = 1, local_mesh%me
          DO n = 1, 3
             i = local_mesh%jj(n,m)
             IF (.NOT.virgin(i)) CYCLE
             virgin(i) = .FALSE.
             count = count + 1 !===Count P1 nodes on local_mesh
             renumber(i) = count
          END DO
       END DO
       np1 = count

       ALLOCATE(field_out(dim,mesh_3d%np))
       DO k = 1, nb_angle
          DO i = 1, local_mesh%np
             n = (k-1)*(local_mesh%np+np1) + i
             field_out(:,n) = field_out_xyz_FFT(k,:,i)
          END DO
       END DO

       DO i = 1, local_mesh%np
          IF (virgin(i)) CYCLE !===This is a vertex
          DO k = 1, nb_angle - 2
             n = (k-1)*(local_mesh%np+np1) + local_mesh%np + renumber(i)
             field_out(:,n) = (3.d0/8)*field_out_xyz_FFT(k,:,i) &
                  + (3.d0/4)*field_out_xyz_FFT(k+1,:,i) - (1.d0/8)*field_out_xyz_FFT(k+2,:,i)
          END DO
       END DO
       k = nb_angle - 1
       DO i = 1, local_mesh%np
          IF (virgin(i)) CYCLE !===This is a vertex
          n = (k-1)*(local_mesh%np+np1) + local_mesh%np + renumber(i)
          field_out(:,n) = (3.d0/8)*field_out_xyz_FFT(k,:,i) &
               + (3.d0/4)*field_out_xyz_FFT(k+1,:,i) - (1.d0/8)*field_out_xyz_FFT(1,:,i)
       END DO
       k = nb_angle
       DO i = 1, local_mesh%np
          IF (virgin(i)) CYCLE !===This is a vertex
          n = (k-1)*(local_mesh%np+np1) + local_mesh%np + renumber(i)
          field_out(:,n) = (3.d0/8)*field_out_xyz_FFT(k,:,i) &
               + (3.d0/4)*field_out_xyz_FFT(1,:,i) - (1.d0/8)*field_out_xyz_FFT(2,:,i)
       END DO
       DEALLOCATE(virgin, renumber)
    ELSE
       CALL error_petsc('Bug in vtu_3d: mesh%gauss%n_w is not correct')
    END IF

    IF (PRESENT(opt_it)) THEN
       it = opt_it
       CALL make_vtu_file_3D(MPI_COMM_WORLD, mesh_3d, header, &
            field_out, name_of_field, what, opt_it = it)
    ELSE
       CALL make_vtu_file_3D(MPI_COMM_WORLD, mesh_3d, header, &
            field_out, name_of_field, what)
    END IF

    !===Cleanup
    DEALLOCATE(dist_field, field_out, loc_to_glob,field_out_xyz_FFT)
    DEALLOCATE(field_out_FFT)
  END SUBROUTINE vtu_3d_base

  SUBROUTINE transfer_fourier_to_real(communicator, mesh, field_in, &
       header, name_of_field, what, list_mode, opt_it)
    USE dyn_line
    USE def_type_mesh
    USE vtk_viz
    USE sft_parallele
    USE input_data
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type),                INTENT(IN) :: mesh
    CHARACTER(*),                   INTENT(IN) :: header, name_of_field, what
    INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode
    INTEGER, OPTIONAL,              INTENT(IN) :: opt_it
    REAL(KIND=8), DIMENSION(:,:,:)             :: field_in
    TYPE(mesh_type)                            :: mesh_3d
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE  :: field_out_FFT
    REAL(KIND=8),DIMENSION(:,:),  ALLOCATABLE  :: field_out
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE  :: dist_field
    LOGICAL,      DIMENSION(:),     POINTER    :: virgin
    INTEGER,      DIMENSION(:),     POINTER    :: renumber
    TYPE(mesh_type)                            :: local_mesh
    TYPE(dyn_int_line), DIMENSION(:), POINTER  :: loc_to_glob
    INTEGER                                    :: gap, width, mode_min
    INTEGER, DIMENSION(1)                      :: loc
    REAL(KIND=8)                               :: pi, dtheta, theta
    INTEGER                                    :: i, k, m, n, np_b, np_e, np_loc, np_loc_max, &
         type, nb_angle, count, dim, np1
    PetscErrorCode :: ierr
    PetscMPIInt    :: rank, nb_procs
    MPI_Comm       :: communicator
    CALL MPI_Comm_rank(communicator, rank, ierr)
    CALL MPI_Comm_size(communicator, nb_procs, ierr)

    IF (inputs%select_mode) THEN
       IF (SIZE(inputs%list_mode_lect)==1) THEN
          gap = inputs%list_mode_lect(1) + 1
       ELSE
          gap = inputs%list_mode_lect(2) - inputs%list_mode_lect(1)
       END IF
    ELSE
       gap = 1
    END IF
    width = SIZE(list_mode)*gap

    CALL divide_mesh(communicator, mesh, local_mesh, loc_to_glob)
    np_loc_max = -1
    DO n = 1, nb_procs
       np_loc_max = MAX(np_loc_max,MAXVAL(loc_to_glob(n)%DIL))
    END DO
    ALLOCATE(dist_field(nb_procs*np_loc_max,SIZE(field_in,2),width))
    dist_field = 0.d0
    mode_min = rank*width
    DO i = 1, width
       IF (MINVAL(ABS(list_mode - (mode_min+i-1)))/=0) CYCLE
       loc = MINLOC(ABS(list_mode - (mode_min+i-1)))
       DO type = 1, SIZE(field_in,2)
          DO n = 1, nb_procs
             np_b = (n-1)*np_loc_max + 1
             np_loc = SIZE(loc_to_glob(n)%DIL)
             np_e = np_b + np_loc - 1
             dist_field(np_b:np_e,type,i) = field_in(loc_to_glob(n)%DIL,type,loc(1))
          END DO
       END DO
    END DO

    CALL FFT_PAR_REAL(communicator, dist_field, field_out_FFT, opt_nb_plane=50)

    IF (SIZE(field_in,2)==2) THEN
       dim = 1
    ELSE IF (SIZE(field_in,2)==6) THEN
       dim = 3
    ELSE
       CALL error_petsc('Bug in transfer_fourier_to_real: SIZE(field_in,2) not correct')
    END IF
    nb_angle = SIZE(field_out_FFT,1)
    pi = ACOS(-1.d0)
    dtheta = 2*pi/nb_angle
    IF (mesh%gauss%n_w==3) THEN
       mesh_3d%np = local_mesh%np*nb_angle
       mesh_3d%me = local_mesh%me*nb_angle
       ALLOCATE(field_out(dim,mesh_3d%np))
       count = 0
       DO k = 1, nb_angle
          DO i = 1, local_mesh%np
             count = count + 1
             field_out(:,count) = field_out_FFT(k,:,i)
          ENDDO
       END DO
       !===Create 3d mesh
       ALLOCATE(mesh_3d%rr(dim,mesh_3d%np))
       count = 0
       DO k = 1, nb_angle
          theta = (k-1)*dtheta
          DO i = 1, local_mesh%np
             count = count + 1
             mesh_3d%rr(1,count) = local_mesh%rr(1,i)*COS(theta)
             mesh_3d%rr(2,count) = local_mesh%rr(1,i)*SIN(theta)
             mesh_3d%rr(3,count) = local_mesh%rr(2,i)
          ENDDO
       END DO
       ALLOCATE(mesh_3d%jj(6,mesh_3d%me)) !===P1 3d VTK wedge 13
       mesh_3d%gauss%n_w = 6
       count = 0
       DO k = 1, nb_angle-1
          DO m = 1, local_mesh%me
             count = count + 1
             mesh_3d%jj(1:3,count) = local_mesh%jj(1:3,m)+(k-1)*local_mesh%np
             mesh_3d%jj(4:6,count) = local_mesh%jj(1:3,m)+k*local_mesh%np
          END DO
       END DO
       DO m = 1, local_mesh%me
          count = count + 1
          mesh_3d%jj(1:3,count) = local_mesh%jj(1:3,m)+(k-1)*local_mesh%np
          mesh_3d%jj(4:6,count) = local_mesh%jj(1:3,m)
       END DO
    ELSE IF (mesh%gauss%n_w==6) THEN
       ALLOCATE(virgin(local_mesh%np), renumber(local_mesh%np))
       virgin = .TRUE.
       renumber = -1
       count = 0
       DO m = 1, local_mesh%me
          DO n = 1, 3
             i = local_mesh%jj(n,m)
             IF (.NOT.virgin(i)) CYCLE
             virgin(i) = .FALSE.
             count = count + 1 !===Count P1 nodes on local_mesh
             renumber(i) = count
          END DO
       END DO
       np1 = count
       mesh_3d%np = (local_mesh%np+np1)*nb_angle
       mesh_3d%me = local_mesh%me*nb_angle
       ALLOCATE(mesh_3d%rr(dim,mesh_3d%np))
       ALLOCATE(mesh_3d%jj(15,mesh_3d%me)) !===P2 3d VTK wedge 26
       mesh_3d%gauss%n_w = 15
       count = 0
       DO k = 1, nb_angle-1
          DO m = 1, local_mesh%me
             count = count + 1
             mesh_3d%jj(1:3,count) = local_mesh%jj(1:3,m) + (k-1)*(local_mesh%np+np1)
             mesh_3d%jj(4:6,count) = local_mesh%jj(1:3,m) + k*(local_mesh%np+np1)
             mesh_3d%jj(7,count)   = local_mesh%jj(6,m)   + (k-1)*(local_mesh%np+np1)
             mesh_3d%jj(8,count)   = local_mesh%jj(4,m)   + (k-1)*(local_mesh%np+np1)
             mesh_3d%jj(9,count)   = local_mesh%jj(5,m)   + (k-1)*(local_mesh%np+np1)
             mesh_3d%jj(10,count)  = local_mesh%jj(6,m)   + k*(local_mesh%np+np1)
             mesh_3d%jj(11,count)  = local_mesh%jj(4,m)   + k*(local_mesh%np+np1)
             mesh_3d%jj(12,count)  = local_mesh%jj(5,m)   + k*(local_mesh%np+np1)
             mesh_3d%jj(13:15,count) = renumber(local_mesh%jj(1:3,m)) &
                  + (k-1)*(local_mesh%np+np1) + local_mesh%np
          END DO
       END DO
       k = nb_angle
       DO m = 1, local_mesh%me
          count = count + 1
          mesh_3d%jj(1:3,count) = local_mesh%jj(1:3,m) + (k-1)*(local_mesh%np+np1)
          mesh_3d%jj(4:6,count) = local_mesh%jj(1:3,m)
          mesh_3d%jj(7,count)   = local_mesh%jj(6,m)   + (k-1)*(local_mesh%np+np1)
          mesh_3d%jj(8,count)   = local_mesh%jj(4,m)   + (k-1)*(local_mesh%np+np1)
          mesh_3d%jj(9,count)   = local_mesh%jj(5,m)   + (k-1)*(local_mesh%np+np1)
          mesh_3d%jj(10,count)  = local_mesh%jj(6,m)
          mesh_3d%jj(11,count)  = local_mesh%jj(4,m)
          mesh_3d%jj(12,count)  = local_mesh%jj(5,m)
          mesh_3d%jj(13:15,count) = renumber(local_mesh%jj(1:3,m)) &
               + (k-1)*(local_mesh%np+np1) + local_mesh%np
       END DO
       !===Create 3d mesh
       DO k = 1, nb_angle
          theta = (k-1)*dtheta
          DO i = 1, local_mesh%np
             n = (k-1)*(local_mesh%np+np1) + i
             mesh_3d%rr(1,n) = local_mesh%rr(1,i)*COS(theta)
             mesh_3d%rr(2,n) = local_mesh%rr(1,i)*SIN(theta)
             mesh_3d%rr(3,n) = local_mesh%rr(2,i)
             IF (virgin(i)) CYCLE !===This is a vertex
             n = (k-1)*(local_mesh%np+np1) + local_mesh%np + renumber(i)
             mesh_3d%rr(1,n) = local_mesh%rr(1,i)*COS(theta+dtheta/2)
             mesh_3d%rr(2,n) = local_mesh%rr(1,i)*SIN(theta+dtheta/2)
             mesh_3d%rr(3,n) = local_mesh%rr(2,i)
          ENDDO
       END DO

       ALLOCATE(field_out(dim,mesh_3d%np))
       DO k = 1, nb_angle
          DO i = 1, local_mesh%np
             n = (k-1)*(local_mesh%np+np1) + i
             field_out(:,n) = field_out_FFT(k,:,i)
          END DO
       END DO

       DO i = 1, local_mesh%np
          IF (virgin(i)) CYCLE !===This is a vertex
          DO k = 1, nb_angle - 2
             n = (k-1)*(local_mesh%np+np1) + local_mesh%np + renumber(i)
             field_out(:,n) = (3.d0/8)*field_out_FFT(k,:,i) &
                  + (3.d0/4)*field_out_FFT(k+1,:,i) - (1.d0/8)*field_out_FFT(k+2,:,i)
          END DO
       END DO
       k = nb_angle - 1
       DO i = 1, local_mesh%np
          IF (virgin(i)) CYCLE !===This is a vertex
          n = (k-1)*(local_mesh%np+np1) + local_mesh%np + renumber(i)
          field_out(:,n) = (3.d0/8)*field_out_FFT(k,:,i) &
               + (3.d0/4)*field_out_FFT(k+1,:,i) - (1.d0/8)*field_out_FFT(1,:,i)
       END DO
       k = nb_angle
       DO i = 1, local_mesh%np
          IF (virgin(i)) CYCLE !===This is a vertex
          n = (k-1)*(local_mesh%np+np1) + local_mesh%np + renumber(i)
          field_out(:,n) = (3.d0/8)*field_out_FFT(k,:,i) &
               + (3.d0/4)*field_out_FFT(1,:,i) - (1.d0/8)*field_out_FFT(2,:,i)
       END DO
       DEALLOCATE(virgin, renumber)

    ELSE
       CALL error_petsc('Bug in transfer_fourier_to_real: mesh%gauss%n_w is not correct')
    END IF

    IF (PRESENT(opt_it)) THEN
       CALL make_vtu_file_3D(MPI_COMM_WORLD, mesh_3d, header, &
            field_out, name_of_field, what, opt_it=opt_it)
    ELSE
       CALL make_vtu_file_3D(MPI_COMM_WORLD, mesh_3d, header, &
            field_out, name_of_field, what)
    END IF

    !===Cleanup
    DO n = 1, nb_procs
       IF (ASSOCIATED(loc_to_glob(n)%DIL)) DEALLOCATE(loc_to_glob(n)%DIL)
    END DO
    DEALLOCATE(loc_to_glob)
    DEALLOCATE(mesh_3d%rr, mesh_3d%jj)
    DEALLOCATE(local_mesh%rr, local_mesh%jj)
    DEALLOCATE(dist_field, field_out)
    DEALLOCATE(field_out_FFT)
  END SUBROUTINE transfer_fourier_to_real

  SUBROUTINE divide_mesh(communicator, mesh, local_mesh, loc_to_glob)
    USE def_type_mesh
    USE vtk_viz
    IMPLICIT NONE
    TYPE(mesh_type),              INTENT(IN)  :: mesh
    TYPE(mesh_type),              INTENT(OUT) :: local_mesh
    INTEGER                                   :: me_b, me_e, me_loc, np_loc
    TYPE(dyn_int_line), DIMENSION(:), POINTER :: loc_to_glob
    INTEGER, POINTER, DIMENSION(:)            :: glob_to_loc
    LOGICAL, DIMENSION(mesh%np)               :: virgin
    INTEGER                                   :: count, m, m_glob, n, i, q, r, loop_rank
    PetscErrorCode :: ierr
    PetscMPIInt    :: rank, nb_procs
    MPI_Comm       :: communicator

    !===Deal with empty mesh
    IF (mesh%me==0) THEN
       local_mesh%me = 0
       local_mesh%np = 0
       local_mesh%gauss%n_w = 0
       ALLOCATE(local_mesh%jj(0,0))
       ALLOCATE(local_mesh%rr(0,0))
       ALLOCATE(loc_to_glob(0))
       RETURN
    END IF

    !===Continue with non empty mesh
    CALL MPI_Comm_rank(communicator,rank,ierr)
    CALL MPI_Comm_size(communicator,nb_procs,ierr)
    ALLOCATE(loc_to_glob(nb_procs))
    ALLOCATE(glob_to_loc(mesh%np))

    !===Quotient and remainder
    q = mesh%me/nb_procs
    r = mesh%me - q*nb_procs

    DO loop_rank = 0, nb_procs-1
       !===Compute me_b and me_b
       IF (loop_rank+1.LE.r) THEN
          me_b = loop_rank*(q+1) + 1
          me_e = me_b + q
       ELSE
          me_b = loop_rank*q + r + 1
          me_e = me_b + q - 1
       END IF
       me_loc = me_e - me_b + 1
       !===Create local mesh
       count = 0
       virgin = .TRUE.
       DO m = me_b, me_e
          DO n = 1, mesh%gauss%n_w
             i = mesh%jj(n,m)
             IF (.NOT.virgin(i)) CYCLE
             virgin(i)=.FALSE.
             count = count + 1
          END DO
       END DO
       np_loc = count
       ALLOCATE(loc_to_glob(loop_rank+1)%DIL(np_loc))
       count = 0
       virgin = .TRUE.
       DO m = me_b, me_e
          DO n = 1, mesh%gauss%n_w
             i = mesh%jj(n,m)
             IF (.NOT.virgin(i)) CYCLE
             virgin(i)=.FALSE.
             count = count + 1
             loc_to_glob(loop_rank+1)%DIL(count) = i
             IF (loop_rank == rank) glob_to_loc(i) = count
          END DO
       END DO

       IF (loop_rank == rank) THEN
          !===Create local mesh
          local_mesh%me = me_loc
          local_mesh%np = np_loc
          local_mesh%gauss%n_w = mesh%gauss%n_w
          ALLOCATE(local_mesh%jj(mesh%gauss%n_w,me_loc))
          DO m = 1, me_loc
             m_glob = me_b + m -1
             local_mesh%jj(:,m) = glob_to_loc(mesh%jj(:,m_glob))
          END DO
          ALLOCATE(local_mesh%rr(2,local_mesh%np))
          local_mesh%rr = mesh%rr(:,loc_to_glob(rank+1)%DIL)
       END IF
    END DO

    DEALLOCATE(glob_to_loc)
  END SUBROUTINE divide_mesh

  SUBROUTINE compute_grad_phi(field_in, field_out)
    USE def_type_mesh
    USE my_util
    USE fem_M_axi
    USE st_matrix
    USE solve_petsc
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)    :: field_in
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT) :: field_out

    INTEGER, POINTER, DIMENSION(:)                :: ifrom
    TYPE(solver_param)                            :: param_grad_phi
    REAL(KIND=8), DIMENSION(SIZE(field_out,1), SIZE(field_out,2)) :: smb_grad_phi
    REAL(KIND=8), DIMENSION(6)                  :: grad_phi_loc
    INTEGER, DIMENSION(mesh_grad_phi%gauss%n_w) :: i_loc
    INTEGER               :: i, mode, m, l, j_loc, k, nj
    REAL(KIND=8)          :: ray, drdthetadz
    PetscErrorCode                              :: ierr


    param_grad_phi%it_max=1
    param_grad_phi%rel_tol=1.d-6
    param_grad_phi%abs_tol=1.d-10
    param_grad_phi%solver='GMRES'
    param_grad_phi%precond='MUMPS'
    param_grad_phi%verbose=.FALSE.

    IF (if_first_grad_phi) THEN
       if_first_grad_phi=.FALSE.

       CALL create_my_ghost(mesh_grad_phi, vizu_grad_phi_LA,ifrom)
       CALL VecCreateGhost(vizu_grad_curl_comm(1), mesh_grad_phi%dom_np, &
            PETSC_DETERMINE, SIZE(ifrom), ifrom, vx_phi, ierr)
       CALL VecGhostGetLocalForm(vx_phi, vx_phi_ghost, ierr)
       CALL VecDuplicate(vx_phi, vb_phi, ierr)

       CALL create_local_petsc_matrix(vizu_grad_curl_comm(1),vizu_grad_phi_LA , mat_grad_phi, clean=.TRUE.)
       CALL qs_diff_mass_scal_M (mesh_grad_phi, vizu_grad_phi_LA, 0.d0, 1.d0, 0.d0, 0, mat_grad_phi)
       CALL init_solver(param_grad_phi,ksp_grad_phi,mat_grad_phi,vizu_grad_curl_comm(1),&
            solver='GMRES',precond='MUMPS')
       CALL MatDestroy(mat_grad_phi,ierr)
    END IF

    smb_grad_phi = 0
    field_out = 0.d0
    DO i = 1, SIZE(vizu_list_mode)
       mode = vizu_list_mode(i)
       smb_grad_phi=0.d0
       DO m = 1, mesh_grad_phi%me
          DO l = 1, mesh_grad_phi%gauss%l_G
             ray = SUM(mesh_grad_phi%rr(1,mesh_grad_phi%jj(:,m))*mesh_grad_phi%gauss%ww(:,l))
             drdthetadz = mesh_grad_phi%gauss%rj(l,m)
             grad_phi_loc = 0.d0
             DO nj = 1,mesh_grad_phi%gauss%n_w
                j_loc = mesh_grad_phi%jj(nj,m)
                grad_phi_loc(1) = grad_phi_loc(1) + field_in(j_loc,1,i)*mesh_grad_phi%gauss%dw(1,nj,l,m)*ray
                grad_phi_loc(2) = grad_phi_loc(2) + field_in(j_loc,2,i)*mesh_grad_phi%gauss%dw(1,nj,l,m)*ray
                grad_phi_loc(3) = grad_phi_loc(3) + field_in(j_loc,2,i)*mesh_grad_phi%gauss%ww(nj,l)*mode
                grad_phi_loc(4) = grad_phi_loc(4) - field_in(j_loc,1,i)*mesh_grad_phi%gauss%ww(nj,l)*mode
                grad_phi_loc(5) = grad_phi_loc(5) + field_in(j_loc,1,i)*mesh_grad_phi%gauss%dw(2,nj,l,m)*ray
                grad_phi_loc(6) = grad_phi_loc(6) + field_in(j_loc,2,i)*mesh_grad_phi%gauss%dw(2,nj,l,m)*ray
             END DO
             i_loc = mesh_grad_phi%jj(:,m)

             DO k = 1, 6
                smb_grad_phi(i_loc,k) = smb_grad_phi(i_loc,k) + grad_phi_loc(k)* &
                     mesh_grad_phi%gauss%ww(:,l)*drdthetadz
             END DO
          END DO
       END DO

       DO k = 1, 6
          CALL VecZeroEntries(vb_phi, ierr)
          CALL VecSetValues(vb_phi, mesh_grad_phi%np, vizu_grad_phi_LA%loc_to_glob(1,:)-1, smb_grad_phi(:,k), ADD_VALUES, ierr)
          CALL VecAssemblyBegin(vb_phi,ierr)
          CALL VecAssemblyEnd(vb_phi,ierr)

          CALL solver(ksp_grad_phi,vb_phi,vx_phi,reinit=.FALSE.,verbose=.FALSE.)

          CALL VecGhostUpdateBegin(vx_phi,INSERT_VALUES,SCATTER_FORWARD,ierr)
          CALL VecGhostUpdateEnd(vx_phi,INSERT_VALUES,SCATTER_FORWARD,ierr)
          IF (mesh_grad_phi%me/=0) THEN
             CALL extract(vx_phi_ghost,1,1,vizu_grad_phi_LA,field_out(:,k,i))
          END IF
       END DO

       IF (mode == 0) THEN
          field_out(:,2,i) = 0.d0
          field_out(:,4,i) = 0.d0
          field_out(:,6,i) = 0.d0
       END IF
    END DO

  END SUBROUTINE compute_grad_phi

  SUBROUTINE compute_rot_h(field_in, field_out)
    USE def_type_mesh
    USE my_util
    USE fem_M_axi
    USE st_matrix
    USE solve_petsc
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)    :: field_in
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT) :: field_out
    INTEGER, POINTER, DIMENSION(:)                :: ifrom
    TYPE(solver_param)                            :: param_rot_h
    REAL(KIND=8), DIMENSION(SIZE(field_out,1), SIZE(field_out,2)) :: smb_rot_h
    REAL(KIND=8), DIMENSION(6)                  :: rot_h_loc
    INTEGER, DIMENSION(mesh_rot_h%gauss%n_w) :: i_loc
    INTEGER               :: i, mode, m, l, j_loc, k, nj
    REAL(KIND=8)          :: ray, drdthetadz
    PetscErrorCode                              :: ierr

    param_rot_h%it_max=1
    param_rot_h%rel_tol=1.d-6
    param_rot_h%abs_tol=1.d-10
    param_rot_h%solver='GMRES'
    param_rot_h%precond='MUMPS'
    param_rot_h%verbose=.FALSE.

    IF (if_first_rot_h) THEN
       if_first_rot_h=.FALSE.

       CALL create_my_ghost(mesh_rot_h, vizu_rot_h_LA, ifrom)
       CALL VecCreateGhost(vizu_grad_curl_comm(1), mesh_rot_h%dom_np, &
            PETSC_DETERMINE, SIZE(ifrom), ifrom, vx_phi, ierr)
       CALL VecGhostGetLocalForm(vx_phi, vx_phi_ghost, ierr)
       CALL VecDuplicate(vx_phi, vb_phi, ierr)

       CALL create_local_petsc_matrix(vizu_grad_curl_comm(1),vizu_rot_h_LA , mat_rot_h, clean=.TRUE.)
       CALL qs_diff_mass_scal_M (mesh_rot_h, vizu_rot_h_LA, 0.d0, 1.d0, 0.d0, 0, mat_rot_h)
       CALL init_solver(param_rot_h,ksp_rot_h,mat_rot_h,vizu_grad_curl_comm(1),&
            solver='GMRES',precond='MUMPS')
       CALL MatDestroy(mat_rot_h,ierr)
    END IF

    smb_rot_h = 0.d0
    field_out = 0.d0
    DO i = 1, SIZE(vizu_list_mode)
       mode = vizu_list_mode(i)
       smb_rot_h=0.d0
       DO m = 1, mesh_rot_h%me
          DO l = 1, mesh_rot_h%gauss%l_G
             ray = SUM(mesh_rot_h%rr(1,mesh_rot_h%jj(:,m))*mesh_rot_h%gauss%ww(:,l))
             drdthetadz = mesh_rot_h%gauss%rj(l,m)
             rot_h_loc = 0.d0
             DO nj = 1,mesh_rot_h%gauss%n_w
                j_loc = mesh_rot_h%jj(nj,m)
                rot_h_loc(1) = rot_h_loc(1) + mode*field_in(j_loc,6,i)*mesh_rot_h%gauss%ww(nj,l) &
                     -field_in(j_loc,3,i)*mesh_rot_h%gauss%dw(2,nj,l,m)*ray
                rot_h_loc(2) = rot_h_loc(2) - mode*field_in(j_loc,5,i)*mesh_rot_h%gauss%ww(nj,l) &
                     -field_in(j_loc,4,i)*mesh_rot_h%gauss%dw(2,nj,l,m)*ray
                rot_h_loc(3) = rot_h_loc(3) + field_in(j_loc,1,i)*mesh_rot_h%gauss%dw(2,nj,l,m)*ray &
                     -field_in(j_loc,5,i)*mesh_rot_h%gauss%dw(1,nj,l,m)*ray
                rot_h_loc(4) = rot_h_loc(4) + field_in(j_loc,2,i)*mesh_rot_h%gauss%dw(2,nj,l,m)*ray &
                     -field_in(j_loc,6,i)*mesh_rot_h%gauss%dw(1,nj,l,m)*ray
                rot_h_loc(5) = rot_h_loc(5) + field_in(j_loc,3,i)*mesh_rot_h%gauss%ww(nj,l) &
                     +field_in(j_loc,3,i)*mesh_rot_h%gauss%dw(1,nj,l,m)*ray &
                     - mode*field_in(j_loc,2,i)*mesh_rot_h%gauss%ww(nj,l)

                rot_h_loc(6) = rot_h_loc(6) + field_in(j_loc,4,i)*mesh_rot_h%gauss%ww(nj,l) &
                     +field_in(j_loc,4,i)*mesh_rot_h%gauss%dw(1,nj,l,m)*ray &
                     + mode*field_in(j_loc,1,i)*mesh_rot_h%gauss%ww(nj,l)
             END DO !=== sum on nj
             i_loc = mesh_rot_h%jj(:,m)

             DO k = 1, 6
                smb_rot_h(i_loc,k) = smb_rot_h(i_loc,k) + rot_h_loc(k)* &
                     mesh_rot_h%gauss%ww(:,l)*drdthetadz
             END DO
          END DO
       END DO

       DO k = 1, 6
          CALL VecZeroEntries(vb_phi, ierr)
          CALL VecSetValues(vb_phi, mesh_rot_h%np, vizu_rot_h_LA%loc_to_glob(1,:)-1, smb_rot_h(:,k), ADD_VALUES, ierr)
          CALL VecAssemblyBegin(vb_phi,ierr)
          CALL VecAssemblyEnd(vb_phi,ierr)

          CALL solver(ksp_rot_h,vb_phi,vx_phi,reinit=.FALSE.,verbose=.FALSE.)

          CALL VecGhostUpdateBegin(vx_phi,INSERT_VALUES,SCATTER_FORWARD,ierr)
          CALL VecGhostUpdateEnd(vx_phi,INSERT_VALUES,SCATTER_FORWARD,ierr)
          IF (mesh_rot_h%me/=0) THEN
             CALL extract(vx_phi_ghost,1,1,vizu_rot_h_LA,field_out(:,k,i))
          END IF
       END DO

       IF (mode == 0) THEN
          field_out(:,2,i) = 0.d0
          field_out(:,4,i) = 0.d0
          field_out(:,6,i) = 0.d0
       END IF
    END DO !===end loop on mode
  END SUBROUTINE compute_rot_h

  SUBROUTINE compute_rot_u(field_in, field_out, mesh_in)
    USE def_type_mesh
    USE my_util
    USE fem_M_axi
    USE st_matrix
    USE solve_petsc
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)    :: field_in
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT) :: field_out
    TYPE(mesh_type), INTENT(IN)                   :: mesh_in
    INTEGER, POINTER, DIMENSION(:)                :: ifrom
    TYPE(solver_param)                            :: param_rot_u
    REAL(KIND=8), DIMENSION(SIZE(field_out,1), SIZE(field_out,2)) :: smb_rot_u
    REAL(KIND=8), DIMENSION(6)                  :: rot_u_loc
    INTEGER, DIMENSION(mesh_in%gauss%n_w) :: i_loc
    INTEGER               :: i, mode, m, l, j_loc, k, nj
    REAL(KIND=8)          :: ray, drdthetadz
    PetscErrorCode                              :: ierr

    param_rot_u%it_max=1
    param_rot_u%rel_tol=1.d-6
    param_rot_u%abs_tol=1.d-10
    param_rot_u%solver='GMRES'
    param_rot_u%precond='MUMPS'
    param_rot_u%verbose=.FALSE.

    IF (if_first_rot_u) THEN
       if_first_rot_u=.FALSE.

       CALL create_my_ghost(mesh_in, vizu_rot_u_LA, ifrom)
       CALL VecCreateGhost(vizu_grad_curl_comm(1), mesh_in%dom_np, &
            PETSC_DETERMINE, SIZE(ifrom), ifrom, vx_phi, ierr)
       CALL VecGhostGetLocalForm(vx_phi, vx_phi_ghost, ierr)
       CALL VecDuplicate(vx_phi, vb_phi, ierr)

       CALL create_local_petsc_matrix(vizu_grad_curl_comm(1),vizu_rot_u_LA , mat_rot_u, clean=.TRUE.)
       CALL qs_diff_mass_scal_M (mesh_in, vizu_rot_u_LA, 0.d0, 1.d0, 0.d0, 0, mat_rot_u)
       CALL init_solver(param_rot_u,ksp_rot_u,mat_rot_u,vizu_grad_curl_comm(1),&
            solver='GMRES',precond='MUMPS')
       CALL MatDestroy(mat_rot_u,ierr)
    END IF

    smb_rot_u = 0
    field_out = 0.d0
    DO i = 1, SIZE(vizu_list_mode)
       mode = vizu_list_mode(i)
       smb_rot_u=0.d0
       DO m = 1, mesh_in%me
          DO l = 1, mesh_in%gauss%l_G
             ray = SUM(mesh_in%rr(1,mesh_in%jj(:,m))*mesh_in%gauss%ww(:,l))
             drdthetadz = mesh_in%gauss%rj(l,m)
             rot_u_loc = 0.d0
             DO nj = 1,mesh_in%gauss%n_w
                j_loc = mesh_in%jj(nj,m)

                rot_u_loc(1) = rot_u_loc(1) + mode*field_in(j_loc,6,i)*mesh_in%gauss%ww(nj,l) &
                     -field_in(j_loc,3,i)*mesh_in%gauss%dw(2,nj,l,m)*ray
                rot_u_loc(2) = rot_u_loc(2) - mode*field_in(j_loc,5,i)*mesh_in%gauss%ww(nj,l) &
                     -field_in(j_loc,4,i)*mesh_in%gauss%dw(2,nj,l,m)*ray

                rot_u_loc(3) = rot_u_loc(3) + field_in(j_loc,1,i)*mesh_in%gauss%dw(2,nj,l,m)*ray &
                     -field_in(j_loc,5,i)*mesh_in%gauss%dw(1,nj,l,m)*ray
                rot_u_loc(4) = rot_u_loc(4) + field_in(j_loc,2,i)*mesh_in%gauss%dw(2,nj,l,m)*ray &
                     -field_in(j_loc,6,i)*mesh_in%gauss%dw(1,nj,l,m)*ray

                rot_u_loc(5) = rot_u_loc(5) + field_in(j_loc,3,i)*mesh_in%gauss%ww(nj,l) &
                     +field_in(j_loc,3,i)*mesh_in%gauss%dw(1,nj,l,m)*ray &
                     - mode*field_in(j_loc,2,i)*mesh_in%gauss%ww(nj,l)
                rot_u_loc(6) = rot_u_loc(6) + field_in(j_loc,4,i)*mesh_in%gauss%ww(nj,l) &
                     +field_in(j_loc,4,i)*mesh_in%gauss%dw(1,nj,l,m)*ray &
                     + mode*field_in(j_loc,1,i)*mesh_in%gauss%ww(nj,l)
             END DO !=== sum on nj
             i_loc = mesh_in%jj(:,m)

             DO k = 1, 6
                smb_rot_u(i_loc,k) = smb_rot_u(i_loc,k) + rot_u_loc(k)* &
                     mesh_in%gauss%ww(:,l)*drdthetadz
             END DO
          END DO
       END DO

       DO k = 1, 6
          CALL VecZeroEntries(vb_phi, ierr)
          CALL VecSetValues(vb_phi, mesh_in%np, vizu_rot_u_LA%loc_to_glob(1,:)-1, smb_rot_u(:,k), ADD_VALUES, ierr)
          CALL VecAssemblyBegin(vb_phi,ierr)
          CALL VecAssemblyEnd(vb_phi,ierr)

          CALL solver(ksp_rot_u,vb_phi,vx_phi,reinit=.FALSE.,verbose=.FALSE.)

          CALL VecGhostUpdateBegin(vx_phi,INSERT_VALUES,SCATTER_FORWARD,ierr)
          CALL VecGhostUpdateEnd(vx_phi,INSERT_VALUES,SCATTER_FORWARD,ierr)
          IF (mesh_in%me/=0) THEN
             CALL extract(vx_phi_ghost,1,1,vizu_rot_u_LA,field_out(:,k,i))
          END IF
       END DO

       IF (mode == 0) THEN
          field_out(:,2,i) = 0.d0
          field_out(:,4,i) = 0.d0
          field_out(:,6,i) = 0.d0
       END IF
    END DO !===end loop on mode
  END SUBROUTINE compute_rot_u

  SUBROUTINE divide_mesh_into_four_subcells(mesh_in,mesh_out)
    !  JLG + HF july 25, 2019
    !  Scheme of the sub_mesh partition
    !     n
    !     |\
    !     |n \
    !     |    \
    !     |      \n1+3
    ! n2+3-        \
    !     |          \
    !     |4*(m-1)+n   \
    !     |              \
    !     |n1     n+3    n2\              this is the element m (large mesh)
    !     |--------|--------|\              divided in 4*(m-1)+1, 4*(m-1)+2,
    !     |\n2    n+3     n1|n \              4*(m-1)+3, 4*m in small mesh
    !     |n\               |    \
    !     |   \        4m   |      \
    !     |     \n1+3       |        \n1+3
    !     |       \     n2+3-n2+3      \
    ! n2+3-    n1+3 \       |            \
    !     |           \     |              \
    !     |             \   |                \
    !     | 4*(m-1)+n1    \n|   4*(m-1)+n2     \
    !     |n1            n2\|n1                n2\
    !     --------|--------------------|-----------
    !   n1       n+3                  n+3         n2
    !      same edge touching (separed for clarity)
    !   n1op     nop+3               nop+3            n2op
    !     ---------|--------------------|----------
    !     |n1op       n2op /|n1op           n2op /
    !     | 4(mop-1)+n1op / |   4(mop-1)+n2op  /
    !     |             /nop|                /
    !     |           /     |              /
    !     |         /       |            /
    !     |       /         |          /
    !     |     /           |        /
    !     |   /   4mop      |      /
    !     | /               |    /
    !     |/                |  /          this is the neigh element mop (large mesh)
    !     |-----------------|/              divided in 4*(mop-1)+1, 4*(mop-1)+2,
    !     |                /                 4*(mop-1)+3, 4*mop in small mesh
    !     |              /
    !     |            /
    !     |          /
    !     |        /
    !     |      /
    !     |    /
    !     |  /
    !     |/
    !     nop
    IMPLICIT NONE
    TYPE(mesh_type), INTENT(IN)      :: mesh_in
    TYPE(mesh_type), INTENT(OUT)     :: mesh_out
    LOGICAL, DIMENSION(mesh_in%np)   :: v_virgin
    LOGICAL, DIMENSION(mesh_in%me)   :: c_virgin
    INTEGER, DIMENSION(mesh_in%np)   :: v_in_to_out
    INTEGER :: edge, nb_vertices, nb_edges, n_dof, &
         m, m_out, j, m_op, m_out_op, n, n1, n2, n_op, n1_op, n2_op, nn, &
         j1, j2
    mesh_out%gauss%k_d = 2
    mesh_out%gauss%n_w = 6
    mesh_out%me = 4*mesh_in%me

    !===Nb of edges
    edge = 0
    DO m = 1, mesh_in%me
       DO n = 1, 3
          IF (mesh_in%neigh(n,m).LE.0) CYCLE  !===ATTENTION, no dg mesh allowed here
          edge = edge + 1
       END DO
    END DO
    edge = edge/2
    DO m = 1, mesh_in%me
       DO n = 1, 3
          IF (mesh_in%neigh(n,m).GE.0) CYCLE  !===ATTENTION, no dg mesh allowed here
          edge = edge + 1
       END DO
    END DO
    nb_edges = edge + mesh_in%mes

    IF (mesh_in%gauss%n_w==10) THEN
       nb_vertices = mesh_in%np - 2*nb_edges - mesh_in%me
    ELSE
       CALL error_petsc('BUG in divide_mesh_into_four_subcells: mesh to be divided is not P3')
    END IF
    mesh_out%np = nb_vertices + 3*nb_edges + 3*mesh_in%me !===Create P2-iso-P4 mesh

    !===Allocation
    ALLOCATE(mesh_out%jj(6,mesh_out%me), mesh_out%rr(2,mesh_out%np))
    v_virgin = .TRUE.
    n_dof = 0
    DO m = 1, mesh_in%me
       DO n = 1, 3
          j = mesh_in%jj(n,m)
          IF (v_virgin(j)) THEN
             v_virgin(j) = .FALSE.
             n_dof = n_dof + 1
             v_in_to_out(j) = n_dof
          END IF
       END DO
    END DO

    IF (n_dof.NE.nb_vertices) THEN
       write(*,*) mesh_in%np, nb_edges, mesh_in%me, nb_vertices, n_dof
       CALL error_petsc('BUG divide_mesh_into_four_subcells: n_dof.NE.nb_vertices')
    END IF

    DO m = 1, mesh_in%me
       DO n = 1, 3 !===take care of vertices
          m_out = 4*(m-1)+n
          j = mesh_in%jj(n,m)
          mesh_out%jj(n,m_out) = v_in_to_out(j)
          mesh_out%rr(:,v_in_to_out(j)) = mesh_in%rr(:,j)
       END DO
    END DO

    c_virgin = .TRUE.
    n_dof = nb_vertices
    DO m = 1, mesh_in%me
       c_virgin(m)=.FALSE.
       DO n = 1, 3 !===We take care of faces and interfaces
          m_op = mesh_in%neigh(n,m)
          IF (m_op.GT.0) THEN
             IF (.NOT.c_virgin(m_op)) CYCLE
          END IF
          IF (m_op.LE.0) THEN
             n_op = n !===Do not create the nodes on the opposite cell
             m_op = m
          ELSE
             DO n_op = 1, 3
                IF (mesh_in%neigh(n_op,m_op)==m) EXIT !===n_op is the local index of the opposite vertex
             END DO
          END IF
          n1_op = MODULO(n_op,3) + 1
          n2_op = MODULO(n_op+1,3) + 1
          n1 = MODULO(n,3) + 1
          n2 = MODULO(n+1,3) + 1
          IF (mesh_in%jj(n1_op,m_op) .NE. mesh_in%jj(n1,m)) THEN
             nn = n1_op
             n1_op = n2_op
             n2_op = nn
          END IF
          IF (mesh_in%jj(n1,m).NE.mesh_in%jj(n1_op,m_op)) THEN
             CALL error_petsc('divide_mesh_into_four_subcells: BUG')
          END IF

          n_dof = n_dof + 1 !===New vertex of P2 sub-mesh
          j1=mesh_in%jj(n1,m)
          j2=mesh_in%jj(n2,m)
          mesh_out%rr(:,n_dof) = (mesh_in%rr(:,j1)+mesh_in%rr(:,j2))/2

          m_out    = 4*(m-1)   +n1
          m_out_op = 4*(m_op-1)+n1_op
          mesh_out%jj(n2,m_out)       = n_dof
          mesh_out%jj(n2_op,m_out_op) = n_dof

          m_out    = 4*(m-1)   +n2
          m_out_op = 4*(m_op-1)+n2_op
          mesh_out%jj(n1,m_out)       = n_dof
          mesh_out%jj(n1_op,m_out_op) = n_dof

          m_out    = 4*m
          m_out_op = 4*m_op
          mesh_out%jj(n,m_out)       = n_dof
          mesh_out%jj(n_op,m_out_op) = n_dof

          n_dof = n_dof +1 !===New midpoint of P2 sub-mesh
          mesh_out%rr(:,n_dof) = mesh_in%rr(:,j1)*3/4 + mesh_in%rr(:,j2)/4

          m_out    = 4*(m-1)   +n1
          m_out_op = 4*(m_op-1)+n1_op
          mesh_out%jj(n+3,m_out)       = n_dof
          mesh_out%jj(n_op+3,m_out_op) = n_dof

          n_dof = n_dof +1 !===New midpoint of P2 sub-mesh
          mesh_out%rr(:,n_dof) = mesh_in%rr(:,j1)/4 + mesh_in%rr(:,j2)*3/4

          m_out    = 4*(m-1)   +n2
          m_out_op = 4*(m_op-1)+n2_op
          mesh_out%jj(n+3,m_out)       = n_dof
          mesh_out%jj(n_op+3,m_out_op) = n_dof

       END DO

       !===Now we take care of the dofs inside the cell m
       DO n = 1, 3
          n1 = MODULO(n,3) + 1
          n2 = MODULO(n+1,3) + 1
          j1=mesh_in%jj(n1,m)
          j2=mesh_in%jj(n2,m)
          j =mesh_in%jj(n ,m)
          n_dof = n_dof + 1
          mesh_out%rr(:,n_dof) = mesh_in%rr(:,j)/2 + mesh_in%rr(:,j1)/4 + mesh_in%rr(:,j2)/4

          m_out    = 4*m
          mesh_out%jj(n+3,m_out) = n_dof
          m_out    = 4*(m-1) + n
          mesh_out%jj(n+3,m_out) = n_dof
       END DO
    END DO
  END SUBROUTINE divide_mesh_into_four_subcells

  SUBROUTINE interpolate(field_in,mesh_in,field_out,mesh_out)
    IMPLICIT NONE
    TYPE(mesh_type),                INTENT(IN)  :: mesh_in, mesh_out
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: field_in
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: field_out
    INTEGER :: m_out, m_in, n1in, n2in, nin, nout
    REAL(KIND=8) :: mesK, mesin, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, x, y
    REAL(KIND=8) :: lambda(3), rout(2), rr(2,3), a(2), b(2)
    INTEGER, DIMENSION(10) :: j10
    f1(x,y) = -0.9d1/0.2d1*x**3 - 0.27d2/0.2d1*x**2*y - 0.27d2/0.2d1*x*y**2 &
         - 0.9d1/0.2d1*y**3 + 0.9d1*x**2 + 0.18d2*x*y + 0.9d1*y**2 &
         - 0.11d2/0.2d1*x - 0.11d2/0.2d1*y + 0.1d1
    f2(x,y) = 0.9d1/0.2d1*x**3 - 0.9d1/0.2d1*x**2 + x
    f3(x,y) = 0.9d1/0.2d1*y**3 - 0.9d1/0.2d1*y**2 + y
    f4(x,y) = 0.27d2/0.2d1*x**2*y - 0.9d1/0.2d1*x*y
    f5(x,y) = 0.27d2/0.2d1*x*y**2 - 0.9d1/0.2d1*x*y
    f6(x,y) = 0.27d2/0.2d1*x**2*y + 0.27d2*x*y**2 + 0.27d2 / 0.2d1*y**3 &
         - 0.45d2/0.2d1*x*y - 0.45d2/0.2d1*y**2 + 0.9d1*y
    f7(x,y) = -0.27d2/0.2d1*x*y**2 - 0.27d2/0.2d1*y**3 + 0.9d1/0.2d1*x*y &
         + 0.18d2*y**2 - 0.9d1/0.2d1*y
    f8(x,y) = 0.27d2/0.2d1*x**3 + 0.27d2*x**2*y + 0.27d2/0.2d1*x*y**2 &
         - 0.45d2/0.2d1*x**2 - 0.45d2/0.2d1*x*y + 0.9d1*x
    f9(x,y) = -0.27d2/0.2d1*x**3 - 0.27d2/0.2d1*x**2*y + 0.18d2*x**2 &
         + 0.9d1/0.2d1*x*y - 0.9d1/0.2d1*x
    f10(x,y) = -27*x**2*y - 27*x*y**2 + 27*x*y

    DO m_out = 1, mesh_out%me
       m_in = (m_out-1)/4 + 1
       rr(:,1:3) = mesh_in%rr(:,mesh_in%jj(1:3,m_in))
       a = rr(:,2)-rr(:,1)
       b = rr(:,3)-rr(:,1)
       mesK = mesureK(a,b)
       DO nout = 1, 6
          rout = mesh_out%rr(:,mesh_out%jj(nout,m_out))
          DO nin = 1, 3
             n1in = MODULO(nin,3) + 1
             n2in = MODULO(nin+1,3) + 1
             a = rr(:,n1in)-rout
             b = rr(:,n2in)-rout
             mesin = mesureK(a,b)
             lambda(nin) = mesin/mesk
          END DO
          IF (ABS(SUM(lambda)-1.d0).GT.1.d-13) THEN
             CALL error_petsc('BUG in interpolate')
          END IF
          j10 = mesh_in%jj(:,m_in)
          field_out(mesh_out%jj(nout,m_out),:) = field_in(j10(1),:)*f1(lambda(2),lambda(3)) &
               + field_in(j10(2),:)*f2(lambda(2),lambda(3)) &
               + field_in(j10(3),:)*f3(lambda(2),lambda(3)) &
               + field_in(j10(4),:)*f4(lambda(2),lambda(3)) &
               + field_in(j10(5),:)*f5(lambda(2),lambda(3)) &
               + field_in(j10(6),:)*f6(lambda(2),lambda(3)) &
               + field_in(j10(7),:)*f7(lambda(2),lambda(3)) &
               + field_in(j10(8),:)*f8(lambda(2),lambda(3)) &
               + field_in(j10(9),:)*f9(lambda(2),lambda(3)) &
               + field_in(j10(10),:)*f10(lambda(2),lambda(3))
       END DO
    END DO
  CONTAINS
    FUNCTION mesureK(a,b) RESULT(mes)
      REAL(KIND=8), DIMENSION(2), INTENT(IN) :: a, b
      REAL(KIND=8)                           :: mes
      mes = ABS(a(1)*b(2)-a(2)*b(1))
    END FUNCTION mesureK
  END SUBROUTINE interpolate

  SUBROUTINE make_vtu_file_2D(communicator, mesh_in, header, field_in, field_name, what, opt_it)
    USE def_type_mesh
    USE my_util
    USE vtk_viz
    IMPLICIT NONE
    TYPE(mesh_type),       INTENT(IN), TARGET :: mesh_in
    TYPE(mesh_type),               POINTER    :: mesh
    CHARACTER(*),                  INTENT(IN) :: header
    CHARACTER(*),                  INTENT(IN) :: field_name, what
    INTEGER, OPTIONAL,             INTENT(IN) :: opt_it
    REAL(KIND=8), DIMENSION(:,:),  INTENT(IN), TARGET :: field_in
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, TARGET :: field_tmp
    REAL(KIND=8), DIMENSION(:,:),  POINTER    :: field
    INTEGER                                   :: j, it
    CHARACTER(LEN=200), DIMENSION(:), POINTER :: file_list
    CHARACTER(LEN=3)                          :: st_rank, st_it
    PetscErrorCode                            :: ierr
    PetscMPIInt                               :: rank, nb_procs
    MPI_Comm                                  :: communicator
    CALL MPI_Comm_rank(communicator, rank, ierr)
    CALL MPI_Comm_size(communicator, nb_procs, ierr)
    ALLOCATE(file_list(nb_procs))

    IF (mesh_in%gauss%n_w==10 .AND. inputs%type_fe_velocity==3) THEN
       ALLOCATE(field_tmp(vv_mesh_p2_iso_p4%np,SIZE(field_in,2)))
       CALL interpolate(field_in,mesh_in,field_tmp,vv_mesh_p2_iso_p4)
       field => field_tmp
       mesh => vv_mesh_p2_iso_p4
    ELSE
       field => field_in
       mesh => mesh_in
    END IF

    IF (PRESENT(opt_it)) THEN
       it = opt_it
       WRITE(st_it,'(I3)') it
       DO j = 1, nb_procs
          WRITE(st_rank,'(I3)') j
          file_list(j) = TRIM(header)//'_proc_'//TRIM(ADJUSTL(st_rank))//&
               '_it_'//TRIM(ADJUSTL(st_it))
       END DO
    ELSE
       DO j = 1, nb_procs
          WRITE(st_rank,'(I3)') j
          file_list(j) = TRIM(header)//'_proc_'//TRIM(ADJUSTL(st_rank))
       END DO
    END IF

    CALL check_list(communicator, file_list, mesh%np)
    IF (rank==0) THEN
       IF (PRESENT(opt_it)) THEN
          it = opt_it
       ELSE
          it = 1
       END IF
       CALL create_pvd_file(file_list, TRIM(header), it, TRIM(what))
    END IF

    IF (SIZE(field,2) == 6) THEN !vector
       CALL create_xml_vtu_vect_file(field, mesh, TRIM(ADJUSTL(file_list(rank+1))), field_name)
    ELSE IF (SIZE(field,2) == 1) THEN
       CALL create_xml_vtu_scal_file(field(:,1), mesh, TRIM(ADJUSTL(file_list(rank+1))), field_name)
    ELSE IF (SIZE(field,2) .GT. 0) THEN
       CALL create_xml_vtu_scal_file(field(:,1), mesh, TRIM(ADJUSTL(file_list(rank+1))), field_name)
    ELSE
       CALL error_Petsc('Bug in make_vtu_file_2D: field needs at least one component')
    END IF

    IF (ALLOCATED(field_tmp)) DEALLOCATE(field_tmp)
  END SUBROUTINE make_vtu_file_2D

END MODULE fourier_to_real_for_vtu
