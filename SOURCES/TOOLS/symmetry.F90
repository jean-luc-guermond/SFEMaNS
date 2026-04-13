MODULE symmetric_field
#include "petsc/finclude/petsc.h"
  USE petsc
  USE def_type_mesh
  USE st_matrix

  IMPLICIT NONE

!   PUBLIC :: symmetric_points, symm_champ, val_ener_sym_centrale, &
  PUBLIC :: symm_champ, val_ener_sym_centrale, &
       val_ener_sym_glob,val_ener_sym_rpi, val_ener_north_south, champ_total_anti_sym

  !Class containg tools necesary for symmetrization --------------------------
  TYPE, PRIVATE ::  symmetric_tools
     INTEGER, DIMENSION(:), ALLOCATABLE :: ltg_LA
     TYPE(petsc_csr_LA), ALLOCATABLE    :: LA
     TYPE(mesh_type), POINTER           :: mesh
     Vec,                ALLOCATABLE :: vec_mz, vec_mz_ghost
   CONTAINS
     PROCEDURE          :: init_symmetry_tool
     PROCEDURE          :: apply_symmetry
  END TYPE symmetric_tools
  !---------------------------------------------------------------------------

  !Symmetrization-------------------------------------------------------------
  TYPE(symmetric_tools), TARGET, ALLOCATABLE, PUBLIC :: vv_sym_tool, H_sym_tool
  TYPE(symmetric_tools), TARGET, ALLOCATABLE, PUBLIC :: T_sym_tool, c_sym_tool, phi_sym_tool
  LOGICAL, PRIVATE                                       :: need_sym=.FALSE.
  !---------------------------------------------------------------------------

  PRIVATE

CONTAINS


   SUBROUTINE init_symmetry_tool(self, mesh, communicator)
      CLASS(symmetric_tools), INTENT(INOUT) :: self
      TYPE(mesh_type), TARGET, INTENT(IN)           :: mesh
      INTEGER,          POINTER, DIMENSION(:)     :: ifrom
      INTEGER                                     :: n
      MPI_Comm       :: communicator
      PetscErrorCode       :: ierr

!=== initialize mapping for symmetric points
      ALLOCATE(self%ltg_LA(mesh%np))
      CALL symmetric_points(mesh, communicator, self%ltg_LA)

!=== build pointers
      self%mesh => mesh

!=== initialize petsc ghost vector for symmetric points
      ALLOCATE(self%LA)
      ALLOCATE(self%vec_mz)
      ALLOCATE(self%vec_mz_ghost)

      ALLOCATE(self%LA%loc_to_glob(1,mesh%np))
      self%LA%loc_to_glob(1,:) = mesh%loc_to_glob
      self%LA%kmax = 1
      CALL create_my_ghost(mesh,self%LA,ifrom)
      n = mesh%dom_np
      CALL VecCreateGhost(communicator, n, &
         PETSC_DETERMINE, SIZE(ifrom), ifrom, self%vec_mz, ierr)
      CALL VecGhostGetLocalForm(self%vec_mz, self%vec_mz_ghost, ierr)

   END SUBROUTINE init_symmetry_tool

   SUBROUTINE apply_symmetry(self, vv_in, vv_out)
      CLASS(symmetric_tools), INTENT(IN)          :: self
      REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: vv_in
      REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: vv_out
      INTEGER                                     :: i, j
      INTEGER, ALLOCATABLE, DIMENSION(:)          :: ix
      PetscScalar, POINTER :: x_loc(:)
      PetscErrorCode       :: ierr

      ALLOCATE(ix(self%mesh%np))
      ix = self%ltg_LA-1
      DO i = 1, SIZE(vv_in,2)
         DO j = 1, SIZE(vv_in,3)
            CALL VecZeroEntries(self%vec_mz, ierr)
            CALL VecSetValues(self%vec_mz, self%mesh%np, ix, vv_in(:,i,j), INSERT_VALUES, ierr)
            CALL VecAssemblyBegin(self%vec_mz, ierr)
            CALL VecAssemblyEnd(self%vec_mz, ierr)
            CALL VecGhostUpdateBegin(self%vec_mz,INSERT_VALUES, SCATTER_FORWARD,ierr)
            CALL VecGhostUpdateEnd(self%vec_mz,INSERT_VALUES, SCATTER_FORWARD, ierr)
            CALL VecGetArrayF90(self%vec_mz_ghost, x_loc, ierr)
            vv_out(:,i,j) = x_loc(:)
            CALL VecRestoreArrayF90(self%vec_mz_ghost, x_loc, ierr)
         END DO
      END DO
      DEALLOCATE(ix)

   END SUBROUTINE apply_symmetry

  ! This routine is called in initialization.f90 s.t. vv_mz_LA(PUBLIC) = ltg_LA

  SUBROUTINE symmetric_points(mesh_loc, communicator, ltg_LA)
      ! USE def_type_mesh
      USE my_util

      IMPLICIT NONE

      TYPE(mesh_type)               :: mesh_loc
      INTEGER, DIMENSION(:)         ::  ltg_LA
     INTEGER, DIMENSION(mesh_loc%np)      :: i_d
      INTEGER             :: i, j, m, dom, n_w
      REAL(KIND=8)        :: xx, zz, eps

      INTEGER    :: code, nb_procs, rank
      MPI_Comm   :: communicator
      INTEGER, DIMENSION(:), ALLOCATABLE :: loc_to_glob, mesh_glob_i_d
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: mesh_glob_rr
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: mesh_glob_jj
      INTEGER :: mesh_glob_me, mesh_glob_np, index_rr
      INTEGER, DIMENSION(:), ALLOCATABLE :: size_per_proc_me, cumulated_size_per_proc_me
      INTEGER, DIMENSION(:), ALLOCATABLE :: size_per_proc_np, cumulated_size_per_proc_np

      
      CALL MPI_COMM_SIZE(communicator,nb_procs,code)
      CALL MPI_COMM_RANK(communicator,rank,code)

      need_sym=.TRUE.
      eps = 1.d-8
      ltg_LA = -1
      n_w = SIZE(mesh_loc%jj,1)

   !============================================ gather me and np from global mesh ============================================
      !=== me
      ALLOCATE(cumulated_size_per_proc_me(0:nb_procs), SOURCE=0)
      ALLOCATE(size_per_proc_me(nb_procs), SOURCE=0)

      CALL MPI_ALLGATHER(mesh_loc%me, 1, MPI_INTEGER, size_per_proc_me, 1, &
                        MPI_INTEGER, communicator, code)
      mesh_glob_me = SUM(size_per_proc_me)
      cumulated_size_per_proc_me = 0
      DO i = 1, nb_procs
         cumulated_size_per_proc_me(i:) = cumulated_size_per_proc_me(i:) + size_per_proc_me(i)
      END DO

      !=== np
      ALLOCATE(cumulated_size_per_proc_np(0:nb_procs), SOURCE=0)
      ALLOCATE(size_per_proc_np(nb_procs), SOURCE=0)

      CALL MPI_ALLGATHER(mesh_loc%np, 1, MPI_INTEGER, size_per_proc_np, 1, &
                        MPI_INTEGER, communicator, code)
      mesh_glob_np = SUM(size_per_proc_np)
      DO i = 1, nb_procs
         cumulated_size_per_proc_np(i:) = cumulated_size_per_proc_np(i:) + size_per_proc_np(i)
      END DO


   !============================================ building mesh_glob_i_d ============================================

      ALLOCATE(mesh_glob_i_d(mesh_glob_me),SOURCE=0)

      CALL MPI_ALLGATHERV(mesh_loc%i_d, mesh_loc%me, MPI_INTEGER, &
                        mesh_glob_i_d, size_per_proc_me, cumulated_size_per_proc_me(0:nb_procs-1), &
                        MPI_INTEGER, communicator, code)

      DO m=1, mesh_loc%me
         i_d(mesh_loc%jj(:,m)) = mesh_loc%i_d(m)
      END DO

   !============================================ send mesh%rr from local to global ============================================

      ALLOCATE(mesh_glob_rr(2, mesh_glob_np), SOURCE=0.d0)

      CALL MPI_ALLGATHERV(mesh_loc%rr, 2*mesh_loc%np, MPI_DOUBLE_PRECISION, &
                        mesh_glob_rr, 2*size_per_proc_np, 2*cumulated_size_per_proc_np(0:nb_procs-1), &
                        MPI_DOUBLE_PRECISION, communicator, code)

   !============================================ send mesh%jj from local to global ============================================

      ALLOCATE(mesh_glob_jj(n_w, mesh_glob_me),SOURCE=0)

      CALL MPI_ALLGATHERV(mesh_loc%jj+cumulated_size_per_proc_np(rank), n_w*mesh_loc%me, MPI_INTEGER, &
                        mesh_glob_jj, n_w*size_per_proc_me, n_w*cumulated_size_per_proc_me(0:nb_procs-1), &
                        MPI_INTEGER, communicator, code)

   !============================================ send mesh%loc_to_glob from local to global ============================================

      ALLOCATE(loc_to_glob(mesh_glob_np),SOURCE=0)

      CALL MPI_ALLGATHERV(mesh_loc%loc_to_glob, mesh_loc%np, MPI_INTEGER, &
                     loc_to_glob, size_per_proc_np, cumulated_size_per_proc_np(0:nb_procs-1), &
                     MPI_INTEGER, communicator, code)

   !=====================================  CHECKS (to be removed) ============================================

      IF (SIZE(loc_to_glob) /= SIZE(mesh_glob_rr, 2)) THEN
         IF (rank == 0) WRITE(*,*) "size loc_to_glob is ", SIZE(loc_to_glob), "size rr global is ", SIZE(mesh_glob_rr, 2)
         CALL error_petsc("BUG in symmetric points: size loc_to_glob and size rr global do not match")
      END IF

      IF (MAXVAL(mesh_glob_jj) /= SIZE(loc_to_glob)) THEN
         IF (rank == 0) WRITE(*,*) "max jj global is ", MAXVAL(mesh_glob_jj), "size loc_to_glob is ", SIZE(loc_to_glob)
         CALL error_petsc("BUG in symmetric points: size loc_to_glob and maxval jj global do not match")
      END IF

      IF (SIZE(mesh_glob_jj, 2) /= mesh_glob_me) THEN
         IF (rank == 0) WRITE(*,*) "size jj global is ", SIZE(mesh_glob_jj, 2), "mesh_glob_me is ", mesh_glob_me
         CALL error_petsc("BUG in symmetric points: size loc_to_glob and maxval jj global do not match")
      END IF

      IF (MAXVAL(loc_to_glob) > mesh_glob_np) THEN
         IF (rank == 0) WRITE(*,*) "max loc_to_glob is ", MAXVAL(loc_to_glob), "glob np is ", mesh_glob_np
         CALL error_petsc("BUG in symmetric points: max loc_to_glob and glob_np do not match")
      END IF

      IF (MAXVAL(loc_to_glob) > SIZE(mesh_glob_rr, 2)) THEN
         IF (rank == 0) WRITE(*,*) "max loc_to_glob is ", MAXVAL(loc_to_glob), "size rr global is ", SIZE(mesh_glob_rr, 2)
         CALL error_petsc("BUG in symmetric points: max loc_to_glob is > size rr global")
      END IF

      IF (MINVAL(loc_to_glob) /= 1) THEN
         IF (rank == 0) WRITE(*,*) "min loc_to_glob is ", MINVAL(loc_to_glob)
         CALL error_petsc("BUG in symmetric points: min loc_to_glob is /= 1")
      END IF

   !=====================================  SYMMETRIZING POINTS about z = 0 ============================================

      DO i = 1, mesh_loc%np
         xx =  mesh_loc%rr(1,i)
         zz = -mesh_loc%rr(2,i)
         DO m = 1, mesh_glob_me
            !=== test necessary for jumps (i.e of mu at interfaces in H)
            IF (mesh_glob_i_d(m) /= i_d(i)) CYCLE
            !=== test necessary for jumps (i.e of mu at interfaces in H)
            DO j = 1, n_w
               index_rr = mesh_glob_jj(j,m)
               IF (ABS(xx-mesh_glob_rr(1,index_rr))+ABS(zz-mesh_glob_rr(2,index_rr)) .LT. eps) THEN
                  ltg_LA(i) = loc_to_glob(index_rr)
                  EXIT
               END IF
            END DO
         END DO
         IF (ltg_LA(i) < 0) THEN
            WRITE(*,*) "rank=", rank, ": BUG in symmetric points: did not match point ", i, " with coordinates ", xx, zz
         END IF
      END DO

      IF (MINVAL(ltg_LA) < 0 .AND. (mesh_loc%me/=0)) THEN
         CALL error_petsc("BUG in symmetric points: some points left unmatched")
      END IF

  END SUBROUTINE symmetric_points

  SUBROUTINE symm_champ(vv_in, vv_out, if_u_h_T_c)
    USE my_util
    IMPLICIT NONE

    TYPE(mesh_type)                             :: mesh
    TYPE(symmetric_tools), POINTER              :: abstract_sym_tool
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: vv_in
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: vv_out
    CHARACTER(len=*), INTENT(IN)                :: if_u_h_T_c != 'u' (velocity) 'h', 'T' (temp), 'c' (conc)
    INTEGER,          POINTER, DIMENSION(:)     :: ifrom
    INTEGER                                     :: n, i, j
    INTEGER, ALLOCATABLE, DIMENSION(:)          :: ix

    IF (if_u_h_T_c=='u') THEN
       abstract_sym_tool => vv_sym_tool
    ELSE IF (if_u_h_T_c=='h') THEN
       abstract_sym_tool => H_sym_tool
    ELSE IF (if_u_h_T_c=='T') THEN
       abstract_sym_tool => T_sym_tool
    ELSE IF (if_u_h_T_c=='c') THEN
       abstract_sym_tool => c_sym_tool
    ELSE IF (if_u_h_T_c=='phi') THEN
       abstract_sym_tool => phi_sym_tool
    ELSE
       CALL error_petsc("BUG in symm_champ: if_u_h_T_c should be 'u', 'h', 'T', 'c', or 'phi', not "&
        // trim(adjustl(if_u_h_T_c)))
    END IF

!=== do not attempt symmetrization here
    IF (abstract_sym_tool%mesh%me == 0) RETURN
    vv_out = 0.d0
    IF (.NOT. need_sym) RETURN
!=== do not attempt symmetrization here

    CALL abstract_sym_tool%apply_symmetry(vv_in, vv_out)
  END SUBROUTINE  symm_champ

  SUBROUTINE val_ener_sym_centrale(communicator, mesh, list_mode, v_in, e_mode, e_mode_sym, e_mode_anti, if_u_h_T_c)
    !type_sym = 1 pour un champ pair
    !type_sym = -1 pour un champ impair

    USE associate_gauss
   !  USE def_type_mesh
    USE tn_axi

    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                :: mesh
    INTEGER, DIMENSION(:),                      INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:)             ,INTENT(IN) :: v_in      !champ
    REAL(KIND=8), DIMENSION(mesh%np, SIZE(v_in,2),SIZE(list_mode)) :: v_sym
    REAL(KIND=8), DIMENSION(SIZE(list_mode))   ,INTENT(OUT):: e_mode, e_mode_sym, e_mode_anti !energie par mode
    CHARACTER(len=*), INTENT(IN)                           :: if_u_h_T_c
    REAL(KIND=8),DIMENSION(3)                              :: type_sym
    !type_sym(r,theta,z)
    REAL(KIND=8), DIMENSION(mesh%np,size(v_in,2),size(list_mode))  :: champ_anti, champ_sym      !champ

    INTEGER,      DIMENSION(:,:), POINTER                  :: jj
    INTEGER,                      POINTER                  :: me
    INTEGER      :: i
    INTEGER      :: m_max_c
    MPI_Comm       :: communicator


    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me
    m_max_c = size(list_mode)
    e_mode      = 0.d0
    e_mode_sym  = 0.d0
    e_mode_anti = 0.d0

    CALL  symm_champ(v_in, v_sym, if_u_h_T_c)

    !===Compute anti-symmetric field
    DO i=1, size(list_mode)
       IF (mod(list_mode(i),2) == 0) THEN  !mode pair
          type_sym(1) = 1.d0
          type_sym(2) = 1.d0
          type_sym(3) = -1.d0
       ELSE
          type_sym(1) = -1.d0
          type_sym(2) = -1.d0
          type_sym(3) = 1.d0
       ENDIF
       champ_anti(:,1:2,i) =  0.5d0*(v_in(:,1:2,i) - type_sym(1)*v_sym(:,1:2,i))
       champ_sym(:,1:2,i)  =  0.5d0*(v_in(:,1:2,i) + type_sym(1)*v_sym(:,1:2,i))
       IF (SIZE(v_in,2) == 6) THEN
          champ_anti(:,3:4,i) =  0.5d0*(v_in(:,3:4,i) - type_sym(2)*v_sym(:,3:4,i))
          champ_sym(:,3:4,i)  =  0.5d0*(v_in(:,3:4,i) + type_sym(2)*v_sym(:,3:4,i))
          champ_anti(:,5:6,i) =  0.5d0*(v_in(:,5:6,i) - type_sym(3)*v_sym(:,5:6,i))
          champ_sym(:,5:6,i)  =  0.5d0*(v_in(:,5:6,i) + type_sym(3)*v_sym(:,5:6,i))
       ENDIF
    ENDDO

    !===Compute energies
    DO i=1, m_max_c
       e_mode(i) = 0.5*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), v_in(:,:,i:i)))**2
       e_mode_anti(i) = 0.5*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), champ_anti(:,:,i:i)))**2
       e_mode_sym(i) = 0.5*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), champ_sym(:,:,i:i)))**2
    ENDDO
!!$    e_mode_sym = e_mode - e_mode_anti

  END SUBROUTINE val_ener_sym_centrale


  SUBROUTINE val_ener_sym_glob(communicator, mesh, list_mode, v_in, e_mode, e_mode_sym, e_mode_anti, type_sym, if_u_h_T_c)
    !type_sym = 1 pour un champ pair
    !type_sym = -1 pour un champ impair

    USE associate_gauss
    USE tn_axi

    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                :: mesh
    INTEGER, DIMENSION(:),                      INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:)             ,INTENT(IN) :: v_in      !champ
    REAL(KIND=8), DIMENSION(mesh%np, SIZE(v_in,2),SIZE(list_mode)) :: v_sym
    REAL(KIND=8), DIMENSION(SIZE(list_mode))   ,INTENT(OUT):: e_mode, e_mode_sym, e_mode_anti !energie par mode
    REAL(KIND=8), DIMENSION(3), INTENT(IN)                 :: type_sym
    CHARACTER(len=*), INTENT(IN)                           :: if_u_h_T_c
    !type_sym(r,theta,z)
    REAL(KIND=8), DIMENSION(mesh%np,size(v_in,2),size(list_mode))  :: champ_anti, champ_sym      !champ

    INTEGER,      DIMENSION(:,:), POINTER                  :: jj
    INTEGER,                      POINTER                  :: me
    INTEGER      :: i, mode
    INTEGER      :: m_max_c
    MPI_Comm       :: communicator


    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me
    m_max_c = size(list_mode)
    e_mode      = 0.d0
    e_mode_sym  = 0.d0
    e_mode_anti = 0.d0

    CALL  symm_champ(v_in, v_sym, if_u_h_T_c)

    !===Compute anti-symmetric field
   champ_anti(:,1:2,:) =  0.5d0*(v_in(:,1:2,:) - type_sym(1)*v_sym(:,1:2,:))
   champ_anti(:,3:4,:) =  0.5d0*(v_in(:,3:4,:) - type_sym(2)*v_sym(:,3:4,:))
   champ_anti(:,5:6,:) =  0.5d0*(v_in(:,5:6,:) - type_sym(3)*v_sym(:,5:6,:))
   champ_sym(:,1:2,:) =  0.5d0*(v_in(:,1:2,:) + type_sym(1)*v_sym(:,1:2,:))
   champ_sym(:,3:4,:) =  0.5d0*(v_in(:,3:4,:) + type_sym(2)*v_sym(:,3:4,:))
   champ_sym(:,5:6,:) =  0.5d0*(v_in(:,5:6,:) + type_sym(3)*v_sym(:,5:6,:))

    !===Compute energies
    DO i=1, m_max_c
       e_mode(i) = 0.5*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), v_in(:,:,i:i)))**2
       e_mode_anti(i) = 0.5*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), champ_anti(:,:,i:i)))**2
       e_mode_sym(i) = 0.5*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), champ_sym(:,:,i:i)))**2
    ENDDO


  END SUBROUTINE val_ener_sym_glob

  SUBROUTINE val_ener_sym_rpi(communicator, mesh, list_mode, v_in, e_mode, e_mode_sym, e_mode_anti, if_u_h_T_c)
    ! SYMETRIE Rpi
    !    type_sym = 1.d0
    !    type_sym =-1.d0
    !    type_sym =-1.d0
    !    type_sym = 1.d0
    !    type_sym =-1.d0
    !    type_sym = 1.d0

    USE associate_gauss
   !  USE def_type_mesh
    USE tn_axi

    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                :: mesh
    INTEGER, DIMENSION(:),                      INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:)             ,INTENT(IN) :: v_in      !champ
    REAL(KIND=8), DIMENSION(mesh%np, SIZE(v_in,2),SIZE(list_mode)) :: v_sym
    REAL(KIND=8), DIMENSION(SIZE(list_mode))   ,INTENT(OUT):: e_mode, e_mode_sym, e_mode_anti !energie par mode
    REAL(KIND=8),DIMENSION(6)                              :: type_sym
    CHARACTER(len=*), INTENT(IN)                           :: if_u_h_T_c
    !type_sym(r,theta,z)
    REAL(KIND=8), DIMENSION(mesh%np,size(v_in,2),size(list_mode))  :: champ_anti, champ_sym      !champ

    INTEGER,      DIMENSION(:,:), POINTER                  :: jj
    INTEGER,                      POINTER                  :: me
    INTEGER      :: i, k
    INTEGER      :: m_max_c
    !#include "petsc/finclude/petsc.h"
    MPI_Comm       :: communicator


    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me
    m_max_c = size(list_mode)
    e_mode      = 0.d0
    e_mode_sym  = 0.d0
    e_mode_anti = 0.d0
    CALL  symm_champ(v_in, v_sym, if_u_h_T_c)
!VB 07/04/2026
    type_sym(1) = 1.d0
    type_sym(2) = -1.d0
    type_sym(3) = -1.d0
    type_sym(4) = 1.d0
    type_sym(5) = -1.d0
    type_sym(6) = 1.d0
!VB 07/04/2026

    !===Compute anti-symmetric field
    DO k = 1, SIZE(v_in, 2)
       champ_anti(:,k,:) =  0.5d0*(v_in(:,k,:) - type_sym(k)*v_sym(:,k,:))
       champ_sym(:,k,:) =  0.5d0*(v_in(:,k,:) + type_sym(k)*v_sym(:,k,:))
    ENDDO

    !===Compute energies
    DO i=1, m_max_c
       e_mode(i) = 0.5*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), v_in(:,:,i:i)))**2
       e_mode_anti(i) = 0.5*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), champ_anti(:,:,i:i)))**2
       e_mode_sym(i) = 0.5*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), champ_sym(:,:,i:i)))**2
    ENDDO

  END SUBROUTINE val_ener_sym_rpi

  SUBROUTINE val_ener_north_south(communicator, mesh, list_mode, v_in, e_north, e_south, e_tot)

    USE associate_gauss
   !  USE def_type_mesh
    USE tn_axi

    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                :: mesh
    INTEGER, DIMENSION(:),                      INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:)             ,INTENT(IN) :: v_in      !champ
    REAL(KIND=8), DIMENSION(mesh%np, SIZE(v_in,2),SIZE(list_mode)) :: vv_n, vv_s
    REAL(KIND=8), DIMENSION(SIZE(list_mode))   ,INTENT(OUT):: e_north, e_south !energie par mode
    REAL(KIND=8), DIMENSION(SIZE(list_mode)), OPTIONAL     :: e_tot

    REAL(KIND=8) :: ray, pi
    REAL(KIND=8), DIMENSION(2) :: e_ns
    REAL(KIND=8) :: epsilon = 1.d-10
    INTEGER      :: i, k, j, l, m
    INTEGER      :: m_max_c
    INTEGER      :: code
    !#include "petsc/finclude/petsc.h"
    MPI_Comm       :: communicator

    m_max_c = size(list_mode)
    vv_n     = 0.d0
    vv_s     = 0.d0
    e_north  = 0.d0
    e_south  = 0.d0
    pi = ACOS(-1.d0)

    DO i = 1, SIZE(list_mode)
       e_ns = 0.d0
       DO m = 1, mesh%me
          IF (MINVAL(mesh%rr(2,mesh%jj(:,m))) > -epsilon) THEN ! l'element est au nord
             j = 1
          ELSE IF (MAXVAL(mesh%rr(2,mesh%jj(:,m))) < epsilon) THEN ! l'element est au sud
             j = 2
          ELSE
             WRITE(*,*) 'Attention, element a cheval entre nord et sud'
             CYCLE
          END IF
          DO l = 1, mesh%gauss%l_G
             ray = SUM(mesh%rr(1,mesh%jj(:,m))*mesh%gauss%ww(:,l))
             DO k = 1, SIZE(v_in,2)
                e_ns(j) = e_ns(j) + ray*SUM(v_in(mesh%jj(:,m),k,i)* mesh%gauss%ww(:,l))**2*mesh%gauss%rj(l,m)
             END DO
          END DO
       END DO

       IF (list_mode(i) /= 0) THEN
          e_ns = 0.5d0*e_ns
       END IF

       CALL MPI_ALLREDUCE(e_ns(1),e_north(i),1,MPI_DOUBLE_PRECISION, MPI_SUM, communicator, code)
       CALL MPI_ALLREDUCE(e_ns(2),e_south(i),1,MPI_DOUBLE_PRECISION, MPI_SUM, communicator, code)
       IF (PRESENT(e_tot)) THEN
          e_tot(i) = 0.5d0*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), v_in(:,:,i:i)))**2
       END IF
    END DO

    e_north = pi*e_north  ! 1/2 * (2pi)
    e_south = pi*e_south

  END SUBROUTINE val_ener_north_south

  SUBROUTINE champ_total_anti_sym(communicator, mesh, list_mode, eps_sa, v_in, v_out, if_u_h_T_c)
    !type_sym = 1 pour un champ pair
    !type_sym = -1 pour un champ impair

    USE associate_gauss
   !  USE def_type_mesh
    USE tn_axi

    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                                 :: mesh
    INTEGER, DIMENSION(:),                                       INTENT(IN) :: list_mode
    REAL(KIND=8),                                                INTENT(IN) :: eps_sa ! eps_sa=+1 (anti), eps_sa=-1 (sym)
    REAL(KIND=8), DIMENSION(:,:,:)                              ,INTENT(IN) :: v_in      !champ
    CHARACTER(len=*), INTENT(IN)                           :: if_u_h_T_c
    REAL(KIND=8), DIMENSION(mesh%np, SIZE(v_in,2),SIZE(list_mode))             :: v_sym
    REAL(KIND=8),DIMENSION(3)                                               :: type_sym
    !type_sym(r,theta,z)
    REAL(KIND=8), DIMENSION(mesh%np,size(v_in,2),size(list_mode)), INTENT(OUT)  :: v_out      !champ antisymetrise

    INTEGER,      DIMENSION(:,:), POINTER                  :: jj
    INTEGER,                      POINTER                  :: me
    INTEGER      :: i
    INTEGER      :: m_max_c
    MPI_Comm       :: communicator


    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me
    m_max_c = size(list_mode)
    v_out = 0.d0

    CALL  symm_champ(v_in, v_sym, if_u_h_T_c)

    !===Compute anti-symmetric field
    DO i=1, size(list_mode)
       IF (mod(list_mode(i),2) == 0) THEN  !mode pair
          type_sym(1) = 1.d0
          type_sym(2) = 1.d0
          type_sym(3) = -1.d0
       ELSE
          type_sym(1) = -1.d0
          type_sym(2) = -1.d0
          type_sym(3) = 1.d0
       ENDIF
       v_out(:,1:2,i) =  0.5d0*(v_in(:,1:2,i) - eps_sa*type_sym(1)*v_sym(:,1:2,i))
       IF (size(v_in,2) == 6) THEN
         v_out(:,3:4,i) =  0.5d0*(v_in(:,3:4,i) - eps_sa*type_sym(2)*v_sym(:,3:4,i))
         v_out(:,5:6,i) =  0.5d0*(v_in(:,5:6,i) - eps_sa*type_sym(3)*v_sym(:,5:6,i))
       ENDIF
    ENDDO

  END SUBROUTINE champ_total_anti_sym
END MODULE symmetric_field
