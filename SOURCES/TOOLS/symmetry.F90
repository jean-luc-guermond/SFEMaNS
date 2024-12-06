MODULE symmetric_field
#include "petsc/finclude/petsc.h"
  USE petsc
  IMPLICIT NONE

  PUBLIC :: symmetric_points, symm_champ, val_ener_sym_centrale, &
       val_ener_sym_glob,val_ener_sym_rpi, val_ener_north_south, champ_total_anti_sym

  !Symmetrization-------------------------------------------------------------
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC             :: vv_mz_LA, H_mz_LA
  LOGICAL, PRIVATE                                       :: need_sym=.FALSE.
  !---------------------------------------------------------------------------
  PRIVATE

CONTAINS

  ! This routine is called in initialization.f90 s.t. vv_mz_LA(PUBLIC) = ltg_LA
  SUBROUTINE symmetric_points(mesh_loc, mesh_glob, ltg_LA)
    USE def_type_mesh
    IMPLICIT NONE

    TYPE(mesh_type)                      :: mesh_loc, mesh_glob
    INTEGER, DIMENSION(:)                ::  ltg_LA
    INTEGER, DIMENSION(mesh_loc%np)      :: i_d
    INTEGER             :: i, j, m, dom, n_w
    REAL(KIND=8)        :: xx, zz, epsilon

    need_sym=.TRUE.
    epsilon = 1.d-8
    ltg_LA = -1
    n_w = SIZE(mesh_loc%jj,1)

    DO m=1, mesh_loc%me
       i_d(mesh_loc%jj(:,m)) = mesh_loc%i_d(m)
    END DO

    DO i = 1, mesh_loc%np
       xx =  mesh_loc%rr(1,i)
       zz = -mesh_loc%rr(2,i)
       dom = i_d(i)
       DO m = 1, mesh_glob%me
          IF (mesh_glob%i_d(m) /= dom) CYCLE
          DO j = 1, n_w
             IF (ABS(xx-mesh_glob%rr(1,mesh_glob%jj(j,m)))+ABS(zz-mesh_glob%rr(2,mesh_glob%jj(j,m))) .LT. epsilon) THEN
                ltg_LA(i) = mesh_glob%jj(j,m)
                EXIT
             END IF
          END DO
       END DO
    END DO

    IF (MINVAL(ltg_LA)<0 .AND. (mesh_loc%me/=0)) THEN
       WRITE(*,*) 'bug in symmetric_points'
    END IF

  END SUBROUTINE symmetric_points

  SUBROUTINE symm_champ(communicator, vv_in, mesh, vv_out, if_u_h)
    USE def_type_mesh
    USE st_matrix

    IMPLICIT NONE

    TYPE(mesh_type)                             :: mesh
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: vv_in
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: vv_out
    CHARACTER(len=1), INTENT(IN)                :: if_u_h
    INTEGER,          POINTER, DIMENSION(:)     :: ifrom
    INTEGER                                     :: n, i, j
    type(petsc_csr_LA), SAVE                    :: LA_u, LA_HH
    INTEGER, ALLOCATABLE, DIMENSION(:)          :: ix

    LOGICAL, SAVE                               :: once_u=.true.
    LOGICAL, SAVE                               :: once_h=.true.
    !#include "petsc/finclude/petscvec.h90"
    Vec, SAVE            :: vv_mz, vv_mz_ghost, H_mz, H_mz_ghost
    PetscErrorCode       :: ierr
    MPI_Comm             :: communicator
    PetscScalar, POINTER :: x_loc(:)


    IF (mesh%me ==0) RETURN

    vv_out = 0.d0
    IF (.NOT. need_sym) RETURN

    IF ((once_u) .AND. (if_u_h=='u')) THEN
       ALLOCATE(LA_u%loc_to_glob(1,mesh%np))
       LA_u%loc_to_glob(1,:) = mesh%loc_to_glob
       LA_u%kmax = 1
       CALL create_my_ghost(mesh,LA_u,ifrom)
       n = mesh%dom_np
       CALL VecCreateGhost(communicator, n, &
            PETSC_DETERMINE, SIZE(ifrom), ifrom, vv_mz, ierr)
       CALL VecGhostGetLocalForm(vv_mz, vv_mz_ghost, ierr)
       once_u = .false.
    END IF

    IF ((once_h) .AND. (if_u_h=='h')) THEN
       ALLOCATE(LA_HH%loc_to_glob(1,mesh%np))
       LA_HH%loc_to_glob(1,:) = mesh%loc_to_glob
       LA_HH%kmax = 1
       CALL create_my_ghost(mesh,LA_HH,ifrom)
       n = mesh%dom_np
       CALL VecCreateGhost(communicator, n, &
            PETSC_DETERMINE, SIZE(ifrom), ifrom, H_mz, ierr)
       CALL VecGhostGetLocalForm(H_mz, H_mz_ghost, ierr)
       once_h = .false.
    END IF

    ! ix begins to 0 because of PETSC and contains the global information for symmetric points
    ALLOCATE(ix(mesh%np))
    IF (if_u_h == 'u') THEN
       ix = vv_mz_LA-1
       DO i = 1, SIZE(vv_in,2)
          DO j = 1, SIZE(vv_in,3)
             CALL VecZeroEntries(vv_mz, ierr)
             CALL VecSetValues(vv_mz, mesh%np, ix, vv_in(:,i,j), INSERT_VALUES, ierr)
             CALL VecAssemblyBegin(vv_mz, ierr)
             CALL VecAssemblyEnd(vv_mz, ierr)
             CALL VecGhostUpdateBegin(vv_mz,INSERT_VALUES, SCATTER_FORWARD,ierr)
             CALL VecGhostUpdateEnd(vv_mz,INSERT_VALUES, SCATTER_FORWARD, ierr)
             CALL VecGetArrayF90(vv_mz_ghost, x_loc, ierr)
             vv_out(:,i,j) = x_loc(:)
             CALL VecRestoreArrayF90(vv_mz_ghost, x_loc, ierr)

          END DO
       END DO
    ELSE IF (if_u_h=='h') THEN
       ix = H_mz_LA-1
       DO i = 1, SIZE(vv_in,2)
          DO j = 1, SIZE(vv_in,3)
             CALL VecZeroEntries(H_mz, ierr)
             CALL VecSetValues(H_mz, mesh%np, ix, vv_in(:,i,j), INSERT_VALUES, ierr)
             CALL VecAssemblyBegin(H_mz, ierr)
             CALL VecAssemblyEnd(H_mz, ierr)
             CALL VecGhostUpdateBegin(H_mz,INSERT_VALUES, SCATTER_FORWARD,ierr)
             CALL VecGhostUpdateEnd(H_mz,INSERT_VALUES, SCATTER_FORWARD, ierr)
             CALL VecGetArrayF90(H_mz_ghost, x_loc, ierr)
             vv_out(:,i,j) = x_loc(:)
             CALL VecRestoreArrayF90(H_mz_ghost, x_loc, ierr)
          END DO
       END DO
    END IF
    DEALLOCATE(ix)
  END SUBROUTINE  symm_champ

  SUBROUTINE val_ener_sym_centrale(communicator, mesh, list_mode, v, e_mode, e_mode_sym, e_mode_anti, if_u_h)
    !type_sym = 1 pour un champ pair
    !type_sym = -1 pour un champ impair

    USE Gauss_points
    USE def_type_mesh
    USE tn_axi

    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                :: mesh
    INTEGER, DIMENSION(:),                      INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:)             ,INTENT(IN) :: v      !champ
    REAL(KIND=8), DIMENSION(mesh%np, SIZE(v,2),SIZE(list_mode)) :: vv
    REAL(KIND=8), DIMENSION(SIZE(list_mode))   ,INTENT(OUT):: e_mode, e_mode_sym, e_mode_anti !energie par mode
    CHARACTER(len=1), INTENT(IN)                           :: if_u_h
    REAL(KIND=8),DIMENSION(3)                              :: type_sym
    !type_sym(r,theta,z)
    REAL(KIND=8), DIMENSION(mesh%np,size(v,2),size(list_mode))  :: champ_anti, champ_sym      !champ

    INTEGER,      DIMENSION(:,:), POINTER                  :: jj
    INTEGER,                      POINTER                  :: me
    INTEGER      :: i
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

    CALL  symm_champ(communicator, v, mesh, vv, if_u_h)

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
       champ_anti(:,1:2,i) =  0.5d0*(v(:,1:2,i) - type_sym(1)*vv(:,1:2,i))
       champ_anti(:,3:4,i) =  0.5d0*(v(:,3:4,i) - type_sym(2)*vv(:,3:4,i))
       champ_anti(:,5:6,i) =  0.5d0*(v(:,5:6,i) - type_sym(3)*vv(:,5:6,i))
       champ_sym(:,1:2,i) =  0.5d0*(v(:,1:2,i) + type_sym(1)*vv(:,1:2,i))
       champ_sym(:,3:4,i) =  0.5d0*(v(:,3:4,i) + type_sym(2)*vv(:,3:4,i))
       champ_sym(:,5:6,i) =  0.5d0*(v(:,5:6,i) + type_sym(3)*vv(:,5:6,i))
    ENDDO

    !===Compute energies
    DO i=1, m_max_c
       e_mode(i) = 0.5*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), v(:,:,i:i)))**2
       e_mode_anti(i) = 0.5*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), champ_anti(:,:,i:i)))**2
       e_mode_sym(i) = 0.5*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), champ_sym(:,:,i:i)))**2
    ENDDO
!!$    e_mode_sym = e_mode - e_mode_anti


  END SUBROUTINE val_ener_sym_centrale

  SUBROUTINE val_ener_sym_glob(communicator, mesh, list_mode, v, e_mode, e_mode_sym, e_mode_anti, type_sym, if_u_h)
    !type_sym = 1 pour un champ pair
    !type_sym = -1 pour un champ impair

    USE Gauss_points
    USE def_type_mesh
    USE tn_axi

    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                :: mesh
    INTEGER, DIMENSION(:),                      INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:)             ,INTENT(IN) :: v      !champ
    REAL(KIND=8), DIMENSION(mesh%np, SIZE(v,2),SIZE(list_mode)) :: vv
    REAL(KIND=8), DIMENSION(SIZE(list_mode))   ,INTENT(OUT):: e_mode, e_mode_sym, e_mode_anti !energie par mode
    REAL(KIND=8),DIMENSION(3),                  INTENT(IN) :: type_sym
    CHARACTER(len=1), INTENT(IN)                           :: if_u_h
    !type_sym(r,theta,z)
    REAL(KIND=8), DIMENSION(mesh%np,size(v,2),size(list_mode))  :: champ_anti, champ_sym      !champ

    INTEGER,      DIMENSION(:,:), POINTER                  :: jj
    INTEGER,                      POINTER                  :: me
    INTEGER      :: i
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

    CALL  symm_champ(communicator, v, mesh, vv, if_u_h)

    !===Compute anti-symmetric field
    champ_anti(:,1:2,:) =  0.5d0*(v(:,1:2,:) - type_sym(1)*vv(:,1:2,:))
    champ_anti(:,3:4,:) =  0.5d0*(v(:,3:4,:) - type_sym(2)*vv(:,3:4,:))
    champ_anti(:,5:6,:) =  0.5d0*(v(:,5:6,:) - type_sym(3)*vv(:,5:6,:))
    champ_sym(:,1:2,:) =  0.5d0*(v(:,1:2,:) + type_sym(1)*vv(:,1:2,:))
    champ_sym(:,3:4,:) =  0.5d0*(v(:,3:4,:) + type_sym(2)*vv(:,3:4,:))
    champ_sym(:,5:6,:) =  0.5d0*(v(:,5:6,:) + type_sym(3)*vv(:,5:6,:))

    !===Compute energies
    DO i=1, m_max_c
       e_mode(i) = 0.5*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), v(:,:,i:i)))**2
       e_mode_anti(i) = 0.5*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), champ_anti(:,:,i:i)))**2
       e_mode_sym(i) = 0.5*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), champ_sym(:,:,i:i)))**2
    ENDDO


  END SUBROUTINE val_ener_sym_glob

  SUBROUTINE val_ener_sym_rpi(communicator, mesh, list_mode, v, e_mode, e_mode_sym, e_mode_anti, type_sym, if_u_h)
    ! SYMETRIE Rpi
    !    type_sym = 1.d0
    !    type_sym =-1.d0
    !    type_sym =-1.d0
    !    type_sym = 1.d0
    !    type_sym =-1.d0
    !    type_sym = 1.d0

    USE Gauss_points
    USE def_type_mesh
    USE tn_axi

    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                :: mesh
    INTEGER, DIMENSION(:),                      INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:)             ,INTENT(IN) :: v      !champ
    REAL(KIND=8), DIMENSION(mesh%np, SIZE(v,2),SIZE(list_mode)) :: vv
    REAL(KIND=8), DIMENSION(SIZE(list_mode))   ,INTENT(OUT):: e_mode, e_mode_sym, e_mode_anti !energie par mode
    REAL(KIND=8),DIMENSION(6),                  INTENT(IN) :: type_sym
    CHARACTER(len=1), INTENT(IN)                           :: if_u_h
    !type_sym(r,theta,z)
    REAL(KIND=8), DIMENSION(mesh%np,size(v,2),size(list_mode))  :: champ_anti, champ_sym      !champ

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

    CALL  symm_champ(communicator, v, mesh, vv, if_u_h)

    !===Compute anti-symmetric field
    DO k = 1,6
       champ_anti(:,k,:) =  0.5d0*(v(:,k,:) - type_sym(k)*vv(:,k,:))
       champ_sym(:,k,:) =  0.5d0*(v(:,k,:) + type_sym(k)*vv(:,k,:))
    ENDDO

    !===Compute energies
    DO i=1, m_max_c
       e_mode(i) = 0.5*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), v(:,:,i:i)))**2
       e_mode_anti(i) = 0.5*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), champ_anti(:,:,i:i)))**2
       e_mode_sym(i) = 0.5*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), champ_sym(:,:,i:i)))**2
    ENDDO


  END SUBROUTINE val_ener_sym_rpi

  SUBROUTINE val_ener_north_south(communicator, mesh, list_mode, v, e_north, e_south, e_tot)

    USE Gauss_points
    USE def_type_mesh
    USE tn_axi

    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                :: mesh
    INTEGER, DIMENSION(:),                      INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:)             ,INTENT(IN) :: v      !champ
    REAL(KIND=8), DIMENSION(mesh%np, SIZE(v,2),SIZE(list_mode)) :: vv_n, vv_s
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
             DO k = 1, SIZE(v,2)
                e_ns(j) = e_ns(j) + ray*SUM(v(mesh%jj(:,m),k,i)* mesh%gauss%ww(:,l))**2*mesh%gauss%rj(l,m)
             END DO
          END DO
       END DO

       IF (list_mode(i) /= 0) THEN
          e_ns = 0.5d0*e_ns
       END IF

       CALL MPI_ALLREDUCE(e_ns(1),e_north(i),1,MPI_DOUBLE_PRECISION, MPI_SUM, communicator, code)
       CALL MPI_ALLREDUCE(e_ns(2),e_south(i),1,MPI_DOUBLE_PRECISION, MPI_SUM, communicator, code)
       IF (PRESENT(e_tot)) THEN
          e_tot(i) = 0.5d0*(norme_L2_champ_par(communicator, mesh, list_mode(i:i), v(:,:,i:i)))**2
       END IF
    END DO

    e_north = pi*e_north  ! 1/2 * (2pi)
    e_south = pi*e_south

  END SUBROUTINE val_ener_north_south

  SUBROUTINE champ_total_anti_sym(communicator, mesh, list_mode, eps_sa, v, v_out, if_u_h)
    !type_sym = 1 pour un champ pair
    !type_sym = -1 pour un champ impair

    USE Gauss_points
    USE def_type_mesh
    USE tn_axi

    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                                 :: mesh
    INTEGER, DIMENSION(:),                                       INTENT(IN) :: list_mode
    REAL(KIND=8),                                                INTENT(IN) :: eps_sa ! eps_sa=+1 (anti), eps_sa=-1 (sym)
    REAL(KIND=8), DIMENSION(:,:,:)                              ,INTENT(IN) :: v      !champ
    CHARACTER(len=1), INTENT(IN)                           :: if_u_h
    REAL(KIND=8), DIMENSION(mesh%np, SIZE(v,2),SIZE(list_mode))             :: vv
    REAL(KIND=8),DIMENSION(3)                                               :: type_sym
    !type_sym(r,theta,z)
    REAL(KIND=8), DIMENSION(mesh%np,size(v,2),size(list_mode)), INTENT(OUT)  :: v_out      !champ antisymetrise

    INTEGER,      DIMENSION(:,:), POINTER                  :: jj
    INTEGER,                      POINTER                  :: me
    INTEGER      :: i
    INTEGER      :: m_max_c
    !#include "petsc/finclude/petsc.h"
    MPI_Comm       :: communicator


    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me
    m_max_c = size(list_mode)
    v_out = 0.d0

    CALL  symm_champ(communicator, v, mesh, vv, if_u_h)

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
       v_out(:,1:2,i) =  0.5d0*(v(:,1:2,i) - eps_sa*type_sym(1)*vv(:,1:2,i))
       v_out(:,3:4,i) =  0.5d0*(v(:,3:4,i) - eps_sa*type_sym(2)*vv(:,3:4,i))
       v_out(:,5:6,i) =  0.5d0*(v(:,5:6,i) - eps_sa*type_sym(3)*vv(:,5:6,i))
    ENDDO

  END SUBROUTINE champ_total_anti_sym
END MODULE symmetric_field
