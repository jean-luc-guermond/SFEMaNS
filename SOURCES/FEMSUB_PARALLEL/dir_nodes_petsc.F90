!
!Authors: Jean-Luc Guermond, Franky Luddens, Copyright 2011
!
MODULE Dir_nodes_petsc

CONTAINS
  !-------------------------------------------------------------------------------
  SUBROUTINE dir_axis_nodes_parallel(mesh, js_d)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                 :: mesh
    INTEGER, DIMENSION(:), POINTER  :: js_d
    LOGICAL, DIMENSION(mesh%dom_np) :: virgin
    INTEGER:: nn, ms, n, p, n_D, nws
    REAL(kind=8) :: eps=1.d-10

    nws = SIZE(mesh%jjs,1)
    nn=0
    virgin=.TRUE.
    DO ms = 1, mesh%dom_mes
       IF (MAXVAL(ABS(mesh%rr(1,mesh%jjs(:,ms)))).GT.eps) CYCLE
       DO n = 1, nws
          p = mesh%jjs(n,ms)
          IF (p>mesh%dom_np) CYCLE
          IF (virgin(p)) THEN
             virgin(p)=.FALSE.
             nn = nn + 1
          END IF
       END DO
    END DO
    n_D = nn
    ALLOCATE(js_D(n_D))
    IF (n_D==0) RETURN

    nn=0
    virgin=.TRUE.
    DO ms = 1, mesh%dom_mes
       IF (MAXVAL(ABS(mesh%rr(1,mesh%jjs(:,ms)))).GT.eps) CYCLE
       DO n = 1, nws
          p = mesh%jjs(n,ms)
          IF (p>mesh%dom_np) CYCLE
          IF (virgin(p)) THEN
             virgin(p)=.FALSE.
             nn = nn + 1
             js_D(nn) = mesh%jjs(n,ms)
          END IF
       END DO
    END DO

  END SUBROUTINE dir_axis_nodes_parallel

  SUBROUTINE dirichlet_nodes_parallel(mesh, list_dirichlet_sides, js_d)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                     :: mesh
    INTEGER, DIMENSION(:),   INTENT(IN) :: list_dirichlet_sides
    INTEGER, DIMENSION(:),   POINTER    :: js_d
    LOGICAL, DIMENSION(:),   POINTER    :: virgin
    INTEGER:: nn, ms, n, p, n_D, nws

    IF (SIZE(list_dirichlet_sides)==0) THEN
       ALLOCATE(js_d(0))
       RETURN
    END IF

    nws = SIZE(mesh%jjs,1)
    nn=0
    ALLOCATE(virgin(mesh%dom_np))
    virgin=.TRUE.
    DO ms = 1, mesh%dom_mes
       IF (MINVAL(ABS(mesh%sides(ms)-list_dirichlet_sides))/=0) CYCLE
       DO n = 1, nws
          p = mesh%jjs(n,ms)
          IF (p>mesh%dom_np) CYCLE
          IF (virgin(p)) THEN
             virgin(p)=.FALSE.
             nn = nn + 1
          END IF
       END DO
    END DO
    n_D = nn
    ALLOCATE(js_D(n_D))
    nn=0
    virgin=.TRUE.
    DO ms = 1, mesh%dom_mes
       IF (MINVAL(ABS(mesh%sides(ms)-list_dirichlet_sides))/=0) CYCLE
       DO n = 1, nws
          p = mesh%jjs(n,ms)
          IF (p>mesh%dom_np) CYCLE
          IF (virgin(p)) THEN
             virgin(p)=.FALSE.
             nn = nn + 1
             js_D(nn) = mesh%jjs(n,ms)
          END IF
       END DO
    END DO
    DEALLOCATE(virgin)
  END SUBROUTINE dirichlet_nodes_parallel

  SUBROUTINE dirichlet_nodes_local(mesh, list_dirichlet_sides, js_d)
    USE def_type_mesh
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type)                     :: mesh
    INTEGER, DIMENSION(:),   INTENT(IN) :: list_dirichlet_sides
    INTEGER, DIMENSION(:),   POINTER    :: js_d
    LOGICAL, DIMENSION(:),   POINTER    :: virgin
    INTEGER:: nn, ms, n, p, n_D, nws

    IF (SIZE(list_dirichlet_sides)==0) THEN
       ALLOCATE(js_d(0))
       RETURN
    END IF

    nws = SIZE(mesh%jjs,1)
    nn=0
    ALLOCATE(virgin(mesh%np))
    virgin=.TRUE.
    DO ms = 1, mesh%dom_mes
       IF (MINVAL(ABS(mesh%sides(ms)-list_dirichlet_sides))/=0) CYCLE
       DO n = 1, nws
          p = mesh%jjs(n,ms)
          IF (p>mesh%np) CALL error_petsc('BUG in dirichlet_nodes_local')
          IF (virgin(p)) THEN
             virgin(p)=.FALSE.
             nn = nn + 1
          END IF
       END DO
    END DO
    n_D = nn
    ALLOCATE(js_D(n_D))
    nn=0
    virgin=.TRUE.
    DO ms = 1, mesh%dom_mes
       IF (MINVAL(ABS(mesh%sides(ms)-list_dirichlet_sides))/=0) CYCLE
       DO n = 1, nws
          p = mesh%jjs(n,ms)
          IF (p>mesh%np) CALL error_petsc('BUG in dirichlet_nodes_local')
          IF (virgin(p)) THEN
             virgin(p)=.FALSE.
             nn = nn + 1
             js_D(nn) = mesh%jjs(n,ms)
          END IF
       END DO
    END DO
    DEALLOCATE(virgin)
  END SUBROUTINE dirichlet_nodes_local

  SUBROUTINE Dirichlet_M_parallel(matrix,glob_js_D)
#include "petsc/finclude/petsc.h"
    use petsc
    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(IN) :: glob_js_D
    INTEGER :: n_D
    INTEGER, DIMENSION(:), POINTER    :: bubu_test
    Mat            :: matrix
    PetscErrorCode :: ierr
    n_D = SIZE(glob_js_D)
    ALLOCATE(bubu_test(n_D))
    IF (n_D/=0) THEN
       bubu_test = glob_js_D-1
    END IF
!!$  CALL MatZeroRows(matrix, n_D, glob_js_D-1, 1.d0, ierr)
    !CALL MatZeroRows(matrix, n_D, bubu_test, 1.d0, PETSC_NULL_OBJECT, PETSC_NULL_OBJECT, ierr) !petsc.3.7.
    CALL MatZeroRows(matrix, n_D, bubu_test, 1.d0, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr) !petsc.3.8.4
    CALL MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)

    DEALLOCATE(bubu_test)
  END SUBROUTINE Dirichlet_M_parallel

  SUBROUTINE dirichlet_rhs(js_D,bs_D,b)
    USE def_type_mesh
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    INTEGER,         DIMENSION(:)  :: js_D
    REAL(KIND=8),    DIMENSION(:)  :: bs_D
    INTEGER :: n_D
    Vec            :: b
    PetscErrorCode :: ierr
    n_D = SIZE(js_D)
    IF (n_D/=0) THEN
       CALL VecSetValues(b, n_D, js_D, bs_D, INSERT_VALUES, ierr)
    END IF
    CALL VecAssemblyBegin(b,ierr)
    CALL VecAssemblyEnd(b,ierr)

  END SUBROUTINE dirichlet_rhs

  SUBROUTINE vector_glob_js_D(vv_mesh, list_mode, vv_3_LA, vv_list_dirichlet_sides, vv_js_D, vv_mode_global_js_D)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type),                  INTENT(IN) :: vv_mesh
    INTEGER,            DIMENSION(:), INTENT(IN) :: list_mode
    TYPE(petsc_csr_LA),               INTENT(IN) :: vv_3_LA
    TYPE(dyn_int_line), DIMENSION(3), INTENT(IN) :: vv_list_dirichlet_sides
    TYPE(dyn_int_line), DIMENSION(:), POINTER    :: vv_mode_global_js_D
    TYPE(dyn_int_line), DIMENSION(3), INTENT(OUT):: vv_js_D
    INTEGER,            DIMENSION(:), POINTER    :: vv_js_axis_D
    INTEGER :: k, m_max_c, i, n1, n2, n3, n123, nalloc, nx
    m_max_c = SIZE(list_mode)

    DO k = 1, 3
       CALL dirichlet_nodes_parallel(vv_mesh, vv_list_dirichlet_sides(k)%DIL, vv_js_D(k)%DIL)
    END DO
    CALL dir_axis_nodes_parallel(vv_mesh, vv_js_axis_d)

    ALLOCATE(vv_mode_global_js_D(m_max_c))
    DO i = 1, m_max_c
       n1 = SIZE(vv_js_D(1)%DIL)
       n2 = SIZE(vv_js_D(2)%DIL)
       n3 = SIZE(vv_js_D(3)%DIL)
       nx = SIZE(vv_js_axis_d)
       n123 = n1+n2+n3
       IF (list_mode(i)==0) THEN
          nalloc = n123 + 2*nx
       ELSE IF (list_mode(i)==1) THEN
          nalloc = n123 + nx
       ELSE
          nalloc = n123 + 3*nx
       END IF
       ALLOCATE(vv_mode_global_js_D(i)%DIL(nalloc))
       vv_mode_global_js_D(i)%DIL(1:n1)                  = vv_3_LA%loc_to_glob(1,vv_js_D(1)%DIL)
       vv_mode_global_js_D(i)%DIL(n1+1:n1+n2)            = vv_3_LA%loc_to_glob(2,vv_js_D(2)%DIL)
       vv_mode_global_js_D(i)%DIL(n1+n2+1:n123)          = vv_3_LA%loc_to_glob(3,vv_js_D(3)%DIL)

       IF (list_mode(i)==0 .AND. nx>0) THEN
          vv_mode_global_js_D(i)%DIL(n123+1:n123+nx)     = vv_3_LA%loc_to_glob(1,vv_js_axis_D)
          vv_mode_global_js_D(i)%DIL(n123+nx+1:)         = vv_3_LA%loc_to_glob(2,vv_js_axis_D)
       ELSE IF (list_mode(i)==1 .AND. nx>0) THEN
          vv_mode_global_js_D(i)%DIL(n123+1:)            = vv_3_LA%loc_to_glob(3,vv_js_axis_D)
       ELSE IF (list_mode(i).GE.2 .AND. nx>0) THEN
          vv_mode_global_js_D(i)%DIL(n123+1:n123+nx)     = vv_3_LA%loc_to_glob(1,vv_js_axis_D)
          vv_mode_global_js_D(i)%DIL(n123+nx+1:n123+2*nx)= vv_3_LA%loc_to_glob(2,vv_js_axis_D)
          vv_mode_global_js_D(i)%DIL(n123+2*nx+1:)       = vv_3_LA%loc_to_glob(3,vv_js_axis_D)
       END IF
    END DO

  END SUBROUTINE vector_glob_js_D

  SUBROUTINE vector_without_bc_glob_js_D(vv_mesh, list_mode, vv_3_LA, vv_mode_global_js_D)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type),                  INTENT(IN) :: vv_mesh
    INTEGER,            DIMENSION(:), INTENT(IN) :: list_mode
    TYPE(petsc_csr_LA),               INTENT(IN) :: vv_3_LA
    TYPE(dyn_int_line), DIMENSION(:), POINTER    :: vv_mode_global_js_D
    INTEGER,            DIMENSION(:), POINTER    :: vv_js_axis_D
    INTEGER :: m_max_c, i, nalloc, nx

    m_max_c = SIZE(list_mode)
    CALL dir_axis_nodes_parallel(vv_mesh, vv_js_axis_d)
    ALLOCATE(vv_mode_global_js_D(m_max_c))
    DO i = 1, m_max_c
       nx = SIZE(vv_js_axis_d)
       IF (list_mode(i)==0) THEN
          nalloc = 2*nx
       ELSE IF (list_mode(i)==1) THEN
          nalloc = nx
       ELSE
          nalloc = 3*nx
       END IF
       ALLOCATE(vv_mode_global_js_D(i)%DIL(nalloc))

       IF (list_mode(i)==0 .AND. nx>0) THEN
          vv_mode_global_js_D(i)%DIL(1:nx)     = vv_3_LA%loc_to_glob(1,vv_js_axis_D)
          vv_mode_global_js_D(i)%DIL(nx+1:)    = vv_3_LA%loc_to_glob(2,vv_js_axis_D)
       ELSE IF (list_mode(i)==1 .AND. nx>0) THEN
          vv_mode_global_js_D(i)%DIL           = vv_3_LA%loc_to_glob(3,vv_js_axis_D)
       ELSE IF (list_mode(i).GE.2 .AND. nx>0) THEN
          vv_mode_global_js_D(i)%DIL(1:nx)     = vv_3_LA%loc_to_glob(1,vv_js_axis_D)
          vv_mode_global_js_D(i)%DIL(nx+1:2*nx)= vv_3_LA%loc_to_glob(2,vv_js_axis_D)
          vv_mode_global_js_D(i)%DIL(2*nx+1:)  = vv_3_LA%loc_to_glob(3,vv_js_axis_D)
       END IF
    END DO

  END SUBROUTINE vector_without_bc_glob_js_D


  SUBROUTINE scalar_with_bc_glob_js_D(pp_mesh, list_mode, pp_1_LA, pp_js_D, pp_mode_global_js_D)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type),                  INTENT(IN) :: pp_mesh
    INTEGER,            DIMENSION(:), INTENT(IN) :: list_mode
    TYPE(petsc_csr_LA),               INTENT(IN) :: pp_1_LA
    TYPE(dyn_int_line), DIMENSION(:), POINTER    :: pp_mode_global_js_D
    INTEGER,            DIMENSION(:), INTENT(IN) :: pp_js_D
    INTEGER,            DIMENSION(:), POINTER    :: pp_js_axis_D
    INTEGER :: m_max_c, i, n, nalloc, nx

    m_max_c = SIZE(list_mode)
    CALL dir_axis_nodes_parallel(pp_mesh, pp_js_axis_D)

    ALLOCATE(pp_mode_global_js_D(m_max_c))
    DO i = 1, m_max_c
       n = SIZE(pp_js_D)
       nx = SIZE(pp_js_axis_d)
       IF (list_mode(i)==0) THEN
          nalloc = n
       ELSE
          nalloc = n + nx
       END IF

       ALLOCATE(pp_mode_global_js_D(i)%DIL(nalloc))
       pp_mode_global_js_D(i)%DIL(1:n)         = pp_1_LA%loc_to_glob(1,pp_js_D)

       IF (list_mode(i).GE.1 .AND. nx>0) THEN
          pp_mode_global_js_D(i)%DIL(n+1:n+nx) = pp_1_LA%loc_to_glob(1,pp_js_axis_D)
       END IF
    END DO
  END SUBROUTINE scalar_with_bc_glob_js_D


  SUBROUTINE scalar_without_glob_js_D(pp_mesh, list_mode, pp_1_LA, pp_mode_global_js_D)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type),                  INTENT(IN) :: pp_mesh
    INTEGER,            DIMENSION(:), INTENT(IN) :: list_mode
    TYPE(petsc_csr_LA),               INTENT(IN) :: pp_1_LA
    TYPE(dyn_int_line), DIMENSION(:), POINTER    :: pp_mode_global_js_D
    INTEGER,            DIMENSION(:), POINTER    :: pp_js_axis_D
    INTEGER :: m_max_c, i, nalloc, nx

    m_max_c = SIZE(list_mode)
    CALL dir_axis_nodes_parallel(pp_mesh, pp_js_axis_d)
    ALLOCATE(pp_mode_global_js_D(m_max_c))
    DO i = 1, m_max_c
       nx = SIZE(pp_js_axis_d)
       IF (list_mode(i)==0) THEN
          nalloc = 0
       ELSE
          nalloc = nx
       END IF
       ALLOCATE(pp_mode_global_js_D(i)%DIL(nalloc))
       IF (list_mode(i).GE.1 .AND. nx>0) THEN
          pp_mode_global_js_D(i)%DIL = pp_1_LA%loc_to_glob(1,pp_js_axis_D)
       END IF
    END DO
  END SUBROUTINE scalar_without_glob_js_D
END MODULE Dir_nodes_petsc
