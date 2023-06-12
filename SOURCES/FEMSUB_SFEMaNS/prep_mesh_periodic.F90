!
!Authors Jean-Luc Guermond, Copyrights 1996
!
MODULE  periodic


  IMPLICIT NONE

  PUBLIC :: prep_periodic, &
       prep_periodic_bloc, &
       prep_periodic_H_p_phi_bc, &
       periodic_matrix_petsc, &
       periodic_rhs_petsc
  PRIVATE

CONTAINS

  !jan 29 2007
  SUBROUTINE prep_periodic(my_periodic, mesh, periodic)
    !=========================================
    USE def_type_mesh
    USE chaine_caractere
    IMPLICIT NONE
    TYPE(periodic_data), INTENT(IN)     :: my_periodic
    TYPE(mesh_type)                     :: mesh
    TYPE(periodic_type)                 :: periodic
    INTEGER,      DIMENSION(:), POINTER :: list_loc, perlist_loc, list_dom, perlist_dom
    INTEGER                             :: n, side1, side2, n_b, nx, i
    REAL(KIND=8), DIMENSION(:), POINTER :: e

    WRITE (*,*) 'Loading periodic-data file ...'

    IF (mesh%np == 0) THEN
       WRITE(*,*) 'no mesh on this proc'
       RETURN
    END IF

    ALLOCATE(e(SIZE(my_periodic%vect_e,1)))

    periodic%n_bord = my_periodic%nb_periodic_pairs
    IF (periodic%n_bord .GT. 20) THEN
       WRITE(*,*) 'PREP_MESH_PERIODIC: trop de bords periodiques'
       STOP
    END IF

    DO n= 1, periodic%n_bord

       side1 = my_periodic%list_periodic(1,n)
       side2 = my_periodic%list_periodic(2,n)
       e = my_periodic%vect_e(:,n)

       CALL list_periodic(mesh%np, mesh%jjs, mesh%sides, mesh%rr, side1, side2, e, &
            list_loc, perlist_loc)

       n_b = SIZE(list_loc)
       ALLOCATE(list_dom(n_b), perlist_dom(n_b))
       nx = 0
       DO i = 1, n_b
          IF (MAX(list_loc(i),perlist_loc(i)) .LE. mesh%dom_np) THEN
             nx = nx+1
             list_dom(nx) = list_loc(i)
             perlist_dom(nx) = perlist_loc(i)
          ELSE IF (MIN(list_loc(i),perlist_loc(i)) .LE. mesh%dom_np) THEN
             WRITE(*,*) 'BUG in prep_periodic'
             STOP
          END IF
       END DO
       IF (n_b /= nx) WRITE(*,*) 'WARNING, I have removed', n_b-nx, ' periodic pairs in prep_periodic'
       n_b = nx

       ALLOCATE (periodic%list(n)%DIL(n_b), periodic%perlist(n)%DIL(n_b))
       periodic%list(n)%DIL = list_dom(1:n_b)
       periodic%perlist(n)%DIL = perlist_dom(1:n_b)

       DEALLOCATE(list_loc,perlist_loc, list_dom, perlist_dom)
    END DO

    DEALLOCATE(e)

    WRITE (*,*) 'Treatment of periodic-data done'

  END SUBROUTINE prep_periodic

  SUBROUTINE prep_periodic_bloc(my_periodic, mesh, periodic, nb_bloc)
    !=========================================
    USE chaine_caractere
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(periodic_data), INTENT(IN)     :: my_periodic
    TYPE(mesh_type)                     :: mesh
    TYPE(periodic_type)                 :: periodic
    INTEGER,                 INTENT(IN) :: nb_bloc
    INTEGER,      DIMENSION(:), POINTER :: list_loc, perlist_loc, list_dom, perlist_dom
    INTEGER                             :: n, side1, side2, nsize, n_b
    INTEGER                             :: k, k_deb, k_fin, nx, i
    REAL(KIND=8), DIMENSION(2)          :: e

    WRITE (*,*) 'Loading periodic-data file ...'

    IF (mesh%np == 0) THEN
       WRITE(*,*) 'no mesh on this proc'
       RETURN
    END IF

    periodic%n_bord = my_periodic%nb_periodic_pairs
    IF (periodic%n_bord .GT. 20) THEN
       WRITE(*,*) 'PREP_MESH_PERIODIC: trop de bords periodiques'
       STOP
    END IF

    DO n= 1, periodic%n_bord

       side1 = my_periodic%list_periodic(1,n)
       side2 = my_periodic%list_periodic(2,n)
       e = my_periodic%vect_e(:,n)

       CALL list_periodic(mesh%np, mesh%jjs, mesh%sides, mesh%rr, side1, side2, e, &
            list_loc, perlist_loc)

       !n_b = SIZE(list_loc)
       n_b = SIZE(perlist_loc)
       ALLOCATE(list_dom(n_b), perlist_dom(n_b))
       nx = 0
       DO i = 1, n_b
          IF (MAX(list_loc(i),perlist_loc(i)) .LE. mesh%dom_np) THEN
             nx = nx+1
             list_dom(nx) = list_loc(i)
             perlist_dom(nx) = perlist_loc(i)
          ELSE IF (MIN(list_loc(i),perlist_loc(i)) .LE. mesh%dom_np) THEN
             WRITE(*,*) 'BUG in prep_periodic_bloc'
             STOP
          END IF
       END DO
       IF (n_b /= nx) WRITE(*,*) 'WARNING, I have removed', n_b-nx, ' periodic pairs in prep_periodic_bloc'
       n_b = nx

       nsize = nb_bloc*n_b !SIZE(list_loc) !n_b

       ALLOCATE(periodic%list(n)%DIL(nsize), periodic%perlist(n)%DIL(nsize))

       DO k = 1, nb_bloc
          k_deb=(k-1)*n_b+1
          k_fin=k*n_b
          periodic%list(n)%DIL(k_deb:k_fin)    = list_dom(1:n_b)    + (k-1)*mesh%dom_np ! First bloc
          periodic%perlist(n)%DIL(k_deb:k_fin) = perlist_dom(1:n_b) + (k-1)*mesh%dom_np ! First bloc
       END DO

       DEALLOCATE(list_loc,perlist_loc, list_dom, perlist_dom)

    END DO

    WRITE (*,*) 'Treatment of periodic-data done'

  END SUBROUTINE prep_periodic_bloc
  !jan 29 2007

  !JLG+FL/Feb 2 2010
  SUBROUTINE prep_periodic_H_p_phi_bc(my_periodic, H_mesh, pmag_mesh, phi_mesh, H_p_phi_per)
    USE chaine_caractere
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(periodic_data), INTENT(IN)     :: my_periodic
    TYPE(mesh_type)                     :: H_mesh, pmag_mesh, phi_mesh
    TYPE(periodic_type)                 :: H_p_phi_per
    INTEGER,      DIMENSION(:), POINTER :: b_list_loc, b_perlist_loc, &
         e_list_loc, e_perlist_loc, p_list_loc, p_perlist_loc
    INTEGER,      DIMENSION(:), POINTER :: b_list_dom, b_perlist_dom, &
         e_list_dom, e_perlist_dom, p_list_dom, p_perlist_dom
    INTEGER                             :: n, side1, side2, nsize, n_b, n_e, n_p, nx, i
    REAL(KIND=8), DIMENSION(2)          :: e

    WRITE (*,*) 'Loading periodic-data file ...'

    H_p_phi_per%n_bord = my_periodic%nb_periodic_pairs
    IF (H_p_phi_per%n_bord .GT. 20) THEN
       WRITE(*,*) 'PREP_MESH_PERIODIC: trop de bords periodiques'
       STOP
    END IF

    DO n= 1, H_p_phi_per%n_bord

       side1 = my_periodic%list_periodic(1,n)
       side2 = my_periodic%list_periodic(2,n)
       e = my_periodic%vect_e(:,n)


       IF (H_mesh%np == 0) THEN
          n_b = 0
       ELSE
          CALL list_periodic(H_mesh%np, H_mesh%jjs, H_mesh%sides, H_mesh%rr, side1, side2, e, &
               b_list_loc, b_perlist_loc)
          n_b = SIZE(b_list_loc)
          ALLOCATE(b_list_dom(n_b), b_perlist_dom(n_b))
          nx = 0
          DO i = 1, n_b
             IF (MAX(b_list_loc(i),b_perlist_loc(i)) .LE. H_mesh%dom_np) THEN
                nx = nx+1
                b_list_dom(nx) = b_list_loc(i)
                b_perlist_dom(nx) = b_perlist_loc(i)
             ELSE IF (MIN(b_list_loc(i),b_perlist_loc(i)) .LE. H_mesh%dom_np) THEN
                WRITE(*,*) 'BUG in prep_periodic_H_p_phi_bc (H)'
                STOP
             END IF
          END DO
          IF (n_b /= nx) WRITE(*,*) 'WARNING, I have removed', n_b-nx, ' periodic pairs on H'
          n_b = nx
       END IF

       IF (pmag_mesh%np == 0) THEN
          n_p = 0
       ELSE
          CALL list_periodic(pmag_mesh%np, pmag_mesh%jjs, pmag_mesh%sides, pmag_mesh%rr, side1, side2, e, &
               p_list_loc, p_perlist_loc)
          n_p = SIZE(p_list_loc)
          ALLOCATE(p_list_dom(n_p), p_perlist_dom(n_p))
          nx = 0
          DO i = 1, n_p
             IF (MAX(p_list_loc(i),p_perlist_loc(i)) .LE. pmag_mesh%dom_np) THEN
                nx = nx+1
                p_list_dom(nx) = p_list_loc(i)
                p_perlist_dom(nx) = p_perlist_loc(i)
             ELSE IF (MIN(p_list_loc(i),p_perlist_loc(i)) .LE. pmag_mesh%dom_np) THEN
                WRITE(*,*) 'BUG in prep_periodic_H_p_phi_bc (pmag) '
                STOP
             END IF
          END DO
          IF (n_p /= nx) WRITE(*,*) 'WARNING, I have removed', n_p-nx, ' periodic pairs on pmag'
          n_p = nx
       END IF

       IF (phi_mesh%np==0) THEN
          n_e = 0
       ELSE
          CALL list_periodic(phi_mesh%np, phi_mesh%jjs, phi_mesh%sides, phi_mesh%rr, side1, side2, e, &
               e_list_loc, e_perlist_loc)
          n_e = SIZE(e_list_loc)
          ALLOCATE(e_list_dom(n_e), e_perlist_dom(n_e))
          nx = 0
          DO i = 1, n_e
             IF (MAX(e_list_loc(i),e_perlist_loc(i)) .LE. phi_mesh%dom_np) THEN
                nx = nx+1
                e_list_dom(nx) = e_list_loc(i)
                e_perlist_dom(nx) = e_perlist_loc(i)
             ELSE IF (MIN(e_list_loc(i),e_perlist_loc(i)) .LE. phi_mesh%dom_np) THEN
                WRITE(*,*) 'BUG in prep_periodic_H_p_phi_bc (phi) '
                STOP
             END IF
          END DO
          IF (n_e /= nx) WRITE(*,*) 'WARNING, I have removed', n_e-nx, ' periodic pairs on phi'
          n_e = nx
       END IF

       !n_b = SIZE(b_list_loc)
       !n_p = SIZE(p_list_loc)
       !n_e = SIZE(e_list_loc)
       nsize = 3*n_b + n_p + n_e

       ALLOCATE(H_p_phi_per%list(n)%DIL(nsize), H_p_phi_per%perlist(n)%DIL(nsize))

       IF (n_b /=0) THEN

          H_p_phi_per%list(n)%DIL(1:n_b)    = b_list_dom(1:n_b)    ! First block
          H_p_phi_per%perlist(n)%DIL(1:n_b) = b_perlist_dom(1:n_b) ! First block

          H_p_phi_per%list(n)%DIL(n_b+1:2*n_b)    = b_list_dom(1:n_b)    + H_mesh%dom_np ! Second block
          H_p_phi_per%perlist(n)%DIL(n_b+1:2*n_b) = b_perlist_dom(1:n_b) + H_mesh%dom_np ! Second block

          H_p_phi_per%list(n)%DIL(2*n_b+1:3*n_b)    = b_list_dom(1:n_b)    + 2*H_mesh%dom_np ! Third block
          H_p_phi_per%perlist(n)%DIL(2*n_b+1:3*n_b) = b_perlist_dom(1:n_b) + 2*H_mesh%dom_np ! Third block

       END IF

       IF (n_p /=0) THEN
          H_p_phi_per%list(n)%DIL(3*n_b+1:3*n_b+n_p)    = p_list_dom(1:n_p)    + 3*H_mesh%dom_np ! Forth block
          H_p_phi_per%perlist(n)%DIL(3*n_b+1:3*n_b+n_p) = p_perlist_dom(1:n_p) + 3*H_mesh%dom_np ! Forth block
       END IF

       IF (n_e/=0) THEN
          H_p_phi_per%list(n)%DIL(3*n_b+n_p+1:)    = e_list_dom(1:n_e)    + 3*H_mesh%dom_np + pmag_mesh%dom_np ! Fourth block
          H_p_phi_per%perlist(n)%DIL(3*n_b+n_p+1:) = e_perlist_dom(1:n_e) + 3*H_mesh%dom_np + pmag_mesh%dom_np ! Fourth block
       END IF

       IF (ASSOCIATED(b_list_loc)) NULLIFY(b_list_loc, b_perlist_loc, b_list_dom, b_perlist_dom)
       IF (ASSOCIATED(p_list_loc)) NULLIFY(p_list_loc, p_perlist_loc, p_list_dom, p_perlist_dom)
       IF (ASSOCIATED(e_list_loc)) NULLIFY(e_list_loc, e_perlist_loc, e_list_dom, e_perlist_dom)

    END DO

    WRITE (*,*) 'Treatment of periodic-data done'

  END SUBROUTINE Prep_periodic_H_p_phi_bc

  SUBROUTINE list_periodic(np, jjs, sides, rr, side1, side2, e, list_out, perlist_out)
    !============================================================================
    IMPLICIT NONE
    INTEGER,                      INTENT(IN)  :: np
    INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jjs
    INTEGER,      DIMENSION(:),   INTENT(IN)  :: sides
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: rr
    INTEGER,                      INTENT(IN)  :: side1, side2
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: e
    INTEGER,      DIMENSION(:),   POINTER     :: list_out, perlist_out
    INTEGER,      DIMENSION(:),   ALLOCATABLE :: list, perlist
    LOGICAL,      DIMENSION(np)               :: virgin
    REAL(KIND=8), DIMENSION(SIZE(rr,1))       :: ri
    INTEGER :: ms, ns, i, j, long, inter
    REAL(KIND=8) :: r, epsilon = 1.d-9
    LOGICAL :: verif

    IF (ALLOCATED(list))    DEALLOCATE(list)
    IF (ALLOCATED(perlist)) DEALLOCATE(perlist)

    ALLOCATE (list(np), perlist(np))
    virgin = .TRUE.

    i = 0; j=0
    DO ms = 1, SIZE(sides)

       IF (sides(ms) .EQ. side1) THEN
          DO ns = 1, SIZE(jjs,1)
             IF (virgin(jjs(ns,ms))) THEN
                i = i + 1
                list(i) = jjs(ns,ms)
                virgin(jjs(ns,ms)) = .FALSE.
             END IF
          END DO
       ELSE IF (sides(ms) .EQ. side2) THEN
          DO ns = 1, SIZE(jjs,1)
             IF (virgin(jjs(ns,ms))) THEN
                j = j + 1
                perlist(j) = jjs(ns,ms)
                virgin(jjs(ns,ms)) = .FALSE.
             END IF
          END DO

       END IF

    END DO

    IF (i .NE. j) THEN
       WRITE(*,*) ' FEM_PERIODIC: side1 and side2 have', &
            ' different numbers of points'
       STOP
    END IF
    long = i

    DO i = 1, long
       ri = rr(:,list(i))+e(:)
       verif = .FALSE.
       !if (i==2) stop
       DO j = i, long
          r = SUM(ABS(ri - rr(:,perlist(j))))
          !if (i==1) write(*,*) ' r',r,'j',  j
          IF (r .LE. epsilon ) THEN
             inter = perlist(i)
             perlist(i) = perlist(j)
             perlist(j) = inter
             verif = .TRUE.
             EXIT
          END IF
       END DO
       IF (.NOT.verif) THEN
          WRITE(*,*) ' BUG dans  data_periodic ou le maillage:', &
               ' side1 + e /= side2'
          WRITE(*,*) ' i = ', i
          !         STOP
       END IF
    END DO

    ALLOCATE (list_out(long))
    list_out(1:long) = list(1:long)
    ALLOCATE (perlist_out(long))
    perlist_out(1:long) = perlist(1:long)

  END SUBROUTINE list_periodic

  SUBROUTINE periodic_matrix_petsc(n_bord, list, perlist, matrix, LA)
    USE dyn_line
    USE def_type_mesh
    USE my_util
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    INTEGER                         , INTENT(IN) :: n_bord
    TYPE(dyn_int_line), DIMENSION(:), INTENT(IN) :: list, perlist
    TYPE(petsc_csr_la)              , INTENT(IN) :: LA
    INTEGER, PARAMETER                           :: nmaxcols = 300
    INTEGER                                      :: ncols
    INTEGER, DIMENSION(nmaxcols)                 :: cols
    REAL(KIND=8), DIMENSION(nmaxcols)            :: vals
    INTEGER, DIMENSION(:), ALLOCATABLE           :: n_cols_i
    INTEGER, DIMENSION(1)                        :: idxn
    INTEGER, DIMENSION(:,:), ALLOCATABLE         :: jdxn
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE    :: vals_pi
    INTEGER                                      :: n, l, i, pi, n_D, k
    Mat                                          :: matrix
    PetscErrorCode                               :: ierr

    WRITE(*,*) 'Entering periodic_matrix_petsc'

    CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
    CALL MatSetOption (matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr)

    DO k = 1, SIZE(LA%loc_to_glob,1)
       DO n = 1, n_bord
          n_D = SIZE(list(n)%DIL)
          IF (n_D /=0) THEN
             ALLOCATE(jdxn(n_D,nmaxcols), vals_pi(n_D,nmaxcols), n_cols_i(n_D))
             jdxn = 0
             vals_pi = 0.d0
             n_cols_i = 0

             DO l = 1, SIZE(list(n)%DIL)
                idxn(1) = LA%loc_to_glob(k,list(n)%DIL(l))
                CALL MatGetRow(matrix, idxn(1)-1,ncols,cols,vals,ierr)
                n_cols_i(l) = ncols
                jdxn(l,1:ncols) = cols(1:ncols)
                vals_pi(l,1:ncols) = vals(1:ncols)
                CALL MatRestoreRow(matrix, idxn(1)-1, ncols,cols, vals, ierr)
             END DO

             DO l= 1, n_D
                idxn(1) = LA%loc_to_glob(k,perlist(n)%DIL(l))-1
                CALL MatSetValues(matrix, 1,idxn, n_cols_i(l), jdxn(l,1:n_cols_i(l)), &
                     vals_pi(l:l,1:n_cols_i(l)), ADD_VALUES, ierr)
             END DO
             DEALLOCATE(jdxn, vals_pi, n_cols_i)

          END IF
       END DO
       CALL MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
       CALL MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)

       DO n = 1, n_bord
          n_D = SIZE(list(n)%DIL)
          !CALL MatZeroRows(matrix, n_D, LA%loc_to_glob(k,list(n)%DIL(:))-1, 1.d0, &
          !     PETSC_NULL_OBJECT, PETSC_NULL_OBJECT, ierr) !petsc.3.7.4
          CALL MatZeroRows(matrix, n_D, LA%loc_to_glob(k,list(n)%DIL(:))-1, 1.d0, &
               PETSC_NULL_VEC, PETSC_NULL_VEC, ierr) !(JLG) Feb 20, 2019, petsc.3.8.4
       END DO
!!$       CALL MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
!!$       CALL MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)

       DO n = 1, n_bord
          DO l = 1, SIZE(list(n)%DIL)
             i = LA%loc_to_glob(k,list(n)%DIL(l))
             pi = LA%loc_to_glob(k,perlist(n)%DIL(l))
             CALL MatSetValue(matrix, i-1, pi-1, -1.d0, INSERT_VALUES, ierr)
          END DO
       END DO
       CALL MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
       CALL MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)

    END DO

  END SUBROUTINE periodic_matrix_petsc

  SUBROUTINE periodic_rhs_petsc(n_bord, list, perlist, v_rhs, LA)
    USE dyn_line
    USE def_type_mesh
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    INTEGER                         , INTENT(IN) :: n_bord
    TYPE(dyn_int_line), DIMENSION(:), INTENT(IN) :: list, perlist
    TYPE(petsc_csr_la)              , INTENT(IN) :: LA
    INTEGER, DIMENSION(:), ALLOCATABLE           :: idxn, jdxn
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE      :: vals, bs
    INTEGER                                      :: n, k, n_D
    Vec                                          :: v_rhs
    PetscErrorCode                               :: ierr


    DO k = 1, SIZE(LA%loc_to_glob,1)
       DO n = 1, n_bord
          n_D = SIZE(list(n)%DIL)
          ALLOCATE(idxn(n_D), vals(n_D), jdxn(n_D), bs(n_D))
          idxn = LA%loc_to_glob(k,list(n)%DIL(:))-1
          jdxn = LA%loc_to_glob(k,perlist(n)%DIL(:))-1
          CALL VecGetValues(v_rhs, n_D, idxn, vals, ierr)
          CALL VecAssemblyBegin(v_rhs,ierr)
          CALL VecAssemblyEnd(v_rhs,ierr)

          bs = 0.d0
          CALL VecSetValues(v_rhs, n_D, jdxn, vals, ADD_VALUES, ierr)
          CALL VecAssemblyBegin(v_rhs,ierr)
          CALL VecAssemblyEnd(v_rhs,ierr)
          CALL VecSetValues(v_rhs, n_D, idxn, bs, INSERT_VALUES, ierr)
          CALL VecAssemblyBegin(v_rhs,ierr)
          CALL VecAssemblyEnd(v_rhs,ierr)

          IF (ALLOCATED(idxn)) DEALLOCATE(idxn, jdxn, vals, bs)
       END DO
    END DO

  END SUBROUTINE periodic_rhs_petsc


END MODULE periodic
