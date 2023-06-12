MODULE solve_petsc
  USE my_util
  TYPE solver_param
     INTEGER:: it_max
     REAL(KIND=8) :: rel_tol, abs_tol
     CHARACTER(LEN=20) :: solver, precond
     LOGICAL           :: verbose
  END TYPE solver_param

CONTAINS
  SUBROUTINE init_solver(my_par,my_ksp,matrix,communicator,solver,precond, opt_re_init)
    USE chaine_caractere
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    LOGICAL, INTENT(IN),    OPTIONAL :: opt_re_init
    TYPE(solver_param)               :: my_par
    CHARACTER(*),           OPTIONAL :: solver, precond
    LOGICAL                          :: re_init
    INTEGER  :: deb, fin

    Mat            :: matrix
    KSP            :: my_ksp
    PC             :: prec
    PetscErrorCode :: ierr
    MPI_Comm       :: communicator

    IF (.NOT.PRESENT(opt_re_init)) THEN
       re_init=.FALSE.
    ELSE
       re_init=opt_re_init
    END IF

    IF (my_par%it_max.LE.0) THEN
       my_par%it_max = 100
    END IF
    IF (my_par%rel_tol.LE.0.d0) THEN
       my_par%rel_tol = 1.d-8
    END IF
    IF (my_par%abs_tol.LE.0.d0) THEN
       my_par%abs_tol = 1.d-14
    END IF

    IF (.NOT.re_init) CALL KSPCreate(communicator,my_ksp,ierr)
    !CALL KSPCreate(communicator,my_ksp,ierr)

    !CALL KSPSetOperators(my_ksp,matrix,matrix,DIFFERENT_NONZERO_PATTERN,ierr)
    !CALL KSPSetOperators(my_ksp,matrix,matrix,SAME_NONZERO_PATTERN ,ierr) !Petsc 3.4.2
    CALL KSPSetOperators(my_ksp,matrix,matrix,ierr) !Petsc 3.7.2

    IF (PRESENT(solver)) THEN
       deb = start_of_string (solver)
       fin = last_of_string (solver)
       IF (solver(deb:fin)=='BCGS') THEN
          CALL KSPSetType(my_ksp, KSPBCGS, ierr)
       ELSE IF (solver(deb:fin)=='GMRES') THEN
          CALL KSPSetType(my_ksp, KSPGMRES, ierr)
       ELSE IF (solver(deb:fin)=='FGMRES') THEN
          CALL KSPSetType(my_ksp, KSPFGMRES, ierr)
       ELSE IF (solver(deb:fin)=='PCR') THEN
          CALL KSPSetType(my_ksp, KSPCR, ierr)
       ELSE IF (solver(deb:fin)=='CHEBYCHEV') THEN
          CALL KSPSetType(my_ksp, KSPCHEBYSHEV, ierr)
       ELSE IF (solver(deb:fin)=='CG') THEN
          CALL KSPSetType(my_ksp, KSPCG, ierr)
       ELSE
          CALL KSPSetType(my_ksp, KSPFGMRES, ierr)
       END IF
    ELSE
       CALL KSPSetType(my_ksp, KSPFGMRES, ierr)
    END IF
    CALL KSPSetTolerances(my_ksp, my_par%rel_tol, my_par%abs_tol, &
         PETSC_DEFAULT_REAL, my_par%it_max, ierr)
    CALL KSPGetPC(my_ksp, prec, ierr)
    IF (PRESENT(precond)) THEN
       deb = start_of_string (precond)
       fin = last_of_string (precond)
       IF (precond(deb:fin)=='JACOBI') THEN
          CALL PCSetType(prec, PCBJACOBI, ierr)
       ELSE IF  (precond(deb:fin)=='HYPRE') THEN
          CALL PCSetType(prec, PCHYPRE, ierr)
       ELSE IF  (precond(deb:fin)=='SSOR') THEN
          CALL PCSetType(prec, PCSOR, ierr)
       ELSE IF  (precond(deb:fin)=='MUMPS') THEN
          CALL PCSetType(prec, PCLU, ierr)
          CALL KSPSetType(my_ksp, KSPPREONLY, ierr)
          !CALL PCFactorSetMatSolverPackage(prec, MATSOLVERMUMPS, ierr) !(JLG) Feb 20, 2019, petsc.3.8.4
          CALL PCFactorSetMatSolverType(prec, MATSOLVERMUMPS, ierr) !
       ELSE
          CALL PCSetType(prec, PCHYPRE, ierr)
       END IF
    ELSE
       CALL PCSetType(prec, PCHYPRE, ierr)
    END IF
    CALL KSPSetFromOptions(my_ksp, ierr)
  END SUBROUTINE init_solver

  SUBROUTINE solver(my_ksp,b,x,reinit,verbose)
#include "petsc/finclude/petsc.h"
    use petsc
    IMPLICIT NONE
    LOGICAL,  OPTIONAL :: reinit, verbose
    INTEGER            :: its
    KSP            :: my_ksp
    PetscErrorCode :: ierr
    Vec            :: x, b
    KSPConvergedReason :: reason
    IF (.NOT.PRESENT(reinit)) reinit=.TRUE.
    IF (.NOT.PRESENT(verbose)) verbose=.FALSE.

    IF (reinit) CALL VecZeroEntries (x,ierr)
    CALL KSPSolve(my_ksp,b,x,ierr)
    IF (verbose) THEN
       CALL KSPGetIterationNumber(my_ksp, its, ierr)
       CALL KSPGetConvergedReason(my_ksp, reason, ierr)
       SELECT CASE(reason)
       CASE(2)
          WRITE(*,*) "KSP_CONVERGED_RTOL, Nb of iterations", its
       CASE(3)
          WRITE(*,*) "KSP_CONVERGED_ATOL, Nb of iterations", its
       CASE(4)
          WRITE(*,*) "Converged after one single iteration of the preconditioner is applied"
       CASE(5,6,7,8)
          WRITE(*,*) "Converge for strange reason:", reason
       CASE(-2)
          WRITE(*,*) "KSP_DIVERGED_NULL"
       CASE(-3)
          WRITE(*,*) "Not converged after it_max", its
       CASE(-4)
          WRITE(*,*) "Not converged: explosion"
       CASE(-5,-6,-7)
          WRITE(*,*) "Not converged for strange reasons", reason
       CASE(-8)
          WRITE(*,*) "Not converged: Indefinite preconditioner"
       CASE(-9)
          WRITE(*,*) "Not converged: NAN"
       CASE(-10)
          WRITE(*,*) "Not converged: Indefinite matrix"
       CASE DEFAULT
          WRITE(*,*) "Something strange happened", reason
       END SELECT
    END IF

  END SUBROUTINE solver

  SUBROUTINE create_local_petsc_matrix(communicator, LA, matrix, clean)
    USE def_type_mesh
#include "petsc/finclude/petsc.h"
    use petsc
    IMPLICIT NONE
    TYPE(petsc_csr_LA)                 :: LA
    LOGICAL, OPTIONAL                  :: clean
    REAL(KIND=8),DIMENSION(:), POINTER :: aa
    INTEGER :: nnzm1, dom_np
    LOGICAL :: test_clean
!!$  INTEGER, DIMENSION(:), POINTER :: ia, ja
    MPI_Comm       :: communicator
    Mat            :: matrix
    PetscErrorCode :: ierr
    !------------------------------------------------------------------------------
    dom_np = SIZE(LA%ia)-1
    nnzm1=LA%ia(dom_np)-LA%ia(0)-1
    ALLOCATE(aa(0:nnzm1))
    aa =0.d0


!!$  ALLOCATE(ia(0:dom_np),ja(0:nnzm1))
!!$  ia = LA%ia
!!$  ja = LA%ja
!!$  CALL MatCreateMPIAIJWithArrays(communicator,dom_np,dom_np,PETSC_DECIDE, &
!!$       PETSC_DECIDE, ia, ja, aa, matrix, ierr)
!!$DEALLOCATE(ia,ja)

    CALL MatCreateMPIAIJWithArrays(communicator,dom_np,dom_np,PETSC_DECIDE, &
         PETSC_DECIDE, LA%ia, LA%ja, aa, matrix, ierr)

    DEALLOCATE(aa)
    IF (PRESENT(clean)) THEN
       test_clean=clean
    ELSE
       test_clean=.TRUE.
    END IF
    IF (test_clean) THEN
       IF (ASSOCIATED(LA%ia)) DEALLOCATE(LA%ia)
       IF (ASSOCIATED(LA%ja)) DEALLOCATE(LA%ja)
    END IF
  END SUBROUTINE create_local_petsc_matrix

  SUBROUTINE create_local_petsc_matrix_a_detruire(communicator, aij, i_loc, matrix)
    USE def_type_mesh
#include "petsc/finclude/petsc.h"
    use petsc
    IMPLICIT NONE
    TYPE(aij_type),         INTENT(IN) :: aij
    INTEGER,     DIMENSION(2)          :: i_loc
    INTEGER,     DIMENSION(:), POINTER :: ia, ja
    REAL(KIND=8),DIMENSION(:), POINTER :: aa
    INTEGER :: nnzm1, dom_np, p, i, n

    MPI_Comm       :: communicator
    Mat            :: matrix
    PetscErrorCode :: ierr
    !------------------------------------------------------------------------------
    dom_np = i_loc(2) - i_loc(1) + 1
    nnzm1 = aij%ia(i_loc(2)+1)-aij%ia(i_loc(1))-1
    ALLOCATE(ia(0:dom_np),ja(0:nnzm1))
    ia(0)=0
    DO i = 1, dom_np
       n = i_loc(1) + i -1
       ia(i) = aij%ia(n+1)-aij%ia(i_loc(1))
       DO p=aij%ia(n), aij%ia(n+1)-1
          ja(p-aij%ia(i_loc(1)))= aij%ja(p)-1
       END DO
    END DO
    !------------------------------------------------------------------------------
    ALLOCATE(aa(0:nnzm1))
    aa =0
    CALL MatCreateMPIAIJWithArrays(communicator,dom_np,dom_np,PETSC_DECIDE, &
         PETSC_DECIDE, ia, ja, aa, matrix, ierr)

    DEALLOCATE(ia,ja,aa)
  END SUBROUTINE create_local_petsc_matrix_a_detruire

  SUBROUTINE create_local_petsc_block_matrix(communicator, n_b, aij, i_loc, matrix)
    USE def_type_mesh
#include "petsc/finclude/petsc.h"
    use petsc
    IMPLICIT NONE
    TYPE(aij_type),           INTENT(IN) :: aij
    INTEGER,     DIMENSION(2)            :: i_loc
    INTEGER                              :: n_b
    INTEGER,     DIMENSION(:),   POINTER :: ia, ja
    REAL(KIND=8),DIMENSION(:),   POINTER :: aa
    INTEGER :: nnzm1, dom_np, p, i, n, ib, k

    MPI_Comm       :: communicator
    Mat            :: matrix
    PetscErrorCode :: ierr

    dom_np = i_loc(2) - i_loc(1) + 1
    nnzm1 = n_b*(aij%ia(i_loc(2)+1)-aij%ia(i_loc(1))-1)
    ALLOCATE(ia(0:n_b*dom_np),ja(0:nnzm1))
    ia(0)=0
    DO k = 1, n_b
       DO i = 1, dom_np
          ib = i + (k-1)*dom_np
          n = i_loc(1) + i - 1
          ia(i) = n_b*(aij%ia(n+1)-aij%ia(i_loc(1)))
          DO p=aij%ia(n), aij%ia(n+1)-1
             ja(p-aij%ia(i_loc(1)))= aij%ja(p)-1
          END DO
       END DO
    END DO
    !------------------------------------------------------------------------------
    ALLOCATE(aa(0:nnzm1))
    aa =0
    CALL MatCreateMPIAIJWithArrays(communicator,dom_np,dom_np,PETSC_DECIDE, &
         PETSC_DECIDE, ia, ja, aa, matrix, ierr)

    DEALLOCATE(ia,ja,aa)
  END SUBROUTINE create_local_petsc_block_matrix


END MODULE solve_petsc
