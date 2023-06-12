MODULE update
#include "petsc/finclude/petsc.h"
  USE petsc
  USE solve_petsc
  USE mesh_handling
  PUBLIC                   :: euler
  PRIVATE
  Mat                      :: mat
  Vec                      :: xx, xghost, bb
  KSP                      :: my_ksp
  PetscMPIInt              :: rank
  PetscErrorCode           :: ierr
  TYPE(solver_param)       :: my_par

CONTAINS

  SUBROUTINE contruct_matrices
    USE st_matrix
    USE fem_M
    USE st_matrix
    USE my_util
    USE Dir_nodes_petsc
    IMPLICIT NONE
    INTEGER, POINTER, DIMENSION(:)  :: ifrom  !===for ghost structure

    !===My rank
    CALL MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)

    !===Create ghost structure
    CALL create_my_ghost(mesh,LA,ifrom)
    CALL VecCreateGhost(PETSC_COMM_WORLD, mesh%dom_np, &
         PETSC_DETERMINE, SIZE(ifrom), ifrom, xx, ierr)
    CALL VecGhostGetLocalForm(xx, xghost, ierr)

    !===Duplicate
    CALL VecDuplicate(xx, bb, ierr)

    !===Create matrices
    CALL create_local_petsc_matrix(PETSC_COMM_WORLD, LA, mat, clean=.FALSE.)

    !===mat
    CALL qs_mass_diff_M (mesh, 1.d0, 1.d0, LA, mat)

    !===Dirichlet Bcs for matrix
    CALL Dirichlet_M_parallel(mat,LA%loc_to_glob(1,js_D_loc))

    !===Initialize solver
    my_par%it_max=100
    my_par%rel_tol=1.d-7
    my_par%abs_tol=1.d-10
    my_par%verbose=.FALSE.
    my_par%solver='GMRES'
    my_par%precond= 'MUMPS'
    CALL init_solver(my_par,my_ksp,mat,PETSC_COMM_WORLD,&
         solver=my_par%solver,precond=my_par%precond)

  END SUBROUTINE contruct_matrices

  SUBROUTINE euler
    USE condlim
    USE fem_rhs
    USE Dir_nodes_petsc
    USE fem_tn
    IMPLICIT NONE
    LOGICAL, SAVE :: once=.true.
    !===
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)  :: phi
    REAL(KIND=8)  :: err1, err2

    IF (once) THEN
       CALL contruct_matrices
       once=.FALSE.
    END IF

    !===Rhs
    CALL qs_00 (mesh, LA, source(mesh%rr), bb)

    !===Bcs for Rhs
    CALL dirichlet_rhs(LA%loc_to_glob(1,js_D_loc)-1,uexact(mesh%rr(:,js_D_loc)),bb)

    !===Solve linear system
    CALL solver(my_ksp,bb,xx,reinit=.FALSE.,verbose=my_par%verbose)

    !===Check errors
    CALL VecGhostUpdateBegin(xx,INSERT_VALUES,SCATTER_FORWARD,ierr)
    CALL VecGhostUpdateEnd(xx,INSERT_VALUES,SCATTER_FORWARD,ierr)
    ALLOCATE(phi(mesh%np))
    phi = uexact(mesh%rr)
    CALL norme_Lp(err1,2,PETSC_COMM_WORLD,mesh,xghost)
    CALL norme_Lp(err2,2,PETSC_COMM_WORLD,mesh,xghost,phi)
    IF (rank==0) WRITE(*,*) ' NORME L^2 ',  err2/err1
    CALL norme_Lp(err1,0,PETSC_COMM_WORLD,mesh,xghost)
    CALL norme_Lp(err2,0,PETSC_COMM_WORLD,mesh,xghost,phi)
    IF (rank==0) WRITE(*,*) ' NORME L^infty ',  err2/err1

  END SUBROUTINE euler

END MODULE update
