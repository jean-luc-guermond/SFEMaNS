PROGRAM scal_cons
#include "petsc/finclude/petsc.h"
  USE petsc
  USE input_data
  USE mesh_handling
  USE two_dim_vtu_xml
  USE sub_plot
  USE my_util
  USE update
  IMPLICIT NONE
  REAL(KIND=8) :: t1
  !---PETSC declaration are at the end of the declarations
  PetscErrorCode :: ierr
  PetscMPIInt    :: rank, nb_proc

  !==Start PETSC
  CALL PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  CALL MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
  CALL MPI_Comm_size(PETSC_COMM_WORLD,nb_proc,ierr)

  !===Read data
  CALL read_my_data('data')

  !===Construct mesh
  IF (inputs%type_fe.NE.1) THEN
     CALL error_Petsc('type_fe>1')
  END IF
  CALL construct_mesh

  !===Run code
  CALL euler
  t1 = user_time()
  CALL euler
  IF (rank==0) WRITE(*,*) ' Time for solution', (user_time()-t1)

  CALL PetscFinalize(ierr)
END PROGRAM scal_cons
