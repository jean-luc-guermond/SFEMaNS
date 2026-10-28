MODULE my_util
CONTAINS
  !
  !Authors: Jean-Luc Guermond, Lugi Quartapelle, Copyright 1994
  !

  FUNCTION user_time() RESULT(time)
    IMPLICIT NONE
    REAL(KIND=8) :: time
    INTEGER :: count, count_rate, count_max
    CALL SYSTEM_CLOCK(COUNT, COUNT_RATE, COUNT_MAX)
    time = (1.d0*count)/count_rate
  END FUNCTION user_time

  SUBROUTINE error_Petsc(string)
    USE petscmpi, ONLY: PETSC_COMM_WORLD
    IMPLICIT NONE
    CHARACTER(LEN=*),       INTENT(IN) :: string
    INTEGER                            :: rank, ierr
    CALL MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
    IF (rank==0) WRITE(*,*) string
    CALL PetscFinalize(ierr)
    STOP
  END SUBROUTINE error_Petsc

END MODULE my_util
