module type_verbose
  TYPE my_verbose
     REAL(KIND=8) :: CFL
     REAL(KIND=8) :: time
     REAL(KIND=8) :: div_L2
     REAL(KIND=8) :: weak_div_L2
     REAL(KIND=8) :: div_B_L2
     REAL(KIND=8) :: total_cpu_time
     REAL(KIND=8) :: total_cpu_time_minus_init
     !CONTAINS
     !PROCEDURE, PUBLIC :: write_verbose
  END type my_verbose
  !CONTAINS
  !  SUBROUTINE write_verbose(a)
  !    USE input_data
  !    CLASS(my_verbose), INTENT(INOUT) :: a
  !    IF (inputs%verbose_CFL) THEN
  !       WRITE(*,'(2(A,e10.3))') ' Time = ', time, ', CFL = ', a%CFL
  !    END IF
  !  END SUBROUTINE write_verbose
END module type_verbose

MODULE verbose
  USE type_verbose
  IMPLICIT NONE
  PUBLIC :: write_verbose
  TYPE(my_verbose), PUBLIC  :: talk_to_me
  PRIVATE

CONTAINS
  SUBROUTINE write_verbose(rank,opt_tps,opt_tploc_max)
    USE input_data
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    REAL(KIND=8), OPTIONAL, INTENT(IN) :: opt_tps, opt_tploc_max
    PetscErrorCode :: code
    PetscMPIInt    :: rank
    IF (inputs%verbose_timing) THEN
       IF (present(opt_tps).AND.present(opt_tploc_max)) THEN
          CALL MPI_ALLREDUCE(opt_tps,talk_to_me%total_cpu_time,1,MPI_DOUBLE_PRECISION,&
               MPI_MAX, PETSC_COMM_WORLD, code)
          IF(inputs%nb_iteration>1) THEN
             CALL MPI_ALLREDUCE(opt_tploc_max,talk_to_me%total_cpu_time_minus_init,1,&
                  MPI_DOUBLE_PRECISION, MPI_MAX, PETSC_COMM_WORLD, code)
          END IF
          IF (rank==0) WRITE(*,'(A,F12.5)') ' Total elapse time ', talk_to_me%total_cpu_time
          IF(inputs%nb_iteration>1) THEN
             IF (rank==0) WRITE(*,'(A,F12.5)') 'Average time in loop (minus initialization) ', &
                  talk_to_me%total_cpu_time_minus_init/(inputs%nb_iteration-1)
          END IF
          RETURN
       END IF
    END IF

    IF (inputs%verbose_CFL) THEN
       IF (rank==0) WRITE(*,'(2(A,e10.3))') ' Time = ', talk_to_me%time, ', CFL = ', talk_to_me%CFL
    END IF
    IF (inputs%verbose_divergence) THEN
       IF (inputs%type_pb=='nst' .OR. inputs%type_pb=='mhd' .OR. inputs%type_pb=='fhd') THEN
          IF (rank==0) WRITE(*,'(2(A,e10.3))') ' Time = ', talk_to_me%time, &
               ', ||div(un)||_L2/||un||_H1 = ', talk_to_me%div_L2
          IF (rank==0) WRITE(*,'(2(A,e10.3))') ' Time = ', talk_to_me%time, &
               ', ||weak_div(un)||_L2/||un||_H1 = ', talk_to_me%weak_div_L2
       END IF
       IF (inputs%type_pb/='nst') THEN
          IF (rank==0) WRITE(*,'(2(A,e10.3))') ' Time = ', talk_to_me%time, &
               ', ||div(Bn)||_L2/||Bn||_L2 = ', talk_to_me%div_B_L2
       END IF
    END IF

  END SUBROUTINE write_verbose
END MODULE verbose
