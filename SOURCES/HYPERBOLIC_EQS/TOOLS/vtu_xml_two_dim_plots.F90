MODULE two_dim_vtu_xml
#include "petsc/finclude/petsc.h"
  USE petsc
  USE vtk_viz
  PUBLIC :: make_vtu_file_scalar_2D
CONTAINS

  SUBROUTINE make_vtu_file_scalar_2D(comm, mesh, header, field, field_name, what, opt_it)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type),               INTENT(IN) :: mesh
    CHARACTER(*),                  INTENT(IN) :: header
    CHARACTER(*),                  INTENT(IN) :: field_name, what
    INTEGER, OPTIONAL,             INTENT(IN) :: opt_it
    REAL(KIND=8), DIMENSION(:),    INTENT(IN) :: field
    INTEGER                                   :: j, it
    CHARACTER(LEN=200), DIMENSION(:), POINTER :: file_list
    CHARACTER(LEN=3)                          :: st_rank, st_it
    PetscErrorCode                            :: ierr
    PetscMPIInt                               :: rank, nb_procs
    MPI_Comm                                  :: comm
    CALL MPI_Comm_rank(comm, rank, ierr)
    CALL MPI_Comm_size(comm, nb_procs, ierr)
    ALLOCATE(file_list(nb_procs))
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

    CALL check_list(comm, file_list, mesh%np)
    IF (rank==0) THEN
       IF (PRESENT(opt_it)) THEN
          it = opt_it
       ELSE
          it = 1
       END IF
       CALL create_pvd_file(file_list, TRIM(header), it, TRIM(what))
    END IF
    CALL create_xml_vtu_scal_file(field, mesh, TRIM(ADJUSTL(file_list(rank+1))), TRIM(ADJUSTL(field_name)))
  END SUBROUTINE make_vtu_file_scalar_2D

END MODULE two_dim_vtu_xml
