MODULE eigenvalue_problem
   PUBLIC :: compute_spectrum

   CONTAINS

    SUBROUTINE compute_spectrum(comm_one_d, H_mesh, phi_mesh, lambda, residuals, X, list_mode)
        USE input_data
        USE exponential_propagator
        USE def_type_mesh
        USE def_type_field
        USE restart
        USE my_util
#include "petsc/finclude/petsc.h"
        USE petsc
        USE fourier_to_real_for_vtu
#ifdef USE_LIGHTKRYLOV
        USE LightKrylov
        USE LightKrylov_Constants
        USE LightKrylov_Logger
        USE stdlib_logger, only: information_level, warning_level, debug_level, error_level, all_level, success
        USE stdlib_io_npy, only: save_npy

        IMPLICIT NONE
        INTEGER, DIMENSION(:), INTENT(IN)               :: list_mode
        COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE,      INTENT(OUT) :: lambda
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE,         INTENT(OUT) :: residuals
        TYPE(mag_field_type), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: X
        INTEGER                                         :: info, it
        CHARACTER(LEN=3)                                :: what
        TYPE(exptA_linop)                               :: exptA
        TYPE(mesh_type), POINTER                        :: H_mesh, phi_mesh
        MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d

        ALLOCATE(X(inputs%LK%nev))
        CALL logger_setup(nio=0, log_level=information_level, log_stdout=.TRUE., log_timestamp=.false.)
        CALL zero_basis(X)
        CALL eigs(exptA, X, lambda, residuals, info, &
                kdim = inputs%LK%kdim, tolerance=inputs%LK%abs_tol)
        CALL check_info(info, 'eigs')

        lambda = log(lambda)/(inputs%dt*inputs%nb_iteration)
        
        CALL save_eigenspectrum(lambda, residuals, "./eigenspectrum.npy")

        what = 'new'
        DO it=1, inputs%LK%nev
            CALL write_restart_maxwell(comm_one_d, H_mesh, phi_mesh, &
                0.d0, list_mode, X(it), inputs%file_name, it, 1)
            IF (inputs%LK%if_vtu_3d) THEN
                CALL vtu_3d(comm_one_d, X(it)%Hn, 'H_mesh', 'Hn', 'Hn', what, opt_it=it)
                what = 'old'
            END IF
        END DO
#else
        CALL error_petsc("LightKrylov not compiled with SFEMaNS, check CMakeLists.txt/variable.cmake")
#endif

      END SUBROUTINE compute_spectrum
END MODULE eigenvalue_problem