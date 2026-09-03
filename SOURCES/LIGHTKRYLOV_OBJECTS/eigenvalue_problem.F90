MODULE eigenvalue_problem
   PUBLIC :: compute_spectrum

   CONTAINS

    SUBROUTINE compute_spectrum(lambda, residuals, X, list_mode)
        USE input_data
        USE exponential_propagator
        USE def_type_mesh
        USE restart
        USE my_util
#include "petsc/finclude/petsc.h"
        USE petsc
        USE fourier_to_real_for_vtu
        USE COMM_LK_SFEM
#ifdef USE_LIGHTKRYLOV
        USE LightKrylov
        USE LightKrylov_Constants
        USE LightKrylov_Logger
        USE stdlib_logger, only: information_level, warning_level, debug_level, error_level, all_level, success
        USE stdlib_io_npy, only: save_npy
#endif

        IMPLICIT NONE
        INTEGER, DIMENSION(:), INTENT(IN)               :: list_mode
        COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE,      INTENT(OUT) :: lambda
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE,         INTENT(OUT) :: residuals
        TYPE(LK_field_type), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: X
        INTEGER                                         :: info, it
        CHARACTER(LEN=3)                                :: what
#ifdef USE_LIGHTKRYLOV
        TYPE(exptA_linop)                               :: exptA
#endif
        TYPE(mesh_type), POINTER                        :: H_mesh, phi_mesh
        MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d

#ifdef USE_LIGHTKRYLOV
        ALLOCATE(X(inputs%LK%nev))
        CALL logger_setup(nio=0, log_level=information_level, log_stdout=.TRUE., log_timestamp=.false.)
        CALL zero_basis(X)
        CALL eigs(exptA, X, lambda, residuals, info, &
                kdim = inputs%LK%kdim, tolerance=inputs%LK%abs_tol)
        CALL check_info(info, 'eigs')

        lambda = log(lambda)/(inputs%dt*inputs%nb_iteration)
        
        CALL post_pro_eigs_LK(X, lambda, residuals)

#else
        !2026/06/05: LC added dummy variables to avoid warning
        ALLOCATE(lambda(0), residuals(0), X(0))
        info=H_mesh%np; info=phi_mesh%np; it=SIZE(comm_one_d)
        it=SIZE(list_mode);
        what = "old"
        !2026/06/05: LC added dummy variables to avoid warning
        CALL error_petsc("LightKrylov not compiled with SFEMaNS, check CMakeLists.txt/variable.cmake")
#endif

      END SUBROUTINE compute_spectrum
END MODULE eigenvalue_problem
