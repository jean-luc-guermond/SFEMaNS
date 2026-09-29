MODULE newton_krylov_module
    PUBLIC :: compute_base_flow

    ABSTRACT INTERFACE
        SUBROUTINE template_newton_initial_guess(bf)
            USE COMM_LK_SFEM, ONLY: LK_field_type
            USE LightKrylov
            IMPLICIT NONE
            CLASS(abstract_vector_rdp), INTENT(INOUT) :: bf
        END SUBROUTINE template_newton_initial_guess
    END INTERFACE

CONTAINS

    SUBROUTINE compute_base_flow(initial_guess)
        USE input_data
        USE jacobian_type_module, ONLY: jacobian_type, sys_SFEMaNS_jacobian, post_pro_newton_LK
        USE def_type_mesh
        USE restart
        USE my_util
        USE fourier_to_real_for_vtu
        USE COMM_LK_SFEM
#ifdef USE_LIGHTKRYLOV
        USE LightKrylov
        USE LightKrylov_Constants
        USE LightKrylov_Logger
        USE LightKrylov_Timing, ONLY: timer => global_lightkrylov_timer
        USE LightKrylov_Utils
        USE lightkrylov_IterativeSolvers, ONLY: gmres_rdp
#endif

        IMPLICIT NONE
        INTEGER                                         :: info
        PROCEDURE(template_newton_initial_guess)        :: initial_guess
#ifdef USE_LIGHTKRYLOV
        TYPE(sys_SFEMaNS_jacobian)                      :: sys_jac
        ! TYPE(jacobian_type)                             :: jacobian
        TYPE(LK_field_type)     :: bf, residual, some_vec
        TYPE(newton_dp_opts)    :: opts
#endif

#ifdef USE_LIGHTKRYLOV
        CALL timer%initialize()
        CALL timer%add_timer('Newton iteration (const. tol)', start=.true.)

        ! Set Jacobian and initial guess for baseflow
        CALL bf%zero(); CALL residual%zero(); CALL some_vec%zero()
        sys_jac%jacobian = jacobian_type()!jacobian
        CALL initial_guess(bf)
        sys_jac%jacobian%X = bf
        CALL sys_jac%eval(bf, residual, inputs%LK%rtol_newton)

write(*,*) "norm ini = ", residual%norm(), bf%norm(), residual%norm()/bf%norm()
CALL residual%rand()
CALL sys_jac%jacobian%matvec(residual, some_vec)
write(*,*) "norm some vec = ", residual%norm(), some_vec%norm(), sys_jac%jacobian%X%norm()
! call error_petsc("end code")

        opts = newton_dp_opts(maxiter=30, ifbisect=.false.)
        CALL newton(sys_jac, bf, gmres_rdp, info, rtol=inputs%LK%rtol_newton, atol=inputs%LK%atol_newton, options=opts, scheduler=constant_tol_dp)
        CALL sys_jac%eval(bf, residual, inputs%LK%rtol_newton)
        CALL timer%stop('Newton iteration (const. tol)')
write(*,*) "Final residual norm = ", residual%norm(), bf%norm(), residual%norm() / bf%norm()

        CALL post_pro_newton_LK(bf, residual)
#else
        !2026/06/05: LC added dummy variables to avoid warning
        info = 0
        !2026/06/05: LC added dummy variables to avoid warning
        CALL error_petsc("LightKrylov not compiled with SFEMaNS, check CMakeLists.txt/variable.cmake")
#endif

      END SUBROUTINE compute_base_flow
END MODULE newton_krylov_module
