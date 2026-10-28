MODULE newton_krylov_module
    PUBLIC :: compute_base_flow

    ABSTRACT INTERFACE
        SUBROUTINE template_newton_initial_guess(bf)
            USE COMM_LK_SFEM, ONLY: LK_field_type
#ifdef USE_LIGHTKRYLOV
            USE LightKrylov
            IMPLICIT NONE
            CLASS(abstract_vector_rdp), INTENT(INOUT) :: bf
#else
            IMPLICIT NONE
            TYPE(LK_FIELD_TYPE), INTENT(INOUT) :: bf
#endif
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
        TYPE(LK_field_type)     :: bf_in, bf_out, residual
        TYPE(newton_dp_opts)    :: opts
        INTEGER, SAVE           :: it=0
        REAL(KIND=8) :: t_dummy
#endif

#ifdef USE_LIGHTKRYLOV
        CALL logger_setup(nio=0, log_level=information_level, log_stdout=.TRUE., log_timestamp=.false.)
        CALL timer%initialize()
        CALL timer%add_timer('Newton iteration (const. tol)', start=.true.)

        ! Set Jacobian and initial guess for baseflow
        CALL bf_in%zero(); CALL bf_out%zero(); CALL residual%zero()
        sys_jac%jacobian = jacobian_type()!jacobian
        IF (inputs%LK%read_suite_initial_guess) THEN
                CALL read_restart_ns(bf_in%vel_n%comm_one_d_ns, t_dummy, bf_in%vel_n%list_mode, bf_in%vel_n%u, bf_in%vel_n1%u, bf_in%vel_n%p, bf_in%vel_n1%p, &
                        bf_in%vel_n%incp, bf_in%vel_n1%incp, inputs%file_name, opt_name="suite_bf_ns")
        ELSE
                CALL initial_guess(bf_in)
        ENDIF
        CALL bf_out%axpby(1.d0, bf_in, 0.d0)
write(*,*) "Initial norm = ", bf_out%norm()
        sys_jac%jacobian%X = bf_out
        CALL sys_jac%eval(bf_out, residual, inputs%LK%rtol_newton)

        opts = newton_dp_opts(maxiter=inputs%LK%maxiter_newton, ifbisect=.false.)
        CALL newton(sys_jac, bf_out, gmres_rdp, info, rtol=inputs%LK%rtol_newton, atol=inputs%LK%atol_newton, options=opts, scheduler=constant_tol_dp)
        CALL sys_jac%eval(bf_out, residual, inputs%LK%rtol_newton)
        CALL timer%stop('Newton iteration (const. tol)')
write(*,*) "Final residual norm = ", residual%norm(), bf_out%norm(), residual%norm() / bf_out%norm()

        CALL post_pro_newton_LK(bf_in, bf_out)
#else
        CALL error_petsc("LightKrylov not compiled with SFEMaNS, check CMakeLists.txt and/or variable.cmake")
        !2026/06/05: LC added dummy variables to avoid warning
        write(*,*) info
        !2026/06/05: LC added dummy variables to avoid warning
#endif

      END SUBROUTINE compute_base_flow
END MODULE newton_krylov_module
