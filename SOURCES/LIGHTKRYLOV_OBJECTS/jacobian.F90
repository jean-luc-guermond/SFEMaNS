MODULE jacobian_type_module
#ifdef USE_LIGHTKRYLOV
    USE LightKrylov
    USE LightKrylov_Logger
    USE LightKrylov_Constants
#endif
    USE COMM_LK_SFEM

#ifdef USE_LIGHTKRYLOV
    TYPE, EXTENDS(abstract_system_rdp), PUBLIC :: sys_SFEMaNS_jacobian
#else
    TYPE, PUBLIC :: sys_SFEMaNS_jacobian
#endif
    CONTAINS
        PRIVATE
        PROCEDURE, PASS(self), PUBLIC :: response => nonlinear_SFEMaNS_solver
    END TYPE sys_SFEMaNS_jacobian

#ifdef USE_LIGHTKRYLOV
    TYPE, EXTENDS(abstract_jacobian_linop_rdp), PUBLIC :: jacobian_type
#else
    TYPE, PUBLIC :: jacobian_type
        CLASS(LK_FIELD_TYPE), ALLOCATABLE :: X
#endif
    CONTAINS
        PRIVATE
        PROCEDURE, PASS(self), PUBLIC :: matvec  => linear_SFEMaNS_solver
        PROCEDURE, PASS(self), PUBLIC :: rmatvec => linear_SFEMaNS_solver
    END TYPE jacobian_type
    
    PUBLIC :: post_pro_newton_LK
    CHARACTER(len=*), PARAMETER, PRIVATE :: this_module = 'jacobian_type_module'
CONTAINS

    SUBROUTINE nonlinear_SFEMaNS_solver(self, vec_in, vec_out, atol)
        USE input_data, ONLY: inputs
        USE exponential_propagator
        IMPLICIT NONE
        REAL(KIND=8), INTENT(IN)  :: atol
        CLASS(sys_SFEMaNS_jacobian),     INTENT(INOUT)  :: self
#ifdef USE_LIGHTKRYLOV
        CLASS(abstract_vector_rdp), INTENT(IN)  :: vec_in
        CLASS(abstract_vector_rdp), INTENT(OUT) :: vec_out
        ! Solver tolerances if needed (forced by LightKrylov, unnecessary inside SFEMaNS)

        SELECT TYPE (vec_in)
        TYPE is (LK_field_type)
            SELECT TYPE (vec_out)
            TYPE is (LK_field_type)
#else
        TYPE(LK_FIELD_TYPE), INTENT(IN)  :: vec_in
        TYPE(LK_FIELD_TYPE), INTENT(OUT) :: vec_out
#endif

                !=== Switch SFEMaNS to nonlinear mode
                inputs%if_linear_navier_stokes_with_u = .FALSE.
                !=== SFEMaNS integration
                CALL SFEMaNS_advance(vec_in, vec_out)
                !=== Evaluate residual
                CALL vec_out%axpby(-1.d0, vec_in, 1.d0)
                !CALL vec_out%sub(vec_in)

#ifdef USE_LIGHTKRYLOV
            CLASS DEFAULT
                CALL TYPE_error('vec_out','LK_field_type','OUT',this_module,'nonlinear_map')
            END SELECT
        CLASS DEFAULT
            CALL TYPE_error('vec_in','LK_field_type','IN',this_module,'nonlinear_map')
        END SELECT
#endif

        RETURN
        !=== dummy to avoid warning
        WRITE(*,*) atol
        !=== dummy to avoid warning
    END SUBROUTINE nonlinear_SFEMaNS_solver

    SUBROUTINE linear_SFEMaNS_solver(self, vec_in, vec_out)
        USE input_data, ONLY: inputs
        USE base_field_module, ONLY: reset_base_flow
        USE exponential_propagator
        IMPLICIT NONE
        CLASS(jacobian_type),     INTENT(INOUT)  :: self
#ifdef USE_LIGHTKRYLOV
        CLASS(abstract_vector_rdp), INTENT(IN)  :: vec_in
        CLASS(abstract_vector_rdp), INTENT(OUT) :: vec_out
        ! Solver tolerances if needed (forced by LightKrylov, unnecessary inside SFEMaNS)

        SELECT TYPE (vec_in)
        TYPE is (LK_field_type)
            SELECT TYPE (vec_out)
            TYPE is (LK_field_type)
#else
        TYPE(LK_field_type), INTENT(IN)  :: vec_in
        TYPE(LK_field_type), INTENT(OUT) :: vec_out
#endif
#ifdef USE_LIGHTKRYLOV
                !=== Switch SFEMaNS to linear mode
                inputs%if_linear_navier_stokes_with_u = .TRUE.
                SELECT TYPE (bf => self%X)
                TYPE is (LK_field_type)
                    CALL reset_base_flow(bf%vel_n%u)
                CLASS DEFAULT
                    CALL TYPE_error('base flow','LK_field_type','OUT',this_module,'linear_map')
                END SELECT
                !=== SFEMaNS integration
                CALL SFEMaNS_advance(vec_in, vec_out)
                !=== Evaluate residual
                !CALL vec_out%sub(vec_in)
                CALL vec_out%axpby(-1.d0, vec_in, 1.d0)

            CLASS DEFAULT
                CALL TYPE_error('vec_out','LK_field_type','OUT',this_module,'linear_map')
            END SELECT
        CLASS DEFAULT
            CALL TYPE_error('vec_in','LK_field_type','IN',this_module,'linear_map')
        END SELECT
#else
                !=== Switch SFEMaNS to linear mode
                inputs%if_linear_navier_stokes_with_u = .TRUE.
                CALL reset_base_flow(self%X%vel_n%u)
                !=== SFEMaNS integration
                CALL SFEMaNS_advance(vec_in, vec_out)
                !=== Evaluate residual
                !CALL vec_out%sub(vec_in)
                CALL vec_out%axpby(-1.d0, vec_in, 1.d0)

#endif
        RETURN
    END SUBROUTINE linear_SFEMaNS_solver

    !=========================================================
    !====== Post processing for Newton method  ===============
    !=========================================================

    SUBROUTINE post_pro_newton_LK(bf_initial_guess, bf_output)
        USE input_data
        USE restart
        USE fourier_to_real_for_vtu
        USE SFEMaNS_logicals_module
#ifdef USE_LIGHTKRYLOV
        USE LightKrylov
        USE LightKrylov_Constants
        USE LightKrylov_Logger
        USE stdlib_logger, only: information_level, warning_level, debug_level, error_level, all_level, success
        USE stdlib_io_npy, only: save_npy

        IMPLICIT NONE
        TYPE(LK_FIELD_TYPE), INTENT(IN) :: bf_initial_guess, bf_output
        INTEGER            :: it, i
        CHARACTER(LEN=3)   :: what
        CHARACTER(LEN=3)   :: st_mode
        CHARACTER(LEN=200) :: header

        it = 0

        IF (if_momentum) THEN
            ASSOCIATE(vn => bf_initial_guess%vel_n, vn1 => bf_initial_guess%vel_n1)
            CALL write_restart_ns(vn%comm_one_d_ns, vn%vv_mesh, vn%pp_mesh, 0.d0, &
                vn%list_mode, vn%u, vn1%u, vn%p, vn1%p, &
                vn%incp, vn1%incp, inputs%file_name, 0, 1, opt_name="suite_bf_ns")
            END ASSOCIATE

            ASSOCIATE(vn => bf_output%vel_n, vn1 => bf_output%vel_n1)
            CALL write_restart_ns(vn%comm_one_d_ns, vn%vv_mesh, vn%pp_mesh, 0.d0, &
                vn%list_mode, vn%u, vn1%u, vn%p, vn1%p, &
                vn%incp, vn1%incp, inputs%file_name, 1, 1, opt_name="suite_bf_ns")
            END ASSOCIATE
            ! ASSOCIATE(vn => X(it)%vel)
            ! CALL write_restart_ns(vn%comm_one_d_ns, vn%vv_mesh, vn%pp_mesh, 0.d0, &
                ! vn%list_mode, vn%u, vn%u, vn%p, vn%p, &
                ! vn%incp, vn%incp, inputs%file_name, it, 1, opt_name="suite_eigs_ns")
            ! END ASSOCIATE
        ENDIF

        IF (if_induction) THEN
            ASSOCIATE(mag_n => bf_output%mag_n, mag_n1 => bf_output%mag_n)
            CALL write_restart_maxwell(mag_n%comm_one_d, mag_n%H_mesh, mag_n%phi_mesh, &
                    0.d0, mag_n%list_mode, mag_n, mag_n1, inputs%file_name, it, 1, opt_name="suite_bf_maxwell")
            END ASSOCIATE
            ! ASSOCIATE(mag => X(it)%mag)
            ! CALL write_restart_maxwell(mag%comm_one_d, mag%H_mesh, mag%phi_mesh, &
                    ! 0.d0, mag%list_mode, mag, mag, inputs%file_name, it, 1, opt_name="suite_eigs_maxwell")
            ! END ASSOCIATE
        ENDIF

        IF (inputs%LK%if_vtu_3d) THEN
            what = 'new'
            IF (if_induction) THEN
                ASSOCIATE(mag_n => bf_output%mag_n)
                CALL vtu_3d(mag_n%comm_one_d, mag_n%H, 'H_mesh', 'Hn_bf', 'Hn_bf', what, opt_it=it)
                END ASSOCIATE
                ! ASSOCIATE(mag => X(it)%mag)
                ! CALL vtu_3d(mag%comm_one_d, mag%H, 'H_mesh', 'Hn_eig', 'Hn_eig', what, opt_it=it)
                ! what = 'old'
                ! END ASSOCIATE
            ENDIF

            what = 'new'
            IF (if_momentum) THEN
                ASSOCIATE(vn => bf_output%vel_n, vn1 => bf_output%vel_n1)
                ! ASSOCIATE(vn => X(it)%vel)
                    CALL vtu_3d(vn%comm_one_d_ns, vn%u, 'vv_mesh', 'vn_bf', 'vn_bf', what, opt_it=it)
                    CALL vtu_3d(vn%comm_one_d_ns, vn%p, 'pp_mesh', 'pn_bf', 'pn_bf', what, opt_it=it)
                    IF (inputs%if_plot_2D) THEN
                        DO i=1, SIZE(vn%list_mode)
                            WRITE(st_mode, '(I3)') vn%list_mode(i)
                            header = 'vn_bf_mode_' // TRIM(ADJUSTL(st_mode))
                            CALL make_vtu_file_2D(vn%comm_one_d_ns(1), vn%vv_mesh, header, vn%u(:, :, i), &
                                            'vn_bf', what, opt_it = it)

                            WRITE(st_mode, '(I3)') vn%list_mode(i)
                            header = 'pn_bf_mode_' // TRIM(ADJUSTL(st_mode))
                            CALL make_vtu_file_2D(vn%comm_one_d_ns(1), vn%pp_mesh, header, vn%p(:, :, i), &
                                            'pn_bf', what, opt_it = it)
                        END DO
                    ENDIF
                END ASSOCIATE
            ENDIF
        ENDIF
#else
        USE my_util
        IMPLICIT NONE
        TYPE(LK_FIELD_TYPE), INTENT(INOUT) :: bf_initial_guess, bf_output
        CALL error_petsc("BUG in post_pro_newton_LK: SFEMaNS compiled without LightKrylov => this subroutine should not be called!!!")
        !== Dummy to avoid warnings
        RETURN
        CALL bf_initial_guess%zero()
        CALL bf_output%zero()
        !== Dummy to avoid warnings
#endif            

    END SUBROUTINE post_pro_newton_LK

END MODULE jacobian_type_module
