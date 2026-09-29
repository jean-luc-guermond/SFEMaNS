MODULE base_field_module
    USE def_type_mesh
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: U_bar, B_bar
    TYPE(mesh_type), POINTER :: H_mesh
    PUBLIC :: init_base_flow, init_imposed_mag_field
    PUBLIC :: get_base_flow, get_project_base_flow, get_imposed_mag_field, get_project_imposed_mag_field
    PUBLIC :: reset_base_flow
CONTAINS

!===============================================
!======= Initialization subroutines ============
!===============================================

    SUBROUTINE init_base_flow(vv_mesh, pp_mesh, comm_one_d_ns, list_mode)
        USE input_data, ONLY: inputs
        USE restart
        USE boundary,   ONLY: base_flow
        IMPLICIT NONE
        TYPE(mesh_type),       INTENT(IN) :: vv_mesh, pp_mesh
        INTEGER, DIMENSION(2), INTENT(IN) :: comm_one_d_ns
        INTEGER, DIMENSION(:), INTENT(IN) :: list_mode
        REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dummy_u, dummy_p
        REAL(KIND=8) :: tps
        INTEGER :: k, i, mode

        ALLOCATE(U_bar(vv_mesh%np, 6, SIZE(list_mode)))
        IF (inputs%if_read_start_linear_navier_stokes_with_u) THEN
            ALLOCATE(dummy_u(vv_mesh%np, 6, SIZE(list_mode)))
            ALLOCATE(dummy_p(pp_mesh%np, 6, SIZE(list_mode)))
            CALL read_restart_ns(comm_one_d_ns, tps, list_mode, U_bar, dummy_u, &
                                dummy_p, dummy_p, dummy_p, dummy_p, inputs%file_name, opt_name="suite_jacob_ns")
            DEALLOCATE(dummy_u, dummy_p)
        ELSE
            DO k = 1, 6
                DO i=1, SIZE(list_mode)
                    mode = list_mode(i)
                    U_bar(:, k, i) = base_flow(k, vv_mesh%rr, mode)
                ENDDO
            ENDDO
        ENDIF

    END SUBROUTINE init_base_flow

    SUBROUTINE init_imposed_mag_field(H_mesh_in, list_mode)
        USE input_data, ONLY: inputs
        USE boundary,   ONLY: H_B_quasi_static
        IMPLICIT NONE
        TYPE(mesh_type), TARGET, INTENT(IN) :: H_mesh_in
        INTEGER, DIMENSION(:),   INTENT(IN) :: list_mode
        INTEGER :: i, mode
        ALLOCATE(B_bar(H_mesh_in%np, 6, SIZE(list_mode)))
        DO i=1, SIZE(list_mode)
            mode = list_mode(i)
            B_bar(:,:,i) = H_B_quasi_static('F', H_mesh_in%rr, mode)
        ENDDO
        H_mesh => H_mesh_in
    END SUBROUTINE init_imposed_mag_field
    
!===================================================
!======= Getting base field subroutines ============
!===================================================

    SUBROUTINE get_base_flow(u)
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT) :: u
        u = U_bar
    END SUBROUTINE get_base_flow

    SUBROUTINE get_project_base_flow(mesh,  which_projection, if_restriction, list_mode, time, u)
        USE def_type_mesh
        USE projection_module, ONLY: projection_velocity
        IMPLICIT NONE
        TYPE(mesh_type),       INTENT(IN) :: mesh
        CHARACTER(LEN=*),      INTENT(IN) :: which_projection
        LOGICAL,               INTENT(IN) :: if_restriction
        INTEGER, DIMENSION(:), INTENT(IN) :: list_mode
        REAL(KIND=8),          INTENT(IN) :: time
        REAL(KIND=8), DIMENSION(mesh%np, 6, SIZE(U_bar, 3)), INTENT(INOUT) :: u

        CALL projection_velocity(mesh, U_bar, which_projection, if_restriction, u, list_mode, time, if_opt_zero_extension=.FALSE.)
    END SUBROUTINE get_project_base_flow

    SUBROUTINE get_imposed_mag_field(b)
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT) :: b
        b = B_bar
    END SUBROUTINE get_imposed_mag_field

    SUBROUTINE get_project_imposed_mag_field(mesh,  which_projection, if_restriction, B)
        USE def_type_mesh
        USE projection_module, ONLY: projection_mag_field
        IMPLICIT NONE
        TYPE(mesh_type),  INTENT(IN) :: mesh
        CHARACTER(LEN=*), INTENT(IN) :: which_projection
        LOGICAL,          INTENT(IN) :: if_restriction
        REAL(KIND=8), DIMENSION(mesh%np, 6, SIZE(B_bar,3)), INTENT(INOUT) :: B

        CALL projection_mag_field(mesh, H_mesh, B_bar, which_projection, if_restriction, B)

    END SUBROUTINE get_project_imposed_mag_field

!===========================================================================================
!======= Update subroutines (mandatory for Newton or time-dependant base field) ============
!===========================================================================================

    SUBROUTINE reset_base_flow(new_U_bar)
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(size(U_bar,1), 6, size(U_bar,3)), INTENT(IN) :: new_U_bar

        U_bar = new_U_bar
    END SUBROUTINE reset_base_flow
    
END MODULE base_field_module