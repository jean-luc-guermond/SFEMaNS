MODULE COMM_LK_SFEM
        USE def_type_field_vel
        USE def_type_field
        IMPLICIT NONE

        TYPE, EXTENDS(abstract_vector_rdp), PUBLIC :: LK_FIELD_TYPE
            TYPE(mag_field_type) :: mag_field
            TYPE(vel_field_type) :: vel_n, vel_n1
            LOGICAL :: if_vel, if_mag
        CONTAINS
            PROCEDURE, PUBLIC, PASS(self) :: LK_2_sfem, sfem_2_LK
            PROCEDURE, PUBLIC, PASS(self) :: my_final_post_pro
    !======     Procedures for LightKrylov
            PROCEDURE, PUBLIC, PASS(self) :: zero     => zero_LK 
            PROCEDURE, PUBLIC, PASS(self) :: dot      => dot_LK 
            PROCEDURE, PUBLIC, PASS(self) :: axpby    => axpby_LK 
            PROCEDURE, PUBLIC, PASS(self) :: rand     => rand_LK 
            PROCEDURE, PUBLIC, PASS(self) :: scal     => scal_LK 
            PROCEDURE, PUBLIC, PASS(self) :: get_size => get_size_LK 
        END TYPE LK_FIELD_TYPE

        LOGICAL, PRIVATE :: if_momentum, if_induction
        PUBLIC :: init_LK

CONTAINS

!=================================
!================================= 
!==== SFEMaNS PROCEDURES =========
!================================= 
!=================================

        SUBROUTINE init_LK
            USE input_data
            IMPLICIT NONE
            
            if_momentum = inputs%type_pb=='nst' .OR. inputs%type_pb=='mhd' .OR. inputs%type_pb=='fhd' &
         .OR. inputs%type_pb=='mhs'
            if_induction = inputs%type_pb=='mxw' .OR. inputs%type_pb=='mhd' .OR. inputs%type_pb=='mxx' &
         .OR. inputs%type_pb=='fhd' .OR. inputs%type_pb=='mhs'

        END SUBROUTINE init_LK

!=====================================
!===================================== 
!==== LIGHTKRYLOV PROCEDURES =========
!===================================== 
!=====================================

        SUBROUTINE zero_LK(self)
            IMPLICIT NONE
            CLASS(LK_field_type), INTENT(INOUT) :: self
            
            IF (if_induction) THEN
                CALL self%mag_field%zero
            ENDIF
            IF (if_momentum) THEN
                CALL self%vel_n%zero
                CALL self%vel_n1%zero
            ENDIF
        END SUBROUTINE zero_LK

        FUNCTION dot_LK(self, vec) RESULT(alpha)
            USE my_util
            IMPLICIT NONE
            CLASS(LK_field_type),       INTENT(IN) :: self
            CLASS(abstract_vector_rdp), INTENT(IN) :: vec
            REAL(KIND=8) :: alpha

            SELECT TYPE(vec)
            TYPE is(LK_field_type)
                alpha = 0.d0
                IF (if_induction) THEN
                    alpha = alpha + self%mag_field%dot(vec%mag_field)
                ENDIF
                IF (if_momentum) THEN
                    alpha = alpha + self%vel_n%dot(vec%vel_n)
                    alpha = alpha + self%vel_n1%dot(vec%vel_n1)
                ENDIF
            CLASS DEFAULT
                CALL error_petsc("BUG in zero_LK => wrong class for LK_field")
            END SELECT

        END FUNCTION dot_LK

        SUBROUTINE axpby_LK(alpha, vec, beta, self)
            USE my_util
            IMPLICIT NONE
            CLASS(LK_field_type),          INTENT(INOUT) :: self
            CLASS(abstract_vector_rdp),    INTENT(IN)    :: vec
            REAL(KIND=8),                  INTENT(IN)    :: alpha, beta

            SELECT TYPE(vec)
            TYPE is(LK_field_type)
                IF (if_induction) THEN
                    CALL self%mag_field%axpby(alpha, vec%mag_field, beta)
                ENDIF
                IF (if_momentum) THEN
                    CALL self%vel_n%axpby(alpha, vec%vel_n, beta)
                    CALL self%vel_n1%axpby(alpha, vec%vel_n1, beta)
                ENDIF
            CLASS DEFAULT
                CALL error_petsc("BUG in axpby_LK => wrong class for LK_field")
            END SELECT

        END SUBROUTINE axpby_LK

        SUBROUTINE rand_LK(self, ifnorm)
            IMPLICIT NONE
            CLASS(LK_field_type), INTENT(INOUT)   :: self
            LOGICAL, OPTIONAL,     INTENT(IN)      :: ifnorm

            IF (if_induction) THEN
                CALL self%mag_field%rand(ifnorm)
            ENDIF
            IF (if_momentum) THEN
                CALL self%vel_n%rand(ifnorm)
                CALL self%vel_n%copy_to(self%vel_n1)
            ENDIF
        END SUBROUTINE rand_LK

        SUBROUTINE scal_LK(self, alpha)
            IMPLICIT NONE
            CLASS(LK_field_type), INTENT(INOUT)       :: self
            REAL(KIND=8),          INTENT(IN)          :: alpha
            
            IF (if_induction) THEN
                CALL self%mag_field%scal(alpha)
            ENDIF
            IF (if_momentum) THEN
                CALL self%vel_n%scal(alpha)
                CALL self%vel_n1%scal(alpha)
            ENDIF
        END SUBROUTINE scal_LK

        INTEGER FUNCTION get_size_LK(self) RESULT(N)
            IMPLICIT NONE
            CLASS(LK_field_type), INTENT(IN) :: self
            N = 0
            IF (if_induction) N = N + self%mag_field%get_size()
            IF (if_momentum) N = N + self%vel_n%get_size() + self%vel_n1%get_size()
        END FUNCTION get_size_LK

!=======================================
!======================================= 
!==== COMMUNICATION PROCEDURES =========
!======================================= 
!=======================================

        SUBROUTINE LK_2_SFEM(self, un, pn, incpn, un_m1, pn_m1, incpn_m1, mag_field)
            USE def_type_field
            IMPLICIT NONE

            CLASS(LK_field_type), INTENT(IN)       :: self
            TYPE(mag_field_type), INTENT(OUT)      :: mag_field
            REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT) :: un, pn, incpn, un_m1, pn_m1, incpn_m1

            IF (if_induction) THEN
                CALL self%mag_field%copy_to(mag_field)
            ENDIF

            IF (if_momentum) THEN
                un = self%vel_n%u
                pn = self%vel_n%p
                incpn = self%vel_n%incp

                un_m1 = self%vel_n1%u
                pn_m1 = self%vel_n1%p
                incpn_m1 = self%vel_n1%incp
            ENDIF

        END SUBROUTINE LK_2_SFEM


        SUBROUTINE SFEM_2_LK(self, un, pn, incpn, un_m1, pn_m1, incpn_m1, mag_field)
            USE def_type_field
            IMPLICIT NONE

            TYPE(mag_field_type), INTENT(IN)   :: mag_field
            CLASS(LK_field_type), INTENT(OUT)  :: self
            REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN) :: un, pn, incpn, un_m1, pn_m1, incpn_m1
            
            IF (if_induction) THEN
                CALL mag_field%copy_to(self%mag_field)
            ENDIF

            IF (if_momentum) THEN
                IF (.NOT. ALLOCATED(self%vel_n%time))  ALLOCATE(self%vel_n%time)
                self%vel_n%time = 0.d0!mag_field%time
                CALL self%vel_n%allocate_velocity_fields()
                self%vel_n%u = un
                self%vel_n%p = pn
                self%vel_n%incp = incpn
                self%vel_n%initialized = .TRUE.

                IF (.NOT. ALLOCATED(self%vel_n1%time))  ALLOCATE(self%vel_n1%time)
                self%vel_n1%time = 0.d0!mag_field%time
                CALL self%vel_n1%allocate_velocity_fields()
                self%vel_n1%u = un_m1
                self%vel_n1%p = pn_m1
                self%vel_n1%incp = incpn_m1
                self%vel_n1%initialized = .TRUE.
            ENDIF

        END SUBROUTINE SFEM_2_LK

        SUBROUTINE my_final_post_pro(self, it)
            USE input_data
            USE restart
            USE fourier_to_real_for_vtu
            IMPLICIT NONE
            CLASS(LK_field_type), INTENT(IN) :: self
            INTEGER,              INTENT(IN) :: it

            !!!!TODOOOOOOOOOOO

            ! CALL write_restart_maxwell(comm_one_d, H_mesh, phi_mesh, &
            !     0.d0, list_mode, self%mag_field, inputs%file_name, it, 1)
            ! IF (inputs%LK%if_vtu_3d) THEN
            !     CALL vtu_3d(comm_one_d, self%mag_field%Hn, 'H_mesh', 'Hn', 'Hn', what, opt_it=it)
            ! END IF

            !!!!TODOOOOOOOOOOO


        END SUBROUTINE my_final_post_pro

END MODULE COMM_LK_SFEM
