MODULE COMM_LK_SFEM
    USE def_type_field_vel
    USE def_type_field_mag
#ifdef USE_LIGHTKRYLOV
    USE LightKrylov
    USE LightKrylov_Constants
    USE LightKrylov_Logger
    USE stdlib_logger, only: information_level, warning_level, debug_level, error_level, all_level, success
    USE stdlib_io_npy, only: save_npy
#endif
    IMPLICIT NONE

#ifdef USE_LIGHTKRYLOV
    TYPE, EXTENDS(abstract_vector_rdp), PUBLIC :: LK_FIELD_TYPE
#else
    TYPE, PUBLIC :: LK_FIELD_TYPE
#endif
        TYPE(mag_field_type) :: mag_n, mag_n1
        TYPE(vel_field_type) :: vel_n, vel_n1
        LOGICAL :: if_vel, if_mag
    CONTAINS
    !======     SFEMaNS <--> LightKrylov Communication Procedures 
        PROCEDURE, PUBLIC, PASS(self) :: LK_2_sfem, sfem_2_LK
    !======     Procedures for LightKrylov
        PROCEDURE, PUBLIC, PASS(self) :: zero     => zero_LK 
        PROCEDURE, PUBLIC, PASS(self) :: dot      => dot_LK 
        PROCEDURE, PUBLIC, PASS(self) :: axpby    => axpby_LK 
        PROCEDURE, PUBLIC, PASS(self) :: rand     => rand_LK 
        PROCEDURE, PUBLIC, PASS(self) :: scal     => scal_LK 
        PROCEDURE, PUBLIC, PASS(self) :: get_size => get_size_LK 
    END TYPE LK_FIELD_TYPE

    LOGICAL, PRIVATE :: if_momentum, if_induction
    PUBLIC :: init_LK, post_pro_eigs_LK

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

    SUBROUTINE post_pro_eigs_LK(X, lambda, residuals)
        USE input_data
        USE restart
        USE fourier_to_real_for_vtu
#ifdef USE_LIGHTKRYLOV
        USE LightKrylov
        USE LightKrylov_Constants
        USE LightKrylov_Logger
        USE stdlib_logger, only: information_level, warning_level, debug_level, error_level, all_level, success
        USE stdlib_io_npy, only: save_npy

        IMPLICIT NONE
        TYPE(LK_FIELD_TYPE), DIMENSION(:), INTENT(IN) :: X
        COMPLEX(KIND=8),     DIMENSION(:), INTENT(IN) :: lambda
        REAL   (KIND=8),     DIMENSION(:), INTENT(IN) :: residuals
        INTEGER          :: it
        CHARACTER(LEN=3) :: what

        CALL save_eigenspectrum(lambda, residuals, "./eigenspectrum.npy")

        IF (inputs%LK%if_write_eigs) THEN
            DO it=1, inputs%LK%nev
                IF (if_momentum) THEN
                    ASSOCIATE(vn => X(it)%vel_n, vn1 => X(it)%vel_n1)
                    CALL write_restart_ns(vn%comm_one_d_ns, vn%vv_mesh, vn%pp_mesh, 0.d0, &
                        vn%list_mode, vn%u, vn1%u, vn%p, vn1%p, &
                        vn%incp, vn1%incp, inputs%file_name, it, 1, opt_name="suite_eigs_ns")
                    END ASSOCIATE
                ENDIF

                IF (if_induction) THEN
                    ASSOCIATE(mag_n => X(it)%mag_n, mag_n1 => X(it)%mag_n)
                    CALL write_restart_maxwell(mag_n%comm_one_d, mag_n%H_mesh, mag_n%phi_mesh, &
                            0.d0, mag_n%list_mode, mag_n, mag_n1, inputs%file_name, it, 1, opt_name="suite_eigs_maxwell")
                    END ASSOCIATE
                ENDIF
            END DO
        ENDIF

        IF (inputs%LK%if_vtu_3d) THEN
            what = 'new'
            IF (if_induction) THEN
                DO it=1, inputs%LK%nev
                    ASSOCIATE(mag_n => X(it)%mag_n)
                    CALL vtu_3d(mag_n%comm_one_d, mag_n%H, 'H_mesh', 'Hn_eig', 'Hn_eig', what, opt_it=it)
                    what = 'old'
                    END ASSOCIATE
                ENDDO
            ENDIF

            what = 'new'
            IF (if_momentum) THEN
                DO it=1, inputs%LK%nev
                    ASSOCIATE(vn => X(it)%vel_n, vn1 => X(it)%vel_n1)
                    CALL vtu_3d(vn%comm_one_d_ns, vn%u, 'vv_mesh', 'vn_eig', 'vn_eig', what, opt_it=it)
                    what = 'old'
                    END ASSOCIATE
                ENDDO
            ENDIF
        ENDIF
#else
        USE my_util
        IMPLICIT NONE
        TYPE(LK_FIELD_TYPE), DIMENSION(:), INTENT(IN) :: X
        COMPLEX(KIND=8),     DIMENSION(:), INTENT(IN) :: lambda
        REAL   (KIND=8),     DIMENSION(:), INTENT(IN) :: residuals
        CALL error_petsc("BUG in post_pro_eigs_LK: SFEMaNS compiled without LightKrylov => this subroutine should not be called!!!")
        !== Dummy to avoid warnings
        RETURN
        write(*,*) size(X), size(lambda), size(residuals)
        !== Dummy to avoid warnings
#endif            

    END SUBROUTINE post_pro_eigs_LK

!=====================================
!===================================== 
!==== LIGHTKRYLOV PROCEDURES =========
!===================================== 
!=====================================

    SUBROUTINE zero_LK(self)
        IMPLICIT NONE
        CLASS(LK_field_type), INTENT(INOUT) :: self
        
        IF (if_induction) THEN
            CALL self%mag_n%zero
            CALL self%mag_n1%zero
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
#ifdef USE_LIGHTKRYLOV
        CLASS(abstract_vector_rdp), INTENT(IN) :: vec
#else
        CLASS(LK_field_type),       INTENT(IN) :: vec
#endif
        REAL(KIND=8) :: alpha

#ifdef USE_LIGHTKRYLOV
        SELECT TYPE(vec)
        TYPE is(LK_field_type)
#endif
            alpha = 0.d0
            IF (if_induction) THEN
                alpha = alpha + self%mag_n%dot(vec%mag_n)
                alpha = alpha + self%mag_n1%dot(vec%mag_n1)
            ENDIF
            IF (if_momentum) THEN
                alpha = alpha + self%vel_n%dot(vec%vel_n)
                alpha = alpha + self%vel_n1%dot(vec%vel_n1)
            ENDIF
#ifdef USE_LIGHTKRYLOV
        CLASS DEFAULT
            CALL error_petsc("BUG in zero_LK => wrong class for LK_field")
        END SELECT
#endif
    END FUNCTION dot_LK

    SUBROUTINE axpby_LK(alpha, vec, beta, self)
        USE my_util
        IMPLICIT NONE
        CLASS(LK_field_type),          INTENT(INOUT) :: self
#ifdef USE_LIGHTKRYLOV
        CLASS(abstract_vector_rdp),    INTENT(IN)    :: vec
#else
        CLASS(LK_field_type),          INTENT(IN)    :: vec
#endif
        REAL(KIND=8),                  INTENT(IN)    :: alpha, beta

#ifdef USE_LIGHTKRYLOV
        SELECT TYPE(vec)
        TYPE is(LK_field_type)
#endif
            IF (if_induction) THEN
                CALL self%mag_n%axpby(alpha, vec%mag_n, beta)
                CALL self%mag_n1%axpby(alpha, vec%mag_n1, beta)
            ENDIF
            IF (if_momentum) THEN
                CALL self%vel_n%axpby(alpha, vec%vel_n, beta)
                CALL self%vel_n1%axpby(alpha, vec%vel_n1, beta)
            ENDIF
#ifdef USE_LIGHTKRYLOV
        CLASS DEFAULT
            CALL error_petsc("BUG in axpby_LK => wrong class for LK_field")
        END SELECT
#endif
    END SUBROUTINE axpby_LK

    SUBROUTINE rand_LK(self, ifnorm)
        IMPLICIT NONE
        CLASS(LK_field_type), INTENT(INOUT)   :: self
        LOGICAL, OPTIONAL,     INTENT(IN)      :: ifnorm

        IF (if_induction) THEN
            CALL self%mag_n%rand(ifnorm)
            CALL self%mag_n%copy_to(self%mag_n1)
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
            CALL self%mag_n%scal(alpha)
            CALL self%mag_n1%scal(alpha)
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
        IF (if_induction) N = N + self%mag_n%get_size() + self%mag_n1%get_size()
        IF (if_momentum) N = N + self%vel_n%get_size() + self%vel_n1%get_size()
    END FUNCTION get_size_LK

!=======================================
!======================================= 
!==== COMMUNICATION PROCEDURES =========
!======================================= 
!=======================================

    SUBROUTINE LK_2_SFEM(self, un, pn, incpn, un_m1, pn_m1, incpn_m1, mag_n, mag_n1)
        USE def_type_field_mag
        IMPLICIT NONE

        CLASS(LK_field_type), INTENT(IN)       :: self
        TYPE(mag_field_type), INTENT(OUT)      :: mag_n, mag_n1
        REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT) :: un, pn, incpn, un_m1, pn_m1, incpn_m1

        IF (if_induction) THEN
            CALL self%mag_n%copy_to(mag_n)
            CALL self%mag_n1%copy_to(mag_n1)
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


    SUBROUTINE SFEM_2_LK(self, un, pn, incpn, un_m1, pn_m1, incpn_m1, mag_n, mag_n1)
        USE def_type_field_mag
        IMPLICIT NONE

        TYPE(mag_field_type), INTENT(IN)   :: mag_n, mag_n1
        CLASS(LK_field_type), INTENT(OUT)  :: self
        REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN) :: un, pn, incpn, un_m1, pn_m1, incpn_m1
            
        IF (if_induction) THEN
            CALL mag_n%copy_to(self%mag_n)
            CALL mag_n1%copy_to(self%mag_n1)
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

END MODULE COMM_LK_SFEM
