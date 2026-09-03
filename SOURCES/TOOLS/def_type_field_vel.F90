MODULE def_type_field_vel
    USE def_type_mesh
#include "petsc/finclude/petsc.h"
    USE petsc
    USE petscmpi
    IMPLICIT NONE

    PUBLIC             :: build_vel_field_params, build_pointers_vel_field, build_global_pointers_vel_field
   
!========================================================
!========================================================
!========    Elements relevant to velocity field comp ===
!========================================================
!========================================================

    TYPE(mesh_type), PRIVATE,               POINTER          :: vv_mesh, pp_mesh, vv_mesh_glob, pp_mesh_glob 
    MPI_Comm,        PRIVATE, DIMENSION(:), POINTER          :: communicator_ns
    INTEGER,         PRIVATE, DIMENSION(:), POINTER          :: list_mode
   

!========================================================
!========================================================
!========    definition of vel_field_type used for nst ==
!========================================================
!========================================================

    TYPE,    PUBLIC ::  vel_field_type
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :, :)   :: u
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :, :)   :: p, incp
        REAL(KIND=8), ALLOCATABLE                       :: time 
        INTEGER, DIMENSION(:), POINTER                  :: list_mode
        TYPE(mesh_type), POINTER                        :: vv_mesh, pp_mesh
        LOGICAL                                         :: initialized = .FALSE. 
        MPI_Comm, DIMENSION(:), POINTER                 :: comm_one_d_ns
    CONTAINS
   !=======    Procedures for SFEMaNS
        PROCEDURE, PUBLIC, PASS(self) :: init_vel_field, allocate_velocity_fields
        PROCEDURE, PUBLIC, PASS(self) :: set_time => set_time_vel
        PROCEDURE, PUBLIC, PASS(self) :: copy_to  => copy_to_vel
   !======     Procedures for LightKrylov
        PROCEDURE, PUBLIC, PASS(self) :: zero     => zero_vel 
        PROCEDURE, PUBLIC, PASS(self) :: dot      => dot_vel 
        PROCEDURE, PUBLIC, PASS(self) :: axpby    => axpby_vel 
        PROCEDURE, PUBLIC, PASS(self) :: rand     => rand_vel 
        PROCEDURE, PUBLIC, PASS(self) :: scal     => scal_vel 
        PROCEDURE, PUBLIC, PASS(self) :: get_size => get_size_vel 
    END TYPE vel_field_type

CONTAINS
!==========================================
!==========================================
!========    Procedures for SFEMaNS =======
!==========================================
!==========================================

    SUBROUTINE build_global_pointers_vel_field(vv_mesh_in, pp_mesh_in)
        USE def_type_mesh

        IMPLICIT NONE
        TYPE(mesh_type), TARGET       ,INTENT(IN)            :: vv_mesh_in, pp_mesh_in

        vv_mesh_glob => vv_mesh_in
        pp_mesh_glob => pp_mesh_in

    END SUBROUTINE build_global_pointers_vel_field

    SUBROUTINE build_vel_field_params()!mu_phi_in, mu_H_field_in, sigma_field_in)
        IMPLICIT NONE
        ! REAL(KIND=8), TARGET, DIMENSION(:), INTENT(IN)       :: mu_H_field_in, sigma_field_in
        ! REAL(KIND=8) :: mu_phi_in
        ! mu_H_field => mu_H_field_in
        ! sigma_field => sigma_field_in
        ! mu_phi = mu_phi_in

    END SUBROUTINE build_vel_field_params
        
    SUBROUTINE build_pointers_vel_field(vv_mesh_in, pp_mesh_in, list_mode_in, communicator_in, if_momentum_in)
        USE def_type_mesh

        IMPLICIT NONE
        TYPE(mesh_type), TARGET       ,INTENT(IN)            :: vv_mesh_in, pp_mesh_in
        MPI_Comm, TARGET, DIMENSION(:),INTENT(IN)            :: communicator_in
        INTEGER, TARGET, DIMENSION(:), INTENT(IN)            :: list_mode_in
        LOGICAL,                       INTENT(IN)            :: if_momentum_in
        vv_mesh           => vv_mesh_in
        pp_mesh           => pp_mesh_in
        communicator_ns   => communicator_in 
        list_mode         => list_mode_in

    END SUBROUTINE build_pointers_vel_field

    SUBROUTINE allocate_velocity_fields(self, if_glob)
        USE def_type_mesh
        
        IMPLICIT NONE

        CLASS(vel_field_type),       INTENT(INOUT)           :: self
        LOGICAL, OPTIONAL,           INTENT(IN)              :: if_glob
        LOGICAL                                              :: if_glob_local, test
        if_glob_local = .FALSE.
        IF (PRESENT(if_glob)) if_glob_local = if_glob
        test = (ALLOCATED(self%u) .AND. ALLOCATED(self%p) .AND. ALLOCATED(self%incp))
        IF (.NOT. test) THEN
            self%comm_one_d_ns  => communicator_ns
            IF (if_glob_local) THEN
                self%vv_mesh => vv_mesh_glob
                self%pp_mesh => pp_mesh_glob
            ELSE
                self%vv_mesh => vv_mesh
                self%pp_mesh => pp_mesh
            END IF
            self%list_mode => list_mode
            ALLOCATE(self%u    (self%vv_mesh%np, 6, SIZE(self%list_mode))) 
            ALLOCATE(self%p    (self%pp_mesh%np, 2, SIZE(self%list_mode)))
            ALLOCATE(self%incp (self%pp_mesh%np, 2, SIZE(self%list_mode)))
        END IF
    END SUBROUTINE allocate_velocity_fields


    SUBROUTINE set_time_vel(self, time)
        USE def_type_mesh
        IMPLICIT NONE
        CLASS(vel_field_type),     INTENT(INOUT)     :: self
        REAL(KIND=8),              INTENT(IN)        :: time
        
        IF (.NOT. ALLOCATED(self%time)) ALLOCATE(self%time)
        self%time = time
    END SUBROUTINE set_time_vel

    SUBROUTINE copy_to_vel(self, vec)
        IMPLICIT NONE
        CLASS(vel_field_type), INTENT(IN)  :: self
        CLASS(vel_field_type), INTENT(OUT) :: vec
        CALL vec%init_vel_field
        vec%u    = self%u
        vec%p    = self%p
        vec%incp = self%incp
    END SUBROUTINE copy_to_vel

!==========================================
!==========================================
!========    Procedures for LightKrylov ===
!==========================================
!==========================================

    SUBROUTINE init_vel_field(self)
        IMPLICIT NONE
        CLASS(vel_field_type),       intent(inout) :: self
        LOGICAL                                    :: test
        test = (ALLOCATED(self%u) .AND. ALLOCATED(self%p) .AND. ALLOCATED(self%incp))
        IF ((.NOT. self%initialized) .OR. (.NOT. test)) THEN
            IF ((.NOT. test) .AND. (self%initialized)) THEN
                WRITE(*,*) "WARNING: somehow vel_field is initialized but with non-allocated stuff"
                WRITE(*,*) 'u', ALLOCATED(self%u) 
                WRITE(*,*) 'p', ALLOCATED(self%p) 
                WRITE(*,*) 'incp', ALLOCATED(self%incp) 
            END IF
            CALL self%zero()
            self%initialized = .TRUE.
        END IF
    END SUBROUTINE init_vel_field

    SUBROUTINE zero_vel(self)
        IMPLICIT NONE
        CLASS(vel_field_type),       INTENT(inout) :: self
        CALL self%allocate_velocity_fields()
        IF (.NOT. ALLOCATED(self%time)) THEN
            ALLOCATE(self%time)
            self%time = 0.d0
        END IF
        self%u(:, :, :)     = 0.d0
        self%p(:, :, :)     = 0.d0
        self%incp(:, :, :)  = 0.d0

    END SUBROUTINE  zero_vel

    FUNCTION dot_vel(self, vec) RESULT(alpha)
        USE tn_axi
        IMPLICIT NONE
        CLASS(vel_field_type),      INTENT(IN)          :: self
        CLASS(vel_field_type),      INTENT(IN)          :: vec
        REAL(KIND=8) :: alpha

        alpha = 0.d0
        alpha = alpha + dot_product_SF(self%comm_one_d_ns, self%vv_mesh, self%list_mode, self%u, vec%u)

    END FUNCTION dot_vel

    SUBROUTINE axpby_vel(alpha, vec, beta, self)
        IMPLICIT NONE
        CLASS(vel_field_type),          INTENT(INOUT) :: self
        CLASS(vel_field_type),          INTENT(IN)    :: vec
        REAL(KIND=8),                   intent(IN)    :: alpha, beta

        CALL self%init_vel_field()
        self%u(:, :, :)    = beta*self%u(:,:,:) + alpha*vec%u(:,:,:)
        self%p(:, :, :)    = beta*self%p(:,:,:) + alpha*vec%p(:,:,:)
        self%incp(:, :, :) = beta*self%incp(:,:,:) + alpha*vec%incp(:,:,:)

    END SUBROUTINE axpby_vel


    SUBROUTINE rand_vel(self, ifnorm)
        IMPLICIT NONE
        CLASS(vel_field_type), intent(inout)         :: self
        logical, optional, intent(in)                :: ifnorm

        REAL(KIND=8), ALLOCATABLE, DIMENSION(:)      :: x
        LOGICAL :: normalize
        REAL(KIND=8) :: alpha
        INTEGER      :: k, i

        ! CALL self%allocate_velocity_fields()
        IF (PRESENT(ifnorm)) THEN
            normalize = ifnorm
        ELSE
            normalize = .TRUE.
        END IF
        CALL self%init_vel_field()
        ALLOCATE(x(self%vv_mesh%np))
        DO k=1,6
            DO i=1, SIZE(self%u, 3)
                IF (MOD(k,2)==0 .AND. list_mode(i)==0) THEN
                    x = 0.d0
                ELSE
                    CALL RANDOM_SEED
                    CALL RANDOM_NUMBER(x)
                ENDIF
                self%u(:, k, i) = x
            END DO
        END DO

        DEALLOCATE(x)
        self%p = 0
        self%incp = 0

        IF (normalize) THEN
            alpha = self%dot(self)
            !alpha = self%norm()
            call self%scal(1.d0/alpha)
        END IF

    END SUBROUTINE rand_vel

    SUBROUTINE scal_vel(self, alpha)
        IMPLICIT NONE
        CLASS(vel_field_type), INTENT(INOUT)       :: self
        REAL(KIND=8),          INTENT(IN)          :: alpha

        self%u    = alpha*self%u
        self%p    = alpha*self%p
        self%incp = alpha*self%incp
    END SUBROUTINE scal_vel

    FUNCTION get_size_vel(self) RESULT(N)
        IMPLICIT NONE
        CLASS(vel_field_type), INTENT(IN) :: self
        INTEGER :: N

        N = SIZE(self%u) + SIZE(self%p)*2
    END FUNCTION get_size_vel

END MODULE def_type_field_vel
