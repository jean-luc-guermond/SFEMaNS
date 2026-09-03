MODULE def_type_field_mag

   USE def_type_mesh

#include "petsc/finclude/petsc.h"
   USE petsc
   USE petscmpi
   IMPLICIT NONE
   CHARACTER(len=*), PARAMETER, PRIVATE :: this_module = 'def_type_field_mag'
    
   PUBLIC             :: build_mag_field_params, build_pointers_mag_field, build_global_pointers_mag_field
   
!========================================================
!========================================================
!========    Elements relevant to magnetic field comp ===
!========================================================
!========================================================

   TYPE(mesh_type), PRIVATE,               POINTER          :: H_mesh, phi_mesh, H_mesh_glob, phi_mesh_glob 
   MPI_Comm,        PRIVATE, DIMENSION(:), POINTER          :: communicator_mxw
   INTEGER,         PRIVATE, DIMENSION(:), POINTER          :: list_mode
   LOGICAL,         PRIVATE                                 :: if_phi_glob, if_induction_glob 
   REAL(KIND=8),    PRIVATE, DIMENSION(:), POINTER          :: mu_H_field   
   LOGICAL, PRIVATE :: if_init_pointers=.FALSE.
   
!========================================================
!========================================================
!========    definition of mag_field_type used for mhd ==
!========================================================
!========================================================

    TYPE,  PUBLIC ::  mag_field_type
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :, :)   :: B, H
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :, :)   :: phi
        REAL(KIND=8), ALLOCATABLE                       :: time 
        LOGICAL                                         :: initialized = .FALSE. 
        TYPE(mesh_type),        POINTER                 :: H_mesh, phi_mesh
        INTEGER, DIMENSION(:),  POINTER                 :: list_mode
        MPI_Comm, DIMENSION(:), POINTER                 :: comm_one_d
    CONTAINS
    !=======    Procedures for SFEMaNS
        PROCEDURE, PUBLIC, PASS(self) :: init_mag_field, allocate_induction_fields
        PROCEDURE, PUBLIC, PASS(self) :: set_time => set_time_mag
        PROCEDURE, PUBLIC, PASS(self) :: copy_to  => copy_to_mag
    !======     Procedures for LightKrylov
        PROCEDURE, PUBLIC, PASS(self) :: zero     => zero_mag
        PROCEDURE, PUBLIC, PASS(self) :: dot      => dot_mag 
        PROCEDURE, PUBLIC, PASS(self) :: axpby    => axpby_mag 
        PROCEDURE, PUBLIC, PASS(self) :: rand     => rand_mag 
        PROCEDURE, PUBLIC, PASS(self) :: scal     => scal_mag 
        PROCEDURE, PUBLIC, PASS(self) :: get_size => get_size_mag 
    END TYPE mag_field_type

CONTAINS
!========================================================
!========================================================
!========    Subroutines to build global pointers =======
!========================================================
!========================================================

    SUBROUTINE build_global_pointers_mag_field(H_mesh_in, phi_mesh_in)
        USE def_type_mesh

        IMPLICIT NONE
        TYPE(mesh_type), TARGET       ,INTENT(IN)            :: H_mesh_in, phi_mesh_in

        H_mesh_glob => H_mesh_in
        phi_mesh_glob => phi_mesh_in

    END SUBROUTINE build_global_pointers_mag_field

    SUBROUTINE build_mag_field_params(mu_H_field_in)
        REAL(KIND=8), TARGET, DIMENSION(:), INTENT(IN)       :: mu_H_field_in
        mu_H_field => mu_H_field_in

    END SUBROUTINE build_mag_field_params
        
    SUBROUTINE build_pointers_mag_field(H_mesh_in, phi_mesh_in, list_mode_in, communicator_in, if_induction_in)
        USE def_type_mesh

        IMPLICIT NONE
        TYPE(mesh_type), TARGET       ,INTENT(IN)            :: H_mesh_in, phi_mesh_in
        MPI_Comm, TARGET, DIMENSION(:),INTENT(IN)            :: communicator_in
        INTEGER, TARGET, DIMENSION(:), INTENT(IN)            :: list_mode_in
        LOGICAL,                       INTENT(IN)            :: if_induction_in

        H_mesh           => H_mesh_in
        phi_mesh         => phi_mesh_in
        communicator_mxw => communicator_in 
        list_mode        => list_mode_in
        
        if_phi_glob       = (phi_mesh%me /= 0)
        if_induction_glob = if_induction_in

    END SUBROUTINE build_pointers_mag_field

!==========================================
!==========================================
!========  Class Procedures for SFEMaNS =======
!==========================================
!==========================================

    SUBROUTINE init_mag_field(self)
        IMPLICIT NONE
        CLASS(mag_field_type),       intent(inout) :: self
        LOGICAL                                    :: test
        test = (ALLOCATED(self%H)) .AND. ALLOCATED(self%B)
        IF ((.NOT. self%initialized) .OR. (.NOT. test)) THEN
            !====== WARNING MESSAGE, could get rid of?
            IF ((.NOT. test) .AND. (self%initialized)) THEN
                WRITE(*,*) "WARNING: somehow mag_field is initialized but with non-allocated stuff"
                WRITE(*,*) 'B', ALLOCATED(self%B) 
                WRITE(*,*) 'H', ALLOCATED(self%H) 
            END IF
            !======= Building useful pointers

            !======= Zeroeing magnetic fields
            CALL self%zero()
            self%initialized = .TRUE.
        END IF
    END SUBROUTINE init_mag_field

    SUBROUTINE set_time_mag(self, time)
        USE def_type_mesh
        IMPLICIT NONE
        CLASS(mag_field_type),     INTENT(INOUT)     :: self
        REAL(KIND=8),              INTENT(IN)        :: time
        
        IF (.NOT. ALLOCATED(self%time)) ALLOCATE(self%time)
        self%time = time
    END SUBROUTINE set_time_mag

    SUBROUTINE allocate_induction_fields(self, if_glob)!, H_mesh, phi_mesh, list_mode, if_induction)!, communicator_mxw)
        USE def_type_mesh
        
        IMPLICIT NONE

        CLASS(mag_field_type),       INTENT(INOUT)           :: self
        LOGICAL, OPTIONAL,           INTENT(IN)              :: if_glob
        LOGICAL                                              :: if_glob_local, test

!          IF (.NOT. if_induction_glob)  RETURN

        if_glob_local = .FALSE.
        IF (PRESENT(if_glob)) if_glob_local = if_glob
        test = (ALLOCATED(self%H) .AND. ALLOCATED(self%B))
        IF (.NOT. test) THEN
            IF (if_glob_local) THEN
                self%H_mesh   => H_mesh_glob
                self%phi_mesh => phi_mesh_glob
            ELSE
                self%H_mesh   => H_mesh
                self%phi_mesh => phi_mesh
            END IF
            self%list_mode        => list_mode
            self%comm_one_d => communicator_mxw
            ALLOCATE(self%H   (self%H_mesh%np, 6, size(self%list_mode))) 
            ALLOCATE(self%B   (self%H_mesh%np, 6, size(self%list_mode)))
            ALLOCATE(self%phi (self%phi_mesh%np, 2, size(self%list_mode)))
        END IF
    END SUBROUTINE allocate_induction_fields


    SUBROUTINE copy_to_mag(self, vec)
        IMPLICIT NONE
        CLASS(mag_field_type),     INTENT(IN)     :: self
        CLASS(mag_field_type),     INTENT(INOUT)  :: vec

        CALL vec%allocate_induction_fields
        
        vec%H = self%H
        vec%B = self%B
        vec%phi = self%phi
        IF (.NOT. ALLOCATED(vec%time)) ALLOCATE(vec%time)
        vec%time = self%time
        vec%initialized = .TRUE.

    END SUBROUTINE copy_to_mag

!==========================================
!==========================================
!========    Procedures for LightKrylov ===
!==========================================
!==========================================

    SUBROUTINE zero_mag(self)
        IMPLICIT NONE
        CLASS(mag_field_type),       intent(inout) :: self
        CALL self%allocate_induction_fields()
        IF (.NOT. ALLOCATED(self%time)) THEN
            ALLOCATE(self%time)
            self%time = 0.d0
        END IF
        self%H(:, :, :) = 0.d0
        self%B(:, :, :) = 0.d0
        IF (if_phi_glob) THEN
            self%phi(:, :, :) = 0.d0
        END IF
    END SUBROUTINE  zero_mag

    FUNCTION dot_mag(self, vec) RESULT(alpha)
        USE tn_axi
        CLASS(mag_field_type),      INTENT(IN)         :: self
        CLASS(mag_field_type),      INTENT(IN)         :: vec
        REAL(KIND=8) :: alpha
        alpha = 0.d0
        alpha = alpha + dot_product_SF(self%comm_one_d, self%H_mesh, self%list_mode, self%H, vec%H)
            !     alpha = alpha + dot_product_SF(self%communicator_mxw, self%H_mesh, self%list_mode, self%B, vec%B)
    END FUNCTION dot_mag

    SUBROUTINE axpby_mag(alpha, vec, beta, self)
        IMPLICIT NONE
        CLASS(mag_field_type),          INTENT(INOUT) :: self
        CLASS(mag_field_type),          INTENT(IN) :: vec
        REAL(KIND=8),                   INTENT(IN)    :: alpha, beta

        CALL self%init_mag_field()
        self%H(:, :, :) = beta*self%H(:,:,:) + alpha*vec%H(:,:,:)
        self%B(:, :, :) = beta*self%B(:,:,:) + alpha*vec%B(:,:,:)
        IF (if_phi_glob) THEN
            self%phi(:, :, :) = beta*self%phi(:,:,:) + alpha*vec%phi(:,:,:)
        END IF

    END SUBROUTINE axpby_mag

    SUBROUTINE rand_mag(self, ifnorm)
        CLASS(mag_field_type), INTENT(INOUT)  :: self
        LOGICAL, OPTIONAL,     INTENT(IN)     :: ifnorm

        REAL(KIND=8), ALLOCATABLE, DIMENSION(:)      :: x
        LOGICAL :: normalize
        REAL(KIND=8) :: alpha
        INTEGER      :: k, i

        IF (PRESENT(ifnorm)) THEN
            normalize = ifnorm
        ELSE
            normalize = .TRUE.
        END IF
        
        CALL self%init_mag_field()
        ALLOCATE(x(self%H_mesh%np))
        DO k=1,6
            DO i=1, SIZE(self%B, 3)
                CALL RANDOM_SEED
                CALL RANDOM_NUMBER(x)
                self%B(:, k, i) = x
                self%H(:, k, i) = x/mu_H_field
            END DO
        END DO

        DEALLOCATE(x)
        IF (if_phi_glob) THEN
            ALLOCATE(x(phi_mesh%np))
            DO k=1,2
                DO i=1, SIZE(self%phi, 3)
                    CALL RANDOM_SEED
                    CALL RANDOM_NUMBER(x)
                    self%phi(:, k, i) = x
                END DO
            END DO
        END IF

        IF (normalize) then
            alpha = self%dot(self)
            call self%scal(1.d0/alpha)
        END IF
    END SUBROUTINE rand_mag

    SUBROUTINE scal_mag(self, alpha)
        CLASS(mag_field_type), INTENT(INOUT)       :: self
        REAL(KIND=8),          INTENT(IN)          :: alpha
        self%H = alpha*self%H
        self%B = alpha*self%B
        IF (if_phi_glob) THEN
            self%phi = alpha*self%phi
        END IF
    END SUBROUTINE scal_mag

    FUNCTION get_size_mag(self) RESULT(N)
        CLASS(mag_field_type), INTENT(IN) :: self
        INTEGER :: N
        N = SIZE(self%B)*2
        IF (if_phi_glob) N = N + SIZE(self%phi)
    END FUNCTION get_size_mag

END MODULE def_type_field_mag
