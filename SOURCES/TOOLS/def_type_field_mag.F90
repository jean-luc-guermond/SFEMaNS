MODULE def_type_field_mag

   USE def_type_mesh

#ifdef USE_LIGHTKRYLOV
   USE LightKrylov
   USE LightKrylov_Logger
   USE stdlib_math, ONLY: linspace
   USE stdlib_optval, ONLY: optval
#endif

#include "petsc/finclude/petsc.h"
   USE petsc
   USE petscmpi
   IMPLICIT NONE
   CHARACTER(len=*), PARAMETER, PRIVATE :: this_module = 'def_type_field_mag'
    
   PUBLIC             :: build_mag_field_params, build_pointers_mag_field, build_global_pointers_mag_field
   
   TYPE(mesh_type), PRIVATE,               POINTER          :: H_mesh, phi_mesh, H_mesh_glob, phi_mesh_glob 
   MPI_Comm,        PRIVATE, DIMENSION(:), POINTER          :: communicator_mxw
   INTEGER,         PRIVATE, DIMENSION(:), POINTER          :: list_mode
   LOGICAL,         PRIVATE                                 :: if_phi_glob, if_induction_glob 
   REAL(KIND=8),    PRIVATE, DIMENSION(:), POINTER          :: mu_H_field   
   LOGICAL, PRIVATE :: if_init_pointers=.FALSE.
   
#ifdef USE_LIGHTKRYLOV
   TYPE, EXTENDS(abstract_vector_rdp), PUBLIC ::  mag_field_type
#else
   TYPE,                               PUBLIC ::  mag_field_type
#endif
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :, :)   :: B, H
        ! REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :, :)   :: Bn, Bn1, Hn, Hn1
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :, :)   :: phi
        ! REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :, :)   :: phin, phin1
        REAL(KIND=8), ALLOCATABLE                       :: time 
        LOGICAL                                         :: initialized = .FALSE. 
        TYPE(mesh_type),        POINTER                 :: H_mesh, phi_mesh
        INTEGER, DIMENSION(:),  POINTER                 :: list_mode
        MPI_Comm, DIMENSION(:), POINTER                 :: comm_one_d
   CONTAINS
   !=======    Procedures for SFEMaNS
        PROCEDURE, PUBLIC, PASS(self) :: init_mag_field, allocate_induction_fields
        PROCEDURE, PUBLIC, PASS(self) :: set_time
        PROCEDURE, PUBLIC, PASS(self) :: copy_to
   !======     Procedures for LightKrylov
        PROCEDURE, PUBLIC, PASS(self) :: zero 
        PROCEDURE, PUBLIC, PASS(self) :: dot 
        PROCEDURE, PUBLIC, PASS(self) :: axpby 
        PROCEDURE, PUBLIC, PASS(self) :: rand 
        PROCEDURE, PUBLIC, PASS(self) :: scal 
        PROCEDURE, PUBLIC, PASS(self) :: get_size 
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

        SUBROUTINE set_time(self, time)
            USE def_type_mesh
            IMPLICIT NONE
            CLASS(mag_field_type),     INTENT(INOUT)     :: self
            REAL(KIND=8),              INTENT(IN)        :: time
            
            IF (.NOT. ALLOCATED(self%time)) ALLOCATE(self%time)
            self%time = time
        END SUBROUTINE set_time

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


      SUBROUTINE copy_to(self, vec)
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

      END SUBROUTINE copy_to

!==========================================
!==========================================
!========    Procedures for LightKrylov ===
!==========================================
!==========================================

   SUBROUTINE zero(self)
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
   END SUBROUTINE  zero

   REAL(KIND=8) function dot(self, vec) result(alpha)
      USE tn_axi
      CLASS(mag_field_type),      intent(in)         :: self
#ifdef USE_LIGHTKRYLOV
      CLASS(abstract_vector_rdp), intent(in)         :: vec
#else
      CLASS(mag_field_type),      intent(in)         :: vec
#endif
#ifdef USE_LIGHTKRYLOV
      SELECT TYPE(vec)
      type is(mag_field_type)
#endif
              alpha = 0.d0
              alpha = alpha + dot_product_SF(self%comm_one_d, self%H_mesh, self%list_mode, self%H, vec%H)
         !     alpha = alpha + dot_product_SF(self%communicator_mxw, self%H_mesh, self%list_mode, self%B, vec%B)
#ifdef USE_LIGHTKRYLOV
          class default
             CALL type_error('vec', 'mag_field_type', 'IN', this_module, 'dot')
      END SELECT 
#endif
   END function dot

   SUBROUTINE axpby(alpha, vec, beta, self)
      IMPLICIT NONE
      CLASS(mag_field_type),          intent(inout) :: self
#ifdef USE_LIGHTKRYLOV
      CLASS(abstract_vector_rdp),     intent(in)    :: vec
#else
      CLASS(mag_field_type),          intent(inout) :: vec
#endif
      REAL(KIND=8),                   intent(in)    :: alpha, beta

#ifdef USE_LIGHTKRYLOV
      SELECT TYPE(vec)
      type is(mag_field_type)
#endif
         CALL self%init_mag_field()
          self%H(:, :, :) = beta*self%H(:,:,:) + alpha*vec%H(:,:,:)
          self%B(:, :, :) = beta*self%B(:,:,:) + alpha*vec%B(:,:,:)
          IF (if_phi_glob) THEN
              self%phi(:, :, :) = beta*self%phi(:,:,:) + alpha*vec%phi(:,:,:)
          END IF
#ifdef USE_LIGHTKRYLOV
      CLASS default
          call type_error('vec', 'state_vector', 'IN', this_module, 'axpby')
      END SELECT
#endif

   END SUBROUTINE axpby


   SUBROUTINE rand(self, ifnorm)
      CLASS(mag_field_type), intent(inout)         :: self
      logical, optional, intent(in)                :: ifnorm

      REAL(KIND=8), ALLOCATABLE, DIMENSION(:)      :: x
      LOGICAL :: normalize
      REAL(KIND=8) :: alpha
      INTEGER      :: k, i

#ifdef USE_LIGHTKRYLOV
      normalize = optval(ifnorm, .true.)
#else
      IF (PRESENT(ifnorm)) THEN
          normalize = ifnorm
      ELSE
          normalize = .TRUE.
      END IF
#endif
      CALL self%init_mag_field()
      ALLOCATE(x(self%H_mesh%np))
    !   ALLOCATE(x(H_mesh%np))
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

      if (normalize) then
         alpha = self%dot(self)
         !alpha = self%norm()
         call self%scal(1.d0/alpha)
      END if
   END SUBROUTINE rand

   SUBROUTINE scal(self, alpha)
      CLASS(mag_field_type), INTENT(INOUT)       :: self
      REAL(KIND=8),          INTENT(IN)          :: alpha
      self%H = alpha*self%H
      self%B = alpha*self%B
      IF (if_phi_glob) THEN
          self%phi = alpha*self%phi
      END IF
    END SUBROUTINE

   integer FUNCTION get_size(self) result(N)
      CLASS(mag_field_type), INTENT(IN) :: self
      N = SIZE(self%B)*2
      IF (if_phi_glob) N = N + SIZE(self%phi)
   END FUNCTION get_size

END MODULE def_type_field_mag
