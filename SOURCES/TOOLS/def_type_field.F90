MODULE def_type_field

   USE def_type_mesh
   ! LightKrylov for linear algebra.
   !USE LightKrylov, ONLY: abstract_vector_rdp

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
    CHARACTER(len=*), PARAMETER, PRIVATE :: this_module = 'def_type_field'
    REAL(KIND=8),                PRIVATE :: code
    
   PUBLIC             :: build_mag_field_params, build_pointers_mag_field, build_global_pointers_mag_field
   
   TYPE(mesh_type), PRIVATE,               POINTER          :: H_mesh, phi_mesh, H_mesh_glob, phi_mesh_glob 
   MPI_Comm,        PRIVATE, DIMENSION(:), POINTER          :: communicator_mxw
   INTEGER,         PRIVATE, DIMENSION(:), POINTER          :: list_mode
   LOGICAL,         PRIVATE                                 :: if_phi_glob, if_induction_glob 
   REAL(KIND=8),    PRIVATE, DIMENSION(:), POINTER          :: mu_H_field   
   
#ifdef USE_LIGHTKRYLOV
   TYPE, EXTENDS(abstract_vector_rdp), PUBLIC ::  mag_field_type
#else
   TYPE,                               PUBLIC ::  mag_field_type
#endif
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :, :)   :: Bn, Bn1, Hn, Hn1
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :, :)   :: phin, phin1
        REAL(KIND=8), ALLOCATABLE                       :: time 
        LOGICAL                                         :: initialized = .FALSE. 
   CONTAINS
   !=======    Procedures for SFEMaNS
        PROCEDURE, PUBLIC, PASS(self) :: init_mag_field, allocate_induction_fields
        PROCEDURE, PUBLIC, PASS(self) :: set_time
   !======     Procedures for LightKrylov
        PROCEDURE, PUBLIC, PASS(self) :: zero 
        PROCEDURE, PUBLIC, PASS(self) :: dot 
        PROCEDURE, PUBLIC, PASS(self) :: axpby 
        PROCEDURE, PUBLIC, PASS(self) :: rand 
        PROCEDURE, PUBLIC, PASS(self) :: scal 
        PROCEDURE, PUBLIC, PASS(self) :: get_size 
   END TYPE mag_field_type

   CONTAINS
!==========================================
!==========================================
!========    Procedures for SFEMaNS =======
!==========================================
!==========================================

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
          H_mesh => H_mesh_in
          phi_mesh => phi_mesh_in
          communicator_mxw => communicator_in 
          list_mode => list_mode_in
          if_phi_glob = (phi_mesh%me /= 0)
          if_induction_glob = if_induction_in

      END SUBROUTINE build_pointers_mag_field

      SUBROUTINE allocate_induction_fields(self, if_glob)!, H_mesh, phi_mesh, list_mode, if_induction)!, communicator_mxw)
          USE def_type_mesh
         
          IMPLICIT NONE

          CLASS(mag_field_type),       INTENT(INOUT)           :: self
          LOGICAL, OPTIONAL,           INTENT(IN)              :: if_glob
          LOGICAL                                              :: if_glob_local, test

!          IF (.NOT. if_induction_glob)  RETURN

          if_glob_local = .FALSE.
          IF (PRESENT(if_glob)) if_glob_local = if_glob
          test = (ALLOCATED(self%Hn) .AND. ALLOCATED(self%Hn1) .AND. ALLOCATED(self%Bn1) .AND. ALLOCATED(self%Bn))
          IF (.NOT. test) THEN
              IF (if_glob_local) THEN
                  ALLOCATE(self%Hn1  (H_mesh_glob%np, 6, size(list_mode))) 
                  ALLOCATE(self%Hn   (H_mesh_glob%np, 6, size(list_mode))) 
                  ALLOCATE(self%Bn1  (H_mesh_glob%np, 6, size(list_mode)))
                  ALLOCATE(self%Bn   (H_mesh_glob%np, 6, size(list_mode)))
                  ALLOCATE(self%phin1(phi_mesh_glob%np, 2, size(list_mode)))
                  ALLOCATE(self%phin (phi_mesh_glob%np, 2, size(list_mode)))
              ELSE
                  ALLOCATE(self%Hn1  (H_mesh%np, 6, size(list_mode))) 
                  ALLOCATE(self%Hn   (H_mesh%np, 6, size(list_mode))) 
                  ALLOCATE(self%Bn1  (H_mesh%np, 6, size(list_mode)))
                  ALLOCATE(self%Bn   (H_mesh%np, 6, size(list_mode)))
                  ALLOCATE(self%phin1(phi_mesh%np, 2, size(list_mode)))
                  ALLOCATE(self%phin (phi_mesh%np, 2, size(list_mode)))
              END IF
          END IF
      END SUBROUTINE allocate_induction_fields


      SUBROUTINE set_time(self, time)
          USE def_type_mesh
          IMPLICIT NONE
          CLASS(mag_field_type),     INTENT(INOUT)     :: self
          REAL(KIND=8),              INTENT(IN)        :: time
          
          IF (.NOT. ALLOCATED(self%time)) ALLOCATE(self%time)
          self%time = time
      END SUBROUTINE set_time
!==========================================
!==========================================
!========    Procedures for LightKrylov ===
!==========================================
!==========================================

   SUBROUTINE init_mag_field(self)
      IMPLICIT NONE
      CLASS(mag_field_type),       intent(inout) :: self
      LOGICAL                                    :: test
      test = (ALLOCATED(self%Hn) .AND. ALLOCATED(self%Hn1) .AND. ALLOCATED(self%Bn1) .AND. ALLOCATED(self%Bn))
      IF ((.NOT. self%initialized) .OR. (.NOT. test)) THEN
         IF ((.NOT. test) .AND. (self%initialized)) THEN
             WRITE(*,*) "WARNING: somehow mag_field is initialized but with non-allocated stuff"
             WRITE(*,*) 'Bn', ALLOCATED(self%Bn) 
             WRITE(*,*) 'Bn1', ALLOCATED(self%Bn1) 
             WRITE(*,*) 'Hn', ALLOCATED(self%Hn) 
             WRITE(*,*) 'Hn1', ALLOCATED(self%Hn1) 
         END IF
         CALL self%zero()
         self%initialized = .TRUE.
      END IF
   END SUBROUTINE init_mag_field

   SUBROUTINE zero(self)
      IMPLICIT NONE
      CLASS(mag_field_type),       intent(inout) :: self
      CALL self%allocate_induction_fields()
      IF (.NOT. ALLOCATED(self%time)) THEN
          ALLOCATE(self%time)
          self%time = 0.d0
      END IF
      self%Hn1(:, :, :) = 0.d0
      self%Hn(:, :, :) = 0.d0
      self%Bn1(:, :, :) = 0.d0
      self%Bn(:, :, :) = 0.d0
      IF (if_phi_glob) THEN
          self%phin1(:, :, :) = 0.d0
          self%phin(:, :, :) = 0.d0
      END IF
   END SUBROUTINE  zero

   REAL(KIND=8) function dot(self, vec) result(alpha)
      USE tn_axi
      CLASS(mag_field_type),      intent(in)          :: self
#ifdef USE_LIGHTKRYLOV
      CLASS(abstract_vector_rdp), intent(in)          :: vec
#else
      CLASS(mag_field_type),      intent(in)          :: vec
#endif
      REAL(KIND=8)                                   :: alpha_loc, factor_mode, H_r_dr_dth_dz           
      INTEGER                                        :: i, index, m, l, type_vec
      INTEGER, DIMENSION(H_mesh%gauss%n_w) :: j_loc
#ifdef USE_LIGHTKRYLOV
      SELECT TYPE(vec)
      type is(mag_field_type)
#endif
              alpha = 0.d0
              alpha = alpha + dot_product_SF(communicator_mxw, H_mesh, list_mode, self%Hn, vec%Hn)
         !     alpha = alpha + dot_product_SF(communicator_mxw, H_mesh, list_mode, self%Bn, vec%Bn)
         !     alpha = alpha + dot_product_SF(communicator_mxw, H_mesh, list_mode, self%Hn1, vec%Hn1)
         !     alpha = alpha + dot_product_SF(communicator_mxw, H_mesh, list_mode, self%Bn1, vec%Bn1)
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
          self%Hn(:, :, :) = beta*self%Hn(:,:,:) + alpha*vec%Hn(:,:,:)
          self%Hn1(:, :, :) = beta*self%Hn1(:,:,:) + alpha*vec%Hn1(:,:,:)
          self%Bn(:, :, :) = beta*self%Bn(:,:,:) + alpha*vec%Bn(:,:,:)
          self%Bn1(:, :, :) = beta*self%Bn1(:,:,:) + alpha*vec%Bn1(:,:,:)
          IF (if_phi_glob) THEN
              self%phin(:, :, :) = beta*self%phin(:,:,:) + alpha*vec%phin(:,:,:)
              self%phin1(:, :, :) = beta*self%phin1(:,:,:) + alpha*vec%phin1(:,:,:)
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
      ALLOCATE(x(H_mesh%np))
      DO k=1,6
          DO i=1, SIZE(self%Bn, 3)
              CALL RANDOM_SEED
              CALL RANDOM_NUMBER(x)
              self%Bn1(:, k, i) = x
              self%Hn1(:, k, i) = x/mu_H_field
              self%Bn(:, k, i) = x
              self%Hn(:, k, i) = x/mu_H_field
          END DO
      END DO

      DEALLOCATE(x)
      IF (if_phi_glob) THEN
          ALLOCATE(x(phi_mesh%np))
          DO k=1,2
              DO i=1, SIZE(self%phin, 3)
                  CALL RANDOM_SEED
                  CALL RANDOM_NUMBER(x)
                  self%phin(:, k, i) = x
                  self%phin1(:, k, i) = x
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
      self%Hn1 = alpha*self%Hn1
      self%Hn = alpha*self%Hn
      self%Bn1 = alpha*self%Bn1
      self%Bn = alpha*self%Bn
      IF (if_phi_glob) THEN
          self%phin = alpha*self%phin
      self%phin1 = alpha*self%phin1
      END IF
    END SUBROUTINE

   integer FUNCTION get_size(self) result(N)
      CLASS(mag_field_type), INTENT(IN) :: self
      N = SIZE(self%Bn)*4
      IF (if_phi_glob) N = N + SIZE(self%phin)*2
   END FUNCTION get_size

END MODULE def_type_field
