MODULE subroutine_mass
  USE my_util
  USE boundary

  PUBLIC :: three_level_mass, reconstruct_variable, total_mass, inject_P1_P2, project_P2_P1
  PRIVATE
CONTAINS

  SUBROUTINE three_level_mass(comm_one_d, time, level_set_LA_P1, level_set_LA_P2, list_mode, &
       mesh_P1, mesh_P2, chmp_vit_P2, max_vel, level_set_per, density_m2, density_m1, density, &
       level_set_m1, level_set, visc_entro_level, level_set_reg, visc_LES_level)
    USE def_type_mesh
    USE input_data
    USE subroutine_level_set
    USE my_util
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    REAL(KIND=8)                                            :: time
    INTEGER,      DIMENSION(:),       INTENT(IN)            :: list_mode
    type(petsc_csr_LA)                                      :: level_set_LA_P1
    type(petsc_csr_LA)                                      :: level_set_LA_P2
    TYPE(periodic_type),              INTENT(IN)            :: level_set_per
    TYPE(mesh_type),                  INTENT(IN)            :: mesh_P1
    TYPE(mesh_type),                  INTENT(IN)            :: mesh_P2
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)            :: chmp_vit_P2
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(INOUT)         :: density, density_m1, density_m2
    REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(INOUT)         :: level_set, level_set_m1
    REAL(KIND=8),                     INTENT(INOUT)         :: max_vel
    REAL(KIND=8), DIMENSION(:,:),     INTENT(IN)            :: visc_entro_level
    REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(OUT)           :: level_set_reg
    REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(IN)            :: visc_LES_level
    REAL(KIND=8), DIMENSION(mesh_P1%np,6,SIZE(list_mode))   :: chmp_vit_P1
    INTEGER                                                 :: n, i, k
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d

    IF (inputs%if_level_set) THEN

       IF (inputs%if_level_set_P2) THEN
          IF (inputs%if_level_set_fixed) THEN
             DO n = 1, inputs%nb_fluid-1
                DO k = 1, 2
                   DO i = 1, SIZE(list_mode)
                      level_set_m1(n,:,k,i)=level_set_exact(n,k,mesh_P2%rr,list_mode(i),time-inputs%dt)
                      level_set(n,:,k,i)=level_set_exact(n,k,mesh_P2%rr,list_mode(i),time)
                   END DO
                END DO
             END DO
          ELSE
             DO n = 1, inputs%nb_fluid-1
                CALL three_level_level_set(comm_one_d, time, level_set_LA_P2, inputs%dt, list_mode, &
                     mesh_P2, level_set_m1(n,:,:,:), level_set(n,:,:,:), chmp_vit_P2, max_vel, &
                     inputs%my_par_level_set, inputs%level_set_list_dirichlet_sides, level_set_per, n, &
                     visc_entro_level, level_set_reg(n,:,:,:), visc_LES_level(n,:,:,:))
             END DO
          END IF
       ELSE
          DO i = 1, SIZE(list_mode)
             DO k = 1, 6
                CALL project_P2_P1(mesh_P2%jj, mesh_P1%jj, chmp_vit_P2(:,k,i), chmp_vit_P1(:,k,i))
             END DO
          END DO
          IF (inputs%if_level_set_fixed) THEN
             DO n = 1, inputs%nb_fluid-1
                DO k = 1, 2
                   DO i = 1, SIZE(list_mode)
                      level_set_m1(n,:,k,i)=level_set_exact(n,k,mesh_P1%rr,list_mode(i),time-inputs%dt)
                      level_set(n,:,k,i)=level_set_exact(n,k,mesh_P1%rr,list_mode(i),time)
                   END DO
                END DO
             END DO
          ELSE
             DO n = 1, inputs%nb_fluid-1
                CALL three_level_level_set(comm_one_d, time, level_set_LA_P1, inputs%dt, list_mode, &
                     mesh_P1, level_set_m1(n,:,:,:), level_set(n,:,:,:), chmp_vit_P1, max_vel, &
                     inputs%my_par_level_set, inputs%level_set_list_dirichlet_sides, level_set_per, n, &
                     visc_entro_level, level_set_reg(n,:,:,:), visc_LES_level(n,:,:,:))
             END DO
          END IF
       END IF

       !===Update densities
       density_m2 = density_m1
       density_m1 = density
       CALL reconstruct_variable(comm_one_d, list_mode, mesh_P1, mesh_P2, level_set, &
            inputs%density_fluid, density)
       RETURN
    ELSE
       RETURN
    END IF

  END SUBROUTINE three_level_mass

  SUBROUTINE reconstruct_variable(comm_one_d, list_mode, mesh_P1, mesh_P2, level_set, values, variable)
    !==============================
    USE def_type_mesh
    USE sft_parallele
    USE input_data
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    INTEGER,      DIMENSION(:),                 INTENT(IN)    :: list_mode
    TYPE(mesh_type),                            INTENT(IN)    :: mesh_P1
    TYPE(mesh_type),                            INTENT(IN)    :: mesh_P2
    REAL(KIND=8), DIMENSION(:,:,:,:),           INTENT(IN)    :: level_set
    REAL(KIND=8), DIMENSION(:),                 INTENT(IN)    :: values
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(INOUT) :: variable
    LOGICAL,                                             SAVE :: once = .TRUE.
    INTEGER,                                             SAVE :: m_max_c
    INTEGER,                                             SAVE :: m_max_pad
    INTEGER,                                             SAVE :: bloc_size
    INTEGER,                                             SAVE :: nb_procs
    INTEGER                                                   :: i, code, k, nb_inter
    REAL(KIND=8), DIMENSION(mesh_P2%np,2,SIZE(list_mode))     :: rho_phi
    REAL(KIND=8), DIMENSION(inputs%nb_fluid-1,mesh_P2%np,2,SIZE(list_mode))   :: level_set_P2
    !Communicators for Petsc, in space and Fourier------------------------------
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d
    !------------------------------END OF DECLARATION--------------------------------------

    IF (once) THEN
       once = .FALSE.
       !-------------DIMENSIONS-------------------------------------------------------
       m_max_c = SIZE(list_mode)
       !-------------FFT VARIABLES----------------------------------------------------
       CALL MPI_COMM_SIZE(comm_one_d(2), nb_procs, code)
       bloc_size = mesh_P2%np/nb_procs+1
       m_max_pad = 3*m_max_c*nb_procs/2
    END IF

    IF (.NOT.inputs%if_level_set) THEN
       ! we don't reconstruct density, viscosity if no level set
       RETURN
    ELSE
       IF (inputs%if_level_set_P2) THEN
          level_set_P2=level_set
       ELSE
          DO nb_inter = 1, inputs%nb_fluid-1
             DO i = 1, SIZE(list_mode)
                DO k = 1, 2
                   CALL inject_P1_P2(mesh_P1%jj, mesh_P2%jj, level_set(nb_inter,:,k,i), level_set_P2(nb_inter,:,k,i))
                END DO
             END DO
          END DO
       END IF
       IF (MAXVAL(ABS(values(1)-values(:))) .LE. 1.d-10*MAXVAL(ABS(values(:)))) THEN
          variable = 0.d0
          DO i = 1, m_max_c
             IF (list_mode(i)==0) THEN
                variable(:,1,i) = values(1)
             END IF
          END DO
       ELSE IF (inputs%level_set_reconstruction_type == 'lin') THEN
          IF (inputs%if_kill_overshoot) THEN
             IF (nb_procs==1.AND.SIZE(list_mode)==1.AND.list_mode(1)==0) THEN !case axisym
                level_set_P2 = MIN(1.d0, level_set_P2)
                level_set_P2 = MAX(0.d0, level_set_P2)
             ELSE !level set depends of theta
                DO k = 1, inputs%nb_fluid-1
                   CALL FFT_NO_OVERSHOOT_LEVEL_SET(comm_one_d(2), level_set_P2(k,:,:,:), &
                        nb_procs, bloc_size, m_max_pad)
                END DO
             END IF
          END IF
          variable = 0.d0
          DO i = 1, m_max_c
             IF (list_mode(i)==0) THEN
                variable(:,1,i) = values(1)
             END IF
          END DO
          variable = variable + (values(2)-values(1))*level_set_P2(1,:,:,:)
          IF (inputs%nb_fluid.GE.3) THEN
             DO i = 1, inputs%nb_fluid-2
                CALL FFT_PAR_PROD_DCL(comm_one_d(2), variable, level_set_P2(i+1,:,:,:), rho_phi, &
                     nb_procs, bloc_size, m_max_pad)
                variable = variable -rho_phi + values(i+2)*level_set_P2(i+1,:,:,:)
             END DO
          END IF
       ELSE
          CALL FFT_HEAVISIDE_DCL(comm_one_d(2), level_set_P2, values, variable, &
               nb_procs, bloc_size, m_max_pad, inputs%level_set_tanh_coeff_reconstruction)
          DO i = 1, m_max_c
             IF (list_mode(i)==0) THEN
                variable(:,2,i) = 0.d0
             END IF
          END DO
       END IF
    END IF
  END SUBROUTINE reconstruct_variable

  SUBROUTINE total_mass(comm_one_d, list_mode, mass_mesh, level_set, mass_tot)
    !==============================
    USE def_type_mesh
    USE sft_parallele
    USE input_data
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    INTEGER,      DIMENSION(:),                 INTENT(IN)    :: list_mode
    TYPE(mesh_type),                            INTENT(IN)    :: mass_mesh
    REAL(KIND=8), DIMENSION(:,:,:,:),           INTENT(IN)    :: level_set
    REAL(KIND=8),                               INTENT(OUT)   :: mass_tot
    REAL(KIND=8), DIMENSION(SIZE(level_set,2),SIZE(level_set,3), &
         SIZE(level_set,4))                   :: density_loc
    INTEGER                                                   :: m_max_pad, bloc_size,  nb_procs
    INTEGER                                                   :: i, code, my_petscworld_rank, m, l
    REAL(KIND=8)                                              :: mass_loc, mass_F, ray
    REAL(KIND=8), DIMENSION(mass_mesh%np,2,SIZE(list_mode))   :: rho_phi
    INTEGER, DIMENSION(mass_mesh%gauss%n_w)                   :: j_loc
    REAL(KIND=8)                                              :: pi= 3.14159265358979323846d0
    !Communicators for Petsc, in space and Fourier------------------------------
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d
    !------------------------------END OF DECLARATION--------------------------------------

    CALL MPI_COMM_RANK(PETSC_COMM_WORLD,my_petscworld_rank,code)
    CALL MPI_COMM_SIZE(comm_one_d(2), nb_procs, code)

    IF(.NOT.inputs%if_level_set) THEN
       CALL error_petsc('BUG in sub_mass : you should not compute any mass')
    ELSE
       IF (inputs%level_set_reconstruction_type == 'lin') THEN
          density_loc = 0.d0
          DO i = 1, SIZE(list_mode)
             IF (list_mode(i)==0) THEN
                density_loc(:,1,i) = inputs%density_fluid(1)
             END IF
          END DO
          density_loc = density_loc + (inputs%density_fluid(2)-inputs%density_fluid(1))*level_set(i,:,:,:)

          bloc_size = SIZE(level_set,2)/nb_procs+1
          m_max_pad = 3*SIZE(list_mode)*nb_procs/2
          IF (inputs%nb_fluid.GE.3) THEN
             DO i = 1, inputs%nb_fluid-2
                CALL FFT_PAR_PROD_DCL(comm_one_d(2), density_loc, level_set(i+1,:,:,:), rho_phi, &
                     nb_procs, bloc_size, m_max_pad)
                density_loc = density_loc -rho_phi + inputs%density_fluid(i+2)*level_set(i+1,:,:,:)
             END DO
          END IF
       ELSE
          bloc_size = SIZE(level_set,2)/nb_procs+1
          m_max_pad = 3*SIZE(list_mode)*nb_procs/2
          CALL FFT_HEAVISIDE_DCL(comm_one_d(2), level_set, inputs%density_fluid, &
               density_loc, nb_procs, bloc_size, m_max_pad, inputs%level_set_tanh_coeff_reconstruction)
       END IF

       mass_loc = 0.d0
       DO i = 1, SIZE(list_mode)
          IF (list_mode(i)==0) THEN
             DO m = 1, mass_mesh%me
                j_loc = mass_mesh%jj(:,m)
                DO l = 1, mass_mesh%gauss%l_G
                   !===Compute radius of Gauss point
                   ray = SUM(mass_mesh%rr(1,j_loc)*mass_mesh%gauss%ww(:,l))
                   mass_loc = mass_loc + SUM(density_loc(j_loc,1,i)* &
                        mass_mesh%gauss%ww(:,l))*ray*mass_mesh%gauss%rj(l,m)
                END DO
             END DO
          END IF
       END DO
       mass_loc = mass_loc*2*pi
       CALL  MPI_ALLREDUCE(mass_loc, mass_F, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
            comm_one_d(2), code)
       CALL MPI_ALLREDUCE(mass_F, mass_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
            comm_one_d(1), code)
    END IF

  END SUBROUTINE total_mass

  SUBROUTINE inject_P1_P2(jj_c, jj_f, pp_c, pp_f)
    USE my_util
    IMPLICIT NONE
    INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj_c, jj_f
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: pp_c
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: pp_f
    REAL(KIND=8) :: half = 0.5
    INTEGER:: m
    IF (SIZE(jj_c,1)==3) THEN
       DO m = 1, SIZE(jj_f,2)
          pp_f(jj_f(1:3,m)) =  pp_c(jj_c(:,m))
          pp_f(jj_f(4,m)) = (pp_c(jj_c(2,m)) + pp_c(jj_c(3,m)))*half
          pp_f(jj_f(5,m)) = (pp_c(jj_c(3,m)) + pp_c(jj_c(1,m)))*half
          pp_f(jj_f(6,m)) = (pp_c(jj_c(1,m)) + pp_c(jj_c(2,m)))*half
       END DO
    ELSE
       CALL error_petsc('BUG in inject_P1_P2: finite element not yet programmed')
    END IF

  END SUBROUTINE inject_P1_P2

  SUBROUTINE project_P2_P1(jj_P2, jj_P1, pp_P2, pp_P1)
    USE my_util
    IMPLICIT NONE
    INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj_P2, jj_P1
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: pp_P2
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: pp_P1
    INTEGER:: m

    IF (SIZE(jj_P1,1)==3) THEN
       DO m = 1, SIZE(jj_P1,2)
          pp_P1(jj_P1(:,m)) =  pp_P2(jj_P2(1:3,m))
       END DO
    ELSE
       CALL error_petsc('BUG in inject_P2_P1: finite element not yet programmed')
    END IF

  END SUBROUTINE project_P2_P1

END MODULE subroutine_mass
