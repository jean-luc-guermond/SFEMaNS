MODULE boundary_test_33
  USE def_type_mesh
  USE input_data
  USE my_util
!!$ATTENTION
!!$Some subroutines have been commented to avoid warning messages when compiling executable.
!!$It can not be done in the module boundary_generic that expects all subroutines to be present.
!!$END ATTENTION
  PUBLIC :: init_velocity_pressure
  PUBLIC :: init_temperature
!!$  PUBLIC :: init_level_set
  PUBLIC :: source_in_NS_momentum
  PUBLIC :: source_in_temperature
!!$  PUBLIC :: source_in_level_set
  PUBLIC :: vv_exact
!!$  PUBLIC :: imposed_velocity_by_penalty
  PUBLIC :: pp_exact
  PUBLIC :: temperature_exact
!!$  PUBLIC :: level_set_exact
!!$  PUBLIC :: penal_in_real_space
  PUBLIC :: extension_velocity
!!$  PUBLIC :: Vexact
!!$  PUBLIC :: H_B_quasi_static
  PUBLIC :: Hexact
  PUBLIC :: Phiexact
  PUBLIC :: Jexact_gauss
  PUBLIC :: Eexact_gauss
  PUBLIC :: init_maxwell
!!$  PUBLIC :: mu_bar_in_fourier_space
!!$  PUBLIC :: grad_mu_bar_in_fourier_space
!!$  PUBLIC :: mu_in_real_space
  PUBLIC :: chi_coeff_law
  PUBLIC :: T_dchi_dT_coeff_law
  PRIVATE

CONTAINS
  !===============================================================================
  !                       Boundary conditions for Navier-Stokes
  !===============================================================================

  !===Initialize velocity, pressure
  SUBROUTINE init_velocity_pressure(mesh_f, mesh_c, time, dt, list_mode, &
       un_m1, un, pn_m1, pn, phin_m1, phin)
    IMPLICIT NONE
    TYPE(mesh_type)                            :: mesh_f, mesh_c
    REAL(KIND=8),                   INTENT(OUT):: time
    REAL(KIND=8),                   INTENT(IN) :: dt
    INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: un_m1, un
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: pn_m1, pn, phin_m1, phin
    INTEGER                                    :: mode, i, j 
    REAL(KIND=8), DIMENSION(mesh_c%np)         :: pn_m2

    time = 0.d0
    DO i= 1, SIZE(list_mode)
       mode = list_mode(i) 
       DO j = 1, 6 
          !===velocity
          un_m1(:,j,i) = vv_exact(j,mesh_f%rr,mode,time-dt)  
          un   (:,j,i) = vv_exact(j,mesh_f%rr,mode,time)
       END DO
       DO j = 1, 2
          !===pressure
          pn_m2(:)       = pp_exact(j,mesh_c%rr,mode,time-2*dt)
          pn_m1  (:,j,i) = pp_exact(j,mesh_c%rr,mode,time-dt)
          pn     (:,j,i) = pp_exact(j,mesh_c%rr,mode,time)
          phin_m1(:,j,i) = pn_m1(:,j,i) - pn_m2(:)
          phin   (:,j,i) = Pn   (:,j,i) - pn_m1(:,j,i)
       ENDDO
    ENDDO
  END SUBROUTINE init_velocity_pressure

  !===Initialize temperature
  SUBROUTINE init_temperature(mesh, time, dt, list_mode, tempn_m1, tempn)
    IMPLICIT NONE
    TYPE(mesh_type)                            :: mesh
    REAL(KIND=8),                   INTENT(OUT):: time
    REAL(KIND=8),                   INTENT(IN) :: dt
    INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: tempn_m1, tempn
    INTEGER                                    :: mode, i, j 

    time = 0.d0
    DO i= 1, SIZE(list_mode)
       mode = list_mode(i) 
       DO j = 1, 2 
          tempn_m1(:,j,i) = temperature_exact(j, mesh%rr, mode, time-dt)
          tempn   (:,j,i) = temperature_exact(j, mesh%rr, mode, time)
       ENDDO
    ENDDO
  END SUBROUTINE init_temperature

!!$  !===Initialize level_set
!!$  SUBROUTINE init_level_set(vv_mesh, time, &
!!$       dt, list_mode, level_set_m1, level_set)
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type)                              :: vv_mesh 
!!$    REAL(KIND=8),                     INTENT(OUT):: time
!!$    REAL(KIND=8),                     INTENT(IN) :: dt
!!$    INTEGER,      DIMENSION(:),       INTENT(IN) :: list_mode
!!$    REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(OUT):: level_set, level_set_m1
!!$    INTEGER                                      :: mode, i, j, n 
!!$    time = 0.d0
!!$    DO i= 1, SIZE(list_mode)
!!$       mode = list_mode(i) 
!!$       DO j = 1, 2
!!$          !===level_set
!!$          DO n = 1, inputs%nb_fluid -1
!!$             level_set_m1(n,:,j,i) = level_set_exact(n,j,vv_mesh%rr,mode,time-dt)  
!!$             level_set   (n,:,j,i) = level_set_exact(n,j,vv_mesh%rr,mode,time)
!!$          END DO
!!$       END DO
!!$    END DO
!!$  END SUBROUTINE init_level_set

  !===Source in momemtum equation. Always called.
  FUNCTION source_in_NS_momentum(TYPE, rr, mode, i, time, Re, ty, opt_density, opt_tempn) RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                             INTENT(IN) :: TYPE
    REAL(KIND=8), DIMENSION(:,:),             INTENT(IN) :: rr
    INTEGER     ,                             INTENT(IN) :: mode, i
    REAL(KIND=8),                             INTENT(IN) :: time   
    REAL(KIND=8),                             INTENT(IN) :: Re
    CHARACTER(LEN=2),                         INTENT(IN) :: ty
    REAL(KIND=8), DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: opt_density
    REAL(KIND=8), DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: opt_tempn 
    REAL(KIND=8), DIMENSION(SIZE(rr,2))                  :: vv, r, z
    REAL(KIND=8)                                         :: alpha, r0 = 0.5d0
    CHARACTER(LEN=2)  :: np

    IF (PRESENT(opt_density)) CALL error_petsc('density should not be present for test 33') 

    alpha = inputs%gravity_coefficient
    r = rr(1,:)
    z = rr(2,:)

    IF (TYPE==5) THEN 
       vv = alpha*(opt_tempn(:,1,i) - temperature_exact(1,rr,mode,time))
    ELSE IF (TYPE==6) THEN
       vv = alpha*(opt_tempn(:,2,i) - temperature_exact(2,rr,mode,time))
    ELSE 
       vv = 0.d0
    END IF

    IF (TYPE==1) THEN
       IF (mode==0) THEN
          vv = vv -(2*r*(r**4 - 2*r**3*r0 + r0**2 + r**2*(-3 + r0**2))*Cos(time)*Cos(z) + &
               (r - r0)*Re*Cos(time)**2*(14*r**3 - 5*r**2*r0 - 5*r*r0**2 + 2*r0**3 + &
               (9*r**5 - 21*r**4*r0 + 5*r*r0**2 - 2*r0**3 + r**3*(-14 + 15*r0**2) + r**2*(5*r0 - 3*r0**3))*Cos(2*z)) - &
               2*r**3*(r - r0)**2*Re*Cos(z)*Sin(time))/(2.*r**3*Re)
       ELSE IF (mode==1) THEN
          vv = vv +(-(r*(r**4 - 2*r*r0 - 2*r**3*r0 + 2*r0**2 + r**2*(-2 + r0**2))*Cos(time)*Cos(z)) - &
               (3*r**2 - 4*r*r0 + r0**2)*Re*Cos(time)**2*(3*r**2 - r0**2 + (2*r**4 - 4*r**3*r0 + r0**2 + &
               r**2*(-3 + 2*r0**2))*Cos(2*z)) + r**3*(r - r0)**2*Re*Cos(z)*Sin(time))/ (r**3*Re)
       ELSE IF (mode==2) THEN
          vv = vv -((r - r0)*Cos(time)**2*(4*r**2 - r*r0 - r0**2 + (3*r**4 - 7*r**3*r0 + r0**2 + &
               r**2*(-4 + 5*r0**2) + r*(r0 - r0**3))*Cos(2*z)))/(2.*r**2)
       END IF
    ELSE IF (TYPE==2) THEN
       IF (mode==1) THEN
          vv = vv + ((r - r0)**2*Cos(time)*(-2*r*Cos(z) + Re*Cos(time)*(r**4 - r*r0 - 2*r**3*r0 + &
               r0**2 + r**2*(-3 + r0**2) + (3*r**2 + r*r0 - r0**2)*Cos(2*z))))/(r**3*Re)
       ELSE IF (mode==2) THEN
          vv = vv + ((r - r0)**2*Cos(time)**2*(r**4 - r*r0 - 2*r**3*r0 + r0**2 + r**2*(-3 + r0**2) + &
               (3*r**2 + r*r0 - r0**2)*Cos(2*z)))/(2.*r**3)
       END IF
    ELSE IF (TYPE==3) THEN
       IF (mode==0) THEN
          vv = vv + (-3*r*(r - r0)**3*(3*r - r0)*Re*Cos(time)**2 + 2*(r**4 - 2*r**3*r0 + r0**2 + &
               r**2*(-3 + r0**2))*Cos(time)*Cos(z) - 2*r**2*(r - r0)**2*Re*Cos(z)*Sin(time))/(2.*r**2*Re)
       ELSE IF (mode==1) THEN
          vv = vv -(-2*r*(r**4 - 2*r*r0 - 2*r**3*r0 + 2*r0**2 + r**2*(-2 + r0**2))*Cos(time)*Cos(z) + &
               (r - r0)**3*(3*r - r0)*Re*Cos(time)**2*(1 + 4*r**2 - Cos(2*z)) + &
               2*r**3*(r - r0)**2*Re*Cos(z)*Sin(time))/(2.*r**3*Re)
       ELSE IF (mode==2) THEN
          vv = vv + ((r - r0)**3*(3*r - r0)*Cos(time)**2*(-1 - r**2 + Cos(2*z)))/(2.*r**3)
       END IF
    ELSE IF (TYPE==4) THEN
       IF (mode==1) THEN
          vv = vv + ((r - r0)**2*Cos(time)*(-4*r*Cos(z) + Re*Cos(time)*((-3*r + r0)**2 + &
               (2*r**4 + 6*r*r0 - 4*r**3*r0 - r0**2 + r**2*(-9 + 2*r0**2))*Cos(2*z))))/(2.*r**3*Re)
       ELSE IF (mode==2) THEN
          vv = vv + ((r - r0)**2*Cos(time)**2*(4*r - 2*r0 + (r*(-4 + (r - r0)**2) + 2*r0)*Cos(2*z)))/(2.*r**2)
       END IF
    ELSE IF (TYPE==5) THEN
       IF (mode==0) THEN
          vv = vv + (((3*r**4 - 4*r**3*r0 - r0**2 + r**2*(-3 + r0**2))*Cos(time) + &
               r*(r - r0)**2*(3*r**4 - r*r0 - 6*r**3*r0 + 2*r0**2 + r**2*(-4 + 3*r0**2))*Re*Cos(time)**2*Cos(z) - &
               r**2*(3*r**2 - 4*r*r0 + r0**2)*Re*Sin(time))*Sin(z))/(r**3*Re)
       ELSE IF (mode==1) THEN
          vv = vv + ((3*r**3 - 4*r0 - 4*r**2*r0 + r*r0**2)*Cos(time)*Sin(z) - r*(3*r**2 - 4*r*r0 + r0**2)*Re*Sin(time)*Sin(z) + &
               ((r - r0)**2*(4*r**4 - 2*r*r0 - 8*r**3*r0 + 3*r0**2 + r**2*(-5 + 4*r0**2))*Re*Cos(time)**2*Sin(2*z))/2.)/(r**2*Re)
       ELSE IF (mode==2) THEN
          vv = vv + ((r - r0)**2*(r**4 - r*r0 - 2*r**3*r0 + r0**2 + r**2*(-1 + r0**2))*Cos(time)**2*Sin(2*z))/(2.*r**2)
       END IF
    ELSE IF (TYPE==6) THEN
       IF (mode==1) THEN
          vv = vv -((((-r**3 + 2*r0 + 2*r**2*r0 - r*r0**2)*Cos(time) + 4*r*(r - r0)**3*Re*Cos(time)**2*Cos(z) + &
               r*(r - r0)**2*Re*Sin(time))*Sin(z))/(r**2*Re))
       ELSE IF (mode==2) THEN
          vv = vv -(((r - r0)**3*Cos(time)**2*Sin(2*z))/r)
       END IF
    END IF
    
    IF (TYPE == 1) THEN
       IF (mode == 0) THEN
          vv = vv - 0.25 * exp(2*z) * r**7 * (215.d0 + 51*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 1) THEN
          vv = vv - 5 * exp(2*z) * r**7 * (11.d0 + 3*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 2) THEN
          vv = vv + 14 * exp(2*z) * r**7 * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 3) THEN
          vv = vv + exp(2*z) * r**7 * (19.d0 + 3*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 4) THEN
          vv = vv + 0.75 * exp(2*z) * r**7 * (5.d0 + r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 2) THEN
       IF (mode == 1) THEN
          vv = vv - 3 * exp(2*z) * r**7 * (24.d0 + 5*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 2) THEN
          vv = vv - 4 * exp(2*z) * r**7 * (13.d0 + 3*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 3) THEN
          vv = vv - exp(2*z) * r**7 * (8.d0 + 3*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 4) THEN
          vv = vv + 2 * exp(2*z) * r**7 * (r - r0)**4 * cos(time)**4 * sin(z)**2 
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 3) THEN
       IF (mode == 0) THEN
          vv = vv - exp(2*z) * r**7 * (15.d0 + 2*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 1) THEN
          vv = vv - 0.5 * exp(2*z) * r**7 * (48.d0 + 7*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 2) THEN
          vv = vv - 2 * exp(2*z) * r**7 * (5.d0 + r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 3) THEN
          vv = vv - 0.5 * exp(2*z) * r**9 * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 4) THEN
          vv = vv + exp(2*z) * r**7 * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 4) THEN
       IF (mode == 1) THEN
          vv = vv - 0.5 * exp(2*z) * r**7 * (25.d0 + 2*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 2) THEN
          vv = vv - 0.25 * exp(2*z) * r**7 * (61.d0 + 6*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2  
       ELSE IF (mode == 3) THEN
          vv = vv - 0.5 * exp(2*z) * r**7 * (17.d0 + 2*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2  
       ELSE IF (mode == 4) THEN
          vv = vv - 0.125 * exp(2*z) * r**7 * (15.d0 + 2*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 5) THEN
       IF (mode == 0) THEN
          vv = vv - 0.125 * exp(2*z) * r**8 * (215.d0 + 34*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 1) THEN
          vv = vv - 2.5 * exp(2*z) * r**8 * (11.d0 + 2*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 2) THEN
          vv = vv + 7 * exp(2*z) * r**8 * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 3) THEN
          vv = vv + 0.5 * exp(2*z) * r**8 * (19.d0 + 2*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 4) THEN
          vv = vv + 0.125 * exp(2*z) * r**8 * (15.d0 + 2*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 6) THEN
       IF (mode == 1) THEN
          vv = vv - exp(2*z) * r**8 * (36.d0 + 5*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 2) THEN
          vv = vv - 2 * exp(2*z) * r**8 * (13.d0 + 2*r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 3) THEN
          vv = vv - exp(2*z) * r**8 * (4.d0 + r**2) * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE IF (mode == 4) THEN
          vv = vv + exp(2*z) * r**8 * (r - r0)**4 * cos(time)**4 * sin(z)**2
       ELSE
          vv = 0.d0
       END IF
    END IF

    RETURN

    !===Dummies variables to avoid warning
    np=ty
    !===Dummies variables to avoid warning
  END FUNCTION source_in_NS_momentum

  !===Extra source in temperature equation. Always called.
  FUNCTION source_in_temperature(TYPE, rr, m, t)RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t   
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv, r, z, c, lambda
    INTEGER                                           :: i
    REAL(KIND=8)                                      :: r0 = 0.5d0

    r = rr(1,:)
    z = rr(2,:)

    DO i=1,SIZE(rr,2)
       IF (r(i).LE.r0) THEN
          c(i) = inputs%vol_heat_capacity(1)
       ELSE
          c(i) = inputs%vol_heat_capacity(2)
       END IF
    END DO

    DO i=1,SIZE(rr,2)
       IF (r(i).LE.r0) THEN
          lambda(i) = inputs%temperature_diffusivity(1)
       ELSE
          lambda(i) = inputs%temperature_diffusivity(2)
       END IF
    END DO

    IF (TYPE==1) THEN
       IF (m==0) THEN
          vv = - (-(r**4 + 18*r*r0 - 2*r**3*r0 - 4*r0**2 + r**2*(-16.d0 + r0**2))*lambda*Cos(t) &
               + c*r**2*(r - r0)**2*Sin(t))*Sin(z)
       ELSE IF (m==1) THEN
          vv = -(-(r**4 + 16*r*r0 - 2*r**3*r0 - 3*r0**2 + r**2*(-15.d0 + r0**2))*lambda*Cos(t) &
               + c*r**2*(r - r0)**2*Sin(t))*Sin(z)
       ELSE
          vv = 0.d0
       END IF
    ELSE
       vv = 0.d0
    END IF

    IF (TYPE==1) THEN
       IF (m==0) THEN
          DO i=1,size(rr,2)
             IF (r(i)>r0) THEN
                vv(i) = vv(i) + c(i)*(-3*r(i)*(r(i) - r0)**4*Cos(t)**2*Sin(2*z(i)))/4.d0
             END IF
          END DO
       ELSE IF (m==1) THEN
          DO i=1,size(rr,2)
             IF (r(i)>r0) THEN
                vv(i) = vv(i) - c(i)*r(i)*(r(i) - r0)**4*Cos(t)**2*Sin(2*z(i))
             END IF
          END DO
       ELSE IF (m==2) THEN
          DO i=1,size(rr,2)
             IF (r(i)>r0) THEN
                vv(i) = vv(i) - c(i)/lambda(i)*(r(i)*(r(i) - r0)**4*Cos(t)**2*Sin(2*z(i)))/4.d0
             END IF
          END DO
       END IF
    END IF
    RETURN
  END FUNCTION source_in_temperature

!!$  FUNCTION source_in_level_set(interface_nb,TYPE, rr, m, t)RESULT(vv)
!!$    IMPLICIT NONE
!!$    INTEGER     ,                        INTENT(IN)   :: TYPE
!!$    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
!!$    INTEGER     ,                        INTENT(IN)   :: m, interface_nb
!!$    REAL(KIND=8),                        INTENT(IN)   :: t   
!!$    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
!!$
!!$    vv=0.d0
!!$    RETURN
!!$  END FUNCTION source_in_level_set

  !===Velocity for boundary conditions in Navier-Stokes.
  !===Can be used also to initialize velocity in: init_velocity_pressure_temperature 
  FUNCTION vv_exact(TYPE,rr,m,t) RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv, r, z
    REAL(KIND=8)                                      :: r0 = 0.5d0

    r = rr(1,:)
    z = rr(2,:)
  
    IF (TYPE==1) THEN
       IF ((m==0).OR.(m==1)) THEN
          vv = -(r-r0)**2*cos(t)*cos(z)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE==3) THEN
       IF ((m==0).OR.(m==1)) THEN
          vv = (r-r0)**2*cos(t)*cos(z)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE==5) THEN
       IF ((m==0).OR.(m==1)) THEN
          vv = (r-r0)*cos(t)*sin(z)/r * (3*r-r0)
       ELSE
          vv = 0.d0
      END IF
    ELSE IF (TYPE==6) THEN
       IF (m==1) THEN
          vv = (r-r0)*cos(t)*sin(z)/r * (r-r0)
       ELSE
          vv = 0.d0
       END IF
    ELSE
       vv = 0.d0
    END IF
    RETURN
  END FUNCTION vv_exact

!!$  !===Solid velocity imposed when using penalty technique
!!$  !===Defined in Fourier space on mode 0 only.
!!$  FUNCTION imposed_velocity_by_penalty(rr,t) RESULT(vv)
!!$    IMPLICIT NONE 
!!$    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
!!$    REAL(KIND=8),                        INTENT(IN)   :: t
!!$    REAL(KIND=8), DIMENSION(SIZE(rr,2),6)             :: vv
!!$
!!$    vv=0.d0
!!$    RETURN
!!$  END FUNCTION imposed_velocity_by_penalty

  !===Pressure for boundary conditions in Navier-Stokes.
  !===Can be used also to initialize pressure in the subroutine init_velocity_pressure.
  !===Use this routine for outflow BCs only.
  !===CAUTION: Do not enfore BCs on pressure where normal component 
  !            of velocity is prescribed.
  FUNCTION pp_exact(TYPE,rr,m,t) RESULT (vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv 
    REAL(KIND=8) :: r
    INTEGER      :: n

    vv=0.d0
    RETURN

    !===Dummies variables to avoid warning
    n=TYPE; n=SIZE(rr,1); n=m; r=t
    !===Dummies variables to avoid warning
  END FUNCTION pp_exact

  !===Temperature for boundary conditions in temperature equation.
  FUNCTION temperature_exact(TYPE,rr,m,t) RESULT (vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv, r, z
    REAL(KIND=8)                                      :: r0=0.5d0

    r = rr(1,:)
    z = rr(2,:)

    IF ((TYPE==1).AND.((m==0).OR.(m==1))) THEN
       vv = r**2*(r-r0)**2*sin(z)*cos(t)
    ELSE
       vv = 0.d0
    END IF
    RETURN
  END FUNCTION temperature_exact

!!$  !===Can be used to initialize level set in the subroutine init_level_set.
!!$  FUNCTION level_set_exact(interface_nb,TYPE,rr,m,t)  RESULT (vv)
!!$    IMPLICIT NONE
!!$    INTEGER     ,                        INTENT(IN)   :: TYPE
!!$    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
!!$    INTEGER     ,                        INTENT(IN)   :: m, interface_nb
!!$    REAL(KIND=8),                        INTENT(IN)   :: t
!!$    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
!!$
!!$    vv = 0.d0
!!$    RETURN
!!$  END FUNCTION level_set_exact

!!$  !===Penalty coefficient (if needed)
!!$  !===This coefficient is equal to zero in subdomain
!!$  !===where penalty is applied
!!$  FUNCTION penal_in_real_space(mesh,rr_gauss,angles,nb_angles,nb,ne,time) RESULT(vv)
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type)                            :: mesh
!!$    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)   :: rr_gauss
!!$    REAL(KIND=8), DIMENSION(:), INTENT(IN)     :: angles
!!$    INTEGER,                    INTENT(IN)     :: nb_angles
!!$    INTEGER,                    INTENT(IN)     :: nb, ne
!!$    REAL(KIND=8),               INTENT(IN)     :: time
!!$    REAL(KIND=8), DIMENSION(nb_angles,ne-nb+1) :: vv
!!$
!!$    vv = 1.d0
!!$    RETURN
!!$  END FUNCTION penal_in_real_space

  !===Extension of the velocity field in the solid.
  !===Used when temperature or Maxwell equations are solved.
  !===It extends the velocity field on the Navier-Stokes domain to a
  !===velocity field on the temperature and the Maxwell domain.
  !===It is also used if problem type=mxw and restart velocity
  !===is set to true in data (type problem denoted mxx in the code).
  FUNCTION extension_velocity(TYPE, H_mesh, mode, t, n_start) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type),                     INTENT(IN)   :: H_mesh     
    INTEGER     ,                        INTENT(IN)   :: TYPE, n_start
    INTEGER,                             INTENT(IN)   :: mode 
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(H_Mesh%np)                :: vv
    REAL(KIND=8) :: r
    INTEGER      :: n

    vv = 0.d0
    RETURN

    !===Dummies variables to avoid warning
    n=H_mesh%np; r=t; n=TYPE; n=mode; n=n_start
    !===Dummies variables to avoid warning
  END FUNCTION extension_velocity

  !===============================================================================
  !                       Boundary conditions for Maxwell
  !===============================================================================
!!$  !===Velocity used in the induction equation.
!!$  !===Used only if problem type is mxw and restart velocity is false
!!$  FUNCTION Vexact(m, H_mesh) RESULT(vv)  !Set uniquement a l'induction
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type),                       INTENT(IN) :: H_mesh 
!!$    INTEGER,                               INTENT(IN) :: m
!!$    REAL(KIND=8), DIMENSION(H_mesh%np,6)              :: vv
!!$
!!$    vv = 0.d0
!!$    RETURN
!!$  END FUNCTION Vexact

!!$  FUNCTION H_B_quasi_static(char_h_b, rr, m) RESULT(vv) 
!!$    IMPLICIT NONE
!!$    CHARACTER(LEN=1),                    INTENT(IN)   :: char_h_b
!!$    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
!!$    INTEGER,                             INTENT(IN)   :: m
!!$    REAL(KIND=8), DIMENSION(SIZE(rr,2),6)             :: vv
!!$
!!$    IF (inputs%if_quasi_static_approx) THEN
!!$       vv = 0.d0
!!$    ELSE 
!!$       CALL error_petsc('H_B_quasi_static should not be called') 
!!$    END IF
!!$    RETURN
!!$  END FUNCTION H_B_quasi_static

  !===Magnetic field for boundary conditions in the Maxwell equations.
  FUNCTION Hexact(H_mesh, TYPE, rr, m, mu_H_field, t) RESULT(vv)  
    IMPLICIT NONE
    TYPE(mesh_type),                     INTENT(IN)   :: H_mesh
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t 
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: mu_H_field
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv, r, z
    INTEGER      :: n

    r = rr(1,:)
    z = rr(2,:)

    IF (TYPE == 1) THEN
       IF (m == 0) THEN
          vv = - exp(z) * r**3 * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 2) THEN
       IF (m == 1) THEN
          vv = - exp(z) * r**3 * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 3) THEN
       IF (m == 0) THEN
          vv = exp(z) * r**3 * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 4) THEN
       IF (m == 1) THEN
          vv = exp(z) * r**3 * cos(t)
       ELSE 
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 5) THEN
       IF (m == 0) THEN
          vv = 4 * exp(z) * r**2 * cos(t)
       ELSE IF (m == 1) THEN
          vv = - exp(z) * r**2 * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE
       IF (m == 1) THEN
          vv = 4 * exp(z) * r**2 * cos(t)
       ELSE
          vv = 0.d0
       END IF
    END IF
    RETURN

    !===Dummies variables to avoid warning
    n=H_mesh%np; n=SIZE(mu_H_field);
    !===Dummies variables to avoid warning
  END FUNCTION Hexact

  !===Scalar potential for boundary conditions in the Maxwell equations.
  FUNCTION Phiexact(TYPE, rr, m, mu_phi,t) RESULT(vv) 
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: mu_phi, t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv   
    REAL(KIND=8) :: r
    INTEGER      :: n

    vv=0.d0
    RETURN

    !===Dummies variables to avoid warning
    n=TYPE; n=SIZE(rr,1); n=m; r=mu_phi; r=t
    !===Dummies variables to avoid warning
  END FUNCTION Phiexact

  !===Current in Ohm's law. Curl(H) = sigma(E + uxB) + current
  FUNCTION Jexact_gauss(TYPE, rr, m, mu_phi, sigma, mu_H, t, mesh_id, opt_B_ext) RESULT(vv) 
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m 
    REAL(KIND=8),                        INTENT(IN)   :: mu_phi, sigma, mu_H, t 
    INTEGER     ,                        INTENT(IN)   :: mesh_id
    REAL(KIND=8), DIMENSION(6), OPTIONAL,INTENT(IN)   :: opt_B_ext 
    REAL(KIND=8)                                      :: vv
    REAL(KIND=8)                                      :: r, z
    INTEGER      :: n

    r = rr(1)
    z = rr(2)

    IF (TYPE == 1) THEN
       IF (m == 0) THEN
          vv = - exp(z) * r**3 * cos(t)
       ELSE IF (m == 1) THEN
          vv = 4 * exp(z) * r * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 2) THEN
       IF (m == 1) THEN
          vv = - exp(z) * r * (-1.d0 + r**2) * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 3) THEN
       IF (m == 0) THEN
          vv = - exp(z) * r * (8.d0 + r**2) * cos(t)
       ELSE IF (m == 1) THEN
          vv = 2 * exp(z) * r * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 4) THEN
       IF (m == 1) THEN
          vv = - exp(z) * r * (8.d0 + r**2) * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 5) THEN
       IF (m == 0) THEN
          vv = 4 * exp(z) *r**2 * cos(t)
       ELSE IF (m == 1) THEN
          vv = exp(z) * r**2 * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE 
       IF (m == 1) THEN
          vv = 4 * exp(z) * r**2 * cos(t)
       ELSE
          vv = 0.d0
       END IF
    END IF
    RETURN

    !===Dummies variables to avoid warning
    r=mu_phi; r=sigma; r=mu_H; n=mesh_id
    IF (PRESENT(opt_B_ext)) r=opt_B_ext(1)
    !===Dummies variables to avoid warning
  END FUNCTION Jexact_gauss

  !===Electric field for Neumann BC (cf. doc)
  FUNCTION Eexact_gauss(TYPE, rr, m, mu_phi, sigma, mu_H, t) RESULT(vv)
    IMPLICIT NONE
    INTEGER,                             INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: mu_phi, sigma, mu_H, t
    REAL(KIND=8)                                      :: vv 
    REAL(KIND=8) :: r
    INTEGER      :: n

    vv = 0.d0
    RETURN

    !===Dummies variables to avoid warning
    r=rr(1); r=mu_phi; r=sigma; r=mu_H; r=t; n=TYPE; n=m
    !===Dummies variables to avoid warning
  END FUNCTION Eexact_gauss

  !===Initialization of magnetic field and scalar potential (if present)
  SUBROUTINE init_maxwell(H_mesh, phi_mesh, time, dt, mu_H_field, mu_phi, &
       list_mode, Hn1, Hn, phin1, phin)
    IMPLICIT NONE
    TYPE(mesh_type)                            :: H_mesh, phi_mesh     
    REAL(KIND=8),                   INTENT(OUT):: time
    REAL(KIND=8),                   INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(:),     INTENT(IN) :: mu_H_field
    REAL(KIND=8),                   INTENT(IN) :: mu_phi
    INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode    
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: Hn, Hn1
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: phin, phin1
    INTEGER                                    :: i, k
    
    time = -dt
    DO k=1,6
       DO i=1, SIZE(list_mode)
          Hn1(:,k,i) = Hexact(H_mesh,k, H_mesh%rr, list_mode(i), mu_H_field, time)
          IF (inputs%nb_dom_phi>0) THEN
             IF (k<3) THEN
                phin1(:,k,i) = Phiexact(k, phi_mesh%rr, list_mode(i) , mu_phi, time)
             ENDIF
          END  IF
       ENDDO
    ENDDO

    time = time + dt
    DO k=1,6
       DO i=1, SIZE(list_mode)
          Hn(:,k,i) = Hexact(H_mesh,k, H_mesh%rr, list_mode(i), mu_H_field, time)
          IF (inputs%nb_dom_phi>0) THEN
             IF (k<3) THEN
                phin(:,k,i) = Phiexact(k, phi_mesh%rr, list_mode(i), mu_phi, time)
             ENDIF
          END  IF
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE init_maxwell

!!$  !===Analytical permeability (if needed)
!!$  !===This function is not needed unless the flag
!!$  !===     ===Use FEM Interpolation for magnetic permeability  (true/false)
!!$  !===is activated and set to .FALSE. in the data data file. Default is .TRUE.
!!$  FUNCTION mu_bar_in_fourier_space(H_mesh,nb,ne,pts,pts_ids) RESULT(vv)
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type), INTENT(IN)                :: H_mesh
!!$    REAL(KIND=8), DIMENSION(ne-nb+1)           :: vv
!!$    INTEGER,     INTENT(IN)                    :: nb, ne
!!$    REAL(KIND=8),DIMENSION(2,ne-nb+1),OPTIONAL :: pts
!!$    INTEGER,     DIMENSION(ne-nb+1), OPTIONAL  :: pts_ids
!!$
!!$    vv = 1.d0
!!$    RETURN
!!$  END FUNCTION mu_bar_in_fourier_space

!!$  !===Analytical mu_in_fourier_space (if needed)
!!$  !===This function is not needed unless the flag
!!$  !===     ===Use FEM Interpolation for magnetic permeability  (true/false)
!!$  !===is activated and set to .FALSE. in the data data file. Default is .TRUE.
!!$  FUNCTION grad_mu_bar_in_fourier_space(pt,pt_id) RESULT(vv)
!!$    IMPLICIT NONE
!!$    REAL(KIND=8),DIMENSION(2), INTENT(in):: pt
!!$    INTEGER,DIMENSION(1), INTENT(in)     :: pt_id
!!$    REAL(KIND=8),DIMENSION(2)            :: vv
!!$
!!$    vv=0.d0
!!$    RETURN
!!$  END FUNCTION grad_mu_bar_in_fourier_space

!!$  !===Analytical permeability, mu in real space (if needed)
!!$  FUNCTION mu_in_real_space(H_mesh,angles,nb_angles,nb,ne,time) RESULT(vv)
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type), INTENT(IN)                :: H_mesh
!!$    REAL(KIND=8), DIMENSION(:), INTENT(IN)     :: angles
!!$    INTEGER, INTENT(IN)                        :: nb_angles
!!$    INTEGER, INTENT(IN)                        :: nb, ne
!!$    REAL(KIND=8), INTENT(IN)                   :: time
!!$    REAL(KIND=8), DIMENSION(nb_angles,ne-nb+1) :: vv
!!$
!!$    vv = 1.d0
!!$    RETURN
!!$  END FUNCTION mu_in_real_space

  !===Coefficient contaning the magnetic susceptibility for magnetic force in ferrofluids: 
  !===F = chi_coeff(T) * grad(H**2/2) (Kelvin form)
  !===or F = -H**2/2 * grad(chi_coeff(T)) (Helmholtz form)
  FUNCTION chi_coeff_law(temp) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: temp
    REAL(KIND=8) :: vv
    
    vv = temp**2
    RETURN
  END FUNCTION chi_coeff_law

  !===Coefficient contaning the temperature dependant factor in the
  !===pyromagnetic coefficient term of the temperature equation for ferrofluids: 
  !===T * dchi/dT(T) * D/Dt(H**2/2)
  FUNCTION T_dchi_dT_coeff_law(temp) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: temp
    REAL(KIND=8) :: vv
    
    vv = 0.d0
    RETURN

    !===Dummies variables to avoid warning
    vv=temp
    !===Dummies variables to avoid warning
  END FUNCTION T_dchi_dT_coeff_law

END MODULE boundary_test_33
