SUBMODULE (boundary_generic_module) boundary_generic
  USE my_util
  USE def_type_mesh
  USE input_data
  USE bessel
  USE user_data

  REAL(KIND=8),  PARAMETER   :: pi = 3.1415926535897932d0
  REAL(KIND=8),  PARAMETER   :: r0 = 0.5d0

  !===Dummy variables to avoid warning
  REAL(KIND=8) :: rd
  INTEGER      :: nd
  CHARACTER(LEN=2)  :: cd2
  CHARACTER(LEN=1)  :: cd1
  !===Dummy variables to avoid warning
CONTAINS
  !===============================================================================
  !                       Boundary conditions for Navier-Stokes
  !===============================================================================

  !===Initialize velocity, pressure
  MODULE SUBROUTINE init_velocity_pressure(mesh_f, mesh_c, time, dt, list_mode, &
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
  MODULE SUBROUTINE init_temperature(mesh, time, dt, list_mode, tempn_m1, tempn)
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

  !===Initialize concentration
  MODULE SUBROUTINE init_concentration(mesh, time, dt, list_mode, concn_m1, concn)
    IMPLICIT NONE
    TYPE(mesh_type)                            :: mesh
    REAL(KIND=8),                   INTENT(OUT):: time
    REAL(KIND=8),                   INTENT(IN) :: dt
    INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: concn_m1, concn
    INTEGER                                    :: mode, i, j

    time = 0.d0
    DO i= 1, SIZE(list_mode)
       mode = list_mode(i)
       DO j = 1, 2
          concn_m1(:,j,i) = concentration_exact(j, mesh%rr, mode, time-dt)
          concn   (:,j,i) = concentration_exact(j, mesh%rr, mode, time)
       ENDDO
    ENDDO
  END SUBROUTINE init_concentration

  !===Initialize level_set
  MODULE SUBROUTINE init_level_set(pp_mesh, time, &
       dt, list_mode, level_set_m1, level_set)
    IMPLICIT NONE
    TYPE(mesh_type)                              :: pp_mesh
    REAL(KIND=8),                     INTENT(OUT):: time
    REAL(KIND=8),                     INTENT(IN) :: dt
    INTEGER,      DIMENSION(:),       INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(OUT):: level_set, level_set_m1
    INTEGER                                      :: mode, i, j, n

    time = 0.d0
    DO i= 1, SIZE(list_mode)
       mode = list_mode(i)
       DO j = 1, 2
          !===level_set
          DO n = 1, inputs%nb_fluid -1
             level_set_m1(n,:,j,i) = level_set_exact(n,j,pp_mesh%rr,mode,time-dt)
             level_set   (n,:,j,i) = level_set_exact(n,j,pp_mesh%rr,mode,time)
          END DO
       END DO
    END DO
  END SUBROUTINE init_level_set

  !===Source in momemtum equation. Always called.
  MODULE FUNCTION source_in_NS_momentum(TYPE, rr, mode, i, time, Re, ty, density, tempn, concn) RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                             INTENT(IN) :: TYPE
    REAL(KIND=8), DIMENSION(:,:),             INTENT(IN) :: rr
    INTEGER     ,                             INTENT(IN) :: mode, i
    REAL(KIND=8),                             INTENT(IN) :: time
    REAL(KIND=8),                             INTENT(IN) :: Re
    CHARACTER(LEN=2),                         INTENT(IN) :: ty
    REAL(KIND=8), DIMENSION(:,:,:),           INTENT(IN) :: density
    REAL(KIND=8), DIMENSION(:,:,:),           INTENT(IN) :: tempn
    REAL(KIND=8), DIMENSION(:,:,:),           INTENT(IN) :: concn
    REAL(KIND=8), DIMENSION(SIZE(rr,2))                  :: vv, r, z, lambda
    REAL(KIND=8)                                         :: alpha
    INTEGER                                              :: j

    alpha = inputs%gravity_coefficient

    r = rr(1,:)
    z = rr(2,:)

    DO j=1,SIZE(rr,2)
       IF (r(j).LE.r0) THEN
          lambda(j) = inputs%temperature_diffusivity(1)
       ELSE
          lambda(j) = inputs%temperature_diffusivity(2)
       END IF
    END DO

    IF (TYPE==5) THEN
       vv = alpha*(tempn(:,1,i) - temperature_exact(1,rr,mode,time))
    ELSE IF (TYPE==6) THEN
       vv = alpha*(tempn(:,2,i) - temperature_exact(2,rr,mode,time))
    ELSE
       vv = 0.d0
    END IF

    IF (TYPE==1) THEN
       IF (mode==0) THEN
          vv = vv -(4*pi*r*(4*pi**2*r**4 - 8*pi**2*r**3*r0 + r0**2 + r**2*(-3 + 4*pi**2*r0**2))*Cos(time)*Cos(2*pi*z) + &
               (r - r0)*Re*Cos(time)**2*(14*r**3 - 5*r**2*r0 - 5*r*r0**2 + 2*r0**3 + &
               (36*pi**2*r**5 - 84*pi**2*r**4*r0 + 5*r*r0**2 - 2*r0**3 + 2*r**3*(-7 + 30*pi**2*r0**2) + &
               r**2*(5*r0 - 12*pi**2*r0**3))*Cos(4*pi*z)) - &
               4*pi*r**3*(r - r0)**2*Re*Cos(2*pi*z)*Sin(time))/(2.*r**3*Re)
       ELSE IF (mode==1) THEN
          vv = vv -2*(2*pi*r*(2*pi**2*r**4 - r*r0 - 4*pi**2*r**3*r0 + r0**2 + r**2*(-1 + 2*pi**2*r0**2))*Cos(time)*Cos(2*pi*z) + &
               ((3*r**2 - 4*r*r0 + r0**2)*Re*Cos(time)**2*(3*r**2 - r0**2 + &
               (8*pi**2*r**4 - 16*pi**2*r**3*r0 + r0**2 + r**2*(-3 + 8*pi**2*r0**2))*Cos(4*pi*z)))/2. - &
               pi*r**3*(r - r0)**2*Re*Cos(2*pi*z)*Sin(time))/(r**3*Re)
       ELSE IF (mode==2) THEN
          vv = vv -((r - r0)*Cos(time)**2*(4*r**2 - r*r0 - r0**2 + (12*pi**2*r**4 - &
               28*pi**2*r**3*r0 + r0**2 + 4*r**2*(-1.d0 + 5*pi**2*r0**2) + &
               r*(r0 - 4*pi**2*r0**3))*Cos(4*pi*z)))/(2.*r**2)
       END IF
    ELSE IF (TYPE==2) THEN
       IF (mode==1) THEN
          vv = vv + ((r - r0)**2*Cos(time)*(-4*pi*r*Cos(2*pi*z) + Re*Cos(time)*(4*pi**2*r**4 - r*r0 - 8*pi**2*r**3*r0 + r0**2 + &
               r**2*(-3.d0 + 4*pi**2*r0**2) + (3*r**2 + r*r0 - r0**2)*Cos(4*pi*z))))/(r**3*Re)
       ELSE IF (mode==2) THEN
          vv = vv + ((r - r0)**2*Cos(time)**2*(4*pi**2*r**4 - r*r0 - 8*pi**2*r**3*r0 + r0**2 + &
               r**2*(-3.d0 + 4*pi**2*r0**2) + (3*r**2 + r*r0 - r0**2)*Cos(4*pi*z)))/(2.*r**3)
       END IF
    ELSE IF (TYPE==3) THEN
       IF (mode==0) THEN
          vv = vv + (2*pi*(-3*pi*r*(r - r0)**3*(3*r - r0)*Re*Cos(time)**2 + &
               (4*pi**2*r**4 - 8*pi**2*r**3*r0 + r0**2 + r**2*(-3.d0 + 4*pi**2*r0**2))*Cos(time)*Cos(2*pi*z) - &
               r**2*(r - r0)**2*Re*Cos(2*pi*z)*Sin(time)))/(r**2*Re)
       ELSE IF (mode==1) THEN
          vv = vv + (4*pi*r*(2*pi**2*r**4 - r*r0 - 4*pi**2*r**3*r0 + r0**2 + &
               r**2*(-1 + 2*pi**2*r0**2))*Cos(time)*Cos(2*pi*z) + &
               ((r - r0)**3*(3*r - r0)*Re*Cos(time)**2*(-1.d0 - 16*pi**2*r**2 + Cos(4*pi*z)))/2. - &
               2*pi*r**3*(r - r0)**2*Re*Cos(2*pi*z)*Sin(time))/(r**3*Re)
       ELSE IF (mode==2) THEN
          vv = vv + ((r - r0)**3*(3*r - r0)*Cos(time)**2*(-1.d0 - 4*pi**2*r**2 + Cos(4*pi*z)))/(2.*r**3)
       END IF
    ELSE IF (TYPE==4) THEN
       IF (mode==1) THEN
          vv = vv + ((r - r0)**2*Cos(time)*(-8*pi*r*Cos(2*pi*z) + Re*Cos(time)*((-3*r + r0)**2 + &
               (8*pi**2*r**4 + 6*r*r0 - 16*pi**2*r**3*r0 - r0**2 + r**2*(-9.d0 + 8*pi**2*r0**2))*Cos(4*pi*z))))/(2.*r**3*Re)
       ELSE IF (mode==2) THEN
          vv = vv + ((r - r0)**2*Cos(time)**2*(2*r - r0 + (2*r*(-1.d0 + pi*(r - r0))*(1 + pi*(r - r0)) + r0)*Cos(4*pi*z)))/r**2
       END IF
    ELSE IF (TYPE==5) THEN
       IF (mode==0) THEN
          vv = vv + ((4*(12*pi**2*r**4 - 16*pi**2*r**3*r0 - r0**2 + &
               r**2*(-3.d0 + 4*pi**2*r0**2))*Cos(time)*Sin(2*pi*z) - &
               4*r**2*(3*r**2 - 4*r*r0 + r0**2)*Re*Sin(time)*Sin(2*pi*z) + &
               pi*r*(r - r0)**2*(12*pi**2*r**4 - r*r0 - 24*pi**2*r**3*r0 + 2*r0**2 + &
               4*r**2*(-1.d0 + 3*pi**2*r0**2))*Re*4*Cos(time)**2*Sin(4*pi*z)))/(4.*r**3*Re)
       ELSE IF (mode==1) THEN
          vv = vv + (4*(-r0 + pi**2*r*(3*r**2 - 4*r*r0 + r0**2))*Cos(time)*Sin(2*pi*z) - &
               r*(3*r**2 - 4*r*r0 + r0**2)*Re*Sin(time)*Sin(2*pi*z) + &
               pi*(r - r0)**2*(16*pi**2*r**4 - 2*r*r0 - 32*pi**2*r**3*r0 + 3*r0**2 + &
               r**2*(-5.d0 + 16*pi**2*r0**2))*Re*Cos(time)**2*Sin(4*pi*z))/(r**2*Re)
       ELSE IF (mode==2) THEN
          vv = vv + (pi*(r - r0)**2*(4*pi**2*r**4 - r*r0 - 8*pi**2*r**3*r0 + r0**2 + &
               r**2*(-1.d0 + 4*pi**2*r0**2))*Cos(time)**2*Sin(4*pi*z))/r**2
       END IF
    ELSE IF (TYPE==6) THEN
       IF (mode==1) THEN
          vv = vv + (2*(2*pi**2*r*(r - r0)**2 - r0)*Cos(time)*Sin(2*pi*z) - r*(r - r0)**2*Re*Sin(time)*Sin(2*pi*z) -&
               4*pi*r*(r - r0)**3*Re*Cos(time)**2*Sin(4*pi*z))/(r**2*Re)
       ELSE IF (mode==2) THEN
          vv = vv -2*pi*(r - r0)**3*Cos(time)**2*Sin(4*pi*z)/r
       END IF
    END IF

    IF ((TYPE==1).AND.(mode==1)) THEN
       vv = vv + 3*r**2*cos(time)*sin(2*pi*z)
    ELSE IF ((TYPE==4).AND.(mode==1)) THEN
       vv = vv - r**2*cos(time)*sin(2*pi*z)
    ELSE IF ((TYPE==5).AND.(mode==1)) THEN
       vv = vv + 2*pi*r**3*cos(time)*cos(2*pi*z)
    END IF

    IF (TYPE == 1) THEN
       IF (mode == 0) THEN
          vv = vv + r**7 * (r - r0)**2 * cos(time)**4 * (215.d0 + 204 * pi**2 * r**2 + &
               (215.d0 - 204 * pi**2 * r**2) * cos(4*pi * z)) * sin(2*pi*z)**2 / &
               (8 * lambda**2)
       ELSE IF (mode == 1) THEN
          vv = vv - 5 * r**7 * (r - r0)**2 * cos(time)**4 * (-11.d0 - 12 * pi**2 * r**2 + &
               (-11.d0 + 12 * pi**2 * r**2) * cos(4*pi*z)) * sin(2*pi*z)**2 / &
               (2 * lambda**2)
       ELSE IF (mode == 2) THEN
          vv = vv - 7 * r**7 * (r - r0)**2 * cos(time)**4 * sin(4*pi*z)**2 / &
               (2 * lambda**2)
       ELSE IF (mode == 3) THEN
          vv = vv + r**7 * (r - r0)**2 * cos(time)**4 * (-19.d0 - 12 * pi**2 *r**2 + &
               (-19.d0 + 12 * pi**2 *r**2) * cos(4*pi*z)) * sin(2*pi*z)**2 / &
               (2 * lambda**2)
       ELSE IF (mode == 4) THEN
          vv = vv + 3 * r**7 * (r - r0)**2 * cos(time)**4 * (-5.d0 - 4 * pi**2 * r**2 + &
               (-5.d0 + 4 * pi**2 * r**2) * cos(4*pi*z)) * sin(2*pi*z)**2 / &
               (8 * lambda**2)
       END IF
    ELSE IF (TYPE == 2) THEN
       IF (mode == 1) THEN
          vv = vv - 6 * r**7 * (r - r0)**2 * cos(time)**4 * (-6.d0 - 5 * pi**2 * r**2 + &
               (-6.d0 + 5.d0 * pi**2 * r**2) * cos(4*pi*z)) * sin(2*pi*z)**2 / &
               lambda**2
       ELSE IF (mode == 2) THEN
          vv = vv + 2 * r**7 * (r - r0)**2 * cos(time)**4 * (13.d0 + 12 * pi**2 * r**2 + &
               (13.d0 - 12 * pi**2 * r**2) * cos(4*pi*z)) * sin(2*pi*z)**2 / &
               lambda**2
       ELSE IF (mode == 3) THEN
          vv = vv + 2 * r**7 * (r - r0)**2 * cos(time)**4 * (2.d0 + 3 * pi**2 * r**2 + &
               (2.d0 - 3 * pi**2 * r**2) * cos(4*pi*z)) * sin(2*pi*z)**2 / &
               lambda**2
       ELSE IF (mode == 4) THEN
          vv = vv - r**7 * (r - r0)**2 * cos(time)**4 * sin(4*pi*z)**2 / &
               (2 * lambda**2)
       END IF
    ELSE IF (TYPE == 3) THEN
       IF (mode == 0) THEN
          vv = vv + r**7 * (r - r0)**2 * cos(time)**4 * (15.d0 + 8 * pi**2 * r**2 + &
               (15.d0 - 8 * pi**2 * r**2) * cos(4*pi*z)) * sin(2*pi*z)**2 / &
               (2 * lambda**2)
       ELSE IF (mode == 1) THEN
          vv = vv + r**7 * (r - r0)**2 * cos(time)**4 * (12.d0 + 7 * pi**2 * r**2 + &
               (12.d0 - 7 * pi**2 * r**2) * cos(4*pi*z)) * sin(2*pi*z)**2 / &
               lambda**2
       ELSE IF (mode == 2) THEN
          vv = vv + r**7 * (r - r0)**2 * cos(time)**4 * (5.d0 + 4 * pi**2 * r**2 + &
               (5.d0 - 4 * pi**2 * r**2) * cos(4*pi*z)) * sin(2*pi*z)**2 / &
               lambda**2
       ELSE IF (mode == 3) THEN
          vv = vv + 2 * pi**2 * r**9 * (r - r0)**2 * cos(time)**4 * sin(2*pi*z)**4 / &
               lambda**2
       ELSE IF (mode == 4) THEN
          vv = vv - r**7 * (r - r0)**2 * cos(time)**4 * sin(4*pi*z)**2 / &
               (4 * lambda**2)
       END IF
    ELSE IF (TYPE == 4) THEN
       IF (mode == 1) THEN
          vv = vv + r**7 * (r - r0)**2 * cos(time)**4 * (25.d0 + 8 * pi**2 * r**2 + &
               (25.d0 - 8 * pi**2 * r**2) * cos(4*pi*z)) * sin(2*pi*z)**2 / &
               (4 * lambda**2)
       ELSE IF (mode == 2) THEN
          vv = vv + r**7 * (r - r0)**2 * cos(time)**4 * (61.d0 + 24 * pi**2 * r**2 + &
               (61.d0 - 24 * pi**2 * r**2) * cos(4*pi*z)) * sin(2*pi*z)**2 / &
               (8 * lambda**2)
       ELSE IF (mode == 3) THEN
          vv = vv + r**7 * (r - r0)**2 * cos(time)**4 * (17.d0 + 8 * pi**2 * r**2 + &
               (17.d0 - 8 * pi**2 * r**2) * cos(4*pi*z)) * sin(2*pi*z)**2 / &
               (4 * lambda**2)
       ELSE IF (mode == 4) THEN
          vv = vv + r**7 * (r - r0)**2 * cos(time)**4 * (15.d0 + 8 * pi**2 * r**2 + &
               (15.d0 - 8 * pi**2 * r**2) * cos(4*pi*z)) * sin(2*pi*z)**2 / &
               (16 * lambda**2)
       END IF
    ELSE IF (TYPE == 5) THEN
       IF (mode == 0) THEN
          vv = vv + pi * r**8 * (-215.d0 + 136 * pi**2 * r**2) * (r - r0)**2 * cos(time)**4 * &
               cos(2*pi*z) * sin(2*pi*z)**3 / (4 * lambda**2)
       ELSE IF (mode == 1) THEN
          vv = vv + 5 * pi * r**8 * (-11.d0 + 8 * pi**2 * r**2) * (r - r0)**2 * cos(time)**4 * &
               cos(2*pi*z) * sin(2*pi*z)**3 / lambda**2
       ELSE IF (mode == 2) THEN
          vv = vv + 14 * pi * r**8 * (r - r0)**2 * cos(time)**4 * &
               cos(2*pi*z) * sin(2*pi*z)**3 / lambda**2
       ELSE IF (mode == 3) THEN
          vv = vv - pi * r**8 * (-19.d0 + 8 * pi**2 * r**2) * (r - r0)**2 * cos(time)**4 * &
               cos(2*pi*z) * sin(2*pi*z)**3 / lambda**2
       ELSE IF (mode == 4) THEN
          vv = vv - pi * r**8 * (-15.d0 + 8 * pi**2 * r**2) * (r - r0)**2 * cos(time)**4 * &
               cos(2*pi*z) * sin(2*pi*z)**3 / (4 * lambda**2)
       END IF
    ELSE IF (TYPE == 6) THEN
       IF (mode == 1) THEN
          vv = vv + 8 * pi * r**8 * (-9.d0 + 5 * pi**2 * r**2) * (r - r0)**2 * cos(time)**4 * &
               cos(2*pi*z) * sin(2*pi*z)**3 / lambda**2
       ELSE IF (mode == 2) THEN
          vv = vv + 4 * pi * r**8 * (-13.d0 + 8 * pi**2 * r**2) * (r - r0)**2 * cos(time)**4 * &
               cos(2*pi*z) * sin(2*pi*z)**3 / lambda**2
       ELSE IF (mode == 3) THEN
          vv = vv + 8 * pi * r**8 * (-1.d0 + pi**2 * r**2) * (r - r0)**2 * cos(time)**4 * &
               cos(2*pi*z) * sin(2*pi*z)**3 / lambda**2
       ELSE IF (mode == 4) THEN
          vv = vv + 2 * pi * r**8 * (r - r0)**2 * cos(time)**4 * &
               cos(2*pi*z) * sin(2*pi*z)**3 / lambda**2
       END IF
    END IF

    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=mode; nd=i; rd=time; rd=Re; cd2=ty
    nd=SIZE(density,1); nd=SIZE(tempn,1); nd=SIZE(concn,1)
    !===Dummy variables to avoid warning
  END FUNCTION source_in_NS_momentum

  !===Extra source in temperature equation. Always called.
  MODULE FUNCTION source_in_temperature(TYPE, rr, m, t)RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv, r, z, c, lambda
    INTEGER                                           :: i

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
          vv = ((-9*r + 4*pi**2*r**3 + 4*r0 - 4*pi**2*r**2*r0) * lambda * Cos(t) + &
               c * r**2 * (-r + r0) * Sin(t)) * Sin(2*pi*z) / lambda
       ELSE IF (m==1) THEN
          vv = - ((-3*r0 + 4*r*(2.d0 + pi**2*r*(-r + r0))) * lambda * Cos(t) + &
               c * r**2 * (r - r0) * Sin(t)) * Sin(2*pi*z) / lambda
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
                vv(i) = vv(i) + 3 * c(i) * pi * r(i) * (r(i) - r0)**2 * r0 * &
                     Cos(t)**2 * Sin(4*pi*z(i)) / (2 * lambda(i))
             END IF
          END DO
       ELSE IF (m==1) THEN
          DO i=1,size(rr,2)
             IF (r(i)>r0) THEN
                vv(i) = vv(i) + 2 * c(i) * pi * r(i) * (r(i) - r0)**2 * r0 * &
                     Cos(t)**2 * Sin(4*pi*z(i)) / lambda(i)
             END IF
          END DO
       ELSE IF (m==2) THEN
          DO i=1,size(rr,2)
             IF (r(i)>r0) THEN
                vv(i) = vv(i) + c(i) * pi * r(i) * (r(i) - r0)**2 * r0 * &
                     Cos(t)**2 * Sin(4*pi*z(i)) / (2 *lambda(i))
             END IF
          END DO
       END IF
    END IF

    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION source_in_temperature

  MODULE FUNCTION source_in_level_set(interface_nb,TYPE, rr, m, t)RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m, interface_nb
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    CALL error_petsc('source_in_level_set: should not be called for this test')
    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; nd=interface_nb; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION source_in_level_set

  !===Velocity for boundary conditions in Navier-Stokes.
  !===Can be used also to initialize velocity in: init_velocity_pressure_temperature
  MODULE FUNCTION vv_exact(TYPE,rr,m,t) RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv, r, z

    r = rr(1,:)
    z = rr(2,:)

    IF (TYPE==1) THEN
       IF ((m==0).OR.(m==1)) THEN
          vv = -2*pi*(r-r0)**2*cos(t)*cos(2*pi*z)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE==3) THEN
       IF ((m==0).OR.(m==1)) THEN
          vv = 2*pi*(r-r0)**2*cos(t)*cos(2*pi*z)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE==5) THEN
       IF ((m==0).OR.(m==1)) THEN
          vv = (r-r0)*cos(t)*sin(2*pi*z)/r * (3*r-r0)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE==6) THEN
       IF (m==1) THEN
          vv = (r-r0)*cos(t)*sin(2*pi*z)/r * (r-r0)
       ELSE
          vv = 0.d0
       END IF
    ELSE
       vv = 0.d0
    END IF

    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION vv_exact

  !===Solid velocity imposed when using penalty technique
  !===Defined in Fourier space on mode 0 only.
  MODULE FUNCTION imposed_velocity_by_penalty(rr,t) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2),6)             :: vv

    vv = 0.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=SIZE(rr,1); rd=t
    !===Dummy variables to avoid warning
  END FUNCTION imposed_velocity_by_penalty

  !===Pressure for boundary conditions in Navier-Stokes.
  !===Can be used also to initialize pressure in the subroutine init_velocity_pressure.
  !===Use this routine for outflow BCs only.
  !===CAUTION: Do not enfore BCs on pressure where normal component
  !            of velocity is prescribed.
  MODULE FUNCTION pp_exact(TYPE,rr,m,t) RESULT (vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv, r, z

    r = rr(1,:)
    z = rr(2,:)

    IF ((TYPE==1).AND.(m==1)) THEN
       vv = r**3*sin(2*pi*z)*cos(t)
    ELSE
       vv = 0.d0
    END IF

    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION pp_exact

  !===Temperature for boundary conditions in temperature equation.
  MODULE FUNCTION temperature_exact(TYPE,rr,m,t) RESULT (vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv, r, z, lambda
    INTEGER                                           :: i

    r = rr(1,:)
    z = rr(2,:)

    DO i=1,SIZE(rr,2)
       IF (r(i).LE.r0) THEN
          lambda(i) = inputs%temperature_diffusivity(1)
       ELSE
          lambda(i) = inputs%temperature_diffusivity(2)
       END IF
    END DO

    IF ((TYPE==1).AND.((m==0).OR.(m==1))) THEN
       vv = r**2*(r-r0)*sin(2*pi*z)*cos(t) / lambda
    ELSE
       vv = 0.d0
    END IF

    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION temperature_exact

  !===Concentration for boundary conditions in concentration equation.
  MODULE FUNCTION concentration_exact(TYPE,rr,m,t) RESULT (vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    CALL error_petsc('concentration_exact: should not be called for this test')
    RETURN
     !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; rd=t
    !===Dummy variables to avoid warning
 END FUNCTION concentration_exact

  !===Can be used to initialize level set in the subroutine init_level_set.
  MODULE FUNCTION level_set_exact(interface_nb,TYPE,rr,m,t)  RESULT (vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m, interface_nb
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    CALL error_petsc('level_set_exact: should not be called for this test')
    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; nd=interface_nb; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION level_set_exact

  !===Penalty coefficient (if needed)
  !===This coefficient is equal to zero in subdomain
  !===where penalty is applied
  MODULE FUNCTION penal_in_real_space(mesh,rr_gauss,angles,nb_angles,nb,ne,time) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type)                            :: mesh
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)   :: rr_gauss
    REAL(KIND=8), DIMENSION(:), INTENT(IN)     :: angles
    INTEGER,                    INTENT(IN)     :: nb_angles
    INTEGER,                    INTENT(IN)     :: nb, ne
    REAL(KIND=8),               INTENT(IN)     :: time
    REAL(KIND=8), DIMENSION(nb_angles,ne-nb+1) :: vv

    vv = 1.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=mesh%np; nd=SIZE(rr_gauss,1); nd=SIZE(angles); nd=nb_angles; nd=nb; nd=ne; rd=time
    !===Dummy variables to avoid warning
  END FUNCTION penal_in_real_space

  !===Extension of the velocity field in the solid.
  !===Used when temperature or Maxwell equations are solved.
  !===It extends the velocity field on the Navier-Stokes domain to a
  !===velocity field on the temperature and the Maxwell domain.
  !===It is also used if problem type=mxw and restart velocity
  !===is set to true in data (type problem denoted mxx in the code).
  MODULE FUNCTION extension_velocity(TYPE, H_mesh, mode, t, n_start) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type),                     INTENT(IN)   :: H_mesh
    INTEGER     ,                        INTENT(IN)   :: TYPE, n_start
    INTEGER,                             INTENT(IN)   :: mode
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(H_Mesh%np)                :: vv

    vv = 0.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=H_mesh%np; nd=TYPE; nd=n_start; nd=mode; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION extension_velocity

  MODULE FUNCTION extension_temperature(TYPE, H_mesh, mode, t, n_start) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type),                     INTENT(IN)   :: H_mesh
    INTEGER     ,                        INTENT(IN)   :: TYPE, n_start
    INTEGER,                             INTENT(IN)   :: mode
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(H_Mesh%np)                :: vv

    vv = 0.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=H_mesh%np; nd=TYPE; nd=n_start; nd=mode; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION extension_temperature

  MODULE FUNCTION extension_concentration(TYPE, vv_mesh, mode, t, n_start) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type),                     INTENT(IN)   :: vv_mesh
    INTEGER     ,                        INTENT(IN)   :: TYPE, n_start
    INTEGER,                             INTENT(IN)   :: mode
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(vv_mesh%np)                :: vv

    vv = 0.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=vv_mesh%np; nd=TYPE; nd=n_start; nd=mode; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION extension_concentration

  !===============================================================================
  !                       Boundary conditions for Maxwell
  !===============================================================================
  !===Velocity used in the induction equation.
  !===Used only if problem type is mxw and restart velocity is false
  MODULE FUNCTION Vexact(m, H_mesh) RESULT(vv)  !Set uniquement a l'induction
    IMPLICIT NONE
    TYPE(mesh_type),                       INTENT(IN) :: H_mesh
    INTEGER,                               INTENT(IN) :: m
    REAL(KIND=8), DIMENSION(H_mesh%np,6)              :: vv

    vv = 0.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=H_mesh%np; nd=m
    !===Dummy variables to avoid warning
  END FUNCTION Vexact

  MODULE FUNCTION H_B_quasi_static(char_h_b, rr, m) RESULT(vv)
    IMPLICIT NONE
    CHARACTER(LEN=1),                    INTENT(IN)   :: char_h_b
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m
    REAL(KIND=8), DIMENSION(SIZE(rr,2),6)             :: vv

    IF (inputs%if_quasi_static_approx) THEN
       vv = 0.d0
    ELSE
       CALL error_petsc('H_B_quasi_static should not be called')
    END IF
    RETURN
    !===Dummy variables to avoid warning
    cd1=char_h_b; nd=SIZE(rr,1); nd=m
    !===Dummy variables to avoid warning
  END FUNCTION H_B_quasi_static

  !===Magnetic field for boundary conditions in the Maxwell equations.
  MODULE FUNCTION Hexact(H_mesh, TYPE, rr, m, mu_H_field, t) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type),                     INTENT(IN)   :: H_mesh
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: mu_H_field
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv, r, z

    r = rr(1,:)
    z = rr(2,:)

    IF (TYPE == 1) THEN
       IF (m == 0) THEN
          vv = 2 * pi * r**3 * sin(2*pi*z) * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 2) THEN
       IF (m == 1) THEN
          vv = 2 * pi * r**3 * sin(2*pi*z) * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 3) THEN
       IF (m == 0) THEN
          vv = - 2 * pi * r**3 * sin(2*pi*z) * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 4) THEN
       IF (m == 1) THEN
          vv = - 2 * pi * r**3 * sin(2*pi*z) * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 5) THEN
       IF (m == 0) THEN
          vv = 4 * r**2 * cos(2*pi*z) * cos(t)
       ELSE IF (m == 1) THEN
          vv = - r**2 * cos(2*pi*z) * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE
       IF (m == 1) THEN
          vv = 4 * r**2 * cos(2*pi*z) * cos(t)
       ELSE
          vv = 0.d0
       END IF
    END IF
    RETURN
    !===Dummy variables to avoid warning
    nd=H_mesh%np; nd=TYPE; nd=SIZE(rr,1); nd=m; rd=t; nd=SIZE(mu_H_field)
    !===Dummy variables to avoid warning
 END FUNCTION Hexact

  !===Scalar potential for boundary conditions in the Maxwell equations.
  MODULE FUNCTION Phiexact(TYPE, rr, m, mu_phi,t) RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: mu_phi, t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; rd=mu_phi; rd=t
    !===Dummy variables to avoid warning
 END FUNCTION Phiexact

  !===Current in Ohm's law. Curl(H) = sigma(E + uxB) + current
  MODULE FUNCTION Jexact_gauss(TYPE, rr, m, mu_phi, sigma, mu_H, t, mesh_id, opt_B_ext) RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: mu_phi, sigma, mu_H, t
    INTEGER     ,                        INTENT(IN)   :: mesh_id
    REAL(KIND=8), DIMENSION(6), OPTIONAL,INTENT(IN)   :: opt_B_ext
    REAL(KIND=8)                                      :: vv, r, z
    REAL(KIND=8)                                      :: pi = acos(-1.d0)

    r = rr(1)
    z = rr(2)

    IF (TYPE == 1) THEN
       IF (m == 0) THEN
          vv = 4 * pi**2 * r**3 * cos(2*pi*z) * cos(t)
       ELSE IF (m == 1) THEN
          vv = 4 * r * cos(2*pi*z) * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 2) THEN
       IF (m == 1) THEN
          vv = r * (1.d0 + 4 * pi**2 * r**2) * cos(2*pi*z) * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 3) THEN
       IF (m == 0) THEN
          vv = 4 * r * (-2.d0 + pi**2 * r**2) * cos(2*pi*z) * cos(t)
       ELSE IF (m == 1) THEN
          vv = 2 * r * cos(2*pi*z) * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 4) THEN
       IF (m == 1) THEN
          vv = 4 * r * (-2.d0 + pi**2 * r**2) * cos(2*pi*z) * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE IF (TYPE == 5) THEN
       IF (m == 0) THEN
          vv = - 8 * pi * r**2 * sin(2*pi*z) * cos(t)
       ELSE IF (m == 1) THEN
          vv = - 2 * pi * r**2 * sin(2*pi*z) * cos(t)
       ELSE
          vv = 0.d0
       END IF
    ELSE
       IF (m == 1) THEN
          vv = - 8 * pi * r**2 * sin(2*pi*z) * cos(t)
       ELSE
          vv = 0.d0
       END IF
    END IF
    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; rd=mu_phi; rd=sigma; rd=mu_H; rd=t; nd=mesh_id
    IF (PRESENT(opt_B_ext)) nd=SIZE(opt_B_ext)
    !===Dummy variables to avoid warning
  END FUNCTION Jexact_gauss

  !===Electric field for Neumann BC (cf. doc)
  MODULE FUNCTION Eexact_gauss(TYPE, rr, m, mu_phi, sigma, mu_H, t) RESULT(vv)
    IMPLICIT NONE
    INTEGER,                             INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: mu_phi, sigma, mu_H, t
    REAL(KIND=8)                                      :: vv

    vv = 0.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; rd=mu_phi; rd=sigma; rd=mu_H; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION Eexact_gauss

  !===Initialization of magnetic field and scalar potential (if present)
  MODULE SUBROUTINE init_maxwell(H_mesh, phi_mesh, time, dt, mu_H_field, mu_phi, &
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
          ENDIF
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
          ENDIF
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE init_maxwell

  !===Analytical permeability (if needed)
  !===This function is not needed unless the flag
  !===     ===Use FEM Interpolation for magnetic permeability  (true/false)
  !===is activated and set to .FALSE. in the data data file. Default is .TRUE.
  MODULE FUNCTION mu_bar_in_fourier_space(H_mesh,nb,ne,pts,pts_ids) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type), INTENT(IN)                :: H_mesh
    REAL(KIND=8), DIMENSION(ne-nb+1)           :: vv
    INTEGER,     INTENT(IN)                    :: nb, ne
    REAL(KIND=8),DIMENSION(2,ne-nb+1),OPTIONAL :: pts
    INTEGER,     DIMENSION(ne-nb+1), OPTIONAL  :: pts_ids

    vv = 1.d0
    CALL error_petsc('mu_bar_in_fourier_space: should not be called for this test')
    RETURN
    !===Dummy variables to avoid warning
    nd=H_mesh%np; nd=nb; nd=ne
    IF (PRESENT(pts)) nd=SIZE(pts,1)
    IF (PRESENT(pts_ids)) nd=SIZE(pts_ids)
    !===Dummy variables to avoid warning
  END FUNCTION mu_bar_in_fourier_space

  !===Analytical mu_in_fourier_space (if needed)
  !===This function is not needed unless the flag
  !===     ===Use FEM Interpolation for magnetic permeability  (true/false)
  !===is activated and set to .FALSE. in the data data file. Default is .TRUE.
  MODULE FUNCTION grad_mu_bar_in_fourier_space(pt,pt_id) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8),DIMENSION(2), INTENT(in):: pt
    INTEGER,DIMENSION(1), INTENT(in)     :: pt_id
    REAL(KIND=8),DIMENSION(2)            :: vv

    vv = 0.d0
    CALL error_petsc('grad_mu_bar_in_fourier_space: should not be called for this test')
    RETURN
    !===Dummy variables to avoid warning
    nd=SIZE(pt,1); nd=SIZE(pt_id)
    !===Dummy variables to avoid warning
  END FUNCTION grad_mu_bar_in_fourier_space

  !===Analytical permeability, mu in real space (if needed)
  MODULE FUNCTION mu_in_real_space(H_mesh,angles,nb_angles,nb,ne,time) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type), INTENT(IN)                :: H_mesh
    REAL(KIND=8), DIMENSION(:), INTENT(IN)     :: angles
    INTEGER, INTENT(IN)                        :: nb_angles
    INTEGER, INTENT(IN)                        :: nb, ne
    REAL(KIND=8), INTENT(IN)                   :: time
    REAL(KIND=8), DIMENSION(nb_angles,ne-nb+1) :: vv

    vv = 1.d0
    CALL error_petsc('mu_in_real_space: should not be called for this test')
    RETURN
    !===Dummy variables to avoid warning
    nd=H_mesh%np; nd=SIZE(angles); rd=time
    !===Dummy variables to avoid warning
  END FUNCTION mu_in_real_space

  MODULE FUNCTION sigma_bar_in_fourier_space(H_mesh) RESULT(vv)
    TYPE(mesh_type), INTENT(IN)                :: H_mesh
    REAL(KIND=8), DIMENSION(SIZE(H_mesh%rr,2)) :: vv

    vv = 1.d0
    CALL error_petsc('sigma_bar_in_fourier_space: should not be called for this test')
    RETURN
  END FUNCTION sigma_bar_in_fourier_space

  !===Coefficient contaning the magnetic susceptibility for magnetic force in ferrofluids:
  !===F = chi_coeff(T) * grad(H**2/2) (Kelvin form)
  !===or F = -H**2/2 * grad(chi_coeff(T)) (Helmholtz form)
  MODULE FUNCTION chi_coeff_law(temp) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: temp
    REAL(KIND=8) :: vv

    vv = - temp**2
    RETURN
  END FUNCTION chi_coeff_law

  !===Coefficient contaning the temperature dependant factor in the
  !===pyromagnetic coefficient term of the temperature equation for ferrofluids:
  !===T * dchi/dT(T) * D/Dt(H**2/2)
  MODULE FUNCTION T_dchi_dT_coeff_law(temp) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: temp
    REAL(KIND=8) :: vv

    vv = 0.d0*temp
    RETURN
  END FUNCTION T_dchi_dT_coeff_law

  MODULE FUNCTION nu_tilde_law(temp) RESULT(vv)
    REAL(KIND=8) :: temp
    REAL(KIND=8) :: vv

    vv = 0.d0*temp
    CALL error_petsc('nu_tilde_law: should not be called for this test')
    RETURN
  END FUNCTION nu_tilde_law

  MODULE FUNCTION rot_H_jump_interface(mesh,rr,list_mode) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type)                                   :: mesh
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,      DIMENSION(:),          INTENT(IN)   :: list_mode
    REAL(KIND=8), DIMENSION(SIZE(rr,2),6,SIZE(list_mode)) :: vv

    vv = 0.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=mesh%np; nd=SIZE(rr,1); nd=SIZE(list_mode)
    !===Dummy variables to avoid warning
  END FUNCTION rot_H_jump_interface

  MODULE FUNCTION Derivative_of_potential_from_rhoLi(delta_rhoLi_phys) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: delta_rhoLi_phys
    REAL(KIND=8) :: vv

    vv = 0.d0*delta_rhoLi_phys
    RETURN
  END FUNCTION Derivative_of_potential_from_rhoLi

  MODULE FUNCTION molar_fraction_from_concentration(delta_rhoLi_phys) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: delta_rhoLi_phys
    REAL(KIND=8) :: vv

    vv = 0.d0*delta_rhoLi_phys
    RETURN
  END FUNCTION molar_fraction_from_concentration

END SUBMODULE BOUNDARY_GENERIC
