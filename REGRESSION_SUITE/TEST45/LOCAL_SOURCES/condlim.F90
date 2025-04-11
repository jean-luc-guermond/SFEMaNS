SUBMODULE (boundary_generic_module) boundary_generic
  USE my_util
  USE def_type_mesh
  USE input_data
  USE bessel
  USE user_data
  REAL(KIND=8),  PARAMETER:: a=2.d0
  REAL(KIND=8),  PARAMETER:: amp=1.d0
  REAL(KIND=8) :: PI=ACOS(-1.d0)
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
    REAL(KIND=8), DIMENSION(SIZE(rr,2))                  :: vv
    REAL(KIND=8), DIMENSION(SIZE(rr,2))                  :: r, z
    INTEGER                                              :: m
    REAL(KIND=8)                                         :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))                  :: ft, fd, fnl, fp, ftemp
    REAL(KIND=8)                                         :: rho1, rho2, eta1, eta2, grav1, grav2

    r = rr(1,:)
    z = rr(2,:)
    m = mode
    t = time

    ft = 0.d0  !source term for time derivative
    fd = 0.d0  !source term for dissipation term
    fnl= 0.d0  !source term for nonlinear term DIV(mxu)
    fp = 0.d0  !source term for pressure gradient
    ftemp = 0.d0  !source term for pressure gradient

    vv = 0.d0*Re  !output=sum all source terms

    IF (inputs%if_level_set) THEN
       rho1=inputs%density_fluid(1)
       rho2=inputs%density_fluid(2)-inputs%density_fluid(1)
       eta1=inputs%dyna_visc_fluid(1)
       eta2=inputs%dyna_visc_fluid(2)-inputs%dyna_visc_fluid(1)
       grav1=inputs%heat_grav_fluid(1)*inputs%density_fluid(1)
       grav2=inputs%heat_grav_fluid(2)*inputs%density_fluid(2)-grav1
    ELSE
       CALL error_petsc('BUG in source_in_NS_momentum, level set should be present for this test')
    END IF

    !Compute ft
    IF (m==0 .AND. TYPE==1) THEN
       ft = 0.25*r**3*Cos(2.*PI*z)*(Cos(t)*(4.*rho1 + 2.*rho2 + rho2*Cos(2.*PI*(t - 1.*z))) &
            - 2.*PI*rho2*Sin(t)*Sin(2.*PI*(t - 1.*z)))
    ELSE IF (m==0 .AND. TYPE==3) THEN
       ft = 0.25*r**2*(4.*rho1*Cos(t) + rho2*Cos(t)*(2. + Cos(2.*PI*(t - 1.*z))) &
            -  2.*PI*rho2*Sin(t)*Sin(2.*PI*(t - 1.*z)))*Sin(2.*PI*z)
    ELSE IF (m==0 .AND. TYPE==5) THEN
       ft = (-0.5*r**2*(Cos(t)*(4.*rho1 + 2.*rho2 + rho2*Cos(2.*PI*(t - 1.*z))) &
            - 2.*PI*rho2*Sin(t)*Sin(2.*PI*(t - 1.*z)))*Sin(2.*PI*z))/PI
    ELSE IF (m==2 .AND. TYPE==2) THEN
       ft = 1.*r**5*rho2*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Cos(2.*PI*z) &
            - 2.*PI*r**5*rho2*Cos(2.*PI*z)*Sin(t)*Sin(0. + 2.*PI*(t - 1.*z))
    ELSE IF (m==2 .AND. TYPE==4) THEN
       ft =  1.*r**4*rho2*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Sin(2.*PI*z) &
            - 2.*PI*r**4*rho2*Sin(t)*Sin(0. + 2.*PI*(t - 1.*z))*Sin(2.*PI*z)
    ELSE IF (m==2 .AND. TYPE==6) THEN
       ft =  (-2.*r**4*rho2*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Sin(2.*PI*z))/PI &
            +  4.*r**4*rho2*Sin(t)*Sin(0. + 2.*PI*(t - 1.*z))*Sin(2.*PI*z)
    ELSE
       ft = 0.d0
    END IF

    !Compute fnl
    IF (m==0 .AND. TYPE==1) THEN
       fnl = 0.125*r**3*Sin(t)**2*(6.*r**2* (4.*rho1 + 2.*rho2 + rho2*Cos(2.*PI*(t - 1.*z)))*Cos(2.*PI*z)**2 &
            + Sin(2.*PI*z)*(-1.*rho2*Sin(2.*PI*t) + (rho2 - 8.*r**2*rho2)*Sin(2.*PI*(t - 2.*z)) &
            +  4.*(-1. + 4.*r**2)*(2.*rho1 + rho2)*Sin(2.*PI*z)))
    ELSE IF (m==2 .AND. TYPE==1) THEN
       fnl = 1.*r**6*rho2*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t)**2*Sin(0. + 4.*PI*z)
    ELSE IF (m==2 .AND. TYPE==2) THEN
       fnl =  - 0.5*r**5*rho2*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t)**2 &
            +  4.5*r**7*rho2*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t)**2 &
            +  0.5*r**5*rho2*Cos(0. + 2.*PI*(t - 1.*z))*Cos(0. + 4.*PI*z)*Sin(t)**2 &
            +  0.5*r**7*rho2*Cos(0. + 2.*PI*(t - 1.*z))*Cos(0. + 4.*PI*z)*Sin(t)**2 &
            -  2.*r**7*rho2*Sin(t)**2*Sin(0. + 2.*PI*(t - 1.*z))*Sin(0. + 4.*PI*z)
    ELSE IF (m==0 .AND. TYPE==3) THEN
       fnl = -0.25*r**4*Sin(t)**2*Sin(2.*PI*z)* ((4.*rho1 + 2.*rho2 + rho2*Cos(2.*PI*(t - 1.*z)))*Cos(2.*PI*z) &
            + 4.*rho2*Sin(2.*PI*(t - 1.*z))*Sin(2.*PI*z))
    ELSE IF (m==2 .AND. TYPE==3) THEN
       fnl = 2.*r**5*rho2*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t)**2*Sin(2.*PI*z)*Sin(0. + 2.*PI*z)
    ELSE IF (m==2 .AND. TYPE==4) THEN
       fnl = 1.*r**6*rho2*Cos(0. + 2.*PI*(t - 1.*z))*Cos(0. + 2.*PI*z)*Sin(t)**2* Sin(2.*PI*z) &
            - 4.*r**6*rho2*Sin(t)**2*Sin(0. + 2.*PI*(t - 1.*z))* Sin(2.*PI*z)*Sin(0. + 2.*PI*z)
    ELSE IF (m==0 .AND. TYPE==5) THEN
       fnl = (r**4*Sin(t)**2*Sin(2.*PI*z)*((4.*rho1 + 2.*rho2 + rho2*Cos(2.*PI*(t - 1.*z)))*Cos(2.*PI*z) &
            + 2.*rho2*Sin(2.*PI*(t - 1.*z))*Sin(2.*PI*z)))/PI
    ELSE IF (m==2 .AND. TYPE==5) THEN
       fnl = (-4.*r**5*rho2*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t)**2*Sin(2.*PI*z)**2)/PI
    ELSE IF (m==2 .AND. TYPE==6) THEN
       fnl = (8.*r**6*rho2*Sin(t)**2*Sin(0. + 2.*PI*(t - 1.*z))*Sin(2.*PI*z)**2)/PI
    ELSE
       fnl = 0.d0
    END IF

    !Compute fd
    IF (m==0 .AND. TYPE==1) THEN
       fd = r*(-2.*eta2*Cos(2.*PI*t) + eta2*PI**2*r**2*Cos(2.*PI*(t - 2.*z)) &
            +2.*(2.*eta1 + eta2)*(-2. + PI**2*r**2)*Cos(2.*PI*z))*Sin(t)
    ELSE IF (m==2 .AND. TYPE==1) THEN
       fd = -2.*eta2*r**2*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t)*Sin(0. + 2.*PI*z)
    ELSE IF (m==2 .AND. TYPE==2) THEN
       fd = -20.*eta2*r**3*Cos(0. + 2.*PI*(t - 1.*z))*Cos(0. + 2.*PI*z)*Sin(t) &
            + 4.*eta2*PI**2*r**5*Cos(0. + 2.*PI*(t - 1.*z))*Cos(0. + 2.*PI*z)*Sin(t) &
            + 8.*eta2*r**3*Sin(t)*Sin(0. + 2.*PI*(t - 1.*z))*Sin(0. + 2.*PI*z) &
            + 4.*eta2*PI**2*r**5*Sin(t)*Sin(0. + 2.*PI*(t - 1.*z))*Sin(0. + 2.*PI*z)
    ELSE IF (m==0 .AND. TYPE==3) THEN
       fd = -0.125*Sin(t)*(3.*eta2*Sin(2.*PI*t) &
            +eta2*(-3. + 8.*PI**2*r**2)*Sin(2.*PI*(t - 2.*z)) &
            -4.*(2.*eta1 + eta2)*(-3. + 4.*PI**2*r**2)*Sin(2.*PI*z))
    ELSE IF (m==2 .AND. TYPE==3) THEN
       fd = -4.*eta2*r**3*Cos(0. + 2.*PI*(t - 1.*z))*Cos(0. + 2.*PI*z)*Sin(t)
    ELSE IF (m==2 .AND. TYPE==4) THEN
       fd = -4.*eta2*PI**2*r**4*Cos(0. + 2.*PI*z)*Sin(t)*Sin(0. + 2.*PI*(t - 1.*z)) &
            -5.*eta2*r**2*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t)*Sin(0. + 2.*PI*z) &
            +4.*eta2*PI**2*r**4*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t)*Sin(0. + 2.*PI*z)
    ELSE IF (m==0 .AND. TYPE==5) THEN
       fd = 4.*eta2*PI*r**2*Cos(2.*PI*z)*Sin(t)*Sin(2.*PI*(t - 1.*z)) &
            -(2.*(-1. + PI**2*r**2)*(4.*eta1 + 2.*eta2 + eta2*Cos(2.*PI*(t - 1.*z)))*Sin(t)*Sin(2.*PI*z))/PI
    ELSE IF (m==2 .AND. TYPE==5) THEN
       fd = -4.*eta2*PI*r**3*Cos(0. + 2.*PI*(t - 1.*z))*Cos(0. + 2.*PI*z)*Sin(t)
    ELSE IF (m==2 .AND. TYPE==6) THEN
       fd = 16.*eta2*PI*r**4*Cos(0. + 2.*PI*z)*Sin(t)*Sin(0. + 2.*PI*(t - 1.*z)) &
            + (16.*eta2*r**2*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t)*Sin(0. + 2.*PI*z))/PI &
            -4.*eta2*PI*r**4*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t)*Sin(0. + 2.*PI*z)
    ELSE
       fd = 0.d0
    END IF
    fd = fd/inputs%Re

    !Compute fp
    IF (m==0 .AND. TYPE==1) THEN
       fp =(3.*Sin(t)**2*(1.*PI**2*r**5*Cos(2*PI*z)**2 &
            +(2.6666666666666665 + 0.6666666666666666*PI**2)*r**3*Sin(2*PI*z)**2))/ PI**2
    ELSE IF (m==0 .AND. TYPE==5) THEN
       fp = (-1.*r**4*(-4. + PI**2*(-1. + r**2))*Sin(t)**2*Sin(4*PI*z))/PI
    ELSE
       fp = 0.d0
    END IF

    !Compute ftemp
    ftemp = 0.d0
    IF (m==0 .AND. TYPE==5) THEN
       ftemp = r**2*Cos(t)*((1.*grav1 + 0.5*grav2)*Sin(2.*PI*z) &
            +Cos(2.*PI*(t - 1.*z))*(grav2*r**2*Cos(2.*PI*z) + 0.25*grav2*Sin(2.*PI*z)))
    ELSE IF (m==1 .AND. TYPE==5) THEN
       ftemp = 1.*grav1*r**2*Cos(t)*Cos(2.*PI*z) + 0.5*grav2*r**2*Cos(t)*Cos(2.*PI*z) &
            + 0.25*grav2*r**2*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Cos(2.*PI*z)
    ELSE IF (m==1 .AND. TYPE==6) THEN
       ftemp = 0.5*grav2*r**4*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Cos(2.*PI*z)
    ELSE IF (m==2 .AND. TYPE==6) THEN
       ftemp = 2.*grav1*r**2*Cos(t)*Cos(0. + 2.*PI*z) &
            + 1.*grav2*r**2*Cos(t)*Cos(0. + 2.*PI*z) &
            + 0.5*grav2*r**2*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Cos(0. + 2.*PI*z) &
            + 1.*grav2*r**4*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Sin(0. + 2.*PI*z)
    ELSE IF (m==3 .AND. TYPE==6) THEN
       ftemp = 0.5*grav2*r**4*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Cos(2.*PI*z)
    ELSE IF (m==4 .AND. TYPE==5) THEN
       ftemp = -1.*grav2*r**4*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Cos(2.*PI*z)
    ELSE
       ftemp = 0.d0
    END IF

    !Sum all source terms
    vv = ft + fd + fnl + fp - ftemp
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
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: r, z
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: ft, fd, fnl
    REAL(KIND=8)                                      :: rho1, rho2, kappa1, kappa2, crho1, crho2

    r = rr(1,:)
    z = rr(2,:)

    ft = 0.d0 ! time derivative of internal energy density
    fd = 0.d0 ! div( k*grad(T) )
    fnl = 0.d0 ! u_grad_e

    vv = 0.d0

    IF (inputs%if_level_set) THEN
       rho1 =inputs%density_fluid(1)
       rho2 =inputs%density_fluid(2)-inputs%density_fluid(1)
       kappa1 =inputs%heat_diffu_fluid(1)
       kappa2 =inputs%heat_diffu_fluid(2)-inputs%heat_diffu_fluid(1)
       crho1=inputs%heat_capacity_fluid(1)*inputs%density_fluid(1)
       crho2=inputs%heat_capacity_fluid(2)*inputs%density_fluid(2) &
            -inputs%heat_capacity_fluid(1)*inputs%density_fluid(1)
    ELSE
       CALL error_petsc('BUG in source_in_temperature, level set should be present for this test')
    END IF

    !Compute ft
    IF (m==0 .AND. TYPE==1) THEN
       ft = r**2*(crho2*PI*Cos(t)*Sin(2.*PI*(t - 1.*z))* (-2.*r**2*Cos(2.*PI*z) - 0.5*Sin(2.*PI*z)) &
            +crho2*Cos(2.*PI*(t - 1.*z))*Sin(t)*(-1.*r**2*Cos(2.*PI*z) - 0.25*Sin(2.*PI*z)) &
            +(-1.*crho1 - 0.5*crho2)*Sin(t)*Sin(2.*PI*z))
    ELSE IF (m==1 .AND. TYPE==1) THEN
       ft = - 1.*crho1*r**2*Cos(2.*PI*z)*Sin(t) - 0.5*crho2*r**2*Cos(2.*PI*z)*Sin(t) &
            - 0.25*crho2*r**2*Cos(0. + 2.*PI*(t - 1.*z))*Cos(2.*PI*z)*Sin(t) &
            -  0.5*crho2*PI*r**2*Cos(t)*Cos(2.*PI*z)*Sin(0. + 2.*PI*(t - 1.*z))
    ELSE IF (m==1 .AND. TYPE==2) THEN
       ft = -0.5*crho2*r**4*Cos(0. + 2.*PI*(t - 1.*z))*Cos(2.*PI*z)*Sin(t) &
            -  1.*crho2*PI*r**4*Cos(t)*Cos(2.*PI*z)*Sin(0. + 2.*PI*(t - 1.*z))
    ELSE IF (m==2 .AND. TYPE==2) THEN
       ft = -2.*crho1*r**2*Cos(0. + 2.*PI*z)*Sin(t) &
            -1.*crho2*r**2*Cos(0. + 2.*PI*z)*Sin(t) &
            -0.5*crho2*r**2*Cos(0. + 2.*PI*(t - 1.*z))*Cos(0. + 2.*PI*z)*Sin(t) &
            -1.*crho2*PI*r**2*Cos(t)*Cos(0. + 2.*PI*z)*Sin(0. + 2.*PI*(t - 1.*z)) &
            -1.*crho2*r**4*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t)*Sin(0. + 2.*PI*z) &
            -2.*crho2*PI*r**4*Cos(t)*Sin(0. + 2.*PI*(t - 1.*z))*Sin(0. + 2.*PI*z)
    ELSE IF (m==3 .AND. TYPE==2) THEN
       ft = -0.5*crho2*r**4*Cos(0. + 2.*PI*(t - 1.*z))*Cos(2.*PI*z)*Sin(t) &
            -1.*crho2*PI*r**4*Cos(t)*Cos(2.*PI*z)*Sin(0. + 2.*PI*(t - 1.*z))
    ELSE IF (m==4 .AND. TYPE==1) THEN
       ft = 1.*crho2*r**4*Cos(0. + 2.*PI*(t - 1.*z))*Cos(2.*PI*z)*Sin(t) &
            +2.*crho2*PI*r**4*Cos(t)*Cos(2.*PI*z)*Sin(0. + 2.*PI*(t - 1.*z))
    ELSE
       ft = 0.d0
    END IF

    !Compute fd
    IF (m==0 .AND. TYPE==1) THEN
       fd = Cos(t)*(kappa2*PI**2*r**2*Cos(2.*PI*z)*Sin(2.*PI*(t - 1.*z)) &
            +2.*(-1.*(2.*kappa1 + kappa2)*(-1. + PI**2*r**2) &
            -2.*kappa2*PI**2*r**4*Sin(2.*PI*(t - 1.*z)))*Sin(2.*PI*z) &
            -1.*kappa2*Cos(2.*PI*(t - 1.*z))* (4.*r**2*(-2. + PI**2*r**2)*Cos(2.*PI*z) &
            +(-1. + PI**2*r**2)*Sin(2.*PI*z)))
    ELSE IF (m==1 .AND. TYPE==1) THEN
       fd =  3.*kappa1*Cos(t)*Cos(0. + 2.*PI*z) + 1.5*kappa2*Cos(t)*Cos(0. + 2.*PI*z) &
            -4.*kappa1*PI**2*r**2*Cos(t)*Cos(0. + 2.*PI*z) &
            -2.*kappa2*PI**2*r**2*Cos(t)*Cos(0. + 2.*PI*z) &
            +0.75*kappa2*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Cos(0. + 2.*PI*z) &
            -1.*kappa2*PI**2*r**2*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Cos(0. + 2.*PI*z) &
            -1.*kappa2*PI**2*r**2*Cos(t)*Sin(0. + 2.*PI*(t - 1.*z))*Sin(0. + 2.*PI*z)
    ELSE IF (m==1 .AND. TYPE==2) THEN
       fd = 4.5*kappa2*r**2*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Cos(0. + 2.*PI*z) &
            -2.*kappa2*PI**2*r**4*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Cos(0. + 2.*PI*z) &
            -2.*kappa2*PI**2*r**4*Cos(t)*Sin(0. + 2.*PI*(t - 1.*z))*Sin(0. + 2.*PI*z)
    ELSE IF (m==2 .AND. TYPE==2) THEN
       fd = -2.*kappa2*PI**2*r**2*Cos(t)*Cos(0. + 2.*PI*(t - 2.*z)) &
            -8.*kappa1*PI**2*r**2*Cos(t)*Cos(0. + 2.*PI*z) &
            -4.*kappa2*PI**2*r**2*Cos(t)*Cos(0. + 2.*PI*z) &
            +4.*kappa2*r**2*Cos(t)*Sin(0. + 2.*PI*t) &
            -4.*kappa2*r**2*Cos(t)*Sin(0. + 2.*PI*(t - 2.*z)) &
            +4.*kappa2*PI**2*r**4*Cos(t)*Sin(0. + 2.*PI*(t - 2.*z))
    ELSE IF (m==3 .AND. TYPE==2) THEN
       fd = 1.25*kappa2*r**2*Cos(t)*Cos(0. + 2.*PI*t) &
            +1.25*kappa2*r**2*Cos(t)*Cos(0. + 2.*PI*(t - 2.*z)) &
            -2.*kappa2*PI**2*r**4*Cos(t)*Cos(0. + 2.*PI*(t - 2.*z))
    ELSE IF (m==4 .AND. TYPE==1) THEN
       fd = 4.*kappa2*PI**2*r**4*Cos(t)*Cos(0. + 2.*PI*(t - 2.*z))
    ELSE
       fd = 0.d0
    END IF

    !Compute fnl
    IF (m==0 .AND. TYPE==1) THEN
       fnl =  r**4*(Sin(2.*PI*z)*(-1.*crho1*Cos(2.*PI*z)*Sin(2.*t) &
            +crho2*Cos(t)*Sin(t)*(Cos(2.*PI*z)*(-1. - 4.*r**2*Sin(2.*PI*(t - 1.*z))) &
            -1.*Sin(2.*PI*(t - 1.*z))*Sin(2.*PI*z))) &
            +crho2*Cos(2.*PI*(t - 1.*z))*Sin(2.*t)*&
            (2.*r**2 - 0.25*Cos((2.0)*PI*z)*Sin((2.0)*PI*z)))
    ELSE IF (m==1 .AND. TYPE==1) THEN
       fnl = 2.*crho1*r**4*Cos(t)*Cos(0. + 2.*PI*z)**2*Sin(t) &
            +1.*crho2*r**4*Cos(t)*Cos(0. + 2.*PI*z)**2*Sin(t) &
            +0.5*crho2*r**4*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Cos(0. + 2.*PI*z)**2*Sin(t) &
            + 0.5*crho2*r**5*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Cos(0. + 2.*PI*z)*Sin(t)*Sin(0. + 2.*PI*z) &
            -1.*crho2*r**4*Cos(t)*Cos(0. + 2.*PI*z)*Sin(t)*Sin(0. + 2.*PI*(t - 1.*z))*Sin(0. + 2.*PI*z) &
            + 4.*crho1*r**4*Cos(t)*Sin(t)*Sin(0. + 2.*PI*z)**2 &
            +2.*crho2*r**4*Cos(t)*Sin(t)*Sin(0. + 2.*PI*z)**2 &
            +1.*crho2*r**4*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t)*Sin(0. + 2.*PI*z)**2
    ELSE IF (m==1 .AND. TYPE==2) THEN
       fnl = 2.*crho2*r**6*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Cos(0. + 2.*PI*z)**2*Sin(t) &
            -1.*crho1*r**3*Cos(t)*Cos(0. + 2.*PI*z)*Sin(t)*Sin(0. + 2.*PI*z) &
            -0.5*crho2*r**3*Cos(t)*Cos(0. + 2.*PI*z)*Sin(t)*Sin(0. + 2.*PI*z) &
            -0.25*crho2*r**3*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Cos(0. + 2.*PI*z)*Sin(t)*Sin(0. + 2.*PI*z) &
            -2.*crho2*r**6*Cos(t)*Cos(0. + 2.*PI*z)*Sin(t)*Sin(0. + 2.*PI*(t - 1.*z))*Sin(0. + 2.*PI*z) &
            +2.*crho2*r**6*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t)*Sin(0. + 2.*PI*z)**2
    ELSE IF (m==2 .AND. TYPE==1) THEN
       fnl = 2.*crho2*r**5*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t)* Sin(0. + 2.*PI*z)**2 &
            +2.*crho1*r**3*Cos(t)*Sin(t)*Sin(0. + 4.*PI*z) &
            +1.*crho2*r**3*Cos(t)*Sin(t)*Sin(0. + 4.*PI*z) &
            +0.5*crho2*r**3*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t)*Sin(0. + 4.*PI*z)
    ELSE IF (m==2 .AND. TYPE==2) THEN
       fnl = 4.*crho1*r**4*Cos(t)*Cos(0. + 2.*PI*z)**2*Sin(t) &
            +2.*crho2*r**4*Cos(t)*Cos(0. + 2.*PI*z)**2*Sin(t) &
            +1.*crho2*r**4*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Cos(0. + 2.*PI*z)**2*Sin(t) &
            +8.*crho1*r**4*Cos(t)*Sin(t)*Sin(0. + 2.*PI*z)**2 &
            +4.*crho2*r**4*Cos(t)*Sin(t)*Sin(0. + 2.*PI*z)**2 &
            +2.*crho2*r**4*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t)*Sin(0. + 2.*PI*z)**2 &
            -4.*crho2*r**6*Cos(t)*Sin(t)*Sin(0. + 2.*PI*(t - 1.*z))*Sin(0. + 2.*PI*z)**2 &
            -1.*crho2*r**4*Cos(t)*Sin(t)*Sin(0. + 2.*PI*(t - 1.*z))*Sin(0. + 4.*PI*z)
    ELSE IF (m==3 .AND. TYPE==1) THEN
       fnl = 0.75*crho2*r**5*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t)*Sin(0. + 4.*PI*z)
    ELSE IF (m==3 .AND. TYPE==2) THEN
       fnl = 2.*crho2*r**6*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t) &
            -1.*crho2*r**6*Cos(t)*Sin(t)*Sin(0. + 2.*PI*(t - 1.*z))*Sin(0. + 4.*PI*z)
    ELSE IF (m==4 .AND. TYPE==1) THEN
       fnl = -4.*crho2*r**6*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t) &
            +2.*crho2*r**6*Cos(t)*Sin(t)*Sin(0. + 2.*PI*(t - 1.*z))*Sin(0. + 4.*PI*z)
    ELSE IF (m==4 .AND. TYPE==2) THEN
       fnl = 2.*crho2*r**5*Cos(t)*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t)*Sin(0. + 4.*PI*z)
    ELSE
       fnl = 0.d0
    END IF

    !Sum up all terms
    vv = ft - fd + fnl
    RETURN

    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION source_in_temperature

  !===Extra source in level set equation. Always called.
  MODULE FUNCTION source_in_level_set(interface_nb,TYPE, rr, m, t)RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m, interface_nb
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: r, z

    r = rr(1,:)
    z = rr(2,:)

    vv = 0.d0

    IF (m==0 .AND. TYPE==1) THEN
       vv = -0.5*Sin(2.*PI*(t - 1.*z))*(PI + 2.*r**2*Sin(t)*Sin(2.*PI*z))
    ELSE IF (m==2 .AND. TYPE==1) THEN
       vv = 2.*r**3*Cos(0. + 2.*PI*(t - 1.*z))*Sin(t)*Sin(0. + 2.*PI*z)
    ELSE IF (m==2 .AND. TYPE==2) THEN
       vv = 2.*r**4*Cos(0. + 2.*PI*(t - 1.*z))*Cos(0. + 2.*PI*z)*Sin(t) &
            - 2.*PI*r**2*Sin(0. + 2.*PI*(t - 1.*z)) &
            - 4.*r**4*Sin(t)*Sin(0. + 2.*PI*(t - 1.*z))*Sin(0. + 2.*PI*z)
    ELSE
       vv = 0.d0
    END IF
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
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: r, z

    r = rr(1,:)
    z = rr(2,:)

    vv = 0.d0
    IF (m==0 .AND. TYPE==1) THEN
       vv = r**3*Cos(2*PI*z)*Sin(t)
    ELSE IF (m==0 .AND. TYPE==3) THEN
       vv= r**2*Sin(t)*Sin(2*PI*z)
    ELSE IF (m==0 .AND. TYPE==5) THEN
       vv = (-2.*r**2*Sin(t)*Sin(2*PI*z))/PI
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
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: r, z

    r = rr(1,:)
    z = rr(2,:)

    IF (m==0.AND.TYPE==1) THEN
       vv = 0.5*Sin(t)**2*(r**6*Cos(2*PI*z)**2 &
            + r**4*Sin(2*PI*z)**2 &
            + (4*r**4*Sin(2*PI*z)**2)/PI**2)
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
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: r, z

    r = rr(1,:)
    z = rr(2,:)

    IF (m==0 .AND. TYPE==1) THEN
       vv = r**2*Cos(t)*Sin(2.*PI*z)
    ELSE IF (m==1 .AND. TYPE==1) THEN
       vv = 1.*r**2*Cos(t)*Cos(2.*PI*z)
    ELSE IF (m==2 .AND. TYPE==2) THEN
       vv = 2.*r**2*Cos(t)*Cos(2.*PI*z)
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
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: r, z

    r = rr(1,:)
    z = rr(2,:)

    IF (interface_nb==1) THEN
       IF (m==0 .AND. TYPE==1) THEN
          vv = 0.25*(2. + Cos(2.*PI*(t - 1.*z)))
       ELSE IF (m==2 .AND. TYPE==2) THEN
          vv = 1.*r**2*Cos(0. + 2.*PI*(t - 1.*z))
       ELSE
          vv = 0.d0
       END IF
    ELSE
       CALL error_petsc(' BUG in level_set_exact, we should compute only 1 level set')
    END IF
    RETURN

    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; nd=interface_nb; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION level_set_exact

  !===Penalty coefficient (if needed)
  !===This coefficient is equal to zero in subdomain
  !===where penalty is applied (penalty is zero in solid)
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
    CALL error_petsc('penal_in_real_space: should not be called for this test')
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
    REAL(KIND=8), DIMENSION(H_mesh%np)                :: vv

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
    CALL error_petsc('Vexact: should not be called for this test')
    RETURN
    !===Dummy variables to avoid warning
    nd=H_mesh%np; nd=m
    !===Dummy variables to avoid warning
  END FUNCTION Vexact

  !===Magnetic field and magnetic induction for quasi-static approximation
  !===if needed
  MODULE FUNCTION H_B_quasi_static(char_h_b, rr, m) RESULT(vv)
    IMPLICIT NONE
    CHARACTER(LEN=1),                    INTENT(IN)   :: char_h_b
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m
    REAL(KIND=8), DIMENSION(SIZE(rr,2),6)             :: vv

    vv = 0.d0
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
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    CALL error_petsc('Hexact: should not be called for this test')
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
    CALL error_petsc('Phiexact: should not be called for this test')
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
    REAL(KIND=8)                                      :: vv

    vv = 0.d0
    CALL error_petsc('Jexact_gauss: should not be called for this test')
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
    CALL error_petsc('Eexact_gauss: should not be called for this test')
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
    IMPLICIT NONE
    TYPE(mesh_type), INTENT(IN)                :: H_mesh
    REAL(KIND=8), DIMENSION(SIZE(H_mesh%rr,2)) :: vv

    vv = 1.d0
    CALL error_petsc('sigma_bar_in_fourier_space: should not be called for this test')
    RETURN
  END FUNCTION sigma_bar_in_fourier_space

  MODULE FUNCTION chi_coeff_law(temp) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: temp
    REAL(KIND=8) :: vv

    vv = 0.d0*temp
    CALL error_petsc('chi_coeff_law: should not be called for this test')
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
    CALL error_petsc('T_dchi_dT_coeff_law: should not be called for this test')
    RETURN
  END FUNCTION T_dchi_dT_coeff_law

  !===Kinematic viscosity's law of temperature
  MODULE FUNCTION nu_tilde_law(temp) RESULT(vv)
    IMPLICIT NONE
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

  MODULE FUNCTION curved_boundary_radius(interface, theta) RESULT(vv)
    INTEGER, INTENT(IN) :: interface
    REAL(KIND = 8), INTENT(IN) :: theta
    REAL(KIND = 8) :: vv

    vv = 0.d0
    RETURN
  END FUNCTION curved_boundary_radius

END SUBMODULE BOUNDARY_GENERIC
