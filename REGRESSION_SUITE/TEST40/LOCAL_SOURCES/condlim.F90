SUBMODULE (boundary_generic_module) boundary_generic  
  USE my_util
  USE def_type_mesh
  USE input_data
  USE bessel
  USE user_data
  REAL(KIND=8),  PARAMETER:: a=2.d0
  REAL(KIND=8),  PARAMETER:: amp=1.d0
  
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

  !===Initialize level_set
  SUBROUTINE init_level_set(pp_mesh, time, &
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
    REAL(KIND=8), DIMENSION(SIZE(rr,2))                  :: vv
    REAL(KIND=8), DIMENSION(SIZE(rr,2))                  :: r, z
    INTEGER                                              :: m
    REAL(KIND=8)                                         :: t
    CHARACTER(LEN=2)  :: np
    REAL(KIND=8), DIMENSION(SIZE(rr,2))                  :: ft, fd, fnl, fp
    REAL(KIND=8)                                         :: rho1, rho2, eta1, eta2

    IF (PRESENT(opt_tempn)) CALL error_petsc('temperature should not be present for test 40')

    r = rr(1,:)
    z = rr(2,:)
    m = mode
    t = time

    ft = 0.d0  !source term for time derivative
    fd = 0.d0  !source term for dissipation term
    fnl= 0.d0  !source term for nonlinear term DIV(mxu)
    fp = 0.d0  !source term for pressure gradient

    vv = 0.d0*Re  !output=sum all source terms
    
    IF (inputs%if_level_set) THEN

       rho1=inputs%density_fluid(1)
       rho2=inputs%density_fluid(2)-inputs%density_fluid(1)
       eta1=inputs%dyna_visc_fluid(1)
       eta2=inputs%dyna_visc_fluid(2)-inputs%dyna_visc_fluid(1)

       !Compute ft
       IF (m==2 .AND. TYPE==1) THEN      !type 1-2
          ft =  amp*r*rho1*z**2*Cos(t - 1.*z) + 0.5d0*amp*r*rho2*z**2*Cos(t - 1.*z) &
                  + 0.25d0*amp*r*rho2*z**2*Cos(0. + 4.*amp*t - 4.*z)*Cos(t - 1.*z) &
                  - amp**2*r*rho2*z**2*Sin(0. + 4.*amp*t - 4.*z)*Sin(t - 1.*z)
       ELSE IF (m==4 .AND. TYPE==2) THEN
          ft =  0.03125d0*a*amp*r**3*rho2*z**2*Cos(0. + t + 4.*amp*t - 5.*z) &
                 + 0.125d0*a*amp**2*r**3*rho2*z**2*Cos(0. + t + 4.*amp*t - 5.*z) &
                 + 0.03125d0*a*amp*r**3*rho2*z**2*Cos(0. + t - 4.*amp*t + 3.*z) &
                 - 0.125d0*a*amp**2*r**3*rho2*z**2*Cos(0. + t - 4.*amp*t + 3.*z)
       ELSE IF (m==0 .AND. TYPE==3) THEN !type 3-4
          ft = 0.25d0*a*amp*r**3*rho2*z**2*(-0.25*Cos(4.*amp*t - 4.*z)*Cos(t - 1.*z) &
                  +amp*Sin(4.*amp*t - 4.*z)*Sin(t - 1.*z))
       ELSE IF (m==2 .AND. TYPE==4) THEN
          ft =  -amp*r*rho1*z**2*Cos(t - 1.*z) - 0.5*amp*r*rho2*z**2*Cos(t - 1.*z) &
                  - 0.25d0*amp*r*rho2*z**2*Cos(0. + 4.*amp*t - 4.*z)*Cos(t - 1.*z) &
                  + amp**2*r*rho2*z**2*Sin(0. + 4.*amp*t - 4.*z)*Sin(t - 1.*z)
       ELSE IF (m==4 .AND. TYPE==3) THEN
          ft =  0.03125d0*a*amp*r**3*rho2*z**2*Cos(0. + t + 4.*amp*t - 5.*z) &
                  + 0.125d0*a*amp**2*r**3*rho2*z**2*Cos(0. + t + 4.*amp*t - 5.*z) &
                  + 0.03125d0*a*amp*r**3*rho2*z**2*Cos(0. + t - 4.*amp*t + 3.*z) &
                  - 0.125d0*a*amp**2*r**3*rho2*z**2*Cos(0. + t - 4.*amp*t + 3.*z)
       ELSE IF (m ==0 .AND. TYPE==5) THEN !type 5-6
          ft = -amp**2*rho2*Sin(4.*amp*t - 4.*z)
       ELSE IF (m==2 .AND. TYPE==6) THEN
          ft = -0.5d0*a*amp**2*r**2*rho2*Sin(4.*amp*t - 4.*z)
       ELSE
          ft = 0.d0
       END IF

       !Compute fnl
        IF (m==0 .AND. TYPE==1) THEN      !type 1-2
          fnl = 0.25d0*amp**2*r*z**4*(4.*rho1 + 2.*rho2 + rho2*Cos(4.*amp*t - 4.*z))*Sin(t - 1.*z)**2
       ELSE IF (m==2 .AND. TYPE==1) THEN
          fnl =  - 1.d0*amp**2*r*rho1*z**2*Cos(t - 1.*z) &
                  - 0.5d0*amp**2*r*rho2*z**2*Cos(t - 1.*z) &
                  - 0.25d0*amp**2*r*rho2*z**2*Cos(0. + 4.*amp*t - 4.*z)*Cos(t - 1.*z) &
                  + 2.d0*amp**2*r*rho1*z*Sin(t - 1.*z) + 1.*amp**2*r*rho2*z*Sin(t - 1.*z) &
                  + 0.5d0*amp**2*r*rho2*z*Cos(0. + 4.*amp*t - 4.*z)*Sin(t - 1.*z) &
                  + amp**2*r*rho2*z**2*Sin(0. + 4.*amp*t - 4.*z)*Sin(t - 1.*z)
       ELSE IF (m==2 .AND. TYPE==2) THEN
          fnl = 0.125d0*a*amp**2*r**3*rho2*z**4*Cos(0. + 4.*amp*t - 4.*z)*Sin(t - 1.*z)**2
       ELSE IF (m==4 .AND. TYPE==2) THEN
          fnl = -0.15625d0*a*amp**2*r**3*rho2*z**2*Cos(0. + t + 4.*amp*t - 5.*z) &
                  + 0.09375d0*a*amp**2*r**3*rho2*z**2*Cos(0. + t - 4.*amp*t + 3.*z) &
                  + 0.125d0*a*amp**2*r**3*rho2*z*Cos(0. + 4.*amp*t - 4.*z)*Sin(t - 1.*z)
       ELSE IF (m==0 .AND. TYPE==3) THEN !type 3-4
          fnl = 0.15625d0*a*amp**2*r**3*rho2*z*(1.*z*Cos(t + 4.*amp*t - 5.*z) &
                  - 0.6*z*Cos(t - 4.*amp*t + 3.*z) - 0.8*Cos(4.*amp*t - 4.*z)*Sin(t - 1.*z))
       ELSE IF (m==2 .AND. TYPE==4) THEN
          fnl = 0.625d0*amp**2*r*rho2*z**2*Cos(0. + t + 4.*amp*t - 5.*z) &
                  + amp**2*r*rho1*z**2*Cos(t - 1.*z) &
                  + 0.5d0*amp**2*r*rho2*z**2*Cos(t - 1.*z) & 
                  - 0.375d0*amp**2*r*rho2*z**2*Cos(0. + t - 4.*amp*t + 3.*z) &
                  - 0.25d0*amp**2*r*rho2*z*Sin(0. + t + 4.*amp*t - 5.*z) &
                  - 2.d0*amp**2*r*rho1*z*Sin(t - 1.*z) - 1.*amp**2*r*rho2*z*Sin(t - 1.*z) &
                  - 0.25d0*amp**2*r*rho2*z*Sin(0. + t - 4.*amp*t + 3.*z)
       ELSE IF (m==4 .AND. TYPE==3) THEN
          fnl =  -0.15625d0*a*amp**2*r**3*rho2*z**2*Cos(0. + t + 4.*amp*t - 5.*z) &
                  + 0.09375d0*a*amp**2*r**3*rho2*z**2*Cos(0. + t - 4.*amp*t + 3.*z) &
                  + 0.125d0*a*amp**2*r**3*rho2*z*Cos(0. + 4.*amp*t - 4.*z)*Sin(t - 1.*z)
       ELSE IF (m==0 .AND. TYPE==5) THEN !type 5-6
          fnl = amp**2*rho2*Sin(4.*amp*t - 4.*z)
       ELSE IF (m==2 .AND. TYPE==6) THEN
          fnl = 0.5d0*a*amp**2*r**2*rho2*Sin(4.*amp*t - 4.*z)
       ELSE
          fnl = 0.d0
       END IF

       !Compute fd
       IF (m==2 .AND. TYPE==1) THEN      !type 1-2
          fd =  1.5*amp*eta2*r*z*Cos(0. + t + 4.*amp*t - 5.*z) &
                  + 4.d0*amp*eta1*r*z*Cos(t - 1.*z) + 2.*amp*eta2*r*z*Cos(t - 1.*z) &
                  - 0.5d0*amp*eta2*r*z*Cos(0. + t - 4.*amp*t + 3.*z) &
                  - 0.25d0*amp*eta2*r*Sin(0. + t + 4.*amp*t - 5.*z) &
                  + 0.625d0*amp*eta2*r*z**2*Sin(0. + t + 4.*amp*t - 5.*z) &
                  - 2.d0*amp*eta1*r*Sin(t - 1.*z) - 1.*amp*eta2*r*Sin(t - 1.*z) &
                  + amp*eta1*r*z**2*Sin(t - 1.*z) + 0.5*amp*eta2*r*z**2*Sin(t - 1.*z) &
                  - 0.25d0*amp*eta2*r*Sin(0. + t - 4.*amp*t + 3.*z) &
                  - 0.375d0*amp*eta2*r*z**2*Sin(0. + t - 4.*amp*t + 3.*z)
       ELSE IF (m==4 .AND. TYPE==2) THEN
          fd =  0.375*a*amp*eta2*r**3*z*Cos(0. + t + 4.*amp*t - 5.*z) &
                  - 0.125*a*amp*eta2*r**3*z*Cos(0. + t - 4.*amp*t + 3.*z) &
                  - 0.0625*a*amp*eta2*r**3*Sin(0. + t + 4.*amp*t - 5.*z) &
                  + 0.15625*a*amp*eta2*r**3*z**2*Sin(0. + t + 4.*amp*t - 5.*z) &
                  - 0.0625*a*amp*eta2*r**3*Sin(0. + t - 4.*amp*t + 3.*z) &
                  - 0.09375*a*amp*eta2*r**3*z**2*Sin(0. + t - 4.*amp*t + 3.*z)
       ELSE IF (m==0 .AND. TYPE==3) THEN !type 3-4
          fd =  -0.0625d0*a*amp*eta2*r*(r**2*z*Sin(4.*amp*t - 4.*z)* &
                  ( 4.d0*z*Cos(t - 1.*z) - 8.d0*Sin(t - 1.*z)) &
                  + Cos(4.*amp*t - 4.*z)*(4.d0*r**2*z*Cos(t - 1.*z) &
                  + (-8.d0*z**2 + r**2*(-2. + z**2))*Sin(t - 1.*z)))
       ELSE IF (m==2 .AND. TYPE==4) THEN
          fd =  -1.5d0*amp*eta2*r*z*Cos(0. + t + 4.*amp*t - 5.*z) &
                  - 4.d0*amp*eta1*r*z*Cos(t - 1.*z) - 2.*amp*eta2*r*z*Cos(t - 1.*z) &
                  + 0.5d0*amp*eta2*r*z*Cos(0. + t - 4.*amp*t + 3.*z) &
                  + 0.25d0*amp*eta2*r*Sin(0. + t + 4.*amp*t - 5.*z) &
                  - 0.625d0*amp*eta2*r*z**2*Sin(0. + t + 4.*amp*t - 5.*z) &
                  + 2.d0*amp*eta1*r*Sin(t - 1.*z) + 1.*amp*eta2*r*Sin(t - 1.*z) &
                  - amp*eta1*r*z**2*Sin(t - 1.*z) - 0.5*amp*eta2*r*z**2*Sin(t - 1.*z) &
                  + 0.25d0*amp*eta2*r*Sin(0. + t - 4.*amp*t + 3.*z) &
                  + 0.375d0*amp*eta2*r*z**2*Sin(0. + t - 4.*amp*t + 3.*z)
       ELSE IF (m==4 .AND. TYPE==3) THEN
          fd =  0.375d0*a*amp*eta2*r**3*z*Cos(0. + t + 4.*amp*t - 5.*z) &
                  - 0.125d0*a*amp*eta2*r**3*z*Cos(0. + t - 4.*amp*t + 3.*z) &
                  - 0.0625d0*a*amp*eta2*r**3*Sin(0. + t + 4.*amp*t - 5.*z) &
                  + 0.15625d0*a*amp*eta2*r**3*z**2*Sin(0. + t + 4.*amp*t - 5.*z) &
                  - 0.0625d0*a*amp*eta2*r**3*Sin(0. + t - 4.*amp*t + 3.*z) &
                  - 0.09375d0*a*amp*eta2*r**3*z**2*Sin(0. + t - 4.*amp*t + 3.*z)
       ELSE
          fd = 0.d0
       END IF
       fd = fd/inputs%Re
       
       !Compute fp
       IF (m==0 .AND. TYPE==1) THEN
          fp = 2.d0*r*z**3*COS(t)
       ELSE IF (m==1 .AND. TYPE==2) THEN
          fp = 2.d0*r*COS(t-z)
       ELSE IF (m==2 .AND. TYPE==1) THEN
          fp = z*SIN(t-r) - r*z*COS(t-r)
       ELSE IF (m==1 .AND. TYPE==3) THEN
          fp = r*COS(t-z)
       ELSE IF (m==2 .AND. TYPE==4) THEN
          fp = -2.d0*z*SIN(t-r)
       ELSE IF (m==0 .AND. TYPE==5) THEN
          fp = 3.d0*r**2*z**2*COS(t)
       ELSE IF (m==1 .AND. TYPE==6) THEN
          fp = r**2*SIN(t-z)
       ELSE IF (m==2 .AND. TYPE==5) THEN
          fp = r*SIN(t-r)
       ELSE
          fp = 0.d0
       END IF
       
       !Sum all source terms
       vv = ft + fd + fnl+ fp
    ELSE
       CALL error_petsc('Error in condlim: if_level_set should be true')
    END IF
    RETURN

    !===Dummies variables to avoid warning
    m=i; m=SIZE(opt_density,1); np=ty
    !===Dummies variables to avoid warning
  END FUNCTION source_in_NS_momentum

  !===Extra source in temperature equation. Always called.
  FUNCTION source_in_temperature(TYPE, rr, m, t)RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t   
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    CALL error_petsc('source_in_temperature: should not be called for this test')
    RETURN
  END FUNCTION source_in_temperature

  !===Extra source in level set equation. Always called.
  FUNCTION source_in_level_set(interface_nb,TYPE, rr, m, t)RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m, interface_nb
    REAL(KIND=8),                        INTENT(IN)   :: t   
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
    REAL(KIND=8) :: r
    INTEGER      :: n

    vv = 0.d0
    RETURN

    !===Dummies variables to avoid warning
    n=TYPE; n=SIZE(rr,1); n=m; n=interface_nb; r=t
    !===Dummies variables to avoid warning
  END FUNCTION source_in_level_set

  !===Velocity for boundary conditions in Navier-Stokes.
  !===Can be used also to initialize velocity in: init_velocity_pressure_temperature 
  FUNCTION vv_exact(TYPE,rr,m,t) RESULT(vv)
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
    IF (m==2 .AND. TYPE==1) THEN
       vv = r*z**2*Sin(t-z)
    ELSE IF (m==2 .AND. TYPE==4) THEN
       vv = -r*z**2*Sin(t-z)
    ELSE IF (m==0 .AND. TYPE==5) THEN
       vv = 1.d0
    ELSE
       vv = 0.d0
    END IF
    vv=amp*vv
    RETURN
  END FUNCTION vv_exact

 !===Solid velocity imposed when using penalty technique
 !===Defined in Fourier space on mode 0 only.
 FUNCTION imposed_velocity_by_penalty(rr,t) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2),6)             :: vv

    vv=0.d0
    RETURN
  END FUNCTION imposed_velocity_by_penalty

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

    IF (m==0.AND.TYPE==1) THEN
       vv(:) = rr(1,:)**2*rr(2,:)**3*COS(t)
    ELSE IF (m==1 .AND. TYPE==2) THEN
       vv = rr(1,:)**2*COS(t-rr(2,:))
    ELSE IF (m==2 .AND. TYPE==1) THEN
       vv = rr(1,:)*rr(2,:)*SIN(t-rr(1,:))
    ELSE
       vv = 0.d0
    END IF
    RETURN
  END FUNCTION pp_exact

  !===Temperature for boundary conditions in temperature equation.
  FUNCTION temperature_exact(TYPE,rr,m,t) RESULT (vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    CALL error_petsc('temperature_exact: should not be called for this test')
    RETURN
  END FUNCTION temperature_exact

  !===Can be used to initialize level set in the subroutine init_level_set.
  FUNCTION level_set_exact(interface_nb,TYPE,rr,m,t)  RESULT (vv)
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
          vv = 0.25d0*(2.d0 + Cos(4.d0*(amp*t-z)))
       ELSE IF (m==2 .AND. TYPE==2) THEN
          vv = 0.125d0*a*r**2*Cos(4.d0*(amp*t-z))
       ELSE
          vv = 0.d0
       END IF
    ELSE 
       CALL error_petsc(' BUG in level_set_exact, we should compute only 1 level set')
    END IF
    RETURN
  END FUNCTION level_set_exact

  !===Penalty coefficient (if needed)
  !===This coefficient is equal to zero in subdomain
  !===where penalty is applied (penalty is zero in solid)
  FUNCTION penal_in_real_space(mesh,rr_gauss,angles,nb_angles,nb,ne,time) RESULT(vv)
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
  END FUNCTION penal_in_real_space

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

    vv = 0.d0
    RETURN

  END FUNCTION extension_velocity

  !===============================================================================
  !                       Boundary conditions for Maxwell
  !===============================================================================
  !===Velocity used in the induction equation.
  !===Used only if problem type is mxw and restart velocity is false
  FUNCTION Vexact(m, H_mesh) RESULT(vv)  !Set uniquement a l'induction
    IMPLICIT NONE
    TYPE(mesh_type),                       INTENT(IN) :: H_mesh 
    INTEGER,                               INTENT(IN) :: m
    REAL(KIND=8), DIMENSION(H_mesh%np,6)              :: vv

    vv = 0.d0
    CALL error_petsc('Vexact: should not be called for this test')
  END FUNCTION Vexact

  !===Magnetic field and magnetic induction for quasi-static approximation
  !===if needed
  FUNCTION H_B_quasi_static(char_h_b, rr, m) RESULT(vv) 
    IMPLICIT NONE
    CHARACTER(LEN=1),                    INTENT(IN)   :: char_h_b
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m
    REAL(KIND=8), DIMENSION(SIZE(rr,2),6)             :: vv

    vv = 0.d0
    RETURN
  END FUNCTION H_B_quasi_static

  !===Magnetic field for boundary conditions in the Maxwell equations.
  FUNCTION Hexact(H_mesh,TYPE, rr, m, mu_H_field, t) RESULT(vv) 
    IMPLICIT NONE
    TYPE(mesh_type),                     INTENT(IN)   :: H_mesh
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t 
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: mu_H_field
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv=0.d0
    CALL error_petsc('Hexact: should not be called for this test')
    RETURN
  END FUNCTION Hexact

  !===Scalar potential for boundary conditions in the Maxwell equations.
 FUNCTION Phiexact(TYPE, rr, m, mu_phi,t) RESULT(vv) 
   IMPLICIT NONE
   INTEGER     ,                        INTENT(IN)   :: TYPE
   REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
   INTEGER     ,                        INTENT(IN)   :: m
   REAL(KIND=8),                        INTENT(IN)   :: mu_phi, t
   REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

   vv=0.d0
   CALL error_petsc('Phiexact: should not be called for this test')
   RETURN
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
   REAL(KIND=8)                                      :: alpha,beta
   
   vv=0.d0
   CALL error_petsc('Jexact_gauss: should not be called for this test')
   RETURN
 END FUNCTION Jexact_gauss

  !===Electric field for Neumann BC (cf. doc)
  FUNCTION Eexact_gauss(TYPE, rr, m, mu_phi, sigma, mu_H, t) RESULT(vv)
    IMPLICIT NONE
    INTEGER,                             INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: mu_phi, sigma, mu_H, t
    REAL(KIND=8)                                      :: vv 

    vv = 0.d0
    CALL error_petsc('Eexact: should not be called for this test')
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

    CALL error_petsc('init_maxwell: should not be called for this test')
  END SUBROUTINE init_maxwell

  !===Analytical permeability (if needed)
  !===This function is not needed unless the flag
  !===     ===Use FEM Interpolation for magnetic permeability  (true/false)
  !===is activated and set to .FALSE. in the data data file. Default is .TRUE.
  FUNCTION mu_bar_in_fourier_space(H_mesh,nb,ne,pts,pts_ids) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type), INTENT(IN)                :: H_mesh
    REAL(KIND=8), DIMENSION(ne-nb+1)           :: vv
    INTEGER,     INTENT(IN)                    :: nb, ne
    REAL(KIND=8),DIMENSION(2,ne-nb+1),OPTIONAL :: pts
    INTEGER,     DIMENSION(ne-nb+1), OPTIONAL  :: pts_ids

    vv = 1.d0
    CALL error_petsc('mu_bar_in_fourier_space: should not be called for this test')
    RETURN
  END FUNCTION mu_bar_in_fourier_space

  !===Analytical mu_in_fourier_space (if needed)
  !===This function is not needed unless the flag
  !===     ===Use FEM Interpolation for magnetic permeability  (true/false)
  !===is activated and set to .FALSE. in the data data file. Default is .TRUE.
  FUNCTION grad_mu_bar_in_fourier_space(pt,pt_id) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8),DIMENSION(2), INTENT(in):: pt
    INTEGER,DIMENSION(1), INTENT(in)     :: pt_id
    REAL(KIND=8),DIMENSION(2)            :: vv
    
    vv=0.d0
    CALL error_petsc('grad_mu_bar_in_fourier_space: should not be called for this test')
    RETURN
  END FUNCTION grad_mu_bar_in_fourier_space

  !===Analytical permeability, mu in real space (if needed)
  FUNCTION mu_in_real_space(H_mesh,angles,nb_angles,nb,ne,time) RESULT(vv)
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
  END FUNCTION mu_in_real_space

  FUNCTION sigma_bar_in_fourier_space(H_mesh) RESULT(vv)
    USE def_type_mesh
    TYPE(mesh_type), INTENT(IN)                :: H_mesh
    REAL(KIND=8), DIMENSION(SIZE(H_mesh%rr,2)) :: vv
  END FUNCTION sigma_bar_in_fourier_space

  FUNCTION chi_coeff_law(temp) RESULT(vv)
    REAL(KIND=8) :: temp
    REAL(KIND=8) :: vv
  END FUNCTION chi_coeff_law

  FUNCTION T_dchi_dT_coeff_law(temp) RESULT(vv)
    REAL(KIND=8) :: temp
    REAL(KIND=8) :: vv
  END FUNCTION T_dchi_dT_coeff_law

  FUNCTION nu_tilde_law(temp) RESULT(vv)
    REAL(KIND=8) :: temp
    REAL(KIND=8) :: vv
  END FUNCTION nu_tilde_law

END SUBMODULE BOUNDARY_GENERIC
