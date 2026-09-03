SUBMODULE (boundary_generic_module) boundary_generic
  USE def_type_field_mag
  USE my_util
  USE def_type_mesh
  USE input_data
  USE bessel
  USE user_data

CONTAINS
  !===============================================================================
  !                       Boundary conditions for Navier-Stokes
  !===============================================================================

  MODULE SUBROUTINE my_assign_boundary
    USE interface_boundary
    IMPLICIT NONE
    base_flow => my_base_flow
  END SUBROUTINE my_assign_boundary

  !===Initialize velocity, pressure
  MODULE SUBROUTINE init_velocity_pressure(mesh_f, mesh_c, time, dt, list_mode, &
       un_m1, un, pn_m1, pn, phin_m1, phin)
    USE user_data
    IMPLICIT NONE
    TYPE(mesh_type)                            :: mesh_f, mesh_c
    REAL(KIND=8),                   INTENT(OUT):: time
    REAL(KIND=8),                   INTENT(IN) :: dt
    INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: un_m1, un
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: pn_m1, pn, phin_m1, phin
    INTEGER                                    :: mode, i, j 
    REAL(KIND=8), DIMENSION(mesh_c%np)         :: pn_m2
    REAL(KIND=8), DIMENSION(mesh_f%np)         :: x, bc_flow, incompressibility_flow
    REAL(KIND=8)                               :: ampl = 5.d-2, radius_ratio, omega_ratio, AA, BB
    REAL(KIND=8), PARAMETER :: amp=1.d-4, H = 4.d0, pi = ACOS(-1.d0)
    time = 0.d0
!    radius_ratio = 2.d0
!    omega_ratio = user%Omega_Out/user%Omega_Int
!    AA = ( omega_Ratio - radius_Ratio**2 ) / ( 1.d0 - radius_Ratio**2 )
!    BB =  ( 1.d0 - omega_Ratio ) / ( 1.d0 - radius_Ratio**2 )
!

    IF (inputs%if_linear_navier_stokes_with_u) THEN
        un = 0.d0
        DO i= 1, SIZE(list_mode)
            mode = list_mode(i)
            ASSOCIATE(Rf => mesh_f%rr(1, :))
            IF (mode == 0) THEN
              bc_flow = (1.d0 - Rf)*(2.d0 - Rf)
              DO j = 1, 6
                  !===velocity
                  IF (j==1) THEN
                    un   (:,j,i) = bc_flow**2 * H/(2.d0*pi)*SIN(2*pi*mesh_f%rr(2,:) / H)
                  ELSEIF (j==3) THEN
                    ! CALL RANDOM_NUMBER(x)
                    ! un   (:,j,i) = amp*(x-0.5d0) * (1.d0 - mesh_f%rr(1, :))*(2.d0 - mesh_f%rr(1, :))
                    un   (:,j,i) = bc_flow
                  ELSEIF (j==5) THEN
                    ! CALL RANDOM_NUMBER(x)
                    ! un   (:,j,i) = (1.d0 - mesh_f%rr(1, :))*(2.d0 - mesh_f%rr(1, :)) * amp*x
                    incompressibility_flow = bc_flow * ((1.d0 - Rf)*(2.d0 - Rf)/Rf + 2.d0*(2.d0*Rf-3.d0))
                    un   (:,j,i) = incompressibility_flow * COS(2*pi*mesh_f%rr(2,:) / H)
                  ENDIF
                  un_m1(:,j,i) = un   (:,j,i) 
              END DO
            END IF
            END ASSOCIATE
            DO j = 1, 2
                !===pressure
                pn_m2(:)       = pp_exact(j,mesh_c%rr,mode,time-2*dt)
                pn_m1  (:,j,i) = pp_exact(j,mesh_c%rr,mode,time-dt)
                pn     (:,j,i) = pp_exact(j,mesh_c%rr,mode,time)
                phin_m1(:,j,i) = pn_m1(:,j,i) - pn_m2(:)
                phin   (:,j,i) = Pn   (:,j,i) - pn_m1(:,j,i)
            ENDDO
        ENDDO
    ELSE
        DO i= 1, SIZE(list_mode)
            mode = list_mode(i) 
            DO j = 1, 6 
              !===velocity
              un(:, j, i) = vv_exact(j, mesh_f%rr, mode, time)
              un_m1(:, j, i) = vv_exact(j, mesh_f%rr, mode, time-dt)
            END DO
            DO j = 1, 2
              !===pressure
              ! IF (mode==0 .AND. j==1) THEN
              !    pn_m2(:)       = 2.d0*AA*BB*log(mesh_c%rr(1,:)) + AA**2/2.d0*(mesh_c%rr(1,:))**2 - BB**2/(mesh_c%rr(1,:))
              ! ELSE
              !    pn_m2 = 0.d0
              ! END IF
              pn_m2(:)       = pp_exact(j,mesh_c%rr,mode,time-2*dt)
              pn_m1  (:,j,i) = pp_exact(j,mesh_c%rr,mode,time-dt)
              pn     (:,j,i) = pp_exact(j,mesh_c%rr,mode,time)
              phin_m1(:,j,i) = pn_m1(:,j,i) - pn_m2(:)
              phin   (:,j,i) = Pn   (:,j,i) - pn_m1(:,j,i)
              !pn_m2(:)       = pp_exact(j,mesh_c%rr,mode,time-2*dt)
              !pn_m1  (:,j,i) = pp_exact(j,mesh_c%rr,mode,time-dt)
              !pn     (:,j,i) = pp_exact(j,mesh_c%rr,mode,time)
              !phin_m1(:,j,i) = pn_m1(:,j,i) - pn_m2(:)
              !phin   (:,j,i) = Pn   (:,j,i) - pn_m1(:,j,i)
            ENDDO
        ENDDO
    END IF

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
  MODULE FUNCTION source_in_NS_momentum(TYPE_VEC, rr, mode, i, time, Re, ty, density, tempn, concn) RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                             INTENT(IN) :: TYPE_VEC
    REAL(KIND=8), DIMENSION(:,:),             INTENT(IN) :: rr
    INTEGER     ,                             INTENT(IN) :: mode, i
    REAL(KIND=8),                             INTENT(IN) :: time
    REAL(KIND=8),                             INTENT(IN) :: Re
    CHARACTER(LEN=2),                         INTENT(IN) :: ty
    REAL(KIND=8), DIMENSION(:,:,:),           INTENT(IN) :: density
    REAL(KIND=8), DIMENSION(:,:,:),           INTENT(IN) :: tempn
    REAL(KIND=8), DIMENSION(:,:,:),           INTENT(IN) :: concn
    REAL(KIND=8), DIMENSION(SIZE(rr,2))                  :: vv

    vv = 0.d0
    RETURN
  END FUNCTION source_in_NS_momentum

  !===Extra source in temperature equation. Always called.
  MODULE FUNCTION source_in_temperature(TYPE_VEC, rr, m, t)RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE_VEC
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t   
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    RETURN
  END FUNCTION source_in_temperature

  MODULE FUNCTION source_in_level_set(interface_nb,TYPE_VEC, rr, m, t)RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE_VEC
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m, interface_nb
    REAL(KIND=8),                        INTENT(IN)   :: t   
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv=0.d0
    RETURN
  END FUNCTION source_in_level_set

  !===Velocity for boundary conditions in Navier-Stokes.
  !===Can be used also to initialize velocity in: init_velocity_pressure_temperature 
  ! MODULE FUNCTION vv_exact(TYPE_VEC,rr,m,t) RESULT(vv)
  !   USE user_data
  !   IMPLICIT NONE
  !   INTEGER     ,                        INTENT(IN)   :: TYPE_VEC
  !   REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
  !   INTEGER,                             INTENT(IN)   :: m
  !   REAL(KIND=8),                        INTENT(IN)   :: t
  !   REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
  !   REAL(KIND=8)                                      :: r,z,radius_ratio, omega_ratio, AA, BB
  !   INTEGER                                           :: n
  !   radius_ratio = .5d0
  !   omega_ratio = user%Omega_Out/user%Omega_Int
  !   AA = ( omega_Ratio - radius_Ratio**2 ) / ( 1.d0 - radius_Ratio**2 )
  !   BB =  ( 1.d0 - omega_Ratio ) / ( 1.d0 - radius_Ratio**2 )
  !   vv=0.0

  !   IF (TYPE_VEC==3 .AND. m==0) THEN
  !      DO n = 1, SIZE(rr,2)
  !         r= rr(1,n)
  !         z= rr(2,n)
  !         vv(n) = r * AA + BB / r
  !      END DO
  !   END IF
  !   RETURN
  ! END FUNCTION vv_exact

  MODULE FUNCTION vv_exact(TYPE_VEC,rr,m,t) RESULT(vv)
    USE user_data
    USE interface_boundary, ONLY: base_flow
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE_VEC
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    IF (inputs%if_linear_navier_stokes_with_u) THEN
      vv = 0.d0    
    ELSE
      vv = base_flow(TYPE_VEC,rr,m)
    END IF
    RETURN
  END FUNCTION vv_exact

  FUNCTION my_base_flow(TYPE_VEC,rr,m) RESULT(vv)
    USE user_data
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE_VEC
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
    REAL(KIND=8)                                      :: r,z,radius_ratio, omega_ratio, AA, BB
    INTEGER                                           :: n
    radius_ratio = .5d0
    omega_ratio = user%Omega_Out/user%Omega_Int
    AA = ( omega_Ratio - radius_Ratio**2 ) / ( 1.d0 - radius_Ratio**2 )
    BB =  ( 1.d0 - omega_Ratio ) / ( 1.d0 - radius_Ratio**2 )
    vv=0.0

    IF (TYPE_VEC==3 .AND. m==0) THEN
       DO n = 1, SIZE(rr,2)
          r= rr(1,n)
          z= rr(2,n)
          vv(n) = r * AA + BB / r
       END DO
    END IF
    RETURN
  END FUNCTION my_base_flow

  !===Solid velocity imposed when using penalty technique
  !===Defined in Fourier space on mode 0 only.
  MODULE FUNCTION imposed_velocity_by_penalty(rr,t) RESULT(vv)
    USE user_data
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2),6)             :: vv
    REAL(KIND=8)                                      :: r, z
    INTEGER :: n 

    vv=0.d0
    RETURN
  END FUNCTION imposed_velocity_by_penalty

  !===Pressure for boundary conditions in Navier-Stokes.
  !===Can be used also to initialize pressure in: init_velocity_pressure_temperature
  !===Use this routine for outflow BCs only.
  !===CAUTION: Do not enfore BCs on pressure where normal component 
  !            of velocity is prescribed.
  MODULE FUNCTION pp_exact(TYPE_VEC,rr,m,t) RESULT (vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE_VEC
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv 

    vv=0.d0
    RETURN
  END FUNCTION pp_exact

  !===Temperature for boundary conditions in temperature equation.
  MODULE FUNCTION temperature_exact(TYPE_VEC,rr,m,t) RESULT (vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE_VEC
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    RETURN
  END FUNCTION temperature_exact

  !===Concentration for boundary conditions in concentration equation.
  MODULE FUNCTION concentration_exact(TYPE_VEC,rr,m,t) RESULT (vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE_VEC
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    RETURN
  END FUNCTION concentration_exact

  MODULE FUNCTION level_set_exact(interface_nb,TYPE_VEC,rr,m,t)  RESULT (vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE_VEC
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m, interface_nb
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    RETURN
  END FUNCTION level_set_exact

  !===Penalty coefficient (if needed)
  !===This coefficient is equal to 1 in fluid subdomain
  !===This coefficient is equal to zero in solid subdomain where penalty is applied


  MODULE FUNCTION penal_in_real_space(mesh,rr_gauss,angles,nb_angles,nb,ne,time) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type)                            :: mesh
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)   :: rr_gauss
    REAL(KIND=8), DIMENSION(:), INTENT(IN)     :: angles
    INTEGER,                    INTENT(IN)     :: nb_angles
    INTEGER,                    INTENT(IN)     :: nb, ne
    REAL(KIND=8),               INTENT(IN)     :: time
    REAL(KIND=8), DIMENSION(nb_angles,ne-nb+1) :: vv
    INTEGER                                    :: n, n_loc
    REAL(KIND=8)                               :: r,z

    vv=0.d0
    RETURN
  END FUNCTION penal_in_real_space


  !===Extension of the velocity field in the solid.
  !===Used only if Maxwell domain for H is larger than
  !===domain for Navier-Stokes velocity, or if there is
  !===no Navier-Stokes domain (i.e. mxw or mxx mode)
  MODULE FUNCTION extension_velocity(TYPE_VEC, H_mesh, mode, t, n_start) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type),                     INTENT(IN)   :: H_mesh     
    INTEGER     ,                        INTENT(IN)   :: TYPE_VEC, n_start
    INTEGER,                             INTENT(IN)   :: mode 
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(H_Mesh%np)                :: vv

    vv = 0.d0
    RETURN
  END FUNCTION extension_velocity

  MODULE FUNCTION extension_temperature(TYPE_VEC, H_mesh, mode, t, n_start) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type),                     INTENT(IN)   :: H_mesh
    INTEGER     ,                        INTENT(IN)   :: TYPE_VEC, n_start
    INTEGER,                             INTENT(IN)   :: mode
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(H_Mesh%np)                :: vv

    vv = 0.d0
    RETURN
  END FUNCTION extension_temperature

  MODULE FUNCTION extension_concentration(TYPE_VEC, vv_mesh, mode, t, n_start) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type),                     INTENT(IN)   :: vv_mesh
    INTEGER     ,                        INTENT(IN)   :: TYPE_VEC, n_start
    INTEGER,                             INTENT(IN)   :: mode
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(vv_mesh%np)                :: vv

    vv = 0.d0
    RETURN
  END FUNCTION extension_concentration
  !===============================================================================
  !                       Boundary conditions for Maxwell
  !===============================================================================
  !===Velocity used in the induction equation.
  !===Used only in mxw mode
  MODULE FUNCTION Vexact(m, H_mesh) RESULT(vv)  !Set uniquement a l'induction
    IMPLICIT NONE
    TYPE(mesh_type),                       INTENT(IN) :: H_mesh 
    INTEGER,                               INTENT(IN) :: m
    REAL(KIND=8), DIMENSION(H_mesh%np,6)              :: vv

    vv = 0.d0
    RETURN
  END FUNCTION Vexact

  MODULE FUNCTION H_B_quasi_static(char_h_b, rr, m) RESULT(vv) 
    USE user_data
    IMPLICIT NONE
    CHARACTER(LEN=1),                    INTENT(IN)   :: char_h_b
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m
    INTEGER :: TYPE_VEC
    REAL(KIND=8), DIMENSION(SIZE(rr,2),6)             :: vv


    IF (inputs%if_quasi_static_approx) THEN
       vv = 0.d0
    ELSE IF (char_h_b == "F") THEN
       vv = 0.d0
       DO TYPE_VEC=1, 6
          IF (TYPE_VEC==5 .AND. m==0) THEN
                  vv(:,TYPE_VEC) = user%Bz_imposed
          END IF
       END DO 
    ELSE
       CALL error_petsc('H_B_quasi_static should not be called') 
    END IF
    RETURN
  END FUNCTION H_B_quasi_static

  !===Magnetic field for boundary conditions in the Maxwell equations.
  MODULE FUNCTION Hexact(H_mesh, TYPE_VEC, rr, m, mu_H_field, t) RESULT(vv)  
    IMPLICIT NONE
    TYPE(mesh_type),                     INTENT(IN)   :: H_mesh
    INTEGER     ,                        INTENT(IN)   :: TYPE_VEC
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t 
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: mu_H_field
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    RETURN
  END FUNCTION Hexact

  FUNCTION Hexact_init(H_mesh, TYPE_VEC, rr, m, mu_H_field, t) RESULT(vv)  
    IMPLICIT NONE
    TYPE(mesh_type),                     INTENT(IN)   :: H_mesh
    INTEGER     ,                        INTENT(IN)   :: TYPE_VEC
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t 
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: mu_H_field
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    !===Constant vertical field Hz
    IF (m == 0) THEN
       IF (TYPE_VEC == 5) THEN
          vv = 1.d-6
       ENDIF
    ENDIF
    !===Constant horizontal field Hx
    IF (m == 1) THEN
       IF (TYPE_VEC == 1) THEN
          vv = 1.d-6
       ELSEIF (TYPE_VEC == 4) THEN
          vv = -1.d-6
       ENDIF
    ENDIF
    RETURN
  END FUNCTION Hexact_init

  !===Scalar potential for boundary conditions in the Maxwell equations.
  MODULE FUNCTION Phiexact(TYPE_VEC, rr, m, mu_phi,t) RESULT(vv) 
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE_VEC
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: mu_phi, t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv   
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: r, z
    INTEGER                                           :: n

    vv=0.d0
    RETURN
  END FUNCTION Phiexact

  !===Current in Ohm's law. Curl(H) = sigma(E + uxB) + current
  MODULE FUNCTION Jexact_gauss(TYPE_VEC, rr, m, mu_phi, sigma, mu_H, t, mesh_id, opt_B_ext) RESULT(vv) 
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE_VEC
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m 
    REAL(KIND=8),                        INTENT(IN)   :: mu_phi, sigma, mu_H, t 
    INTEGER     ,                        INTENT(IN)   :: mesh_id
    REAL(KIND=8), DIMENSION(6), OPTIONAL,INTENT(IN)   :: opt_B_ext 
    REAL(KIND=8)                                      :: vv

    vv = 0.d0
    RETURN
  END FUNCTION Jexact_gauss

  !===Electric field for Neumann BC (cf. doc)
  MODULE FUNCTION Eexact_gauss(TYPE_VEC, rr, m, mu_phi, sigma, mu_H, t) RESULT(vv)
    IMPLICIT NONE
    INTEGER,                             INTENT(IN)   :: TYPE_VEC
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: mu_phi, sigma, mu_H, t
    REAL(KIND=8)                                      :: vv 

    vv = 0.d0
    RETURN
  END FUNCTION Eexact_gauss

  !===Initialization of magnetic field and scalar potential (if present)
  MODULE SUBROUTINE init_maxwell(H_mesh, phi_mesh, time, dt, mu_H_field, mu_phi, &
       list_mode, mag_n, mag_n1)
      !  list_mode, magn, magn1)
    IMPLICIT NONE
    TYPE(mesh_type)                            :: H_mesh, phi_mesh     
    REAL(KIND=8),                   INTENT(OUT):: time
    REAL(KIND=8),                   INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(:),     INTENT(IN) :: mu_H_field
    REAL(KIND=8),                   INTENT(IN) :: mu_phi
    INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode    
    TYPE(mag_field_type), INTENT(OUT)          :: mag_n, mag_n1
    ! TYPE(mag_field_type), INTENT(OUT)        :: magn, magn1
    INTEGER                                    :: i, k, mode, m_max_c
    REAL(KIND=8), DIMENSION(H_mesh%np)         :: x
    REAL(KIND=8)                               :: ampl = 1.d-6

    CALL mag_n1%allocate_induction_fields; CALL mag_n%allocate_induction_fields
    time = 0.d0
    m_max_c = SIZE(list_mode)

     time = -dt

    time = time + dt
    DO k=1,6
       DO i=1, m_max_c
          mode = list_mode(i)
          mag_n%H(:,k,i) = Hexact_init(H_mesh, k, H_mesh%rr, mode, mu_H_field, time)
          IF (k<3) THEN
             mag_n%phi(:,k,i) = Phiexact(k, phi_mesh%rr, mode, mu_phi, time)
          ENDIF
          IF (k==3 .AND. list_mode(i)==0) THEN
             CALL RANDOM_SEED
             CALL RANDOM_NUMBER(x)
             mag_n%H(:,k,i) = ampl*(x - 0.5d0)
          END IF
!          IF (list_mode(i)==1) THEN
!             CALL RANDOM_SEED
!             CALL RANDOM_NUMBER(x)
!             Hn(:,k,i) = Hn(:,k,i) + ampl*(x - 0.5d0)
!          END IF
!          IF (list_mode(i).GE.2) THEN
!             CALL RANDOM_SEED
!             CALL RANDOM_NUMBER(x)
!             Hn(:,k,i) = ampl*(x - 0.5d0)
!          END IF
       ENDDO
    ENDDO
    mag_n1%H = mag_n%H
    mag_n1%phi = mag_n%phi
    RETURN
  END SUBROUTINE init_maxwell

  !===Analytical permeability (if needed)
  !===This function is not needed unless the flag
  !===     ===Use FEM Interpolation for magnetic permeability  (true/false)
  !===is activated and set to .FALSE. in the data data file. Default is .TRUE.
  MODULE FUNCTION mu_bar_in_fourier_space(H_mesh,nb,ne,pts,pts_ids) RESULT(vv)
    USE user_data
    IMPLICIT NONE
    TYPE(mesh_type), INTENT(IN)                :: H_mesh
    REAL(KIND=8), DIMENSION(ne-nb+1)           :: vv
    INTEGER,     INTENT(IN)                    :: nb, ne
    REAL(KIND=8),DIMENSION(2,ne-nb+1),OPTIONAL :: pts
    INTEGER,     DIMENSION(ne-nb+1), OPTIONAL  :: pts_ids

    INTEGER                                    :: n, n_loc
    INTEGER, DIMENSION(H_mesh%np)              :: id
    REAL(KIND=8),DIMENSION(ne-nb+1)            :: r,z

    vv = 0.d0
    RETURN
  END FUNCTION mu_bar_in_fourier_space

  !===Analytical mu_in_fourier_space (if needed)
  !===This function is not needed unless the flag
  !===     ===Use FEM Interpolation for magnetic permeability  (true/false)
  !===is activated and set to .FALSE. in the data data file. Default is .TRUE.
  MODULE FUNCTION grad_mu_bar_in_fourier_space(pt,pt_id) RESULT(vv)
    USE user_data
    IMPLICIT NONE
    REAL(KIND=8),DIMENSION(2), INTENT(in):: pt
    INTEGER,DIMENSION(1), INTENT(in)     :: pt_id
    REAL(KIND=8),DIMENSION(2)            :: vv

    vv = 0.d0
    RETURN
  END FUNCTION grad_mu_bar_in_fourier_space

  !===Analytical permeability, mu in real space (if needed)
  MODULE FUNCTION mu_in_real_space(H_mesh,angles,nb_angles,nb,ne,time) RESULT(vv)
    USE user_data
    IMPLICIT NONE
    TYPE(mesh_type), INTENT(IN)                :: H_mesh
    REAL(KIND=8), DIMENSION(:), INTENT(IN)     :: angles
    INTEGER, INTENT(IN)                        :: nb_angles
    INTEGER, INTENT(IN)                        :: nb, ne
    REAL(KIND=8), INTENT(IN)                   :: time
    REAL(KIND=8), DIMENSION(nb_angles,ne-nb+1) :: vv

    vv = 0.d0
    RETURN
  END FUNCTION mu_in_real_space


  !=== Only used for multiphase flow with variable electrical conductivity
  !===This function is not needed unless the following line is in the data
  !=== ===Conductivity of fluid 0, fluid 1, ...
  !=== The above line should be followed by the value of the conductivity in the different fluids.
  !=== The result vv has to smaller than the effective conductivity on every nodes of the mesh.
  MODULE FUNCTION sigma_bar_in_fourier_space(H_mesh) RESULT(vv)
    USE input_data
    IMPLICIT NONE
    TYPE(mesh_type), INTENT(IN)                :: H_mesh
    REAL(KIND=8), DIMENSION(SIZE(H_mesh%rr,2)) :: vv

    vv=0.9d0*MINVAL(inputs%sigma_fluid)
    RETURN
  END FUNCTION sigma_bar_in_fourier_space

  !===Coefficient contaning the magnetic susceptibility for magnetic force in
  !ferrofluids: 
  !===F = chi_coeff(T) * grad(H**2/2) (Kelvin form)
  !===or F = -H**2/2 * grad(chi_coeff(T)) (Helmholtz form)
  MODULE FUNCTION chi_coeff_law(temp) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: temp
    REAL(KIND=8) :: vv

    vv = 0.d0
    RETURN
  END FUNCTION chi_coeff_law

  !===Coefficient contaning the temperature dependant factor in the
  !===pyromagnetic coefficient term of the temperature equation for ferrofluids: 
  !===T * dchi/dT(T) * D/Dt(H**2/2)
  MODULE FUNCTION T_dchi_dT_coeff_law(temp) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: temp
    REAL(KIND=8) :: vv

    vv = 0.d0
    RETURN
  END FUNCTION T_dchi_dT_coeff_law

  !===Kinematic viscosity's law of temperature
  MODULE FUNCTION nu_tilde_law(temp) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: temp
    REAL(KIND=8) :: vv

    vv = 0.d0
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
    REAL(KIND=8), INTENT(IN) :: theta
    REAL(KIND=8) :: vv

    vv = 0.d0
    RETURN
  END FUNCTION curved_boundary_radius

END SUBMODULE BOUNDARY_GENERIC
