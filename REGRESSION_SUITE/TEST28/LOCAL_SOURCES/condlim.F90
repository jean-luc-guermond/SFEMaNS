SUBMODULE (boundary_generic_module) boundary_generic
  USE my_util
  USE def_type_mesh
  USE input_data
  USE bessel
  USE user_data

  REAL(KIND=8),  PARAMETER   :: pi = 3.14159265358979323846d0
  REAL(KIND=8),  PARAMETER   :: twopi=2*pi
  REAL(KIND=8),  PARAMETER:: mu_disk= 5.0d1  !mu, Permeability of blades and disk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !===TM73
  !Parameters for discontinuous blades
  !some offset to begin at the same vertical axis as the bottom propeller
  REAL(KIND=8),  PARAMETER:: top_propeller_angle_offset =  0.d0
  INTEGER :: nblades = 8 !number of blades
  REAL(KIND=8),  PARAMETER::  lw = 0.0125 !lw= is the half thickness of the blades
  REAL(KIND=8),  PARAMETER::  disk_bot=-1.1d0, top_of_disk_bot=-0.9d0, top_of_blade_bot=-0.7d0
  REAL(KIND=8),  PARAMETER::  disk_top= 1.1d0, bot_of_disk_top= 0.9d0, bot_of_blade_top= 0.7d0
  REAL(KIND=8),  PARAMETER:: disk_radius=0.75d0, hole_radius=0.1d0
  !For straight blades use two_rp=0.d0
  REAL(KIND=8),  PARAMETER::  two_rp = disk_radius/SIN(twopi*24.d0/360.d0)
  !Volume of the cylinder (all domain)
  REAL(KIND=8)         :: omega_Vol= pi*(1.0)*(2.0) !pi*r^2*h

  !Parameters for smooth_blades
  REAL(KIND=8),  PARAMETER::  wjump_hole = 0.06*(1.0),    wjump= 0.04*(1.0)
  REAL(KIND=8),  PARAMETER::  hole_r = hole_radius,       hole_rp=hole_r - wjump_hole
  REAL(KIND=8),  PARAMETER::  disk_r = disk_radius, disk_rp= disk_r- wjump
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Bottom Smooth_propeller
  REAL(KIND=8),  PARAMETER::  cyl_bott = -1.0
  REAL(KIND=8),  PARAMETER::  Bdisk_z =  cyl_bott  + 0.3d0,          bdisk_z_p= Bdisk_z  - wjump
  REAL(KIND=8),  PARAMETER::  zbot =     cyl_bott  + 0.1d0,          zbot_p = zbot - wjump
  REAL(KIND=8),  PARAMETER::  zbot_bar =     zbot  - 0.04d0,         zbot_bar_p = zbot_bar - wjump
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Top Smooth_propeller
  REAL(KIND=8),  PARAMETER::  cyl_top = 1.0
  REAL(KIND=8),  PARAMETER::  Tdisk_z = cyl_top - 0.3d0,           Tdisk_z_p= Tdisk_z  + wjump
  REAL(KIND=8),  PARAMETER::  ztop =    cyl_top - 0.1d0,           ztop_p = ztop + wjump
  REAL(KIND=8),  PARAMETER::  ztop_bar =  ztop  +  0.04d0,         ztop_bar_p = ztop_bar + wjump
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Parameters for  both smooth_blades
  REAL(KIND=8),  PARAMETER::   alpha=200.d0*(1.0), alpha_th = 80.d0*(1.0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Do we want both propellers?
  LOGICAL,PARAMETER::  if_bottom_prop=.TRUE.
  LOGICAL,PARAMETER::  if_top_prop=.TRUE.
  REAL(KIND=8), PARAMETER::  solid_vel=1.0;

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

    vv = 0.d0
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

    vv = 0.d0
    CALL error_petsc('source_in_temperature: should not be called for this test')
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

    vv = 0.d0
    CALL error_petsc('source_in_temperature: should not be called for this test')
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
    REAL(KIND=8)                                      :: r,z
    INTEGER                                           :: n

    vv=0.0

    IF (type==3 .AND. m==0) THEN
       DO n = 1, SIZE(rr,2)
          r= rr(1,n)
          z= rr(2,n)
          ! Are we in the Bottom propeller?
          IF ( if_bottom_prop .AND. r <disk_radius .AND. z <  top_of_blade_bot ) then
             vv(n)=solid_vel*r
          END IF

          !are we in the top Propeller?
          IF ( if_top_prop  .AND. r <disk_radius .AND. z  >   bot_of_blade_top) then
             vv(n)=-solid_vel*r

          END IF
       END DO
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
    REAL(KIND=8)                                      :: r, z
    INTEGER                                           :: n

    vv=0.d0
    DO n = 1, SIZE(rr,2)
       r= rr(1,n)
       z= rr(2,n)
       IF (z<-0.5d0) THEN
          vv(n,3) =  solid_vel*rr(1,n)
       ELSE
          vv(n,3) =  -solid_vel*rr(1,n)
       ENDIF
    END DO
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

    vv = 0.d0
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

    vv = 0.d0
    CALL error_petsc('temperature_exact: should not be called for this test')
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !1) USE Smooth Blades
    vv=smooth_penal_in_real_space(mesh,rr_gauss,angles,nb_angles,nb,ne,time)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  MODULE FUNCTION Hexact(H_mesh,TYPE, rr, m, mu_H_field, t) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type),                     INTENT(IN)   :: H_mesh
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: mu_H_field
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
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
    REAL(KIND=8)                                      :: vv

    vv = 0.d0
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
    CALL error_petsc('Eexact: should not be called for this test')
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
    INTEGER,     INTENT(IN)                    :: nb, ne
    REAL(KIND=8),DIMENSION(2,ne-nb+1),OPTIONAL :: pts
    INTEGER,     DIMENSION(ne-nb+1), OPTIONAL  :: pts_ids
    REAL(KIND=8), DIMENSION(ne-nb+1)           :: vv
    INTEGER                                    :: n
    REAL(KIND=8),DIMENSION(ne-nb+1)            :: r,z

    IF( PRESENT(pts) .AND. PRESENT(pts_ids) ) THEN !Computing mu at  pts
       r=pts(1,nb:ne)
       z=pts(2,nb:ne)
    ELSE
       r=H_mesh%rr(1,nb:ne) !Computing mu  at nodes
       z=H_mesh%rr(2,nb:ne)
    END IF

    DO n = 1, ne - nb + 1
       vv(n)=mu_bar_func(r(n),z(n))
    END DO

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

    vv=grad_mu_bar_func(pt(1),pt(2))
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

    vv = (1.d0 - penal_in_real_space(H_mesh,H_mesh%rr(:,nb:ne),angles,nb_angles,nb,ne,time))*(mu_disk - 1.d0 ) +1.d0
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

  MODULE FUNCTION chi_coeff_law(temp) RESULT(vv)
    REAL(KIND=8) :: temp
    REAL(KIND=8) :: vv

    vv = 0.d0*temp
    CALL error_petsc('chi_coeff_law: should not be called for this test')
    RETURN
  END FUNCTION chi_coeff_law

  MODULE FUNCTION T_dchi_dT_coeff_law(temp) RESULT(vv)
    REAL(KIND=8) :: temp
    REAL(KIND=8) :: vv

    vv = 0.d0*temp
    CALL error_petsc('T_dchi_dT_coeff_law: should not be called for this test')
    RETURN
  END FUNCTION T_dchi_dT_coeff_law

  MODULE FUNCTION nu_tilde_law(temp) RESULT(vv)
    REAL(KIND=8) :: temp
    REAL(KIND=8) :: vv

    vv = 0.d0*temp
    CALL error_petsc('nu_tilde_law: should not be called for this test')
    RETURN
  END FUNCTION nu_tilde_law

  FUNCTION mu_bar_func(r,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8)                               :: r,z,vv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !mu bar for Blades , a disks with no hole

    IF ( if_bottom_prop .AND.  if_top_prop) THEN
       vv=(top_mu_bar_func(r,z) + bottom_mu_bar_func(r,z))*(mu_disk-1.0) + 1.0
    ELSE IF (if_bottom_prop) THEN
       vv= bottom_mu_bar_func(r,z)*(mu_disk-1.0) + 1.0
    ELSE
       vv= top_mu_bar_func(r,z)*(mu_disk-1.0) + 1.0
    END IF
    RETURN
  END FUNCTION mu_bar_func

  FUNCTION bottom_mu_bar_func(r,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8)                               :: r,z,vv
    REAL(KIND=8)                               :: r2,r3,z0,z1
    REAL(KIND=8)                               :: psi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !mu bar for Blades , a disks with no hole
    r2=disk_rp
    r3=disk_r

    z0=zbot_bar_p
    z1=zbot_bar

    psi=0.d0
    IF ( z .GE.  z1  .AND. r .GE. r3 )  THEN
       psi = 0.d0 ;
    ELSE IF (r.LE.r2) THEN
       IF(z.LE.z0) THEN
          psi=1.0;
       ELSE IF (z.GE.z1) THEN
          psi=0.0;
       ELSE
          psi=smooth_jump_down(z,z0,z1);
       END IF
    ELSE IF (r.GE.r2 .AND. r.LE.r3) THEN
       IF(z.GE.z1) THEN
          psi=0.0;
       ELSE IF(z.LE.z0) THEN
          psi=smooth_jump_down(r,r2,r3);
       ELSE
          psi=smooth_jump_down(r,r2,r3)*smooth_jump_down(z,z0,z1);
       END IF
    END IF

    vv = psi
    RETURN
  END FUNCTION bottom_mu_bar_func
  !
  FUNCTION top_mu_bar_func(r,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8)                               :: r,z,vv,psi
    REAL(KIND=8)                               :: r2,r3,z2,z3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !mu bar for Blades , a disks with no hole
    r2=disk_rp
    r3=disk_r

    z2=ztop_bar_p
    z3=ztop_bar

    psi=0.d0
    IF ( z .LE.  z3 .AND. r .GE. r3)  THEN
       psi = 0.d0 ;
    ELSE IF(r.LE.r2) THEN
       IF(z.GE.z2) THEN
          psi=1.0;
       ELSE IF (z.LE.z3) THEN
          psi=0.0;
       ELSE
          psi=smooth_jump_up(z,z3,z2);
       END IF
    ELSE IF (r.GE.r2 .AND. r.LE.r3) THEN
       IF(z.LE.z3) THEN
          psi=0.0;
       ELSE IF(z.GE.z2) THEN
          psi=smooth_jump_down(r,r2,r3);
       ELSE
          psi=smooth_jump_down(r,r2,r3)*smooth_jump_up(z,z3,z2);
       END IF
    END IF

    vv=psi
    RETURN
  END FUNCTION top_mu_bar_func

  FUNCTION grad_mu_bar_func(r,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8)                               :: r,z
    REAL(KIND=8),DIMENSION(2)                  :: vv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !mu bar for Blades , a disks with no hole

    IF ( if_bottom_prop .AND.  if_top_prop) THEN
       vv=( grad_top_mu_bar_func(r,z) + grad_bottom_mu_bar_func(r,z) )*(mu_disk-1.0)
    ELSE IF (if_bottom_prop) THEN
       vv= grad_bottom_mu_bar_func(r,z)*(mu_disk-1.0)
    ELSE
       vv= grad_top_mu_bar_func(r,z)*(mu_disk-1.0)
    END IF
    RETURN
  END FUNCTION grad_mu_bar_func

  FUNCTION grad_bottom_mu_bar_func(r,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8)                               :: r,z
    REAL(KIND=8),DIMENSION(2)                  :: vv
    REAL(KIND=8)                               :: r2,r3,z0,z1
    REAL(KIND=8)                               :: DFr,DFz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !mu bar for Blades , a disks with no hole
    r2=disk_rp
    r3=disk_r

    z0=zbot_bar_p
    z1=zbot_bar

    DFr=0.d0
    DFz=0.d0

    IF ( z .GE.  z1  .AND. r .GE. r3 )  THEN
       DFr=0.d0
       DFz=0.d0
    ELSE IF (r.LE.r2) THEN
       IF(z.LE.z0) THEN
          DFr=0.d0
          DFz=0.d0
       ELSE IF (z.GE.z1) THEN
          DFr=0.d0
          DFz=0.d0
       ELSE
          DFr=0.d0
          DFz=Dsmooth_jump_down(z,z0,z1);
       END IF
    ELSE IF (r.GE.r2 .AND. r.LE.r3) THEN
       IF(z.GE.z1) THEN
          DFr=0.d0
          DFz=0.d0
       ELSE IF(z.LE.z0) THEN
          DFr=Dsmooth_jump_down(r,r2,r3);
          DFz=0.d0
       ELSE
          DFr=Dsmooth_jump_down(r,r2,r3)*smooth_jump_down(z,z0,z1);
          DFz=smooth_jump_down(r,r2,r3)*Dsmooth_jump_down(z,z0,z1);
       END IF
    END IF

    vv(1)=DFr
    vv(2)=DFz
    RETURN
  END FUNCTION grad_bottom_mu_bar_func

  FUNCTION grad_top_mu_bar_func(r,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8)                               :: r,z
    REAL(KIND=8),DIMENSION(2)                  :: vv
    REAL(KIND=8)                               :: r2,r3,z2,z3
    REAL(KIND=8)                               :: DFr,DFz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !mu bar for Blades , a disks with no hole
    r2=disk_rp
    r3=disk_r

    z2=ztop_bar_p
    z3=ztop_bar

    DFr=0.d0
    DFz=0.d0
    IF ( z .LE.  z3 .AND. r .GE. r3)  THEN
       DFr=0.d0
       DFz=0.d0
    ELSE IF(r.LE.r2) THEN
       IF(z.GE.z2) THEN
          DFr=0.d0
          DFz=0.d0
       ELSE IF (z.LE.z3) THEN
          DFr=0.d0
          DFz=0.d0
       ELSE
          DFr=0.d0
          DFz=Dsmooth_jump_up(z,z3,z2);
       END IF
    ELSE IF (r.GE.r2 .AND. r.LE.r3) THEN
       IF(z.LE.z3) THEN
          DFr=0.d0
          DFz=0.d0
       ELSE IF(z.GE.z2) THEN
          DFr=Dsmooth_jump_down(r,r2,r3);
          DFz=0.d0
       ELSE
          DFr=Dsmooth_jump_down(r,r2,r3)*smooth_jump_up(z,z3,z2);
          DFz=smooth_jump_down(r,r2,r3)*Dsmooth_jump_up(z,z3,z2);
       END IF
    END IF

    vv(1)=DFr
    vv(2)=DFz
    RETURN
  END FUNCTION grad_top_mu_bar_func

  !This is 1 in the fluid, and  0 in the solid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION smooth_penal_in_real_space(mesh,rr_gauss,angles,nb_angles,nb,ne,time) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type)                            :: mesh
    REAL(KIND=8), DIMENSION(:,:)               :: rr_gauss
    REAL(KIND=8), DIMENSION(:)                 :: angles
    INTEGER                                    :: nb_angles
    INTEGER                                    :: nb, ne
    REAL(KIND=8), DIMENSION(nb_angles,ne-nb+1) :: vv
    INTEGER                                    :: n, n_loc
    REAL(KIND=8)                               :: r,z,time

    DO n = nb, ne
       n_loc = n - nb + 1

       r = rr_gauss(1,n_loc)
       z = rr_gauss(2,n_loc)

       ! DCQ   14/08/2014
       IF (if_bottom_prop .AND.  if_top_prop) THEN
          vv(:,n_loc) =   smooth_top_propeller(r,z,angles, nb_angles,time)&
               *smooth_bottom_propeller(r,z,angles, nb_angles,time)
       ELSE IF (if_bottom_prop) THEN
          vv(:,n_loc) =   smooth_bottom_propeller(r,z,angles, nb_angles,time)
       ELSE
          vv(:,n_loc) =   smooth_top_propeller(r,z,angles, nb_angles,time)
       END IF
    END DO
    RETURN

    !===Dummy variables to avoid warning
    n=mesh%np
    !===Dummy variables to avoid warning
  END FUNCTION smooth_penal_in_real_space

  ! DCQ   14/08/2014
  ! 1 - Characteristic Function of  top_propeller
  FUNCTION smooth_top_propeller(r,z,angles,nb_angles,time) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8)                               :: r,theta,z,time
    REAL(KIND=8), DIMENSION(:)                 :: angles
    INTEGER                                    :: nb_angles
    REAL(KIND=8), DIMENSION(nb_angles)         :: vv
    INTEGER                                    :: na
    REAL(KIND=8)                               :: g, a,alphaz,alphar
    REAL(KIND=8)                               :: r0,r1,r2,r3,z0,z1,z2,z3
    REAL(KIND=8)                               :: psi, tanhp,tanhm, tanhd, r_theta

    !Get characteristic function of the supporting disk-cylinder
    psi=smooth_top_supporting_disk(r,z)

    !If we are outside of the supporting disk (respect to r)
    IF(ABS(psi) .LE. 1.d-8) THEN
       DO na = 1, nb_angles
          vv(na) = 1-psi
       END DO
       RETURN
    END IF

    !Do Blades stuff
    r0=hole_rp
    r1=hole_r
    r2=disk_rp
    r3=disk_r

    z0=ztop_p - wjump
    z1=ztop   - wjump
    z2=Tdisk_z_p
    z3=Tdisk_z

    ! Parabolic jump
    IF (z .LE. z1 )  THEN
       alphaz = 1.d0;
    ELSE IF(z .LE. z0 .AND. z .GE. z1) THEN
       alphaz= smooth_jump_down(z,z1,z0);
    ELSE
       alphaz=0.0;
    END IF

    If ( r .LE. r0 )  THEN
       alphar = 0.d0;
    ELSE IF(r .GE. r0 .AND. r .LE. r1) THEN
       alphar= smooth_jump_up(r,r0,r1);
    ELSE
       alphar=1.0;
    END IF
    alphaz= alpha_th*alphaz*alphar

    r_theta = ASIN(r/two_rp)

    DO na = 1, nb_angles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !DCQ go backwards and do the test
       !These blades  rotate the other way (-lbd), (+user%solid_vel*time) and they begin at the same angle,

       theta=   angles(na) + solid_vel*time
       theta=   theta -   top_propeller_angle_offset !some offset to begin at the same vertical axis as the bottom propeller

       !a = theta + (-lbd*r) - FLOOR(( theta + (-lbd)*r)/(twopi/nblades))*twopi/nblades &
       !     - twopi/(2*nblades)

       !JL-CN 01/2015
       a = theta - r_theta - FLOOR(( theta - r_theta)/(twopi/nblades))*twopi/nblades &
            - twopi/(2*nblades)

       tanhp = tanh(alphaz*(r*a+lw+lw*r))
       tanhm = tanh(alphaz*(r*a-lw-lw*r))
       tanhd = tanh(alphaz*(lw+lw*r))
       g=(1+tanhp)*(1-tanhm)/(1+tanhd)**2
       vv(na) = 1-g*psi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END DO
    RETURN
  END FUNCTION smooth_top_propeller

  ! DCQ   14/08/2014
  ! 1 - Characteristic Function of  bottom_propeller
  FUNCTION smooth_bottom_propeller(r,z,angles,nb_angles,time) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8)                               :: r,theta,z,time
    REAL(KIND=8), DIMENSION(:)                 :: angles
    INTEGER                                    :: nb_angles
    REAL(KIND=8), DIMENSION(nb_angles)         :: vv
    INTEGER                                    :: na
    REAL(KIND=8)                               :: g, a,alphaz,alphar
    REAL(KIND=8)                               :: r0,r1,r2,r3,z0,z1,z2,z3
    REAL(KIND=8)                               :: psi, tanhp,tanhm, tanhd, r_theta

    !Supporting disk stuff
    !Get characteristic function of the supporting disk-cylinder
    psi=smooth_bottom_supporting_disk(r,z)

    !If we are outside of the supporting disk (respect to r)
    IF(ABS(psi) .LE. 1.d-8) THEN
       DO na = 1, nb_angles
          vv(na) = 1-psi
       END DO
       RETURN
    END IF

    ! Do blades stuff
    !! Blades with no hole in the disk
    r0=hole_rp
    r1=hole_r
    r2=disk_rp
    r3=disk_r

    z0=zbot_p + wjump
    z1=zbot   + wjump
    z2=Bdisk_z_p
    z3=Bdisk_z

    ! Parabolic jump
    IF (z .LE. z0 )  THEN
       alphaz = 0.d0;
    ELSE IF(z .GE. z0 .AND. z .LE. z1) THEN
       alphaz= smooth_jump_up(z,z0,z1);
    ELSE
       alphaz=1.0;
    END IF

    IF ( r .LE. r0 )  THEN
       alphar = 0.d0;
    ELSE IF(r .GE. r0 .AND. r .LE. r1) THEN
       alphar= smooth_jump_up(r,r0,r1);
    ELSE
       alphar=1.0;
    END IF
    alphaz= alpha_th*alphaz*alphar

    r_theta = ASIN(r/two_rp)

    DO na = 1, nb_angles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !DCQ go backwards and do the test
       theta=   angles(na) - solid_vel*time
       !a = theta + lbd*r - FLOOR(( theta +  lbd*r)/(twopi/nblades))*twopi/nblades &
       !     - twopi/(2*nblades)

       !JL-CN 01/2015
       a = theta + r_theta - FLOOR(( theta +  r_theta)/(twopi/nblades))*twopi/nblades &
            - twopi/(2*nblades)

       tanhp = tanh(alphaz*(r*a+lw+lw*r))
       tanhm = tanh(alphaz*(r*a-lw-lw*r))
       tanhd = tanh(alphaz*(lw+lw*r))
       g=(1+tanhp)*(1-tanhm)/(1+tanhd)**2
       vv(na) = 1-g*psi
    END DO
    RETURN
  END FUNCTION smooth_bottom_propeller

  !Characteristic Function of  top_supporting_disk
  FUNCTION smooth_top_supporting_disk(r,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8)                               :: r,z
    REAL(KIND=8)                               :: vv
    REAL(KIND=8)                               :: r0,r1,r2,r3,z0,z1,z2,z3
    REAL(KIND=8)                               :: curve_1,curve_2,psi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Blades with no hole in the disk
    !Supporting disk stuff
    r0=hole_rp
    r1=hole_r
    r2=disk_rp
    r3=disk_r

    z0=ztop_p - wjump
    z1=ztop   - wjump
    z2=Tdisk_z_p
    z3=Tdisk_z

    psi=0.d0
    IF ( z .LE.  z3 .AND. r .GE. r3)  THEN
       psi = 0.d0 ;
    ELSE IF (r.LE.r0) THEN
       IF(z.GE.z0) THEN
          psi=1.0;
       ELSE IF (z.LE.z1) THEN
          psi=0.0;
       ELSE
          psi=smooth_jump_up(z,z1,z0);
       END IF
    ELSE IF(r.GE.r0 .AND. r.LE.r1) THEN
       curve_2= smooth_jump_up(r,r0,r1)*(z3-z1)+z1;
       curve_1= smooth_jump_up(r,r0,r1)*(z2-z0)+z0;
       IF(z.LE.curve_2) THEN
          psi=0.0;
       ELSE IF(z.GE.curve_1) THEN
          psi=1.0;
       ELSE
          psi = 1 - smooth_jump_up(z ,curve_1,curve_2);
       END IF
    ELSE IF(r.GE.r1 .AND. r.LE.r2) THEN
       IF(z.GE.z2) THEN
          psi=1.0;
       ELSE IF (z.LE.z3) THEN
          psi=0.0;
       ELSE
          psi=smooth_jump_up(z,z3,z2);
       END IF
    ELSE IF (r.GE.r2 .AND. r.LE.r3) THEN
       IF(z.LE.z3) THEN
          psi=0.0;
       ELSE IF(z.GE.z2) THEN
          psi=smooth_jump_down(r,r2,r3);
       ELSE
          psi=smooth_jump_down(r,r2,r3)*smooth_jump_up(z,z3,z2);
       END IF
    END IF

    vv=psi

  END FUNCTION smooth_top_supporting_disk

  !Characteristic Function of  bot_supporting_disk
  FUNCTION smooth_bottom_supporting_disk(r,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8)                               :: r,z
    REAL(KIND=8)                               :: vv
    REAL(KIND=8)                               :: r0,r1,r2,r3,z0,z1,z2,z3
    REAL(KIND=8)                               :: curve_1,curve_2,psi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Blades with no hole in the disk
    !Supporting disk stuff

    r0=hole_rp
    r1=hole_r
    r2=disk_rp
    r3=disk_r

    z0=zbot_p + wjump
    z1=zbot   + wjump
    z2=Bdisk_z_p
    z3=Bdisk_z

    psi=0.d0
    IF ( z .GE.  z3 .AND. r .GE. r3)  THEN
       psi = 0.d0 ;
    ELSE IF (r.LE.r0) THEN
       IF(z.LE.z0) THEN
          psi=1.0;
       ELSE IF (z.GE.z1) THEN
          psi=0.0;
       ELSE
          psi=smooth_jump_up(z,z1,z0);
       END IF
    ELSE IF(r.GE.r0 .AND. r.LE.r1) THEN
       curve_2= smooth_jump_up(r,r0,r1)*(z3-z1)+z1;
       curve_1= smooth_jump_up(r,r0,r1)*(z2-z0)+z0;
       IF(z.GE.curve_2) THEN
          psi=0.0;
       ELSE IF(z.LE.curve_1) THEN
          psi=1.0;
       ELSE
          psi = 1 - smooth_jump_up(z ,curve_1,curve_2);
       END IF
    ELSE IF(r.GE.r1 .AND. r.LE.r2) THEN
       IF(z.LE.z2) THEN
          psi=1.0;
       ELSE IF (z.GE.z3) THEN
          psi=0.0;
       ELSE
          psi=smooth_jump_up(z,z3,z2);
       END IF
    ELSE IF (r.GE.r2 .AND. r.LE.r3) THEN
       IF(z.GE.z3) THEN
          psi=0.0;
       ELSE IF(z.LE.z2) THEN
          psi=smooth_jump_down(r,r2,r3);
       ELSE
          psi=smooth_jump_down(r,r2,r3)*smooth_jump_up(z,z3,z2);
       END IF
    END IF

    vv=psi

  END FUNCTION smooth_bottom_supporting_disk

  !A cubic profile, which is 1 at x0 and 0 at x1
  FUNCTION smooth_jump_down(x,x0,x1) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8)                               :: x,x0,x1
    REAL(KIND=8)                               :: vv
    REAL(KIND=8)                               :: a0,a1,a2,a3

    !Cubic
    a0 = x1**2*(3*x0-x1)/(x0-x1)**3;
    a1 = -6.0*x0*x1/(x0-x1)**3;
    a2 = (3.0*(x0+x1))/(x0-x1)**3;
    a3 = -2.0/(x0-x1)**3;
    vv = a0+a1*x+a2*x*x + a3*x*x*x
    RETURN
  END FUNCTION smooth_jump_down

  !A cubic profile, which is 0 at x0 and 1 at x1
  FUNCTION smooth_jump_up(x,x0,x1) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8)                               :: x,x0,x1
    REAL(KIND=8)                               :: vv

    vv = 1.d0 - smooth_jump_down( x,x0,x1 );
    RETURN
  END FUNCTION smooth_jump_up

  !derivative with respect to x
  FUNCTION Dsmooth_jump_down(x,x0,x1) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8)                               :: x,x0,x1
    REAL(KIND=8)                               :: vv
    REAL(KIND=8)                               :: a0,a1,a2,a3

    !Cubic Factorized
    a0 = x1**2*(3*x0-x1)/(x0-x1)**3;
    a1 = -6.0*x0*x1/(x0-x1)**3;
    a2 = (3.0*(x0+x1))/(x0-x1)**3;
    a3 = -2.0/(x0-x1)**3;

    vv = a1+2.d0*a2*x + 3.d0*a3*x*x
    RETURN
  END FUNCTION Dsmooth_jump_down

  !derivative with respect to x
  FUNCTION Dsmooth_jump_up(x,x0,x1) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8)                               :: x,x0,x1
    REAL(KIND=8)                               :: vv

    vv =  - Dsmooth_jump_down( x,x0,x1 );
    RETURN
  END FUNCTION Dsmooth_jump_up

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
