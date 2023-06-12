!
!Authors Jean-Luc Guermond, Raphael Laguerre, Copyrights 2005
!Revised: Jean-Luc Guermond Jan 2011
!
MODULE tn_axi

CONTAINS

  FUNCTION dot_product_SF(communicator, mesh, list_mode, v, w) RESULT(norm)
    USE def_type_mesh
    USE chaine_caractere
    USE fem_tn_NS_MHD
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: mesh !type de maillage
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) :: v, w
    REAL(KIND=8) :: norm_loc, norm_tot, norm
    INTEGER  :: code
    REAL(KIND=8) :: pi
    MPI_Comm, DIMENSION(2)      :: communicator

    norm = 0.d0
    norm_tot = 0.d0
    IF (mesh%me==0) THEN
       norm_loc = 0.d0
    ELSE
       norm_loc = dot_product_champ(mesh, list_mode, v, w)
    END IF

    CALL MPI_ALLREDUCE(norm_loc,norm_tot,1,MPI_DOUBLE_PRECISION, MPI_SUM, communicator(2), code)
    CALL MPI_ALLREDUCE(norm_tot,norm,1,MPI_DOUBLE_PRECISION, MPI_SUM, communicator(1), code)
    pi = ACOS(-1.d0)
    norm = norm*2*pi

  END FUNCTION dot_product_SF

  FUNCTION norm_SF(communicator, norm_type, mesh, list_mode, v) RESULT(norm)
    USE def_type_mesh
    USE chaine_caractere
    USE fem_tn_NS_MHD
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: mesh !type de maillage
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) :: v
    CHARACTER(*),                    INTENT(IN) :: norm_type
    REAL(KIND=8) :: norm_loc, norm_tot, norm
    INTEGER  :: deb, fin, code
    REAL(KIND=8) :: pi
    MPI_Comm, DIMENSION(2)      :: communicator

    norm = 0.d0
    norm_tot = 0.d0
    IF (mesh%me==0) THEN
       norm_loc = 0.d0
    ELSE

       deb = start_of_string (norm_type)
       fin = last_of_string (norm_type)
       IF (norm_type(deb:fin)=='L2') THEN
          norm_loc = norme_L2_champ(mesh, list_mode, v)
       ELSE IF (norm_type(deb:fin)=='H1') THEN
          norm_loc = SQRT(norme_L2_champ(mesh, list_mode, v)**2 &
               + norme_H1_champ(mesh, list_mode, v)**2)
       ELSE IF (norm_type(deb:fin)=='sH1') THEN
          norm_loc = norme_H1_champ(mesh, list_mode, v)
       ELSE IF (norm_type(deb:fin)=='div') THEN
          norm_loc = norme_div(mesh, list_mode, v)
       ELSE IF (norm_type(deb:fin)=='curl') THEN
          norm_loc = norme_curl(mesh, list_mode, v)
       ELSE
          WRITE(*,*) ' BUG in norm, norm_type not programmed yet' , norm_type(deb:fin)
          STOP
       END IF
    END IF
    CALL MPI_ALLREDUCE(norm_loc**2,norm_tot,1,MPI_DOUBLE_PRECISION, MPI_SUM, communicator(2), code)
    CALL MPI_ALLREDUCE(norm_tot,norm,1,MPI_DOUBLE_PRECISION, MPI_SUM, communicator(1), code)
    pi = ACOS(-1.d0)
    norm = SQRT(norm*2*pi)

  END FUNCTION norm_SF


  FUNCTION norm_S(communicator, norm_type, mesh, list_mode, v) RESULT(norm)
    USE def_type_mesh
    USE chaine_caractere
    USE fem_tn_NS_MHD
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: mesh !type de maillage
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) :: v
    CHARACTER(*),                    INTENT(IN) :: norm_type
    REAL(KIND=8) :: norm_loc, norm_tot, norm
    INTEGER  :: deb, fin, code
    REAL(KIND=8) :: pi
    MPI_Comm, DIMENSION(2)      :: communicator

    norm = 0.d0
    norm_tot = 0.d0
    IF (mesh%me==0) THEN
       norm_loc = 0.d0
    ELSE

       deb = start_of_string (norm_type)
       fin = last_of_string (norm_type)
       IF (norm_type(deb:fin)=='L2') THEN
          norm_loc = norme_L2_champ(mesh, list_mode, v)
       ELSE IF (norm_type(deb:fin)=='H1') THEN
          norm_loc = SQRT(norme_L2_champ(mesh, list_mode, v)**2 &
               + norme_H1_champ(mesh, list_mode, v)**2)
       ELSE IF (norm_type(deb:fin)=='sH1') THEN
          norm_loc = norme_H1_champ(mesh, list_mode, v)
       ELSE IF (norm_type(deb:fin)=='div') THEN
          norm_loc = norme_div(mesh, list_mode, v)
       ELSE IF (norm_type(deb:fin)=='curl') THEN
          norm_loc = norme_curl(mesh, list_mode, v)
       ELSE
          WRITE(*,*) ' BUG in norm, norm_type not programmed yet' , norm_type(deb:fin)
          STOP
       END IF
    END IF
    CALL MPI_ALLREDUCE(norm_loc**2,norm,1,MPI_DOUBLE_PRECISION, MPI_SUM, communicator(1), code)
    pi = ACOS(-1.d0)
    norm = SQRT(norm*2*pi)

  END FUNCTION norm_S

  !DCQ, compute S_L1_norm using only the zero mode
  FUNCTION norm_S_L1_zero_mode(communicator,mesh, list_mode, v) RESULT(norm)
    USE def_type_mesh
    USE chaine_caractere
    USE fem_tn_NS_MHD
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: mesh !type de maillage
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) :: v
    REAL(KIND=8) :: norm_loc, norm_tot, norm
    INTEGER  :: code
    INTEGER, DIMENSION(1) :: zero_mode
    REAL(KIND=8) :: pi
    MPI_Comm, DIMENSION(2)      :: communicator

    norm = 0.d0
    norm_tot = 0.d0
    IF (mesh%me==0) THEN
       norm_loc = 0.d0
    ELSE
       IF (MINVAL(list_mode)==0) THEN !Just mode zero
          zero_mode = MINLOC(list_mode)
          norm_loc = norme_L1_one_mode(mesh, zero_mode(1), v)
       ELSE
          norm_loc=0;
       END IF
    END IF
    CALL MPI_ALLREDUCE(norm_loc,norm,1,MPI_DOUBLE_PRECISION, MPI_SUM, communicator(1), code)
    pi = ACOS(-1.d0)
    norm =(norm*2*pi)
  END FUNCTION norm_S_L1_zero_mode

  SUBROUTINE integration_mode(communicator,norm_loc,norm_tot)
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN)  :: norm_loc
    REAL(KIND=8), INTENT(OUT) :: norm_tot
    INTEGER                   :: code
    REAL(KIND=8) :: pi
    MPI_Comm                  :: communicator

    pi = ACOS(-1.d0)
    CALL MPI_ALLREDUCE(norm_loc,norm_tot,1,MPI_DOUBLE_PRECISION, MPI_SUM, communicator, code)
    !CN-AR Tue Jan 13 2009
    norm_tot = norm_tot*(2*pi)
  END SUBROUTINE integration_mode


  FUNCTION norme_L2_champ_par(communicator, mesh, list_mode, v) RESULT(norm)
    USE def_type_mesh
    USE fem_tn_NS_MHD
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: mesh !type de maillage
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) :: v
    REAL(KIND=8) :: norm_loc, norm
    MPI_Comm                  :: communicator

    IF (mesh%me==0) THEN
       norm_loc = 0.d0
    ELSE
       norm_loc = norme_L2_champ(mesh, list_mode, v)
    ENDIF
    CALL integration_mode(communicator,norm_loc**2, norm)
    norm = SQRT(norm)
  END FUNCTION norme_L2_champ_par

  FUNCTION norme_H1_champ_par(communicator, mesh, list_mode, v) RESULT(norm)
    USE def_type_mesh
    USE fem_tn_NS_MHD
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: mesh !type de maillage
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) :: v
    REAL(KIND=8) :: norm_loc, norm
    MPI_Comm                  :: communicator

    norm_loc = norme_H1_champ(mesh, list_mode, v)
    CALL integration_mode(communicator,norm_loc**2, norm)
    norm = SQRT(norm)
  END FUNCTION norme_H1_champ_par

  FUNCTION norme_div_par(communicator, H_mesh, list_mode, v) RESULT(norm)
    USE def_type_mesh
    USE fem_tn_NS_MHD
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: H_mesh !type de maillage
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) :: v
    REAL(KIND=8) :: norm_loc, norm
    MPI_Comm                  :: communicator

    norm_loc = norme_div(H_mesh, list_mode, v)
    CALL integration_mode(communicator,norm_loc**2, norm)
    norm = SQRT(norm)
  END FUNCTION norme_div_par

  FUNCTION norme_curl_par(communicator, H_mesh, list_mode, v) RESULT(norm)
    USE def_type_mesh
    USE fem_tn_NS_MHD
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: H_mesh !type de maillage
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) :: v
    REAL(KIND=8) :: norm_loc, norm
    MPI_Comm                  :: communicator

    norm_loc = norme_curl(H_mesh, list_mode, v)
    CALL integration_mode(communicator,norm_loc**2, norm)
    norm = SQRT(norm)
  END FUNCTION norme_curl_par

  SUBROUTINE angular_momentum(mesh, list_mode, field, moments_out)
    USE def_type_mesh
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN)  :: mesh
    INTEGER, DIMENSION(:),           INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: field
    REAL(KIND=8), DIMENSION(3),      INTENT(OUT) :: moments_out

    REAL(KIND=8), DIMENSION(3) :: moments_loc
    REAL(KIND=8)               :: m_x, m_y, ray, zed, urc, urs, utc, uts, uzc, uzs
    INTEGER                    :: m, l, code, i

    moments_loc=0.d0
    moments_out=0.d0
    DO i = 1, SIZE(list_mode)
       IF (list_mode(i)==0) THEN
          DO m = 1, mesh%me
             DO l = 1, mesh%gauss%l_G
                ray = SUM(mesh%rr(1,mesh%jj(:,m))*mesh%gauss%ww(:,l))
                utc = SUM(field(mesh%jj(:,m),3,i)*mesh%gauss%ww(:,l))

                moments_loc(3) = moments_loc(3) + ray**2*utc*mesh%gauss%rj(l,m)
             END DO
          END DO
       ELSE IF (list_mode(i)==1) THEN
          DO m = 1, mesh%me
             DO l = 1, mesh%gauss%l_G
                ray = SUM(mesh%rr(1,mesh%jj(:,m))*mesh%gauss%ww(:,l))
                zed = SUM(mesh%rr(2,mesh%jj(:,m))*mesh%gauss%ww(:,l))

                urc = SUM(field(mesh%jj(:,m),1,i)*mesh%gauss%ww(:,l))
                urs = SUM(field(mesh%jj(:,m),2,i)*mesh%gauss%ww(:,l))
                utc = SUM(field(mesh%jj(:,m),3,i)*mesh%gauss%ww(:,l))
                uts = SUM(field(mesh%jj(:,m),4,i)*mesh%gauss%ww(:,l))
                uzc = SUM(field(mesh%jj(:,m),5,i)*mesh%gauss%ww(:,l))
                uzs = SUM(field(mesh%jj(:,m),6,i)*mesh%gauss%ww(:,l))

                m_x = ray*(-zed*(urs+utc)+ray*uzs)
                m_y = ray*(zed*(urc-uts)-ray*uzc)

                moments_loc(1) = moments_loc(1) + m_x*mesh%gauss%rj(l,m)
                moments_loc(2) = moments_loc(2) + m_y*mesh%gauss%rj(l,m)
             END DO
          END DO
       END IF
    END DO

    CALL MPI_ALLREDUCE(moments_loc, moments_out, 3, MPI_DOUBLE_PRECISION, MPI_SUM, PETSC_COMM_WORLD, code)
    moments_out(1) = ACOS(-1.d0)*moments_out(1)
    moments_out(2) = ACOS(-1.d0)*moments_out(2)
    moments_out(3) = 2*ACOS(-1.d0)*moments_out(3)

  END SUBROUTINE angular_momentum

  SUBROUTINE moments(communicator, H_mesh, list_mode, v, dipole_out, quadripole_out)
    USE def_type_mesh
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: H_mesh
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) :: v
    REAL(KIND=8), DIMENSION(3),      INTENT(OUT):: dipole_out
    REAL(KIND=8), DIMENSION(3,3),    INTENT(OUT):: quadripole_out
    REAL(KIND=8), DIMENSION(3)                  :: dipole
    REAL(KIND=8), DIMENSION(3,3)                :: quadripole
    REAL(KIND=8)         :: jr, ray, zed, pi
    INTEGER              :: m_max_c, mode, k, m, l, ni, i
    INTEGER              :: code
    REAL(KIND=8), DIMENSION(6) :: c

    MPI_Comm                  :: communicator

    pi = ACOS(-1.d0)
    m_max_c = SIZE(list_mode)
    dipole = 0
    quadripole = 0

    IF (SIZE(v,1)/=H_mesh%np .OR. SIZE(v,2)/=6 .OR. SIZE(v,3)/=m_max_c ) THEN
       WRITE(*,*) ' BUG in MOMENTS', SIZE(v,1), H_mesh%np
       STOP
    END IF
    DO m = 1, H_mesh%me
       DO l = 1, H_mesh%gauss%l_G

          !--------On calcule le rayon et z du point gauss
          ray = 0
          zed = 0
          DO ni = 1, H_mesh%gauss%n_w;  i = H_mesh%jj(ni,m)
             ray = ray + H_mesh%rr(1,i)*H_mesh%gauss%ww(ni,l)
             zed = zed + H_mesh%rr(2,i)*H_mesh%gauss%ww(ni,l)
          END DO
          jr = ray * H_mesh%gauss%rj(l,m)

          DO k=1, m_max_c
             mode = list_mode(k)
             IF (mode /=0 .AND. mode /=1 .AND. mode /=2) CYCLE

             !--------Compute Curl------
             c = 0
             DO ni = 1,H_mesh%gauss%n_w; i = H_mesh%jj(ni,m)
                !--------Composante r------
                c(1) = c(1) + ( mode/ray*v(i,6,k)*H_mesh%gauss%ww(ni,l) &
                     - v(i,3,k)*H_mesh%gauss%dw(2,ni,l,m))
                c(2) = c(2) + (-mode/ray*v(i,5,k)*H_mesh%gauss%ww(ni,l) &
                     - v(i,4,k)*H_mesh%gauss%dw(2,ni,l,m))
                !--------Composante theta------
                c(3) = c(3) + (v(i,1,k)*H_mesh%gauss%dw(2,ni,l,m) &
                     - v(i,5,k)*H_mesh%gauss%dw(1,ni,l,m))
                c(4) = c(4) + (v(i,2,k)*H_mesh%gauss%dw(2,ni,l,m) &
                     - v(i,6,k)*H_mesh%gauss%dw(1,ni,l,m))
                !--------Composante z------
                c(5) = c(5) + (v(i,3,k)*H_mesh%gauss%dw(1,ni,l,m) &
                     + v(i,3,k)*H_mesh%gauss%ww(ni,l)/ray &
                     - mode/ray*v(i,2,k)*H_mesh%gauss%ww(ni,l))
                c(6) = c(6) + (v(i,4,k)*H_mesh%gauss%dw(1,ni,l,m) &
                     + v(i,4,k)*H_mesh%gauss%ww(ni,l)/ray &
                     + mode/ray*v(i,1,k)*H_mesh%gauss%ww(ni,l))
             ENDDO

             !--------Compute dipole and quadripole------
             IF (mode == 0) THEN
                dipole(3) = dipole(3) + 2*pi*ray*c(3)*jr
                quadripole(1,1) = quadripole(1,1) + pi*(-zed*ray*c(3))*jr
                quadripole(1,2) = quadripole(1,2) + pi*(-zed*ray*c(1)+ray*ray*c(5))*jr
                quadripole(2,1) = quadripole(2,1) + pi*( zed*ray*c(1)-ray*ray*c(5))*jr
                quadripole(2,2) = quadripole(2,2) + pi*(-zed*ray*c(3))*jr
                quadripole(3,3) = quadripole(3,3)+2*pi*( zed*ray*c(3))*jr
             ELSE IF (mode == 1) THEN
                dipole(1) = dipole(1) + pi*(-zed*c(3) -zed*c(2) +ray*c(6))*jr
                dipole(2) = dipole(2) + pi*(-zed*c(4) +zed*c(1) -ray*c(5))*jr
                quadripole(1,3) = quadripole(1,3) + pi*(-zed*c(3)-zed*c(2)+ray*c(6))*zed*jr
                quadripole(2,3) = quadripole(2,3) + pi*(-zed*c(4)+zed*c(1)-ray*c(5))*zed*jr
                quadripole(3,1) = quadripole(3,1) + pi*(ray*ray*c(3))*jr
                quadripole(3,2) = quadripole(3,2) + pi*(ray*ray*c(4))*jr
             ELSE IF (mode == 2) THEN
                quadripole(1,1) = quadripole(1,1) + pi*(-zed*c(3)-zed*c(2)+ray*c(6))*ray*jr/2
                quadripole(1,2) = quadripole(1,2) + pi*(-zed*c(4)+zed*c(1)-ray*c(5))*ray*jr/2
                quadripole(2,1) = quadripole(2,1) + pi*(-zed*c(4)+zed*c(1)-ray*c(5))*ray*jr/2
                quadripole(2,2) = quadripole(2,2) + pi*( zed*c(3)+zed*c(2)-ray*c(6))*ray*jr/2
             END IF

          END DO
       END DO
    END DO

    !--------Collect from everybody------
    CALL MPI_ALLREDUCE(dipole,    dipole_out,    3,MPI_DOUBLE_PRECISION, &
         MPI_SUM, communicator, code)
    CALL MPI_ALLREDUCE(quadripole,quadripole_out,9,MPI_DOUBLE_PRECISION, &
         MPI_SUM, communicator, code)
    RETURN
  END SUBROUTINE moments

END MODULE tn_axi
