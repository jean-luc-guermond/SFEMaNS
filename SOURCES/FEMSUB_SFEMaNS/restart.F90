!
!Authors Jean-Luc Guermond, Raphael Laguerre, Caroline Nore, Copyrights 2005
!
MODULE restart
#include "petsc/finclude/petsc.h"
  USE petsc
CONTAINS

  SUBROUTINE write_restart_ns(communicator, vv_mesh, pp_mesh, time, list_mode, &
       un, un_m1, pn, pn_m1, incpn, incpn_m1, filename, it, freq_restart, &
       opt_level_set, opt_level_set_m1, opt_max_vel, opt_mono, opt_dt)

    USE def_type_mesh
    USE chaine_caractere
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                    :: vv_mesh, pp_mesh
    REAL(KIND=8),                                   INTENT(IN) :: time
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: communicator
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(IN) :: un, un_m1
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(IN) :: pn, pn_m1, incpn, incpn_m1
    REAL(KIND=8), DIMENSION(:,:,:,:), OPTIONAL,     INTENT(IN) :: opt_level_set, opt_level_set_m1
    REAL(KIND=8),                     OPTIONAL,     INTENT(IN) :: opt_max_vel, opt_dt
    LOGICAL, OPTIONAL,                              INTENT(IN) :: opt_mono
    CHARACTER(len=200),                             INTENT(IN) :: filename
    INTEGER,                                        INTENT(IN) :: it, freq_restart
    INTEGER                           :: code, n, i, rang_S, rang_F, nb_procs_S, nb_procs_F
    INTEGER                           :: l, lblank
    CHARACTER(len=3)                  :: tit, tit_S
    LOGICAL                           :: mono=.FALSE.
    LOGICAL                           :: skip
    CHARACTER(len=250)                :: out_name

    CALL MPI_COMM_SIZE(communicator(1),nb_procs_S,code)
    CALL MPI_COMM_SIZE(communicator(2),nb_procs_F,code)
    CALL MPI_COMM_RANK(communicator(1),rang_S,code)
    CALL MPI_COMM_RANK(communicator(2),rang_F,code)

    WRITE(tit,'(i3)') it/freq_restart
    lblank = eval_blank(3,tit)
    DO l = 1, lblank - 1
       tit(l:l) = '0'
    END DO
    WRITE(tit_S,'(i3)') rang_S
    lblank = eval_blank(3,tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO

    IF (PRESENT(opt_mono)) THEN
       mono = opt_mono
    END IF

    IF (mono) THEN
       out_name = 'suite_ns_I'//tit//'.'//filename
    ELSE
       out_name = 'suite_ns_S'//tit_S//'_I'//tit//'.'//filename
    END IF

    skip = (mono .AND. rang_S /= 0)

    DO n = 1, nb_procs_F
       IF ( (rang_F == n-1) .AND. (.NOT. skip) ) THEN
          IF (rang_F == 0) THEN
             OPEN(UNIT = 10, FILE = out_name, POSITION='append', &
                  FORM = 'unformatted', STATUS = 'replace')
             IF (PRESENT(opt_dt)) THEN
                IF (mono) THEN
                   WRITE(10) time, vv_mesh%np , pp_mesh%np , nb_procs_F, SIZE(list_mode), opt_dt
                ELSE
                   WRITE(10) time, nb_procs_S, nb_procs_F, SIZE(list_mode), opt_dt
                END IF
             ELSE
                IF (mono) THEN
                   WRITE(10) time, vv_mesh%np , pp_mesh%np , nb_procs_F, SIZE(list_mode)
                ELSE
                   WRITE(10) time, nb_procs_S, nb_procs_F, SIZE(list_mode)
                END IF
             END IF
          ELSE
             OPEN(UNIT = 10, FILE = out_name, POSITION='append', &
                  FORM = 'unformatted', STATUS = 'unknown')
          END IF

          DO i= 1, SIZE(list_mode)
             WRITE(10) list_mode(i)
             WRITE(10) un(:,:,i)
             WRITE(10) un_m1(:,:,i)
             WRITE(10) pn(:,:,i)
             WRITE(10) pn_m1(:,:,i)
             WRITE(10) incpn(:,:,i)
             WRITE(10) incpn_m1(:,:,i)
             IF (PRESENT(opt_level_set) .AND. PRESENT(opt_level_set_m1)) THEN
                WRITE(10) opt_level_set(:,:,:,i)
                WRITE(10) opt_level_set_m1(:,:,:,i)
                WRITE(10) opt_max_vel
             END IF
          END DO
          CLOSE(10)
       END IF
       CALL MPI_BARRIER(communicator(2),code)
    END DO

  END SUBROUTINE write_restart_ns

  SUBROUTINE read_restart_ns(communicator, time, list_mode, &
       un, un_m1, pn, pn_m1, incpn, incpn_m1, filename, val_init, interpol, &
       opt_level_set, opt_level_set_m1, opt_max_vel, opt_mono, &
       opt_it, opt_dt) !===HF may 2020

    USE def_type_mesh
    USE chaine_caractere
    USE my_util
    IMPLICIT NONE
    REAL(KIND=8),                                   INTENT(OUT):: time
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: communicator
    INTEGER,      DIMENSION(:)                                 :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(OUT):: un, un_m1
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(OUT):: pn, pn_m1, incpn, incpn_m1
    REAL(KIND=8), DIMENSION(:,:,:,:), OPTIONAL,     INTENT(OUT):: opt_level_set, opt_level_set_m1
    REAL(KIND=8),                     OPTIONAL,     INTENT(OUT):: opt_max_vel
    CHARACTER(len=200),                             INTENT(IN) :: filename
    REAL(KIND=8), OPTIONAL,                         INTENT(IN) :: val_init
    LOGICAL     , OPTIONAL,                         INTENT(IN) :: interpol
    LOGICAL     , OPTIONAL,                         INTENT(IN) :: opt_mono
    INTEGER     , OPTIONAL,                         INTENT(IN) :: opt_it
    REAL(KIND=8),                     OPTIONAL,     INTENT(IN) :: opt_dt
    REAL(KIND=8)                                               :: max_vel_loc, dt_read, dt_ratio
    INTEGER     :: code, n, i, mode, j, rang_S, nb_procs_S, rang_F, nb_procs_F, nlignes, rank
    INTEGER     :: m_max_cr, nb_procs_r, nb_procs_Sr
    INTEGER     :: m_max_c, nb_mode_r, mode_cherche
    LOGICAL     :: trouve, okay
    INTEGER     :: npv, npp
    INTEGER           :: l, lblank
    CHARACTER(len=3)  :: tit_S
    !===HF may 2020
    CHARACTER(len=3)  :: tit
    !===HF may 2020
    LOGICAL           :: mono=.FALSE.
    CHARACTER(len=250):: in_name
    CALL MPI_COMM_RANK(communicator(2),rang_F,code)
    CALL MPI_COMM_SIZE(communicator(2),nb_procs_F,code)
    CALL MPI_COMM_RANK(communicator(1),rang_S,code)
    CALL MPI_COMM_SIZE(communicator(1),nb_procs_S,code)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,code)

    max_vel_loc = 0.d0

    nlignes = 6
    IF (PRESENT(opt_level_set) .AND. PRESENT(opt_level_set_m1)) THEN
       nlignes = nlignes + 3
    END IF

    !=== HF may 2020
    IF (PRESENT(opt_it)) THEN
       WRITE(tit,'(i3)') opt_it
       lblank = eval_blank(3,tit)
       DO l = 1, lblank - 1
          tit(l:l) = '0'
       END DO
    END IF
    !=== HF may 2020

    WRITE(tit_S,'(i3)') rang_S
    lblank = eval_blank(3,tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO

    IF (PRESENT(opt_mono)) THEN
       mono = opt_mono
    END IF

    IF (mono) THEN
       !=== HF may 2020
       IF (PRESENT(opt_it)) THEN
          in_name = 'suite_ns_I'//tit//'.'//filename
       ELSE
          in_name = 'suite_ns.'//filename
       END IF
    ELSE
       IF (PRESENT(opt_it)) THEN
          in_name = 'suite_ns_S'//tit_S//'_I'//tit//'.'//filename
       ELSE
          in_name = 'suite_ns_S'//tit_S//'.'//filename
       END IF
       !=== HF may 2020
    END IF

    !=== HF may 2020
    IF (PRESENT(opt_it)) THEN
       WRITE(*,*) 'restart Navier-Stokes for it', opt_it
    ELSE
       WRITE(*,*) 'restart Navier-Stokes'
    END IF
    !=== HF may 2020
    OPEN(UNIT = 10, FILE = in_name, FORM = 'unformatted', STATUS = 'unknown')

    IF (PRESENT(opt_dt)) THEN
       IF (mono) THEN
          READ(10) time, npv, npp, nb_procs_r, m_max_cr, dt_read
          nb_procs_Sr = -1
       ELSE
          READ(10) time, nb_procs_Sr, nb_procs_r, m_max_cr, dt_read
       END IF
    ELSE
       IF (mono) THEN
          READ(10) time, npv, npp, nb_procs_r, m_max_cr
          nb_procs_Sr = -1
       ELSE
          READ(10) time, nb_procs_Sr, nb_procs_r, m_max_cr
       END IF
    END IF

    CLOSE(10)

    IF ((nb_procs_Sr /= nb_procs_S) .AND. (.NOT. mono)) THEN
       CALL error_petsc('BUG in read_restart: nb_procs_Sr /= nb_procs_S')
       !STOP
    END IF

    IF (rang_F == 0) THEN
       WRITE(*,*) 'File name ', trim(adjustl(in_name))
       WRITE(*,*) 'Time = ', time
       WRITE(*,*) 'Number of processors from restart file = ',nb_procs_r
       WRITE(*,*) 'Number of modes per processor from restart file = ',m_max_cr
    ENDIF

    m_max_c   = SIZE(list_mode)      !nombre de modes par proc pour le calcul
    nb_mode_r = nb_procs_r*m_max_cr  !nombre total de modes contenus dans le suite

    !June 7 2007, JLG
    IF (nb_procs_F*m_max_c /= nb_mode_r) THEN
       !CALL error_petsc('Bug in read_restart_ns: nb_procs_F*m_max_c /= nb_mode_r')
       WRITE(*,*) 'Warning in read_restart_ns: nb_procs_F*m_max_c /= nb_mode_r'
       !STOP
    END IF

    okay = .FALSE.
    IF (PRESENT(interpol)) THEN
       IF (interpol) THEN
          okay =.TRUE.
       END IF
    END IF
    !June 7 2007, JLG

    IF (rank==0) THEN
       WRITE(*,*) 'Reading Navier-Stokes modes ...'
    END IF
    DO i=1, m_max_c                  !pour tout les modes du processeur courant
       !ouverture du fichier
       OPEN(UNIT = 10, FILE = in_name, FORM = 'unformatted', STATUS = 'unknown')
       !on saute la premiere ligne du fichier qui contient des donnees
       READ(10)
       mode_cherche = list_mode(i)
       !recherche du bon mode
       trouve = .FALSE.
       DO j=1, nb_mode_r             !pour tout les modes ecris dans le suite.
          !lecture du mode
          READ(10) mode
          !June 7 2007, JLG
          IF (okay) THEN
             IF (j/=rang_F*m_max_c+i) THEN
                DO n=1, nlignes
                   READ(10)
                ENDDO
                CYCLE
             ELSE
                list_mode(i) = mode
                mode_cherche = mode
             END IF
          END IF
          !June 7 2007, JLG
          IF (mode == mode_cherche) THEN   !on a trouve le bon mode
             READ(10) un(:,:,i)
             READ(10) un_m1(:,:,i)
             READ(10) pn(:,:,i)
             READ(10) pn_m1(:,:,i)
             READ(10) incpn(:,:,i)
             READ(10) incpn_m1(:,:,i)
             IF (PRESENT(opt_level_set) .AND. PRESENT(opt_level_set_m1)) THEN
                READ(10) opt_level_set(:,:,:,i)
                READ(10) opt_level_set_m1(:,:,:,i)
                READ(10) max_vel_loc
             END IF
             WRITE(*,'(A,i4,A)') 'mode ns ', mode_cherche,' found '
             trouve = .TRUE.
             EXIT                        !car on a trouve le bon mode
          ELSE                             !on passe au mode suivant en sautant 6 lignes
             DO n=1, nlignes
                READ(10)
             ENDDO
          ENDIF
       ENDDO

       IF (.NOT.trouve) THEN               !mode_cherche non trouve
          IF (PRESENT(val_init)) THEN ! not implemented yet
             un(:,:,i)    = val_init  ; un_m1(:,:,i)   = val_init
             pn(:,:,i)    = val_init ; pn_m1(:,:,i)    = val_init
             incpn(:,:,i) = val_init ; incpn_m1(:,:,i) = val_init
             IF (PRESENT(opt_level_set) .AND. PRESENT(opt_level_set_m1)) THEN
                opt_level_set(:,:,:,i)    = val_init
                opt_level_set_m1(:,:,:,i) = val_init
                max_vel_loc = val_init
             END IF
             WRITE(*,'(A,i4,A)') 'mode ns', mode_cherche,' not found'
          ELSE
             un(:,:,i)    = 0.d0 ; un_m1(:,:,i)    = 0.d0
             pn(:,:,i)    = 0.d0 ; pn_m1(:,:,i)    = 0.d0
             incpn(:,:,i) = 0.d0 ; incpn_m1(:,:,i) = 0.d0
             IF (PRESENT(opt_level_set) .AND. PRESENT(opt_level_set_m1)) THEN
                opt_level_set(:,:,:,i) = 0.d0
                opt_level_set_m1(:,:,:,i) = 0.d0
             END IF
             WRITE(*,*) 'mode ns', mode_cherche, ' not found'
          ENDIF
       ENDIF
       CLOSE(10)                          !fermeture du fichier suite
    ENDDO

    IF (PRESENT(opt_max_vel)) THEN
       CALL MPI_ALLREDUCE(max_vel_loc, opt_max_vel, 1, MPI_DOUBLE_PRECISION, &
            MPI_MAX, communicator(2), code)
    END IF

    !
    IF (PRESENT(opt_dt)) THEN
       IF (ABS((opt_dt - dt_read)/opt_dt).GT.1d-4) THEN
          dt_ratio = opt_dt/dt_read
          IF (rank==0) THEN
             WRITE(*,*) 'In Navier-Stokes restart, suite_time_step different from inputs%dt ...'
             WRITE(*,*) ' opt_dt, dt_read =', opt_dt, dt_read
          END IF
          un_m1 = dt_ratio * un_m1 +(1.d0 - dt_ratio)* un
          pn_m1 = dt_ratio * pn_m1 +(1.d0 - dt_ratio)* pn
          incpn_m1 = dt_ratio * incpn_m1 +(1.d0 - dt_ratio)* incpn
          IF (PRESENT(opt_level_set) .AND. PRESENT(opt_level_set_m1)) THEN
             opt_level_set_m1 = dt_ratio * opt_level_set_m1 +(1.d0 - dt_ratio)* opt_level_set
          END IF
       END IF
    END IF

  END SUBROUTINE read_restart_ns

  SUBROUTINE write_restart_LES(communicator, vv_mesh, pp_mesh, time, list_mode, &
       filename, it, freq_restart, opt_LES_NS, opt_LES_level, opt_mono, opt_dt)

    USE def_type_mesh
    USE chaine_caractere
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                    :: vv_mesh, pp_mesh
    REAL(KIND=8),                                   INTENT(IN) :: time
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: communicator
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: list_mode
    REAL(KIND=8),                     OPTIONAL,     INTENT(IN) :: opt_dt
    REAL(KIND=8), DIMENSION(:,:,:,:), OPTIONAL,     INTENT(IN) :: opt_LES_NS
    REAL(KIND=8), DIMENSION(:,:,:,:), OPTIONAL,     INTENT(IN) :: opt_LES_level
    LOGICAL, OPTIONAL,                              INTENT(IN) :: opt_mono
    CHARACTER(len=200),                             INTENT(IN) :: filename
    INTEGER,                                        INTENT(IN) :: it, freq_restart
    INTEGER                           :: code, n, i, rang_S, rang_F, nb_procs_S, nb_procs_F
    INTEGER                           :: l, lblank
    CHARACTER(len=3)                  :: tit, tit_S
    LOGICAL                           :: mono=.FALSE.
    LOGICAL                           :: skip
    CHARACTER(len=250)                :: out_name

    CALL MPI_COMM_SIZE(communicator(1),nb_procs_S,code)
    CALL MPI_COMM_SIZE(communicator(2),nb_procs_F,code)
    CALL MPI_COMM_RANK(communicator(1),rang_S,code)
    CALL MPI_COMM_RANK(communicator(2),rang_F,code)

    WRITE(tit,'(i3)') it/freq_restart
    lblank = eval_blank(3,tit)
    DO l = 1, lblank - 1
       tit(l:l) = '0'
    END DO
    WRITE(tit_S,'(i3)') rang_S
    lblank = eval_blank(3,tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO

    IF (PRESENT(opt_mono)) THEN
       mono = opt_mono
    END IF

    IF (mono) THEN
       out_name = 'suite_LES_I'//tit//'.'//filename
    ELSE
       out_name = 'suite_LES_S'//tit_S//'_I'//tit//'.'//filename
    END IF

    skip = (mono .AND. rang_S /= 0)

    DO n = 1, nb_procs_F
       IF ( (rang_F == n-1) .AND. (.NOT. skip) ) THEN
          IF (rang_F == 0) THEN
             OPEN(UNIT = 10, FILE = out_name, POSITION='append', &
                  FORM = 'unformatted', STATUS = 'replace')
             IF (PRESENT(opt_dt)) THEN
                IF (mono) THEN
                   WRITE(10) time, vv_mesh%np , pp_mesh%np , nb_procs_F, SIZE(list_mode), opt_dt
                ELSE
                   WRITE(10) time, nb_procs_S, nb_procs_F, SIZE(list_mode), opt_dt
                END IF
             ELSE
                IF (mono) THEN
                   WRITE(10) time, vv_mesh%np , pp_mesh%np , nb_procs_F, SIZE(list_mode)
                ELSE
                   WRITE(10) time, nb_procs_S, nb_procs_F, SIZE(list_mode)
                END IF
             END IF
          ELSE
             OPEN(UNIT = 10, FILE = out_name, POSITION='append', &
                  FORM = 'unformatted', STATUS = 'unknown')
          END IF

          DO i= 1, SIZE(list_mode)
             WRITE(10) list_mode(i)
             IF (PRESENT(opt_LES_NS)) THEN
                WRITE(10) opt_LES_NS(:,:,:,i)
             END IF
             IF (PRESENT(opt_LES_level)) THEN
                WRITE(10) opt_LES_level(:,:,:,i)
             END IF
          END DO
          CLOSE(10)
       END IF
       CALL MPI_BARRIER(communicator(2),code)
    END DO

  END SUBROUTINE write_restart_LES

  SUBROUTINE read_restart_LES(communicator, time, list_mode, &
       filename, val_init, interpol, opt_LES_NS, opt_LES_level, &
       opt_mono, opt_it, opt_dt)

    USE def_type_mesh
    USE chaine_caractere
    USE my_util
    IMPLICIT NONE
    REAL(KIND=8),                                   INTENT(OUT):: time
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: communicator
    INTEGER,      DIMENSION(:)                                 :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:,:), OPTIONAL,     INTENT(OUT):: opt_LES_NS
    REAL(KIND=8), DIMENSION(:,:,:,:), OPTIONAL,     INTENT(OUT):: opt_LES_level
    CHARACTER(len=200),                             INTENT(IN) :: filename
    REAL(KIND=8), OPTIONAL,                         INTENT(IN) :: val_init
    LOGICAL     , OPTIONAL,                         INTENT(IN) :: interpol
    LOGICAL     , OPTIONAL,                         INTENT(IN) :: opt_mono
    INTEGER     , OPTIONAL,                         INTENT(IN) :: opt_it
    REAL(KIND=8),                     OPTIONAL,     INTENT(IN) :: opt_dt
    REAL(KIND=8)                                               :: max_vel_loc, dt_read, dt_ratio
    INTEGER     :: code, n, i, mode, j, rang_S, nb_procs_S, rang_F, nb_procs_F, nlignes, rank
    INTEGER     :: m_max_cr, nb_procs_r, nb_procs_Sr
    INTEGER     :: m_max_c, nb_mode_r, mode_cherche
    LOGICAL     :: trouve, okay
    INTEGER     :: npv, npp
    INTEGER           :: l, lblank
    CHARACTER(len=3)  :: tit_S
    !===HF may 2020
    CHARACTER(len=3)  :: tit
    !===HF may 2020
    LOGICAL           :: mono=.FALSE.
    CHARACTER(len=250):: in_name
    CALL MPI_COMM_RANK(communicator(2),rang_F,code)
    CALL MPI_COMM_SIZE(communicator(2),nb_procs_F,code)
    CALL MPI_COMM_RANK(communicator(1),rang_S,code)
    CALL MPI_COMM_SIZE(communicator(1),nb_procs_S,code)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,code)

    max_vel_loc = 0.d0

    nlignes = 0
    IF (PRESENT(opt_LES_NS)) THEN
       nlignes = nlignes + 1
    END IF
    IF (PRESENT(opt_LES_level)) THEN
       nlignes = nlignes + 1
    END IF

    !=== HF may 2020
    IF (PRESENT(opt_it)) THEN
       WRITE(tit,'(i3)') opt_it
       lblank = eval_blank(3,tit)
       DO l = 1, lblank - 1
          tit(l:l) = '0'
       END DO
    END IF
    !=== HF may 2020

    WRITE(tit_S,'(i3)') rang_S
    lblank = eval_blank(3,tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO

    IF (PRESENT(opt_mono)) THEN
       mono = opt_mono
    END IF

    IF (mono) THEN
       !=== HF may 2020
       IF (PRESENT(opt_it)) THEN
          in_name = 'suite_LES_I'//tit//'.'//filename
       ELSE
          in_name = 'suite_LES.'//filename
       END IF
    ELSE
       IF (PRESENT(opt_it)) THEN
          in_name = 'suite_LES_S'//tit_S//'_I'//tit//'.'//filename
       ELSE
          in_name = 'suite_LES_S'//tit_S//'.'//filename
       END IF
       !=== HF may 2020
    END IF

    !=== HF may 2020
    IF (PRESENT(opt_it)) THEN
       WRITE(*,*) 'restart LES for it', opt_it
    ELSE
       WRITE(*,*) 'restart LES'
    END IF
    !=== HF may 2020
    OPEN(UNIT = 10, FILE = in_name, FORM = 'unformatted', STATUS = 'unknown')

    IF (PRESENT(opt_dt)) THEN
       IF (mono) THEN
          READ(10) time, npv, npp, nb_procs_r, m_max_cr, dt_read
          nb_procs_Sr = -1
       ELSE
          READ(10) time, nb_procs_Sr, nb_procs_r, m_max_cr, dt_read
       END IF
    ELSE
       IF (mono) THEN
          READ(10) time, npv, npp, nb_procs_r, m_max_cr
          nb_procs_Sr = -1
       ELSE
          READ(10) time, nb_procs_Sr, nb_procs_r, m_max_cr
       END IF
    END IF

    CLOSE(10)

    IF ((nb_procs_Sr /= nb_procs_S) .AND. (.NOT. mono)) THEN
       CALL error_petsc('BUG in read_restart_LES: nb_procs_Sr /= nb_procs_S')
       !STOP
    END IF

    IF (rang_F == 0) THEN
       WRITE(*,*) 'File name ', trim(adjustl(in_name))
       WRITE(*,*) 'Time = ', time
       WRITE(*,*) 'Number of processors from restart file = ',nb_procs_r
       WRITE(*,*) 'Number of modes per processor from restart file = ',m_max_cr
    ENDIF

    m_max_c   = SIZE(list_mode)      !nombre de modes par proc pour le calcul
    nb_mode_r = nb_procs_r*m_max_cr  !nombre total de modes contenus dans le suite

    !June 7 2007, JLG
    IF (nb_procs_F*m_max_c /= nb_mode_r) THEN
       !CALL error_petsc('Bug in read_restart_LES: nb_procs_F*m_max_c /= nb_mode_r')
       WRITE(*,*) 'Warning in read_restart_LES: nb_procs_F*m_max_c /= nb_mode_r'
       !STOP
    END IF

    okay = .FALSE.
    IF (PRESENT(interpol)) THEN
       IF (interpol) THEN
          okay =.TRUE.
       END IF
    END IF
    !June 7 2007, JLG

    IF (rank==0) THEN
       WRITE(*,*) 'Reading LES modes ...'
    END IF
    DO i=1, m_max_c
       !Open restart file
       OPEN(UNIT = 10, FILE = in_name, FORM = 'unformatted', STATUS = 'unknown')
       !Skip first line (only contains data on mesh)
       READ(10)
       mode_cherche = list_mode(i)
       !Find correct Fourier mode
       trouve = .FALSE.
       DO j=1, nb_mode_r
          !Read Fourier mode
          READ(10) mode
          IF (okay) THEN
             IF (j/=rang_F*m_max_c+i) THEN
                DO n=1, nlignes
                   READ(10)
                ENDDO
                CYCLE
             ELSE
                list_mode(i) = mode
                mode_cherche = mode
             END IF
          END IF
          !June 7 2007, JLG
          IF (mode == mode_cherche) THEN
             IF (PRESENT(opt_LES_NS)) THEN 
                READ(10) opt_LES_NS(:,:,:,i)
             END IF
             IF (PRESENT(opt_LES_level)) THEN 
                READ(10) opt_LES_level(:,:,:,i)
             END IF
             WRITE(*,'(A,i4,A)') 'mode ns ', mode_cherche,' found '
             trouve = .TRUE.
             EXIT
          ELSE
             DO n=1, nlignes !go to next mode by skypping nlignes number of lines
                READ(10)
             ENDDO
          ENDIF
       ENDDO

       IF (.NOT.trouve) THEN               !mode_cherche non trouve
          IF (PRESENT(val_init)) THEN ! not implemented yet
             IF (PRESENT(opt_LES_NS)) THEN
                opt_LES_NS(:,:,:,i) = val_init
             END IF
             IF (PRESENT(opt_LES_level)) THEN
                opt_LES_level(:,:,:,i) = val_init
             END IF
             WRITE(*,'(A,i4,A)') 'mode ns', mode_cherche,' not found'
          ELSE
             IF (PRESENT(opt_LES_NS)) THEN
                opt_LES_NS(:,:,:,i) = 0.d0
             END IF
             IF (PRESENT(opt_LES_level)) THEN
                opt_LES_level(:,:,:,i) = 0.d0
             END IF
             WRITE(*,*) 'mode ns', mode_cherche, ' not found'
          ENDIF
       ENDIF
       CLOSE(10)                          !fermeture du fichier suite
    ENDDO

    IF (PRESENT(opt_dt)) THEN
       IF (ABS((opt_dt - dt_read)/opt_dt).GT.1d-4) THEN
          dt_ratio = opt_dt/dt_read
          IF (rank==0) THEN
             WRITE(*,*) 'In LES restart, suite_time_step different from inputs%dt ...'
             WRITE(*,*) 'opt_dt, dt_read =', opt_dt, dt_read
          END IF
          IF (PRESENT(opt_LES_NS)) THEN
             opt_LES_NS = dt_ratio * opt_LES_NS +(1.d0 - dt_ratio)* opt_LES_NS
          END IF
          IF (PRESENT(opt_LES_level)) THEN
             opt_LES_level = dt_ratio * opt_LES_level +(1.d0 - dt_ratio)* opt_LES_level
          END IF
       END IF
    END IF

  END SUBROUTINE read_restart_LES

  SUBROUTINE write_restart_ns_taylor(communicator, vv_mesh, pp_mesh, time, list_mode, &
       un, der_un, pn, der_pn, filename, it, freq_restart, &
       opt_level_set, opt_level_set_m1, opt_max_vel, opt_mono)
    USE input_data
    USE def_type_mesh
    USE chaine_caractere
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                    :: vv_mesh,pp_mesh
    REAL(KIND=8),                                   INTENT(IN) :: time
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: communicator
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: list_mode
    TYPE(dyn_real_array_three), DIMENSION(:),       INTENT(IN) :: der_un
    TYPE(dyn_real_array_three), DIMENSION(:),       INTENT(IN) :: der_pn
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(IN) :: pn
    REAL(KIND=8), DIMENSION(:,:,:,:), OPTIONAL,     INTENT(IN) :: opt_level_set, opt_level_set_m1
    REAL(KIND=8),                     OPTIONAL,     INTENT(IN) :: opt_max_vel
    LOGICAL, OPTIONAL,                              INTENT(IN) :: opt_mono
    CHARACTER(len=200),                             INTENT(IN) :: filename
    INTEGER,                                        INTENT(IN) :: it, freq_restart
    INTEGER                           :: code, n, i, rang_S, rang_F, nb_procs_S, nb_procs_F
    INTEGER                           :: l, lblank, kp
    CHARACTER(len=3)                  :: tit, tit_S
    LOGICAL                           :: mono=.FALSE.
    LOGICAL                           :: skip
    CHARACTER(len=250)                :: out_name

    CALL MPI_COMM_SIZE(communicator(1),nb_procs_S,code)
    CALL MPI_COMM_SIZE(communicator(2),nb_procs_F,code)
    CALL MPI_COMM_RANK(communicator(1),rang_S,code)
    CALL MPI_COMM_RANK(communicator(2),rang_F,code)

    WRITE(tit,'(i3)') it/freq_restart
    lblank = eval_blank(3,tit)
    DO l = 1, lblank - 1
       tit(l:l) = '0'
    END DO
    WRITE(tit_S,'(i3)') rang_S
    lblank = eval_blank(3,tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO

    IF (PRESENT(opt_mono)) THEN
       mono = opt_mono
    END IF

    IF (mono) THEN
       out_name = 'suite_ns_I'//tit//'.'//filename
    ELSE
       out_name = 'suite_ns_S'//tit_S//'_I'//tit//'.'//filename
    END IF

    skip = (mono .AND. rang_S /= 0)

    DO n = 1, nb_procs_F
       IF ( (rang_F == n-1) .AND. (.NOT. skip) ) THEN
          IF (rang_F == 0) THEN
             OPEN(UNIT = 10, FILE = out_name, POSITION='append', &
                  FORM = 'unformatted', STATUS = 'replace')
             IF (mono) THEN
                WRITE(10) time, vv_mesh%np , pp_mesh%np , nb_procs_F, SIZE(list_mode), inputs%taylor_order
             ELSE
                WRITE(10) time, nb_procs_S, nb_procs_F, SIZE(list_mode), inputs%taylor_order
             END IF
          ELSE
             OPEN(UNIT = 10, FILE = out_name, POSITION='append', &
                  FORM = 'unformatted', STATUS = 'unknown')
          END IF
          DO i= 1, SIZE(list_mode)
             WRITE(10) list_mode(i)
             WRITE(10) un(:,:,i)
             DO kp = 1, inputs%taylor_order-1
                WRITE(10) der_un(kp)%DRT(:,:,i)
             END DO
             WRITE(10) pn(:,:,i)
             DO kp = 1, inputs%taylor_order-1
                WRITE(10) der_pn(kp)%DRT(:,:,i)
             END DO
             IF (PRESENT(opt_level_set) .AND. PRESENT(opt_level_set_m1)) THEN
                WRITE(10) opt_level_set(:,:,:,i)
                WRITE(10) opt_level_set_m1(:,:,:,i)
                WRITE(10) opt_max_vel
             END IF
          END DO
          CLOSE(10)
       END IF
       CALL MPI_BARRIER(communicator(2),code)
    END DO

  END SUBROUTINE write_restart_ns_taylor

  SUBROUTINE read_restart_ns_taylor(communicator, time, list_mode, &
       un, der_un, pn, der_pn, filename, val_init, interpol, &
       opt_level_set, opt_level_set_m1, opt_max_vel, opt_mono, &
       opt_it) !===HF may 2020
    USE input_data
    USE def_type_mesh
    USE chaine_caractere
    USE my_util
    IMPLICIT NONE
    REAL(KIND=8),                                   INTENT(OUT):: time
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: communicator
    INTEGER,      DIMENSION(:)                                 :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(OUT):: un
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(OUT):: pn
    TYPE(dyn_real_array_three), DIMENSION(:),       INTENT(OUT):: der_un
    TYPE(dyn_real_array_three), DIMENSION(:),       INTENT(OUT):: der_pn
    REAL(KIND=8), DIMENSION(:,:,:,:), OPTIONAL,     INTENT(OUT):: opt_level_set, opt_level_set_m1
    REAL(KIND=8),                     OPTIONAL,     INTENT(OUT):: opt_max_vel
    CHARACTER(len=200),                             INTENT(IN) :: filename
    REAL(KIND=8), OPTIONAL,                         INTENT(IN) :: val_init
    LOGICAL     , OPTIONAL,                         INTENT(IN) :: interpol
    LOGICAL     , OPTIONAL,                         INTENT(IN) :: opt_mono
    INTEGER     , OPTIONAL,                         INTENT(IN) :: opt_it
    REAL(KIND=8)                                               :: max_vel_loc
    INTEGER     :: code, n, i, mode, j, rang_S, nb_procs_S, rang_F, nb_procs_F, nlignes, rank
    INTEGER     :: m_max_cr, nb_procs_r, nb_procs_Sr
    INTEGER     :: m_max_c, nb_mode_r, mode_cherche, taylor_order, taylor_order_min
    LOGICAL     :: trouve, okay
    INTEGER     :: npv, npp, kp
    INTEGER           :: l, lblank
    CHARACTER(len=3)  :: tit_S
    !===HF may 2020
    CHARACTER(len=3)  :: tit
    !===HF may 2020
    LOGICAL           :: mono=.FALSE.
    CHARACTER(len=250):: in_name
    CALL MPI_COMM_RANK(communicator(2),rang_F,code)
    CALL MPI_COMM_SIZE(communicator(2),nb_procs_F,code)
    CALL MPI_COMM_RANK(communicator(1),rang_S,code)
    CALL MPI_COMM_SIZE(communicator(1),nb_procs_S,code)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,code)

    max_vel_loc = 0.d0

    WRITE(tit_S,'(i3)') rang_S
    lblank = eval_blank(3,tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO

    !=== HF may 2020
    IF (PRESENT(opt_it)) THEN
       WRITE(tit,'(i3)') opt_it
       lblank = eval_blank(3,tit)
       DO l = 1, lblank - 1
          tit(l:l) = '0'
       END DO
    END IF
    !=== HF may 2020

    IF (PRESENT(opt_mono)) THEN
       mono = opt_mono
    END IF

    IF (mono) THEN
       !===HF may 2020
       IF (PRESENT(opt_it)) THEN
          in_name = 'suite_ns_I'//tit//'.'//filename
       ELSE
          in_name = 'suite_ns.'//filename
       END IF
    ELSE
       IF (PRESENT(opt_it)) THEN
          in_name = 'suite_ns_S'//tit_S//'_I'//tit//'.'//filename
       ELSE
          in_name = 'suite_ns_S'//tit_S//'.'//filename
       END IF
       !===HF may 2020
    END IF

    !=== HF may 2020
    IF (PRESENT(opt_it)) THEN
       WRITE(*,*) 'restart Navier-Stokes for it', opt_it
    ELSE
       WRITE(*,*) 'restart Navier-Stokes'
    END IF
    !=== HF may 2020

    OPEN(UNIT = 10, FILE = in_name, FORM = 'unformatted', STATUS = 'unknown')

    IF (mono) THEN
       READ(10) time, npv, npp, nb_procs_r, m_max_cr, taylor_order
       nb_procs_Sr = -1
    ELSE
       READ(10) time, nb_procs_Sr, nb_procs_r, m_max_cr, taylor_order
    END IF
    CLOSE(10)
    taylor_order_min =  MIN(inputs%taylor_order,taylor_order)

    nlignes = 2*taylor_order
    IF (PRESENT(opt_level_set) .AND. PRESENT(opt_level_set_m1)) THEN
       nlignes = nlignes + 3
    END IF

    IF ((nb_procs_Sr /= nb_procs_S) .AND. (.NOT. mono)) THEN
       CALL error_petsc('BUG in read_restart: nb_procs_Sr /= nb_procs_S')
       !STOP
    END IF

    IF (rang_F == 0) THEN
       WRITE(*,*) 'File name ', trim(adjustl(in_name))
       WRITE(*,*) 'Time = ', time
       WRITE(*,*) 'Number of processors from restart file = ',nb_procs_r
       WRITE(*,*) 'Number of modes per processor from restart file = ',m_max_cr
    ENDIF

    m_max_c   = SIZE(list_mode)      !nombre de modes par proc pour le calcul
    nb_mode_r = nb_procs_r*m_max_cr  !nombre total de modes contenus dans le suite

    !June 7 2007, JLG
    IF (nb_procs_F*m_max_c /= nb_mode_r) THEN
       !CALL error_petsc('Bug in read_restart_ns: nb_procs_F*m_max_c /= nb_mode_r')
       WRITE(*,*) 'Warning in read_restart_ns: nb_procs_F*m_max_c /= nb_mode_r'
       !STOP
    END IF

    okay = .FALSE.
    IF (PRESENT(interpol)) THEN
       IF (interpol) THEN
          okay =.TRUE.
       END IF
    END IF
    !June 7 2007, JLG

    IF (rank==0) THEN
       WRITE(*,*) 'Reading Navier-Stokes modes ...'
    END IF
    DO i=1, m_max_c                  !pour tout les modes du processeur courant
       !ouverture du fichier
       OPEN(UNIT = 10, FILE = in_name, FORM = 'unformatted', STATUS = 'unknown')
       !on saute la premiere ligne du fichier qui contient des donnees
       READ(10)
       mode_cherche = list_mode(i)
       !recherche du bon mode
       trouve = .FALSE.
       DO j=1, nb_mode_r             !pour tout les modes ecris dans le suite.
          !lecture du mode
          READ(10) mode
          !June 7 2007, JLG
          IF (okay) THEN
             IF (j/=rang_F*m_max_c+i) THEN
                DO n=1, nlignes
                   READ(10)
                ENDDO
                CYCLE
             ELSE
                list_mode(i) = mode
                mode_cherche = mode
             END IF
          END IF
          !June 7 2007, JLG
          IF (mode == mode_cherche) THEN   !on a trouve le bon mode
             READ(10) un(:,:,i)
             DO kp = 1, taylor_order_min-1
                READ(10) der_un(kp)%DRT(:,:,i)
             END DO
             DO kp = taylor_order_min, taylor_order-1 !===inputs%taylor_order<taylor_order
                READ(10) !===Read empty stuff
             END DO
             READ(10) pn(:,:,i)
             DO kp = 1, taylor_order_min-1
                READ(10) der_pn(kp)%DRT(:,:,i)
             END DO
             DO kp = taylor_order_min, taylor_order-1!===inputs%taylor_order<taylor_order
                READ(10) !===Read empty stuff
             END DO
             DO kp = taylor_order_min , inputs%taylor_order-1 !===if taylor_order<inputs%taylor_order
                der_un(kp)%DRT(:,:,i) = 0.d0
                der_pn(kp)%DRT(:,:,i) = 0.d0
             END DO
             IF (PRESENT(opt_level_set) .AND. PRESENT(opt_level_set_m1)) THEN
                READ(10) opt_level_set(:,:,:,i)
                READ(10) opt_level_set_m1(:,:,:,i)
                READ(10) max_vel_loc
             END IF
             WRITE(*,'(A,i4,A)') 'mode ns ', mode_cherche,' found '
             trouve = .TRUE.
             EXIT                        !car on a trouve le bon mode
          ELSE                             !on passe au mode suivant en sautant 6 lignes
             DO n=1, nlignes
                READ(10)
             ENDDO
          ENDIF
       ENDDO

       IF (.NOT.trouve) THEN               !mode_cherche non trouve
          DO kp = 1, inputs%taylor_order-1
             der_un(kp)%DRT(:,:,i) = 0.d0
             der_pn(kp)%DRT(:,:,i) = 0.d0
          END DO
          IF (PRESENT(val_init)) THEN ! not implemented yet
             un(:,:,i)    = val_init
             pn(:,:,i)    = val_init
             IF (PRESENT(opt_level_set) .AND. PRESENT(opt_level_set_m1)) THEN
                opt_level_set(:,:,:,i)    = val_init
                opt_level_set_m1(:,:,:,i) = val_init
                max_vel_loc = val_init
             END IF
             WRITE(*,'(A,i4,A)') 'mode ns', mode_cherche,' not found'
          ELSE
             un(:,:,i)    = 0.d0
             pn(:,:,i)    = 0.d0
             IF (PRESENT(opt_level_set) .AND. PRESENT(opt_level_set_m1)) THEN
                opt_level_set(:,:,:,i)=0.d0
                opt_level_set_m1(:,:,:,i)=0.d0
             END IF
             WRITE(*,*) 'mode ns', mode_cherche, ' not found'
          ENDIF
       ENDIF
       CLOSE(10)                          !fermeture du fichier suite
    ENDDO

    IF (PRESENT(opt_max_vel)) THEN
       CALL MPI_ALLREDUCE(max_vel_loc, opt_max_vel, 1, MPI_DOUBLE_PRECISION, &
            MPI_MAX, communicator(2), code)
    END IF

  END SUBROUTINE read_restart_ns_taylor

  SUBROUTINE write_restart_maxwell(communicator, H_mesh, phi_mesh, time, list_mode, Hn, Hn1, Bn, Bn1, phin, phin1, &
       filename, it, freq_restart, opt_mono, opt_dt)

    USE def_type_mesh
    USE chaine_caractere
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                    :: H_mesh,phi_mesh
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: communicator
    REAL(KIND=8),                                   INTENT(IN) :: time
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(IN) :: Hn, Hn1, Bn, Bn1
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(IN) :: phin, phin1
    CHARACTER(len=200),                              INTENT(IN) :: filename
    INTEGER,                                        INTENT(IN) :: it, freq_restart
    LOGICAL, OPTIONAL,                              INTENT(IN) :: opt_mono
    REAL(KIND=8),                     OPTIONAL,     INTENT(IN) :: opt_dt

    INTEGER                           :: rang_S, rang_F, code, nb_procs_S, nb_procs_F, n, i
    INTEGER                           :: l, lblank
    CHARACTER(len=3)                  :: tit, tit_S
    CHARACTER(len=250)                :: out_name
    LOGICAL                           :: mono=.FALSE.
    LOGICAL                           :: skip

    CALL MPI_COMM_SIZE(communicator(1),nb_procs_S,code)
    CALL MPI_COMM_SIZE(communicator(2),nb_procs_F,code)
    CALL MPI_COMM_RANK(communicator(1),rang_S,code)
    CALL MPI_COMM_RANK(communicator(2),rang_F,code)

    WRITE(tit,'(i3)') it/freq_restart
    lblank = eval_blank(3,tit)
    DO l = 1, lblank - 1
       tit(l:l) = '0'
    END DO
    WRITE(tit_S,'(i3)') rang_S
    lblank = eval_blank(3,tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO

    IF (PRESENT(opt_mono)) THEN
       mono = opt_mono
    END IF

    IF (mono) THEN
       out_name = 'suite_maxwell_I'//tit//'.'//filename
    ELSE
       out_name = 'suite_maxwell_S'//tit_S//'_I'//tit//'.'//filename
    END IF
    skip = (mono .AND. rang_S /= 0)

    DO n = 1, nb_procs_F
       IF ( (rang_F == n-1) .AND. (.NOT. skip) ) THEN
          IF (rang_F == 0) THEN
             OPEN(UNIT = 10, FILE = out_name, POSITION='append', &
                  FORM = 'unformatted', STATUS = 'replace')
             IF (PRESENT(opt_dt)) THEN
                IF (mono) THEN
                   WRITE(10) time, H_mesh%np , phi_mesh%np , nb_procs_F, SIZE(list_mode), opt_dt
                ELSE
                   WRITE(10) time, nb_procs_S, nb_procs_F, SIZE(list_mode), opt_dt
                END IF
             ELSE
                IF (mono) THEN
                   WRITE(10) time, H_mesh%np , phi_mesh%np , nb_procs_F, SIZE(list_mode)
                ELSE
                   WRITE(10) time, nb_procs_S, nb_procs_F, SIZE(list_mode)
                END IF
             END IF
          ELSE
             OPEN(UNIT = 10, FILE = out_name, POSITION='append', &
                  FORM = 'unformatted', STATUS = 'unknown')
          END IF
          DO i= 1, SIZE(list_mode)
             WRITE(10) list_mode(i)
             IF (H_mesh%me /=0) THEN
                WRITE(10) Hn(:,:,i)
                WRITE(10) Hn1(:,:,i)
                WRITE(10) Bn(:,:,i)
                WRITE(10) Bn1(:,:,i)
             ELSE
                WRITE(10) 1
                WRITE(10) 1
                WRITE(10) 1
                WRITE(10) 1
             END IF
             IF (phi_mesh%me /=0) THEN
                WRITE(10) phin(:,:,i)
                WRITE(10) phin1(:,:,i)
             ELSE
                WRITE(10) 1
                WRITE(10) 1
             END IF
          END DO
          CLOSE(10)
       END IF
       CALL MPI_BARRIER(communicator(2),code)
    END DO

  END SUBROUTINE write_restart_maxwell


  SUBROUTINE read_restart_maxwell(communicator, H_mesh, phi_mesh, time, list_mode, Hn, Hn1, Bn, Bn1, phin, phin1, &
       filename, val_init, interpol, opt_mono &
       , opt_it, opt_dt) !===HF may 2020

    USE def_type_mesh
    USE chaine_caractere
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                    :: H_mesh,phi_mesh
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: communicator
    REAL(KIND=8),                                   INTENT(OUT):: time
    INTEGER,      DIMENSION(:)                                 :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(OUT):: Hn, Hn1, Bn, Bn1
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(OUT):: phin, phin1
    CHARACTER(len=200),                              INTENT(IN):: filename
    REAL(KIND=8), OPTIONAL,                         INTENT(IN) :: val_init
    LOGICAL     , OPTIONAL,                         INTENT(IN) :: interpol
    LOGICAL     , OPTIONAL,                         INTENT(IN) :: opt_mono
    !===HF may 2020
    INTEGER     , OPTIONAL,                         INTENT(IN) :: opt_it
    !===HF may 2020
    REAL(KIND=8),                     OPTIONAL,     INTENT(IN) :: opt_dt
    REAL(KIND=8)                                               :: dt_read, dt_ratio

    INTEGER     :: code, n, i, mode, j, rang_S, rang_F, nb_procs_F, nb_procs_S, rank
    INTEGER     :: m_max_cr, nb_procs_r, nb_procs_Sr
    INTEGER     :: m_max_c, nb_mode_r, mode_cherche
    LOGICAL     :: trouve, okay
    INTEGER     :: nph, npp
    INTEGER           :: l, lblank
    CHARACTER(len=3)  :: tit_S
    !===HF may 2020
    CHARACTER(len=3)  :: tit
    !===HF may 2020
    CHARACTER(len=250):: in_name
    LOGICAL           :: mono=.FALSE.

    CALL MPI_COMM_RANK(communicator(2),rang_F,code)
    CALL MPI_COMM_SIZE(communicator(2),nb_procs_F,code)
    CALL MPI_COMM_RANK(communicator(1),rang_S,code)
    CALL MPI_COMM_SIZE(communicator(1),nb_procs_S,code)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,code)

    !=== HF may 2020
    IF (PRESENT(opt_it)) THEN
       WRITE(tit,'(i3)') opt_it
       lblank = eval_blank(3,tit)
       DO l = 1, lblank - 1
          tit(l:l) = '0'
       END DO
    END IF
    !=== HF may 2020

    WRITE(tit_S,'(i3)') rang_S
    lblank = eval_blank(3,tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO
    IF (PRESENT(opt_mono)) THEN
       mono = opt_mono
    END IF

    IF (mono) THEN
       !=== HF may 2020
       IF (PRESENT(opt_it)) THEN
          in_name = 'suite_maxwell_I'//tit//'.'//filename
       ELSE
          in_name = 'suite_maxwell.'//filename
       END IF
    ELSE
       IF (PRESENT(opt_it)) THEN
          in_name = 'suite_maxwell_S'//tit_S//'_I'//tit//'.'//filename
       ELSE
          in_name = 'suite_maxwell_S'//tit_S//'.'//filename
       END IF
       !=== HF may 2020
    END IF

    !=== HF may 2020
    IF (PRESENT(opt_it)) THEN
       WRITE(*,*) 'restart Maxwell for it', opt_it
    ELSE
       WRITE(*,*) 'restart Maxwell'
    END IF
    !=== HF may 2020

    OPEN(UNIT = 10, FILE = in_name, FORM = 'unformatted', STATUS = 'unknown')
    IF (PRESENT(opt_dt)) THEN
       IF (mono) THEN
          READ(10) time, nph, npp, nb_procs_r, m_max_cr, dt_read
       ELSE
          READ(10) time, nb_procs_Sr, nb_procs_r, m_max_cr, dt_read
       END IF
    ELSE
       IF (mono) THEN
          READ(10) time, nph, npp, nb_procs_r, m_max_cr
       ELSE
          READ(10) time, nb_procs_Sr, nb_procs_r, m_max_cr
       END IF
    END IF

    IF ((nb_procs_Sr /= nb_procs_S) .AND. (.NOT. mono)) THEN
       CALL error_petsc('BUG in read_restart: nb_procs_Sr /= nb_procs_S')
       !STOP
    END IF

    CLOSE(10)

    IF (rang_F == 0) THEN
       WRITE(*,*) 'File name ', trim(adjustl(in_name))
       WRITE(*,*) 'Time = ', time
       WRITE(*,*) 'Number of processors from restart file = ',nb_procs_r
       WRITE(*,*) 'Number of modes per processor from restart file = ',m_max_cr
    ENDIF

    m_max_c   = SIZE(list_mode)      !nombre de modes par proc pour le calcul
    nb_mode_r = nb_procs_r*m_max_cr  !nombre total de modes contenus dans le suite

    !June 7 2007, JLG
    IF (nb_procs_F*m_max_c /= nb_mode_r) THEN
       WRITE(*,*) 'Warning in read_restart_maxwell: nb_procs_F*m_max_c /= nb_mode_r'
       !STOP
    END IF

    okay = .FALSE.
    IF (PRESENT(interpol)) THEN
       IF (interpol) THEN
          okay =.TRUE.
       END IF
    END IF
    !June 7 2007, JLG

    IF (rank==0) THEN
       WRITE(*,*) 'Reading Maxwell modes ...'
    END IF
    DO i=1, m_max_c                  !pour tout les modes du processeur courant
       !ouverture du fichier
       OPEN(UNIT = 10, FILE = in_name, FORM = 'unformatted', STATUS = 'unknown')
       !on saute la premier ligne du fichier qui contient des donnes
       READ(10)
       mode_cherche = list_mode(i)
       !recherche du bon mode
       trouve = .FALSE.
       DO j=1, nb_mode_r             !pour tout les modes ecris dans le suite.
          !lecture du mode
          READ(10) mode
          !June 7 2007, JLG
          IF (okay) THEN
             IF (j/=rang_F*m_max_c+i) THEN
                DO n=1, 6
                   READ(10)
                ENDDO
                CYCLE
             ELSE
                list_mode(i) = mode
                mode_cherche = mode
             END IF
          END IF
          !June 7 2007, JLG
          IF (mode == mode_cherche) THEN   !on a trouve le bon mode
             IF (H_mesh%me /=0) THEN
                READ(10) Hn(:,:,i)
                READ(10) Hn1(:,:,i)
                READ(10) Bn(:,:,i)
                READ(10) Bn1(:,:,i)
             ELSE
                READ(10)
                READ(10)
                READ(10)
                READ(10)
             END IF
             IF (phi_mesh%me /=0) THEN
                READ(10) phin(:,:,i)
                READ(10) phin1(:,:,i)
             ELSE
                READ(10)
                READ(10)
             END IF
             WRITE(*,*) 'mode maxwell',mode_cherche,' trouve '
             trouve = .TRUE.
             EXIT                        !car on a trouve le bon mode
          ELSE                             !on passe au mode suivant en sautant 4 lignes
             DO n=1, 6
                READ(10)
             ENDDO
          ENDIF
       ENDDO
       IF (.NOT.trouve) THEN               !mode_cherche non trouve
          IF (PRESENT(val_init)) THEN
             Hn(:,:,i)   = val_init ; Hn1(:,:,i)   = val_init
             Bn(:,:,i)   = val_init ; Bn1(:,:,i)   = val_init
             phin(:,:,i) = val_init ; phin1(:,:,i) = val_init
             WRITE(*,*) 'mode maxwell',mode_cherche,' non trouve'
          ELSE
             Hn(:,:,i)   = 0.d0 ; Hn1(:,:,i)   = 0.d0
             Bn(:,:,i)   = 0.d0 ; Bn1(:,:,i)   = 0.d0
             phin(:,:,i) = 0.d0 ; phin1(:,:,i) = 0.d0
             WRITE(*,*) 'mode maxwell',mode_cherche,' non trouve'
          ENDIF
       ENDIF
       CLOSE(10)                          !fermeture du fichier suite
    ENDDO
    !
    IF (PRESENT(opt_dt)) THEN
       IF (ABS((opt_dt - dt_read)/opt_dt).GT.1d-4) THEN
          dt_ratio = opt_dt/dt_read
          IF (rank==0) THEN
             WRITE(*,*) 'In Maxwell restart, suite_time_step different from inputs%dt ...'
             WRITE(*,*) ' opt_dt, dt_read =', opt_dt, dt_read
          END IF
          Hn1 = dt_ratio * Hn1 +(1.d0 - dt_ratio)* Hn
          Bn1 = dt_ratio * Bn1 +(1.d0 - dt_ratio)* Bn
          phin1 = dt_ratio * phin1 +(1.d0 - dt_ratio)* phin
       END IF
    END IF
  END SUBROUTINE read_restart_maxwell


  SUBROUTINE write_restart_temp(communicator, temp_mesh, time, list_mode, &
       tempn, tempn_m1, filename, it, freq_restart, opt_mono)
    USE def_type_mesh
    USE chaine_caractere
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                    :: temp_mesh
    REAL(KIND=8),                                   INTENT(IN) :: time
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: communicator
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(IN) :: tempn, tempn_m1
    LOGICAL, OPTIONAL,                              INTENT(IN) :: opt_mono
    CHARACTER(len=200),                             INTENT(IN) :: filename
    INTEGER,                                        INTENT(IN) :: it, freq_restart
    INTEGER                           :: code, n, i, rang_S, rang_F, nb_procs_S, nb_procs_F
    INTEGER                           :: l, lblank
    CHARACTER(len=3)                  :: tit, tit_S
    LOGICAL                           :: mono=.FALSE.
    LOGICAL                           :: skip
    CHARACTER(len=250)                :: out_name

    CALL MPI_COMM_SIZE(communicator(1),nb_procs_S,code)
    CALL MPI_COMM_SIZE(communicator(2),nb_procs_F,code)
    CALL MPI_COMM_RANK(communicator(1),rang_S,code)
    CALL MPI_COMM_RANK(communicator(2),rang_F,code)

    WRITE(tit,'(i3)') it/freq_restart
    lblank = eval_blank(3,tit)
    DO l = 1, lblank - 1
       tit(l:l) = '0'
    END DO
    WRITE(tit_S,'(i3)') rang_S
    lblank = eval_blank(3,tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO

    IF (PRESENT(opt_mono)) THEN
       mono = opt_mono
    END IF

    IF (mono) THEN
       out_name = 'suite_temp_I'//tit//'.'//filename
    ELSE
       out_name = 'suite_temp_S'//tit_S//'_I'//tit//'.'//filename
    END IF

    skip = (mono .AND. rang_S /= 0)

    DO n = 1, nb_procs_F
       IF ( (rang_F == n-1) .AND. (.NOT. skip) ) THEN
          IF (rang_F == 0) THEN
             OPEN(UNIT = 10, FILE = out_name, POSITION='append', &
                  FORM = 'unformatted', STATUS = 'replace')
             IF (mono) THEN
                WRITE(10) time, temp_mesh%np , nb_procs_F, SIZE(list_mode)
             ELSE
                WRITE(10) time, nb_procs_S, nb_procs_F, SIZE(list_mode)
             END IF
          ELSE
             OPEN(UNIT = 10, FILE = out_name, POSITION='append', &
                  FORM = 'unformatted', STATUS = 'unknown')
          END IF

          DO i= 1, SIZE(list_mode)
             WRITE(10) list_mode(i)
             WRITE(10) tempn(:,:,i)
             WRITE(10) tempn_m1(:,:,i)
          END DO
          CLOSE(10)
       END IF
       CALL MPI_BARRIER(communicator(2),code)
    END DO

  END SUBROUTINE write_restart_temp

  SUBROUTINE read_restart_temp(communicator, time, list_mode, &
       tempn, tempn_m1, filename, val_init, interpol, opt_mono)
    USE def_type_mesh
    USE chaine_caractere
    USE my_util
    IMPLICIT NONE
    REAL(KIND=8),                                   INTENT(OUT):: time
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: communicator
    INTEGER,      DIMENSION(:)                                 :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(OUT):: tempn, tempn_m1
    CHARACTER(len=200),                             INTENT(IN) :: filename
    REAL(KIND=8), OPTIONAL,                         INTENT(IN) :: val_init
    LOGICAL     , OPTIONAL,                         INTENT(IN) :: interpol
    LOGICAL     , OPTIONAL,                         INTENT(IN) :: opt_mono
    INTEGER     :: code, n, i, mode, j, rang_S, nb_procs_S, rang_F, nb_procs_F, nlignes, rank
    INTEGER     :: m_max_cr, nb_procs_r, nb_procs_Sr
    INTEGER     :: m_max_c, nb_mode_r, mode_cherche
    LOGICAL     :: trouve, okay
    INTEGER     :: np
    INTEGER           :: l, lblank
    CHARACTER(len=3)  :: tit_S
    LOGICAL           :: mono=.FALSE.
    CHARACTER(len=250):: in_name
    CALL MPI_COMM_RANK(communicator(2),rang_F,code)
    CALL MPI_COMM_SIZE(communicator(2),nb_procs_F,code)
    CALL MPI_COMM_RANK(communicator(1),rang_S,code)
    CALL MPI_COMM_SIZE(communicator(1),nb_procs_S,code)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,code)

    nlignes = 2

    WRITE(tit_S,'(i3)') rang_S
    lblank = eval_blank(3,tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO

    IF (PRESENT(opt_mono)) THEN
       mono = opt_mono
    END IF

    IF (mono) THEN
       in_name = 'suite_temp.'//filename
    ELSE
       in_name = 'suite_temp_S'//tit_S//'.'//filename
    END IF

    WRITE(*,*) 'restart temperature'
    OPEN(UNIT = 10, FILE = in_name, FORM = 'unformatted', STATUS = 'unknown')

    IF (mono) THEN
       READ(10) time, np, nb_procs_r, m_max_cr
       nb_procs_Sr = -1
    ELSE
       READ(10) time, nb_procs_Sr, nb_procs_r, m_max_cr
    END IF
    CLOSE(10)

    IF ((nb_procs_Sr /= nb_procs_S) .AND. (.NOT. mono)) THEN
       CALL error_petsc('BUG in read_restart: nb_procs_Sr /= nb_procs_S')
       !STOP
    END IF

    IF (rang_F == 0) THEN
       WRITE(*,*) 'File name ', trim(adjustl(in_name))
       WRITE(*,*) 'Time = ', time
       WRITE(*,*) 'Number of processors from restart file = ',nb_procs_r
       WRITE(*,*) 'Number of modes per processor from restart file = ',m_max_cr
    ENDIF

    m_max_c   = SIZE(list_mode)      !nombre de modes par proc pour le calcul
    nb_mode_r = nb_procs_r*m_max_cr  !nombre total de modes contenus dans le suite

    !June 7 2007, JLG
    IF (nb_procs_F*m_max_c /= nb_mode_r) THEN
       !CALL error_petsc('Bug in read_restart_ns: nb_procs_F*m_max_c /= nb_mode_r')
       WRITE(*,*) 'Warning in read_restart_temp: nb_procs_F*m_max_c /= nb_mode_r'
       !STOP
    END IF

    okay = .FALSE.
    IF (PRESENT(interpol)) THEN
       IF (interpol) THEN
          okay =.TRUE.
       END IF
    END IF
    !June 7 2007, JLG

    IF (rank==0) THEN
       WRITE(*,*) 'Reading temperature modes ...'
    END IF
    DO i=1, m_max_c                  !pour tout les modes du processeur courant
       !ouverture du fichier
       OPEN(UNIT = 10, FILE = in_name, FORM = 'unformatted', STATUS = 'unknown')
       !on saute la premiere ligne du fichier qui contient des donnees
       READ(10)
       mode_cherche = list_mode(i)
       !recherche du bon mode
       trouve = .FALSE.
       DO j=1, nb_mode_r             !pour tout les modes ecris dans le suite.
          !lecture du mode
          READ(10) mode
          !June 7 2007, JLG
          IF (okay) THEN
             IF (j/=rang_F*m_max_c+i) THEN
                DO n=1, nlignes
                   READ(10)
                ENDDO
                CYCLE
             ELSE
                list_mode(i) = mode
                mode_cherche = mode
             END IF
          END IF
          !June 7 2007, JLG
          IF (mode == mode_cherche) THEN   !on a trouve le bon mode
             READ(10) tempn(:,:,i)
             READ(10) tempn_m1(:,:,i)
             WRITE(*,'(A,i4,A)') 'mode temp ', mode_cherche,' found '
             trouve = .TRUE.
             EXIT                        !car on a trouve le bon mode
          ELSE                             !on passe au mode suivant en sautant 6 lignes
             DO n=1, nlignes
                READ(10)
             ENDDO
          ENDIF
       ENDDO

       IF (.NOT.trouve) THEN               !mode_cherche non trouve
          IF (PRESENT(val_init)) THEN ! not implemented yet
             tempn(:,:,i)    = val_init  ; tempn_m1(:,:,i)    = val_init
             WRITE(*,'(A,i4,A)') 'mode temp', mode_cherche,' not found'
          ELSE
             tempn(:,:,i)    = 0.d0 ; tempn_m1(:,:,i)    = 0.d0
             WRITE(*,*) 'mode temp', mode_cherche, ' not found'
          ENDIF
       ENDIF
       CLOSE(10)                          !fermeture du fichier suite
    ENDDO

  END SUBROUTINE read_restart_temp

  SUBROUTINE write_restart_conc(communicator, conc_mesh, time, list_mode, &
       concn, concn_m1, filename, it, freq_restart, opt_mono)
    USE def_type_mesh
    USE chaine_caractere
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                    :: conc_mesh
    REAL(KIND=8),                                   INTENT(IN) :: time
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: communicator
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(IN) :: concn, concn_m1
    LOGICAL, OPTIONAL,                              INTENT(IN) :: opt_mono
    CHARACTER(len=200),                             INTENT(IN) :: filename
    INTEGER,                                        INTENT(IN) :: it, freq_restart
    INTEGER                           :: code, n, i, rang_S, rang_F, nb_procs_S, nb_procs_F
    INTEGER                           :: l, lblank
    CHARACTER(len=3)                  :: tit, tit_S
    LOGICAL                           :: mono=.FALSE.
    LOGICAL                           :: skip
    CHARACTER(len=250)                :: out_name

    CALL MPI_COMM_SIZE(communicator(1),nb_procs_S,code)
    CALL MPI_COMM_SIZE(communicator(2),nb_procs_F,code)
    CALL MPI_COMM_RANK(communicator(1),rang_S,code)
    CALL MPI_COMM_RANK(communicator(2),rang_F,code)

    WRITE(tit,'(i3)') it/freq_restart
    lblank = eval_blank(3,tit)
    DO l = 1, lblank - 1
       tit(l:l) = '0'
    END DO
    WRITE(tit_S,'(i3)') rang_S
    lblank = eval_blank(3,tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO

    IF (PRESENT(opt_mono)) THEN
       mono = opt_mono
    END IF

    IF (mono) THEN
       out_name = 'suite_conc_I'//tit//'.'//filename
    ELSE
       out_name = 'suite_conc_S'//tit_S//'_I'//tit//'.'//filename
    END IF

    skip = (mono .AND. rang_S /= 0)

    DO n = 1, nb_procs_F
       IF ( (rang_F == n-1) .AND. (.NOT. skip) ) THEN
          IF (rang_F == 0) THEN
             OPEN(UNIT = 10, FILE = out_name, POSITION='append', &
                  FORM = 'unformatted', STATUS = 'replace')
             IF (mono) THEN
                WRITE(10) time, conc_mesh%np , nb_procs_F, SIZE(list_mode)
             ELSE
                WRITE(10) time, nb_procs_S, nb_procs_F, SIZE(list_mode)
             END IF
          ELSE
             OPEN(UNIT = 10, FILE = out_name, POSITION='append', &
                  FORM = 'unformatted', STATUS = 'unknown')
          END IF

          DO i= 1, SIZE(list_mode)
             WRITE(10) list_mode(i)
             WRITE(10) concn(:,:,i)
             WRITE(10) concn_m1(:,:,i)
          END DO
          CLOSE(10)
       END IF
       CALL MPI_BARRIER(communicator(2),code)
    END DO

  END SUBROUTINE write_restart_conc

  SUBROUTINE read_restart_conc(communicator, time, list_mode, &
       concn, concn_m1, filename, val_init, interpol, opt_mono)
    USE def_type_mesh
    USE chaine_caractere
    USE my_util
    IMPLICIT NONE
    REAL(KIND=8),                                   INTENT(OUT):: time
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: communicator
    INTEGER,      DIMENSION(:)                                 :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(OUT):: concn, concn_m1
    CHARACTER(len=200),                             INTENT(IN) :: filename
    REAL(KIND=8), OPTIONAL,                         INTENT(IN) :: val_init
    LOGICAL     , OPTIONAL,                         INTENT(IN) :: interpol
    LOGICAL     , OPTIONAL,                         INTENT(IN) :: opt_mono
    INTEGER     :: code, n, i, mode, j, rang_S, nb_procs_S, rang_F, nb_procs_F, nlignes, rank
    INTEGER     :: m_max_cr, nb_procs_r, nb_procs_Sr
    INTEGER     :: m_max_c, nb_mode_r, mode_cherche
    LOGICAL     :: trouve, okay
    INTEGER     :: np
    INTEGER           :: l, lblank
    CHARACTER(len=3)  :: tit_S
    LOGICAL           :: mono=.FALSE.
    CHARACTER(len=250):: in_name
    CALL MPI_COMM_RANK(communicator(2),rang_F,code)
    CALL MPI_COMM_SIZE(communicator(2),nb_procs_F,code)
    CALL MPI_COMM_RANK(communicator(1),rang_S,code)
    CALL MPI_COMM_SIZE(communicator(1),nb_procs_S,code)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,code)

    nlignes = 2

    WRITE(tit_S,'(i3)') rang_S
    lblank = eval_blank(3,tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO

    IF (PRESENT(opt_mono)) THEN
       mono = opt_mono
    END IF

    IF (mono) THEN
       in_name = 'suite_conc.'//filename
    ELSE
       in_name = 'suite_conc_S'//tit_S//'.'//filename
    END IF

    WRITE(*,*) 'restart concentration'
    OPEN(UNIT = 10, FILE = in_name, FORM = 'unformatted', STATUS = 'unknown')

    IF (mono) THEN
       READ(10) time, np, nb_procs_r, m_max_cr
       nb_procs_Sr = -1
    ELSE
       READ(10) time, nb_procs_Sr, nb_procs_r, m_max_cr
    END IF
    CLOSE(10)

    IF ((nb_procs_Sr /= nb_procs_S) .AND. (.NOT. mono)) THEN
       CALL error_petsc('BUG in read_restart: nb_procs_Sr /= nb_procs_S')
       !STOP
    END IF

    IF (rang_F == 0) THEN
       WRITE(*,*) 'File name ', trim(adjustl(in_name))
       WRITE(*,*) 'Time = ', time
       WRITE(*,*) 'Number of processors from restart file = ',nb_procs_r
       WRITE(*,*) 'Number of modes per processor from restart file = ',m_max_cr
    ENDIF

    m_max_c   = SIZE(list_mode)      !nombre de modes par proc pour le calcul
    nb_mode_r = nb_procs_r*m_max_cr  !nombre total de modes contenus dans le suite

    !June 7 2007, JLG
    IF (nb_procs_F*m_max_c /= nb_mode_r) THEN
       !CALL error_petsc('Bug in read_restart_ns: nb_procs_F*m_max_c /= nb_mode_r')
       WRITE(*,*) 'Warning in read_restart_conc: nb_procs_F*m_max_c /= nb_mode_r'
       !STOP
    END IF

    okay = .FALSE.
    IF (PRESENT(interpol)) THEN
       IF (interpol) THEN
          okay =.TRUE.
       END IF
    END IF
    !June 7 2007, JLG

    IF (rank==0) THEN
       WRITE(*,*) 'Reading concentration modes ...'
    END IF
    DO i=1, m_max_c                  !pour tout les modes du processeur courant
       !ouverture du fichier
       OPEN(UNIT = 10, FILE = in_name, FORM = 'unformatted', STATUS = 'unknown')
       !on saute la premiere ligne du fichier qui contient des donnees
       READ(10)
       mode_cherche = list_mode(i)
       !recherche du bon mode
       trouve = .FALSE.
       DO j=1, nb_mode_r             !pour tout les modes ecris dans le suite.
          !lecture du mode
          READ(10) mode
          !June 7 2007, JLG
          IF (okay) THEN
             IF (j/=rang_F*m_max_c+i) THEN
                DO n=1, nlignes
                   READ(10)
                ENDDO
                CYCLE
             ELSE
                list_mode(i) = mode
                mode_cherche = mode
             END IF
          END IF
          !June 7 2007, JLG
          IF (mode == mode_cherche) THEN   !on a trouve le bon mode
             READ(10) concn(:,:,i)
             READ(10) concn_m1(:,:,i)
             WRITE(*,'(A,i4,A)') 'mode conc ', mode_cherche,' found '
             trouve = .TRUE.
             EXIT                        !car on a trouve le bon mode
          ELSE                             !on passe au mode suivant en sautant 6 lignes
             DO n=1, nlignes
                READ(10)
             ENDDO
          ENDIF
       ENDDO

       IF (.NOT.trouve) THEN               !mode_cherche non trouve
          IF (PRESENT(val_init)) THEN ! not implemented yet
             concn(:,:,i)    = val_init  ; concn_m1(:,:,i)    = val_init
             WRITE(*,'(A,i4,A)') 'mode conc', mode_cherche,' not found'
          ELSE
             concn(:,:,i)    = 0.d0 ; concn_m1(:,:,i)    = 0.d0
             WRITE(*,*) 'mode conc', mode_cherche, ' not found'
          ENDIF
       ENDIF
       CLOSE(10)                          !fermeture du fichier suite
    ENDDO

  END SUBROUTINE read_restart_conc

END MODULE restart
