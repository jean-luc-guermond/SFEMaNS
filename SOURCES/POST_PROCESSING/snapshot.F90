MODULE snapshot
  USE petsc
#include "petsc/finclude/petsc.h"
  USE input_data
  IMPLICIT NONE
CONTAINS

  SUBROUTINE WRITE_PHYS_FIELD(communicator, mesh, V1_in, nb_procs, bloc_size, &
       type_field, opt_I, opt_dir, opt_nb_plane)
    USE my_util
    USE def_type_mesh
    USE chaine_caractere
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2

    CHARACTER(*), INTENT(IN) :: type_field
    CHARACTER(*), INTENT(IN) :: opt_I, opt_dir
    TYPE(mesh_type), INTENT(IN) :: mesh
    REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: V1_in
    INTEGER, INTENT(IN) :: nb_procs, bloc_size
    INTEGER, OPTIONAL :: opt_nb_plane
    INTEGER :: np, nb_field, &
         m_max, m_max_c, MPID, m_max_pad, N_r_pad, np_tot
    INTEGER(KIND = 8) :: fftw_plan_multi_c2r
    COMPLEX(KIND = 8), ALLOCATABLE, DIMENSION(:, :, :) :: cu
    INTEGER :: nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :, :) :: dist_field, combined_field
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :, :) :: V_out

    INTEGER :: fft_dim, howmany, istride, ostride, idist, odist, l, lblank
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    LOGICAL, SAVE :: once_fft = .TRUE.
    !LOGICAL,  DIMENSION(11),                 SAVE   ::  &
    !once=[.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE.]
    INTEGER :: rang_S, rang_F
    CHARACTER(len=4)                  :: tit_S
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !v5.0!#include "petsc/finclude/petsc.h"
    MPI_Comm, DIMENSION(:), POINTER :: communicator


    CALL MPI_COMM_RANK(communicator(1), rang_S, code)
    CALL MPI_COMM_RANK(communicator(2), rang_F, code)

    np = SIZE(V1_in, 1)
    nb_field = SIZE(V1_in, 2)    ! Number of fields
    m_max_c = SIZE(V1_in, 3)    ! Number of complex (cosines + sines) coefficients per point
    m_max = m_max_c * nb_procs ! Number of complex coefficients per point per processor
    np_tot = nb_procs * bloc_size

    IF (once_fft) THEN
       once_fft = .FALSE.
       !WRITE(*,*) 'np,nb_field,m_max_c,m_max,np_tot ', np,nb_field,m_max_c,m_max,np_tot
    END IF

    IF (MOD(nb_field, 2)/=0 .OR. m_max_c==0) THEN
       CALL error_petsc('Bug in WRITE_PHYS_FIELD: MOD(nb_field,2)/=0 .OR. m_max_c==0')
    END IF

    !===Bloc_size is the number of points that are handled by one processor
    !===once the Fourier modes are all collected on the processor

    IF (PRESENT(opt_nb_plane)) THEN
       IF (opt_nb_plane> 2 * m_max - 1) THEN
          m_max_pad = (opt_nb_plane + 1) / 2
       ELSE
          m_max_pad = m_max
       END IF
    ELSE
       m_max_pad = m_max
    END IF
    N_r_pad = 2 * m_max_pad - 1

    ALLOCATE(cu(m_max_pad, nb_field / 2, bloc_size))
    ALLOCATE(dist_field(m_max_c, nb_field, np_tot))
    ALLOCATE(combined_field(m_max_c, nb_field, np_tot))
    ALLOCATE(V_out(N_r_pad, nb_field / 2, bloc_size))

    DO i = 1, m_max_c
       dist_field(i, :, 1:np) = TRANSPOSE(V1_in(:, :, i))
    END DO

    IF (np/=np_tot) dist_field(:, :, np + 1:np_tot) = 1.d100

    longueur_tranche = bloc_size * m_max_c * nb_field

    MPID = MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field, longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator(2), code)

    cu = 0.d0
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb - 1) * bloc_size
          shiftl = (nb - 1) * m_max_c
          jindex = n + shiftc
          DO nf = 1, nb_field / 2
             !===Put real and imaginary parts in a complex
             !===nf=1,2,3 => V1_in
             !===INPUT ARE COSINE AND SINE COEFFICIENTS
             !===THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl + 1:shiftl + m_max_c, nf, n) = 0.5d0 * CMPLX(combined_field(:, 2 * nf - 1, jindex), &
                  -combined_field(:, 2 * nf, jindex), KIND = 8)
          END DO
       END DO
    END DO
    cu(1, :, :) = 2 * CMPLX(REAL(cu(1, :, :), KIND = 8), 0.d0, KIND = 8)
    !===Padding is done by initialization of cu: cu = 0
    !===This is equivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0

    !===Set the parameters for dfftw
    fft_dim = 1
    istride = 1
    ostride = 1
    idist = N_r_pad
    inembed(1) = N_r_pad
    DIM(1) = N_r_pad
    odist = m_max_pad
    onembed(1) = m_max_pad
    howmany = bloc_size * nb_field / 2

    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, V_out, inembed, istride, idist, FFTW_ESTIMATE)
    CALL dfftw_execute(fftw_plan_multi_c2r)

    WRITE(tit_S, '(i4)') rang_S
    lblank = eval_blank(4, tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO

    CALL WRITE_REAL_GAUSS_PDF(communicator, mesh, bloc_size, V_out, trim(type_field), &
         opt_I, TRIM(ADJUSTL(opt_dir)))

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)
    DEALLOCATE(cu, dist_field, combined_field)

  END SUBROUTINE WRITE_PHYS_FIELD

  SUBROUTINE FFT_INV_REAL_GAUSS_PDF(communicator, mesh, V_out, nb_procs, bloc_size, m_max_pad, &
             type_field, opt_I, opt_nb_plane)
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    USE my_util
    USE def_type_mesh
    USE chaine_caractere
    IMPLICIT NONE
    INCLUDE 'fftw3.f'

    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND = 8), DIMENSION(:, :, :), INTENT(OUT) :: V_out
    CHARACTER(*), INTENT(IN) :: type_field
    CHARACTER(*), INTENT(IN) :: opt_I
    TYPE(mesh_type), INTENT(IN) :: mesh

    INTEGER, OPTIONAL :: opt_nb_plane

    INTEGER, INTENT(IN) :: bloc_size, nb_procs
    INTEGER :: np, np_tot, nb_field, m_max, m_max_c, MPID, N_r_pad, m_max_pad
    INTEGER(KIND = 8) :: fftw_plan_multi_r2c
    COMPLEX(KIND = 8), DIMENSION(SIZE(V_out, 2) / 2, bloc_size) :: intermediate
    COMPLEX(KIND = 8), DIMENSION(SIZE(V_out, 2) / 2, bloc_size, SIZE(V_out, 3) * nb_procs) :: combined_prod_cu
    COMPLEX(KIND = 8), DIMENSION(SIZE(V_out, 2) / 2, bloc_size, SIZE(V_out, 3) * nb_procs) :: dist_prod_cu
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :, :) :: V_read
    COMPLEX(KIND = 8), ALLOCATABLE, DIMENSION(:, :, :) :: prod_cu
    INTEGER :: i_field
    INTEGER :: nb, shiftc, shiftl, longueur_tranche, i, n, code
    REAL(KIND = 8) :: t
    ! FFTW parameters
    INTEGER :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    CHARACTER(len = 1) :: tit_field

    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm, DIMENSION(:), POINTER :: communicator

    np = SIZE(V_out, 1)
    nb_field = SIZE(V_out, 2) ! Number of fields
    m_max_c = SIZE(V_out, 3) ! Number of complex (cosines + sines) coefficients per point
    m_max = m_max_c * nb_procs! Number of comlex coefficients per point per processor
    IF (PRESENT(opt_nb_plane)) THEN
       IF (opt_nb_plane> 2 * m_max - 1) THEN
          m_max_pad = (opt_nb_plane + 1) / 2
       ELSE
          m_max_pad = m_max
       END IF
    ELSE
       m_max_pad = m_max
    END IF
    N_r_pad = 2 * m_max_pad - 1
    np_tot = nb_procs * bloc_size
    ALLOCATE(V_read(N_r_pad, nb_field / 2, bloc_size))
    ALLOCATE(prod_cu(m_max_pad, SIZE(V_out, 2) / 2, bloc_size))
    IF (MOD(nb_field, 2)/=0 .OR. m_max_c==0) THEN
       WRITE(*, *) ' BUG '
       STOP
    END IF
    V_read = 0.d0
    IF (nb_field / 2 > 1) THEN
       DO n = 1, nb_field / 2
          WRITE(tit_field, '(i1)') n
          CALL READ_REAL_GAUSS_PDF(communicator, mesh, V_read(:, n, :), trim(type_field) // tit_field, opt_I)
       END DO
    ELSE
       CALL READ_REAL_GAUSS_PDF(communicator, mesh, V_read(:, 1, :), trim(type_field), opt_I)
    END IF

    fft_dim = 1; istride = 1; ostride = 1;
    idist = N_r_pad;   inembed(1) = N_r_pad; DIM(1) = N_r_pad
    odist = m_max_pad; onembed(1) = m_max_pad

    howmany = bloc_size * nb_field
    howmany = howmany / 2
    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, V_read, &
         inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE)
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    prod_cu = (1.d0 / N_r_pad) * prod_cu !Scaling
    combined_prod_cu(:, :, 1) = prod_cu(1, :, :)
    DO n = 2, m_max
       !combined_prod_cu(:,:,n)=prod_cu(n,:,:)
       combined_prod_cu(:, :, n) = 2 * CONJG(prod_cu(n, :, :))
    END DO
    t = MPI_WTIME()
    longueur_tranche = bloc_size * m_max_c * nb_field
    MPID = MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu, longueur_tranche, MPID, dist_prod_cu, longueur_tranche, &
         MPID, communicator(2), code)

    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb - 1) * bloc_size
          shiftl = (nb - 1) * m_max_c
          intermediate = dist_prod_cu(:, :, shiftl + i)
          DO n = 1, bloc_size
             IF (n + shiftc > np) CYCLE
             DO i_field = 1, nb_field / 2
                V_out(n + shiftc, i_field * 2 - 1, i) = REAL (intermediate(i_field, n), KIND = 8)
                V_out(n + shiftc, i_field * 2, i) = AIMAG(intermediate(i_field, n))
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE FFT_INV_REAL_GAUSS_PDF


  SUBROUTINE nodes_to_gauss(vv_mesh, field_in, field_out)
    USE my_util
    USE def_type_mesh
    USE chaine_caractere
    IMPLICIT NONE
    TYPE(mesh_type), INTENT(IN) :: vv_mesh
    REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: field_in
    REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%me, SIZE(field_in, 2), SIZE(field_in, 3)), INTENT(OUT) :: field_out
    REAL(KIND = 8), DIMENSION(vv_mesh%gauss%n_w, SIZE(field_in, 2)) :: Vs
    INTEGER, DIMENSION(vv_mesh%gauss%n_w) :: j_loc
    INTEGER :: m_max_c, i, index, m, k, l, nb_field
    m_max_c = SIZE(field_in, 3)
    nb_field = SIZE(field_in, 2)
    IF (MOD(nb_field, 2)/=0 .OR. m_max_c==0) THEN
       CALL error_petsc('nodes_to_gauss PROBLEM')
    END IF
    !=== Computation of dield_in on gauss points
    DO i = 1, m_max_c
       index = 0
       DO m = 1, vv_mesh%dom_me
          j_loc = vv_mesh%jj(:, m)

          DO k = 1, nb_field
             Vs(:, k) = field_in(j_loc, k, i)
          END DO

          DO l = 1, vv_mesh%gauss%l_G
             index = index + 1
             DO k = 1, nb_field
                field_out(index, k, i) = SUM(Vs(:, k) * vv_mesh%gauss%ww(:, l))
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE nodes_to_gauss

  SUBROUTINE WRITE_REAL_GAUSS_PDF(communicator, mesh, bloc_size, field, type_field, opt_I, opt_dir)
    USE my_util
    USE def_type_mesh
    USE chaine_caractere
    USE my_data_module
    REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: field
    INTEGER,                            INTENT(IN) :: bloc_size
    CHARACTER(*), INTENT(IN) :: type_field
    CHARACTER(200) :: filename
    CHARACTER(*), INTENT(IN) :: opt_I, opt_dir
    TYPE(mesh_type), INTENT(IN) :: mesh
    MPI_Comm, DIMENSION(:), POINTER :: communicator
    INTEGER :: code, nb_angles, nb_procs_S, nb_procs_F, lblank, k, l, rang_S, rang_F, n
    INTEGER :: ntot, d
    CHARACTER(len = 4) :: tit_S
    CHARACTER(len = 250) :: out_name

    CALL MPI_COMM_RANK(communicator(1), rang_S, code)
    CALL MPI_COMM_RANK(communicator(2), rang_F, code)
    CALL MPI_COMM_SIZE(communicator(1), nb_procs_S, code)
    CALL MPI_COMM_SIZE(communicator(2), nb_procs_F, code)

    WRITE(tit_S, '(i4)') rang_S
    lblank = eval_blank(4, tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO
    filename = inputs%file_name
    out_name = TRIM(ADJUSTL(opt_dir)) // '/phys_' // trim(type_field) // '_S' // tit_S // trim(opt_I) // '.' // filename

    IF (inputs%if_snapshot_gauss) THEN
       IF (rang_F==nb_procs_F - 1) THEN
          ntot = mesh%gauss%l_G * mesh%me
       ELSE
          ntot = rang_F * bloc_size + bloc_size
       END IF
    ELSE
       IF (rang_F==nb_procs_F - 1) THEN
          ntot = mesh%gauss%n_w
       ELSE
          ntot = rang_F * bloc_size + bloc_size
       END IF
    END IF

    ntot = ntot - rang_F * bloc_size
    nb_angles = SIZE(field,1)
    DO d = 1, SIZE(field,2)
       DO k = 1, nb_angles
          DO n = 1, nb_procs_F
             IF (rang_F == n - 1) THEN
                IF (rang_F == 0 .AND. k == 1) THEN
                   OPEN(UNIT = 10, FILE = out_name, POSITION = 'append', &
                        FORM = 'unformatted', ACCESS = 'stream', STATUS = 'replace')
                ELSE
                   OPEN(UNIT = 10, FILE = out_name, POSITION = 'append', &
                        FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown')
                END IF
                IF (rang_F == nb_procs_F - 1) THEN
                   WRITE(10) field(k, d, 1:ntot)
                ELSE
                   WRITE(10) field(k, d, :)
                END IF
                CLOSE(10)

             END IF
             CALL MPI_BARRIER(communicator(2), code)
          END DO
       END DO
    END DO

  END SUBROUTINE WRITE_REAL_GAUSS_PDF

  SUBROUTINE READ_REAL_GAUSS_PDF(communicator, mesh, field, type_field, opt_I)
    USE my_util
    USE def_type_mesh
    USE chaine_caractere
    REAL(KIND = 8), DIMENSION(:, :) :: field
    REAL(KIND = 8), DIMENSION(SIZE(field, 1), SIZE(field, 2)) :: field_trash
    CHARACTER(*), INTENT(IN) :: type_field
    CHARACTER(200) :: filename
    CHARACTER(*), INTENT(IN) :: opt_I
    TYPE(mesh_type), INTENT(IN) :: mesh
    MPI_Comm, DIMENSION(:), POINTER :: communicator
    INTEGER :: code, nb_angles, ntot_gauss, ntot_gauss_last, nb_procs_S, nb_procs_F, lblank
    INTEGER :: k, l, rang_S, rang_F, bloc_size_gauss, n, tot
    CHARACTER(len = 4) :: tit_S
    CHARACTER(len = 250) :: out_name

    CALL MPI_COMM_RANK(communicator(1), rang_S, code)
    CALL MPI_COMM_RANK(communicator(2), rang_F, code)
    CALL MPI_COMM_SIZE(communicator(1), nb_procs_S, code)
    CALL MPI_COMM_SIZE(communicator(2), nb_procs_F, code)

    bloc_size_gauss = mesh%gauss%l_G * mesh%me / nb_procs_F + 1

    WRITE(tit_S, '(i4)') rang_S
    lblank = eval_blank(4, tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO

    filename = inputs%file_name
    out_name = trim(inputs%folder_for_snapshot) // '/phys_' // trim(type_field) // '_S' // tit_S // &
            trim(opt_I) // '.' // filename

    ! added n_tot_last<ntot_gauss for last procF : RB MC 11/06/24
    ntot_gauss_last = mesh%gauss%l_G * mesh%me - (nb_procs_F-1)*bloc_size_gauss
    ntot_gauss = bloc_size_gauss
    nb_angles = SIZE(field, 1)
    OPEN(UNIT = 10, FILE = out_name, &
         FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown')
    tot=0
    DO k = 1, nb_angles
       DO n = 1, nb_procs_F
          ! IF (n == nb_procs_F) THEN
          !     tot= tot + ntot_gauss_last
          !     WRITE(*,*) "RANK S", rang_S, "LAST RANK F ", rang_F, 'ntot_gauss =',ntot_gauss_last
          ! ELSE
          !     tot= tot + ntot_gauss
          !     WRITE(*,*) "RANK S", rang_S, "RANK F ", rang_F, 'ntot_gauss =',ntot_gauss
          ! END IF
          IF (rang_F == n - 1) THEN
             IF (n == nb_procs_F) THEN
                ! IF (rang_F == nb_procs_F - 1) THEN deleted , RB MC 12/06/24
                READ(10) field(k, 1:ntot_gauss_last)
             ELSE
                READ(10) field(k, 1:ntot_gauss)
             END IF
          ELSE
             IF (n == nb_procs_F) THEN
                ! IF (rang_F == nb_procs_F - 1) THEN deleted , RB MC 12/06/24
                READ(10) field_trash(k, 1:ntot_gauss_last)
             ELSE
                READ(10) field_trash(k, 1:ntot_gauss)
             END IF
          END IF
       END DO
    END DO
    CLOSE(10)
    ! WRITE(*,*) 'tot =',tot
  END SUBROUTINE READ_REAL_GAUSS_PDF

  SUBROUTINE WRITE_PLAN_GAUSS(communicator, field, type_field, opt_I, opt_dir)
    USE my_util
    USE def_type_mesh
    USE chaine_caractere
    REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: field
    CHARACTER(*) :: type_field
    CHARACTER(*) :: opt_I, opt_dir
    MPI_Comm, DIMENSION(:), POINTER :: communicator
    INTEGER :: code, nb_procs_S, rang_S, lblank, l
    CHARACTER(len = 4) :: tit_S
    CHARACTER(len = 250) :: out_name
    CHARACTER(len = 200) :: filename

    !===MPI comms
    CALL MPI_COMM_SIZE(communicator(1), nb_procs_S, code)
    CALL MPI_COMM_RANK(communicator(1), rang_S, code)

    WRITE(tit_S, '(i4)') rang_S
    lblank = eval_blank(4, tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO

    CALL system('mkdir -p ' // TRIM(ADJUSTL(opt_dir)))
    filename = inputs%file_name
    out_name = TRIM(ADJUSTL(opt_dir)) // '/' // trim(type_field) // '_S' // tit_S // trim(opt_I) // &
            '.' // filename
    OPEN(UNIT = 10, FILE = out_name, POSITION = 'append', FORM = 'unformatted', ACCESS = 'stream', &
            STATUS = 'replace')

    WRITE(10) field(:)

  END SUBROUTINE WRITE_PLAN_GAUSS

  SUBROUTINE WRITE_MESH(communicator, mesh, opt_dir,mesh_name)
    USE my_util
    USE def_type_mesh
    USE chaine_caractere
    TYPE(mesh_type), INTENT(IN) :: mesh
    CHARACTER(*) :: opt_dir, mesh_name
    MPI_Comm, DIMENSION(:), POINTER :: communicator
    INTEGER :: code, nb_procs_S, rang_S, lblank, l
    CHARACTER(len = 4) :: tit_S
    CHARACTER(len = 250) :: out_name
    CHARACTER(len = 200) :: filename

    !===MPI comms
    CALL MPI_COMM_SIZE(communicator(1), nb_procs_S, code)
    CALL MPI_COMM_RANK(communicator(1), rang_S, code)

    WRITE(tit_S, '(i4)') rang_S
    lblank = eval_blank(4, tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO

    CALL system('mkdir -p ' // TRIM(ADJUSTL(opt_dir)))
    filename = inputs%file_name

    !start writing
    !WRITE mesh%dom_me, mesh%gauss%l_G ...  e, clair dans un fichier séparé
    out_name = TRIM(ADJUSTL(opt_dir)) // '/' // mesh_name // 'mesh_info' // '_S' // tit_S // '.txt'
    OPEN(UNIT = 667, FILE = out_name, STATUS="unknown")
    WRITE(667,*) "n_w : nombre de noeud par maille"
    WRITE(667,*) mesh%gauss%n_w
    WRITE(667,*) "l_G : nombre de point de gauss par maille"
    WRITE(667,*) mesh%gauss%l_G
    WRITE(667,*) "dom_me : nombre de maille dans le domaine S"
    WRITE(667,*) mesh%dom_me
    WRITE(667,*) "np : nombre de noeud dans le domaine S"
    WRITE(667,*) mesh%np

    ! added 15/06/25
    WRITE(667,*) "n_ws : nombre de noeud par maille apppartenant au bord du domaine global"
    WRITE(667,*) mesh%gauss%n_ws
    WRITE(667,*) "l_Gs : nombre de point de gauss par maille sur le bord du domaine global"
    WRITE(667,*) mesh%gauss%l_Gs
    WRITE(667,*) "dom_mes : nombre de maille appartenant au bord du domaine global"
    WRITE(667,*) mesh%dom_mes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    WRITE(667,*) "ww are the gauss point weights :: has shape (n_w,l_G) "
    WRITE(667,*) "wws are the gauss point weights at domain boundaries :: has shape (n_ws,l_Gs) "
    WRITE(667,*) "jj is the connectivity :: has shape (n_w,dom_me) "
    WRITE(667,*) "jjs is the connectivity of boundary nodes written in global index :: has shape (n_w,dom_me) "
    WRITE(667,*) "j_s are the indicies of boundary nodes :: has shape (n_ws,) "
    WRITE(667,*) "gauss_rj is the jacobian weights :: has shape (l_G,dom_me) "
    WRITE(667,*) "gauss_rjs is the jacobian weights at domain boundaries :: has shape (l_Gs,dom_mes) "
    WRITE(667,*) "gauss_dw is the derivative of the test functions :: has shape (dim=2,n_w,l_G,dom_me) "
    WRITE(667,*) "gauss_dw_s is the derivative of the test functions at domain boundaries :: has shape (dim=2,n_w,l_Gs,dom_mes) "

    CLOSE(667)

    out_name = TRIM(ADJUSTL(opt_dir)) // '/' // mesh_name // 'mesh_gauss_ww' // '_S' // tit_S // '.' &
            // filename
    OPEN(UNIT = 10, FILE = out_name, FORM = 'unformatted', ACCESS = 'stream',STATUS = 'unknown')
    WRITE(10) mesh%gauss%ww(:,:)
    CLOSE(10)

    out_name = TRIM(ADJUSTL(opt_dir)) // '/' // mesh_name // 'mesh_gauss_wws' // '_S' // tit_S // '.' &
            // filename
    OPEN(UNIT = 10, FILE = out_name, FORM = 'unformatted', ACCESS = 'stream',STATUS = 'unknown')
    WRITE(10) mesh%gauss%wws(:,:)
    CLOSE(10)

    out_name = TRIM(ADJUSTL(opt_dir)) // '/'// mesh_name // 'mesh_jj' // '_S' // tit_S // '.' // filename
    OPEN(UNIT = 10, FILE = out_name, FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown')
    WRITE(10) mesh%jj(:,:)
    CLOSE(10)

    out_name = TRIM(ADJUSTL(opt_dir)) // '/'// mesh_name // 'mesh_jjs' // '_S' // tit_S // '.' // filename
    OPEN(UNIT = 10, FILE = out_name, FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown')
    WRITE(10) mesh%jjs(:,:)
    CLOSE(10)

    out_name = TRIM(ADJUSTL(opt_dir)) // '/'// mesh_name // 'mesh_j_s' // '_S' // tit_S // '.' // filename
    OPEN(UNIT = 10, FILE = out_name, FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown')
    WRITE(10) mesh%j_s(:)
    CLOSE(10)

    out_name = TRIM(ADJUSTL(opt_dir)) // '/'// mesh_name // 'mesh_gauss_dw' // '_S' // tit_S // &
            '.' // filename
    OPEN(UNIT = 10, FILE = out_name, FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown')
    WRITE(10) mesh%gauss%dw(:,:,:,:) !=== shape (2,
    CLOSE(10)

    out_name = TRIM(ADJUSTL(opt_dir)) // '/'// mesh_name // 'mesh_gauss_dw_s' // '_S' // tit_S // &
            '.' // filename
    OPEN(UNIT = 10, FILE = out_name, FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown')
    WRITE(10) mesh%gauss%dw_s(:,:,:,:) !=== shape (2,
    CLOSE(10)

    out_name = TRIM(ADJUSTL(opt_dir)) // '/'// mesh_name // 'mesh_gauss_rj' // '_S' // tit_S // &
            '.' // filename
    OPEN(UNIT = 10, FILE = out_name, FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown')
    WRITE(10) mesh%gauss%rj(:,:) !=== shape (n_w,dom_me)
    CLOSE(10)

    out_name = TRIM(ADJUSTL(opt_dir)) // '/'// mesh_name // 'mesh_gauss_rjs' // '_S' // tit_S // &
            '.' // filename
    OPEN(UNIT = 10, FILE = out_name, FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown')
    WRITE(10) mesh%gauss%rjs(:,:) !=== shape (n_ws,dom_mes)
    CLOSE(10)

    out_name = TRIM(ADJUSTL(opt_dir)) // '/'// mesh_name // 'mesh_rr_node' // '_S' // tit_S // &
            '.' // filename
    OPEN(UNIT = 10, FILE = out_name, FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown')
    WRITE(10) mesh%rr(1,:)
    CLOSE(10)

    out_name = TRIM(ADJUSTL(opt_dir)) // '/'// mesh_name // 'mesh_zz_node' // '_S' // tit_S // &
            '.' // filename
    OPEN(UNIT = 10, FILE = out_name, FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown')
    WRITE(10) mesh%rr(2,:)
    CLOSE(10)

  END SUBROUTINE WRITE_MESH


  SUBROUTINE WRITE_FOURIER_FIELD(communicator, list_mode, field, name_field, opt_I, opt_dir)
    !!
    !! THIS WRITES IN PARALLEL IN FILE :: WORKS ONLY IF the total number of modes
    !!                                                   is a multiple of nb_procs_F
    USE def_type_mesh
    USE chaine_caractere
    REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: field
    INTEGER, DIMENSION(:), INTENT(IN) :: list_mode
    CHARACTER(*) :: name_field
    CHARACTER(*) :: opt_I, opt_dir
    INTEGER :: ntot, nb_field, nb_modes, rang_S, rang_F, nb_procs_F, nb_procs_S, code, lblank, l
    CHARACTER(len = 4) :: tit_S
    CHARACTER(len = 250) :: out_name
    CHARACTER(len = 200) :: filename
    MPI_Comm, DIMENSION(:), POINTER :: communicator
    INTEGER :: filehandle , arraytype, local_size, ierr
    INTEGER (KIND=MPI_OFFSET_KIND) :: cursor_offset

    ntot = SIZE(field, 1)
    nb_field = SIZE(field, 2)
    nb_modes = SIZE(field, 3)

    CALL MPI_COMM_RANK(communicator(1), rang_S, code)
    CALL MPI_COMM_RANK(communicator(2), rang_F, code)
    CALL MPI_COMM_SIZE(communicator(1), nb_procs_S, code)
    CALL MPI_COMM_SIZE(communicator(2), nb_procs_F, code)

    WRITE(tit_S, '(i4)') rang_S
    lblank = eval_blank(4, tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO

    filename = inputs%file_name

    CALL system('mkdir -p ' // TRIM(ADJUSTL(opt_dir)))
    out_name = TRIM(ADJUSTL(opt_dir)) // '/fourier_' // trim(adjustl(name_field)) // '_S' // &
            tit_S // trim(opt_I) // '.' // filename

    local_size = SIZE(field,1) * SIZE(field,2) * SIZE(list_mode)
    ! off set = taille d'un mode * nb_mode * rang_F
    cursor_offset = local_size * rang_F
    WRITE(*,*) local_size, cursor_offset

    CALL MPI_File_open(communicator(2), out_name, MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, &
           filehandle, ierr)
    CALL MPI_Type_contiguous(local_size,MPI_DOUBLE_PRECISION,arraytype,ierr)
    CALL MPI_Type_commit(arraytype,ierr)
    CALL MPI_File_set_view(filehandle, cursor_offset*8_MPI_OFFSET_KIND, MPI_DOUBLE_PRECISION, &
            arraytype, 'native', MPI_INFO_NULL, ierr)
    CALL MPI_File_write_all(filehandle, field(:,:,:), local_size, MPI_DOUBLE_PRECISION, &
            MPI_STATUS_IGNORE, ierr)
    CALL MPI_File_close(filehandle, ierr)

  END SUBROUTINE WRITE_FOURIER_FIELD

  SUBROUTINE WRITE_FOURIER_MODES(communicator, list_mode, field, name_field, opt_dir)
    USE chaine_caractere
    USE petsc
    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(IN) :: communicator
    CHARACTER(*), INTENT(IN) :: name_field, opt_dir
    REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: field
    INTEGER, DIMENSION(:), INTENT(IN) :: list_mode
    CHARACTER(len = 200) :: filename
    INTEGER :: code, i, rang_S, rang_F, nb_procs_S, nb_procs_F
    INTEGER :: l, lblank, nb_field, unit_to_write
    CHARACTER(len = 4) :: tit_S, tit_F
    CHARACTER(len = 250) :: out_name

    CALL MPI_COMM_SIZE(communicator(1), nb_procs_S, code)
    CALL MPI_COMM_SIZE(communicator(2), nb_procs_F, code)
    CALL MPI_COMM_RANK(communicator(1), rang_S, code)
    CALL MPI_COMM_RANK(communicator(2), rang_F, code)

    ! CN 24/01/2024
    CALL system('mkdir -p ' // TRIM(ADJUSTL(opt_dir)))

    WRITE(tit_S, '(i4)') rang_S
    lblank = eval_blank(4, tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO
    nb_field = SIZE(field, 2)



    !DO n = 1, nb_field
    !IF (nb_field > 2) THEN
    !   IF (MOD(n, 2) == 0) THEN
    !      WRITE(tit_n, '(i1)') n / 2
    !      type_field = trim(name_field) // tit_n // 's'
    !   ELSE
    !      WRITE(tit_n, '(i1)') (n + 1) / 2
    !      type_field = trim(name_field) // tit_n // 'c'
    !   END IF
    !ELSE
    !   IF (n == 1) THEN
    !      type_field = trim(name_field) // 'c'
    !   ELSE
    !      type_field = trim(name_field) // 's'
    !   END IF
    !END IF

    filename = inputs%file_name

    DO i = 1, SIZE(field, 3)  !fourier modes
       !===Creation of the suite file name
       WRITE(tit_F, '(i4)') list_mode(i)
       lblank = eval_blank(4, tit_F)
       DO l = 1, lblank - 1
          tit_F(l:l) = '0'
       END DO

       out_name = trim(adjustl(opt_dir)) // '/fourier_' // trim(adjustl(name_field)) // '_S' &
            // tit_S // '_F' // tit_F // '.' // filename

       OPEN(NEWUNIT = unit_to_write, FILE = out_name, action = 'write', position = 'append', &
               FORM = 'unformatted', ACCESS = 'stream')
       WRITE(unit_to_write) field(:, :, i)
       CLOSE(unit_to_write)
    END DO
    CALL MPI_BARRIER(communicator(2), code)
    !END DO

  END SUBROUTINE WRITE_FOURIER_MODES

  SUBROUTINE R_mrgrnk (XDONT, IRNGT)
    ! __________________________________________________________
    !   MRGRNK = Merge-sort ranking of an array
    !   For performance reasons, the first 2 passes are taken
    !   out of the standard loop, and use dedicated coding.
    ! __________________________________________________________
    ! _________________________________________________________
    Real(KIND = 8), Dimension (:), Intent (In) :: XDONT
    Integer, Dimension (:), Intent (Out) :: IRNGT
    ! __________________________________________________________
    Real(KIND = 8) :: XVALA, XVALB
    !
    Integer, Dimension (SIZE(IRNGT)) :: JWRKT
    Integer :: LMTNA, LMTNC, IRNG1, IRNG2
    Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
    !
    NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
    Select Case (NVAL)
    Case (:0)
       Return
    Case (1)
       IRNGT (1) = 1
       Return
    Case Default
       Continue
    End Select
    !
    !  Fill-in the index array, creating ordered couples
    !
    Do IIND = 2, NVAL, 2
       If (XDONT(IIND - 1) <= XDONT(IIND)) Then
          IRNGT (IIND - 1) = IIND - 1
          IRNGT (IIND) = IIND
       Else
          IRNGT (IIND - 1) = IIND
          IRNGT (IIND) = IIND - 1
       End If
    End Do
    If (Modulo(NVAL, 2) /= 0) Then
       IRNGT (NVAL) = NVAL
    End If
    !
    !  We will now have ordered subsets A - B - A - B - ...
    !  and merge A and B couples into     C   -   C   - ...
    !
    LMTNA = 2
    LMTNC = 4
    !
    !  First iteration. The length of the ordered subsets goes from 2 to 4
    !
    Do
       If (NVAL <= 2) Exit
       !
       !   Loop on merges of A and B into C
       !
       Do IWRKD = 0, NVAL - 1, 4
          If ((IWRKD + 4) > NVAL) Then
             If ((IWRKD + 2) >= NVAL) Exit
             !
             !   1 2 3
             !
             If (XDONT(IRNGT(IWRKD + 2)) <= XDONT(IRNGT(IWRKD + 3))) Exit
             !
             !   1 3 2
             !
             If (XDONT(IRNGT(IWRKD + 1)) <= XDONT(IRNGT(IWRKD + 3))) Then
                IRNG2 = IRNGT (IWRKD + 2)
                IRNGT (IWRKD + 2) = IRNGT (IWRKD + 3)
                IRNGT (IWRKD + 3) = IRNG2
                !
                !   3 1 2
                !
             Else
                IRNG1 = IRNGT (IWRKD + 1)
                IRNGT (IWRKD + 1) = IRNGT (IWRKD + 3)
                IRNGT (IWRKD + 3) = IRNGT (IWRKD + 2)
                IRNGT (IWRKD + 2) = IRNG1
             End If
             Exit
          End If
          !
          !   1 2 3 4
          !
          If (XDONT(IRNGT(IWRKD + 2)) <= XDONT(IRNGT(IWRKD + 3))) Cycle
          !
          !   1 3 x x
          !
          If (XDONT(IRNGT(IWRKD + 1)) <= XDONT(IRNGT(IWRKD + 3))) Then
             IRNG2 = IRNGT (IWRKD + 2)
             IRNGT (IWRKD + 2) = IRNGT (IWRKD + 3)
             If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD + 4))) Then
                !   1 3 2 4
                IRNGT (IWRKD + 3) = IRNG2
             Else
                !   1 3 4 2
                IRNGT (IWRKD + 3) = IRNGT (IWRKD + 4)
                IRNGT (IWRKD + 4) = IRNG2
             End If
             !
             !   3 x x x
             !
          Else
             IRNG1 = IRNGT (IWRKD + 1)
             IRNG2 = IRNGT (IWRKD + 2)
             IRNGT (IWRKD + 1) = IRNGT (IWRKD + 3)
             If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD + 4))) Then
                IRNGT (IWRKD + 2) = IRNG1
                If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD + 4))) Then
                   !   3 1 2 4
                   IRNGT (IWRKD + 3) = IRNG2
                Else
                   !   3 1 4 2
                   IRNGT (IWRKD + 3) = IRNGT (IWRKD + 4)
                   IRNGT (IWRKD + 4) = IRNG2
                End If
             Else
                !   3 4 1 2
                IRNGT (IWRKD + 2) = IRNGT (IWRKD + 4)
                IRNGT (IWRKD + 3) = IRNG1
                IRNGT (IWRKD + 4) = IRNG2
             End If
          End If
       End Do
       !
       !  The Cs become As and Bs
       !
       LMTNA = 4
       Exit
    End Do
    !
    !  Iteration loop. Each time, the length of the ordered subsets
    !  is doubled.
    !
    Do
       If (LMTNA >= NVAL) Exit
       IWRKF = 0
       LMTNC = 2 * LMTNC
       !
       !   Loop on merges of A and B into C
       !
       Do
          IWRK = IWRKF
          IWRKD = IWRKF + 1
          JINDA = IWRKF + LMTNA
          IWRKF = IWRKF + LMTNC
          If (IWRKF >= NVAL) Then
             If (JINDA >= NVAL) Exit
             IWRKF = NVAL
          End If
          IINDA = 1
          IINDB = JINDA + 1
          !
          !   Shortcut for the case when the max of A is smaller
          !   than the min of B. This line may be activated when the
          !   initial set is already close to sorted.
          !
          !          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
          !
          !  One steps in the C subset, that we build in the final rank array
          !
          !  Make a copy of the rank array for the merge iteration
          !
          JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
          !
          XVALA = XDONT (JWRKT(IINDA))
          XVALB = XDONT (IRNGT(IINDB))
          !
          Do
             IWRK = IWRK + 1
             !
             !  We still have unprocessed values in both A and B
             !
             If (XVALA > XVALB) Then
                IRNGT (IWRK) = IRNGT (IINDB)
                IINDB = IINDB + 1
                If (IINDB > IWRKF) Then
                   !  Only A still with unprocessed values
                   IRNGT (IWRK + 1:IWRKF) = JWRKT (IINDA:LMTNA)
                   Exit
                End If
                XVALB = XDONT (IRNGT(IINDB))
             Else
                IRNGT (IWRK) = JWRKT (IINDA)
                IINDA = IINDA + 1
                If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                XVALA = XDONT (JWRKT(IINDA))
             End If
             !
          End Do
       End Do
       !
       !  The Cs become As and Bs
       !
       LMTNA = 2 * LMTNA
    End Do
    !
    Return
    !
  END SUBROUTINE R_mrgrnk

  SUBROUTINE write_snapshot(communicator, mesh, list_mode, field, filename, it, freq_snapshot)
    USE def_type_mesh
#include "petsc/finclude/petsc.h"
    USE petsc
    USE chaine_caractere
    USE restart
    !INTEGER, DIMENSION(:),     INTENT(IN) :: communicator
    MPI_Comm, DIMENSION(:), POINTER :: communicator
    TYPE(mesh_type),           INTENT(IN) :: mesh
    INTEGER, DIMENSION(:),     INTENT(IN) :: list_mode
    REAL(KIND = 8), DIMENSION(:,:,:), INTENT(IN) :: field
    CHARACTER(*),              INTENT(IN) :: filename
    INTEGER,                   INTENT(IN) :: it, freq_snapshot

    REAL(KIND = 8), DIMENSION(mesh%gauss%l_G * mesh%me, SIZE(field,2), SIZE(list_mode)) :: field_gauss
    INTEGER :: nb_procs_F, bloc_size_gauss, bloc_size_node, rang_F, rang_S, nb_procs_S
    INTEGER :: rank, code, m_max_c, lblank, l
    CHARACTER(len = 4) :: tit_I
    CHARACTER(len = 250) :: opt_I, opt_dir
    !    INTEGER :: index, m, k, i
    !    CHARACTER(len = 3) :: tit
    !    INTEGER, DIMENSION(mesh%gauss%n_w) :: j_loc
    !    REAL(KIND = 8), DIMENSION(mesh%gauss%l_G * mesh%me) :: scalar_gauss
    !    REAL(KIND = 8), DIMENSION(2, mesh%gauss%l_G * mesh%dom_me) :: rr_gauss

    CALL MPI_COMM_RANK(PETSC_COMM_WORLD, rank, code)
    CALL MPI_COMM_SIZE(communicator(1), nb_procs_S, code)
    CALL MPI_COMM_RANK(communicator(1), rang_S, code)
    CALL MPI_COMM_SIZE(communicator(2), nb_procs_F, code)
    CALL MPI_COMM_RANK(communicator(2), rang_F, code)

    m_max_c = SIZE(list_mode)
    bloc_size_gauss = mesh%gauss%l_G * mesh%me / nb_procs_F + 1
    bloc_size_node = mesh%np / nb_procs_F + 1

    opt_dir = inputs%folder_for_snapshot


    WRITE(tit_I, '(i4)') it/freq_snapshot
    lblank = eval_blank(4, tit_I)
    DO l = 1, lblank - 1
       tit_I(l:l) = '0'
    END DO
    opt_I = '_I' // tit_I


    IF (inputs%if_snapshot_gauss) THEN
       CALL nodes_to_gauss(mesh, field, field_gauss)
       IF (inputs%if_snapshot_phys) THEN
          CALL WRITE_PHYS_FIELD(communicator, mesh, field_gauss, nb_procs_F, bloc_size_gauss, &
               filename, opt_I, opt_dir)
       END IF
       IF (inputs%if_snapshot_fourier) THEN
          CALL WRITE_FOURIER_FIELD(communicator, list_mode, field_gauss, &
               filename, opt_I, opt_dir)
       END IF
       IF (inputs%if_snapshot_fourier_per_mode) THEN
          CALL WRITE_FOURIER_MODES(communicator, list_mode, field_gauss, &
               filename, opt_dir)
       END IF
    ELSE
       IF (inputs%if_snapshot_phys) THEN
          CALL WRITE_PHYS_FIELD(communicator, mesh, field, nb_procs_F, bloc_size_node, &
               filename, opt_I, opt_dir)
       END IF
       IF (inputs%if_snapshot_fourier) THEN
          CALL WRITE_FOURIER_FIELD(communicator, list_mode, field, &
               filename, opt_I, opt_dir)
       END IF
       IF (inputs%if_snapshot_fourier_per_mode) THEN
          CALL WRITE_FOURIER_MODES(communicator, list_mode, field, &
               filename, opt_dir)
       END IF
    END IF

  END SUBROUTINE write_snapshot

END MODULE snapshot
