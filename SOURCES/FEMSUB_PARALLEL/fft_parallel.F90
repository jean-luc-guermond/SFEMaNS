!
!Authors:  Katarzyna Boronska, Jean-Luc Guermond, Copyrights 2007
!
MODULE sft_parallele
#include "petsc/finclude/petsc.h"
  USE petsc
  IMPLICIT NONE
  PUBLIC :: FFT_PAR_REAL, FFT_PAR_CROSS_PROD_DCL, &
       FFT_PAR_DOT_PROD_DCL, FFT_PAR_PROD_DCL, FFT_PAR_DIV_DCL, FFT_HEAVISIDE_DCL, &
       FFT_SCALAR_VECT_DCL, FFT_MAX_VEL_DCL, FFT_TENSOR_DCL,  &
       FFT_PAR_VAR_ETA_PROD_T_DCL, FFT_PAR_VAR_ETA_PROD_GAUSS_DCL, &
       FFT_CHECK_INTERFACE, FFT_COMPUTE_ENTROPY_VISC, FFT_COMPUTE_ENTROPY_VISC_MOM, &
       FFT_COMPUTE_DIFFU_MOM, FFT_PAR_COMPR_ENTRO_VISC_DCL, &
       FFT_COMPUTE_ENTROPY_VISC_GRAD_MOM, FFT_NO_OVERSHOOT_LEVEL_SET, &
       FFT_COMPRESSION_LEVEL_SET_DCL, regul, FFT_PAR_ENTRO_VISC_DCL, &
       FFT_PAR_SCAL_FUNCT, FFT_SCALAR_VECT_NO_OVERSHOOT, FFT_MAX_MIN_VEL_DCL, &
       FFT_SQUAREROOT_SCALAR
  PRIVATE
  !===
  !Comments about the use of select_mode.
  !The first mode must be 0 and all the other modes must be multiple of one non-zero integer
  !If select_mode is used, the FFT does not care and computes the real values
  !as if al the modes were contiguous. It does not matter since the values
  !in the real spaces are computed at the same fake angles, but they are sent back
  !in the fourier space on the correct fourier modes.
  !===
CONTAINS
  SUBROUTINE FFT_PAR_REAL(communicator, V1_in, V_out, opt_nb_plane)
    USE my_util
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)    :: V1_in
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE    :: V_out
    ! FL-CN possible faute ??? 19/03/2013
    !REAL(KIND=8), DIMENSION(:,:,:), POINTER        :: V_out
    INTEGER, OPTIONAL                              :: opt_nb_plane
    INTEGER                                        :: np, bloc_size, nb_field, &
         m_max, m_max_c, rang, nb_procs, MPID, m_max_pad, N_r_pad
    INTEGER(KIND=8)                                :: fftw_plan_multi_c2r
    COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: cu
    INTEGER :: nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)    :: dist_field, combined_field
    INTEGER               :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator

    CALL MPI_COMM_SIZE(communicator, nb_procs, code)
    CALL MPI_COMM_RANK(communicator, rang, code)

    np       = SIZE(V1_in,1)
    nb_field = SIZE(V1_in,2)    ! Number of fields
    m_max_c  = SIZE(V1_in,3)    ! Number of complex (cosines + sines) coefficients per point
    m_max    = m_max_c*nb_procs ! Number of comlex coefficients per point per processor
    IF (MOD(nb_field,2)/=0 .OR. m_max_c==0) THEN
       CALL error_petsc('Bug in FFT_PAR_REAL: MOD(nb_field,2)/=0 .OR. m_max_c==0')
    END IF

    !===Bloc_size is the number of points that are handled by one processor
    !===once the Fourier modes are all collected on the processor
    IF (MODULO(np,nb_procs)==0) THEN
       bloc_size = np/nb_procs
    ELSE
       bloc_size=1
       CALL error_petsc('Bug in FFT_PAR_REAL: np is not a multiple of nb_procs')
    END IF

    IF (PRESENT(opt_nb_plane)) THEN
       IF (opt_nb_plane> 2*m_max-1) THEN
          m_max_pad = (opt_nb_plane+1)/2
       ELSE
          m_max_pad = m_max
       END IF
    ELSE
       m_max_pad = m_max
    END IF
    N_r_pad=2*m_max_pad-1

    ALLOCATE(cu(m_max_pad,nb_field/2, bloc_size))
    ALLOCATE(dist_field(m_max_c,nb_field,np))
    ALLOCATE(combined_field(m_max_c,nb_field,np))

    DO i = 1, m_max_c
       dist_field(i,:,:) = TRANSPOSE(V1_in(:,:,i))
    END DO

    longueur_tranche=bloc_size*m_max_c*nb_field

    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)

    cu = 0.d0
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, nb_field/2
             !===Put real and imaginary parts in a complex
             !===nf=1,2,3 => V1_in
             !===INPUT ARE COSINE AND SINE COEFFICIENTS
             !===THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !===Padding is done by initialization of cu: cu = 0
    !===This is equivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0

    !===Set the parameters for dfftw
    fft_dim   = 1
    istride   = 1
    ostride   = 1
    idist     = N_r_pad
    inembed(1)= N_r_pad
    DIM(1)    = N_r_pad
    odist     = m_max_pad
    onembed(1)= m_max_pad
    howmany   = bloc_size*nb_field/2

    ! FL-CN possible faute ??? 19/03/2013
    !   IF (ASSOCIATED(V_out)) NULLIFY(V_out)
    !   IF (ASSOCIATED(V_out)) DEALLOCATE(V_out)
    ! pb sur la ligne suivante
    IF (ALLOCATED(V_out)) DEALLOCATE(V_out)
    ALLOCATE(V_out(N_r_pad,nb_field/2,bloc_size))

    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, V_out, inembed, istride, idist, FFTW_ESTIMATE)
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    DEALLOCATE(cu, dist_field, combined_field)

  END SUBROUTINE FFT_PAR_REAL

  SUBROUTINE FFT_PAR_CROSS_PROD_DCL(communicator,V1_in, V2_in, V_out, nb_procs, bloc_size, m_max_pad, temps)
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V1_in, V2_in
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    INTEGER, INTENT(IN)                                 :: bloc_size, m_max_pad, nb_procs
    INTEGER          :: np, np_tot, nb_field, m_max, m_max_c, MPID, N_r_pad
    INTEGER(KIND=8)  :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    COMPLEX(KIND=8), DIMENSION(m_max_pad, SIZE(V1_in,2), bloc_size)  :: cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,SIZE(V1_in,2),bloc_size)   :: ru
    COMPLEX(KIND=8), DIMENSION(m_max_pad,SIZE(V1_in,2)/2,bloc_size)  :: prod_cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,SIZE(V1_in,2)/2,bloc_size) :: prod_ru
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size)            :: intermediate
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),2*SIZE(V1_in,2),bloc_size*nb_procs)    :: dist_field, combined_field
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: dist_prod_cu
    INTEGER :: i_field
    INTEGER :: nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code
    REAL(KIND=8) :: t
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator

    IF (PRESENT(temps)) temps = 0.d0

    np      = SIZE(V1_in,1)
    nb_field= SIZE(V1_in,2) ! Number of fields
    m_max_c = SIZE(V1_in,3) ! Number of complex (cosines + sines) coefficients per point
    m_max = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    N_r_pad=2*m_max_pad-1
    np_tot = nb_procs*bloc_size

    IF (MOD(nb_field,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_PAR_CROSS_PROD_DCL '
       STOP
    END IF

    ! Bloc_size is the number of points that are handled by one processor
    ! once the Fourier modes are all collected
    ! Computation of bloc_size and np_tot
    ! fin de la repartition des points

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()

    DO i = 1, m_max_c
       dist_field(i,1:nb_field,1:np) = TRANSPOSE(V1_in(:,:,i))
       dist_field(i,nb_field+1:2*nb_field,1:np) = TRANSPOSE(V2_in(:,:,i))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100

    longueur_tranche=bloc_size*m_max_c*nb_field*2

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, nb_field
             ! Put real and imaginary parts in a complex
             ! nf=1,2,3 => V1_in
             ! nf=4,5,6 => V2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*nb_field

    t = MPI_WTIME()
    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_CROSS_PROD: fftw_plan_multi_c2r', fftw_plan_multi_c2r
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    ! CROSS PRODUCT
    IF (nb_field==6) THEN
       prod_ru(:,1,:) = ru(:,2,:)*ru(:,6,:) - ru(:,3,:)*ru(:,5,:)
       prod_ru(:,2,:) = ru(:,3,:)*ru(:,4,:) - ru(:,1,:)*ru(:,6,:)
       prod_ru(:,3,:) = ru(:,1,:)*ru(:,5,:) - ru(:,2,:)*ru(:,4,:)
    END IF
    ! CROSS PRODUCT

    howmany = howmany/2
    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
         inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_CROSS_PROD: fftw_plan_multi_r2c', fftw_plan_multi_r2c
    CALL dfftw_execute(fftw_plan_multi_r2c)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    !JLG, FEB 4, 2011
!!$       prod_cu = prod_cu/N_r !Scaling
    prod_cu = (1.d0/N_r_pad)*prod_cu !Scaling
    !JLG, FEB 4, 2011
    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

    !Now we need to redistribute the Fourier coefficients on each processor

    t = MPI_WTIME()
    combined_prod_cu(:,:,1)=prod_cu(1,:,:)
    DO n=2, m_max
       !combined_prod_cu(:,:,n)=prod_cu(n,:,:)
       combined_prod_cu(:,:,n)=2*CONJG(prod_cu(n,:,:))
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    t = MPI_WTIME()
    longueur_tranche=bloc_size*m_max_c*nb_field
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu,longueur_tranche,MPID, dist_prod_cu,longueur_tranche, &
         MPID,communicator,code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t
    ! dimensions:
    t = MPI_WTIME()
    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          intermediate = dist_prod_cu(:,:,shiftl+i)
          DO n = 1, bloc_size
             IF (n+shiftc > np ) CYCLE
             DO i_field = 1, nb_field/2
                v_out(n+shiftc, i_field*2-1, i) = REAL (intermediate(i_field,n),KIND=8)
                v_out(n+shiftc, i_field*2 , i)  = AIMAG(intermediate(i_field,n))
             END DO
          END DO
       END DO
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

  END SUBROUTINE FFT_PAR_CROSS_PROD_DCL

  SUBROUTINE FFT_PAR_DOT_PROD_DCL(communicator,V1_in, V2_in, c_out, nb_procs, bloc_size, m_max_pad, temps)
    !FFT (FFT(-1) V1 . FFT(-1) V2) = c_out
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V1_in, V2_in
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT) :: c_out
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    INTEGER, INTENT(IN)                                 :: nb_procs, bloc_size, m_max_pad
    COMPLEX(KIND=8), DIMENSION(m_max_pad, SIZE(V1_in,2), bloc_size)  :: cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,SIZE(V1_in,2),bloc_size)   :: ru
    COMPLEX(KIND=8), DIMENSION(m_max_pad,bloc_size)  :: prod_cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size) :: prod_ru
    COMPLEX(KIND=8), DIMENSION(bloc_size)            :: intermediate
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),2*SIZE(V1_in,2),bloc_size*nb_procs)    :: dist_field, combined_field
    COMPLEX(KIND=8), DIMENSION(bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu
    COMPLEX(KIND=8), DIMENSION(bloc_size,SIZE(V1_in,3)*nb_procs) :: dist_prod_cu
    INTEGER  :: np, np_tot, nb_field,  m_max, m_max_c, MPID,  N_r_pad
    INTEGER(KIND=8) :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    INTEGER ::   nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code
    REAL(KIND=8) :: t
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator

    IF (PRESENT(temps)) temps = 0.d0

    np      = SIZE(V1_in,1)
    nb_field= SIZE(V1_in,2) ! Number of fields
    m_max_c = SIZE(V1_in,3) ! Number of complex (cosines + sines) coefficients per point
    m_max = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    np_tot = nb_procs*bloc_size
    N_r_pad=2*m_max_pad-1

    IF (MOD(nb_field,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_PAR_DOT_PROD_DCL '
       STOP
    END IF

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()

    DO i = 1, m_max_c
       dist_field(i,1:nb_field,1:np) = TRANSPOSE(V1_in(:,:,i))
       dist_field(i,nb_field+1:2*nb_field,1:np) = TRANSPOSE(V2_in(:,:,i))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100

    longueur_tranche=bloc_size*m_max_c*nb_field*2

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, nb_field
             ! Put real and imaginary parts in a complex
             ! nf=1,2,3 => V1_in
             ! nf=4,5,6 => V2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*nb_field


    t = MPI_WTIME()
    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_DOT_PROD: fftw_plan_multi_c2r', fftw_plan_multi_c2r
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    ! DOT PRODUCT
    IF (nb_field==6) THEN
       prod_ru(:,:) = ru(:,1,:)*ru(:,4,:) + ru(:,2,:)*ru(:,5,:) + ru(:,3,:)*ru(:,6,:)
    END IF
    ! DOT PRODUCT

    howmany = bloc_size*1
    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
         inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_DOT_PROD: fftw_plan_multi_r2c', fftw_plan_multi_r2c
    CALL dfftw_execute(fftw_plan_multi_r2c)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    !JLG, FEB 4, 2011
!!$       prod_cu = prod_cu/N_r !Scaling
    prod_cu = prod_cu*(1.d0/N_r_pad) !Scaling
    !JLG, FEB 4, 2011
    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

    !Now we need to redistribute the Fourier coefficients on each processor
    t = MPI_WTIME()
    combined_prod_cu(:,1)=prod_cu(1,:)
    DO n=2, m_max
       combined_prod_cu(:,n)=2*CONJG(prod_cu(n,:))
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    t = MPI_WTIME()
    longueur_tranche=bloc_size*m_max_c*2
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu,longueur_tranche,MPID, dist_prod_cu,longueur_tranche, &
         MPID,communicator,code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          intermediate = dist_prod_cu(:,shiftl+i)
          DO n = 1, bloc_size
             IF (n+shiftc > np ) CYCLE
             c_out(n+shiftc, 1, i) = REAL (intermediate(n),KIND=8)
             c_out(n+shiftc, 2 , i)  = AIMAG(intermediate(n))
          END DO
       END DO
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

  END SUBROUTINE FFT_PAR_DOT_PROD_DCL

  SUBROUTINE FFT_PAR_PROD_DCL(communicator, c1_in, c2_in, c_out, nb_procs, bloc_size, m_max_pad, temps)
    !FFT (FFT(-1) V1 . FFT(-1) V2) = c_out
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    USE my_util
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: c_1in(1:np,1:2,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: c1_in, c2_in
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT) :: c_out
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    INTEGER, INTENT(IN)                                 :: nb_procs, bloc_size, m_max_pad
    COMPLEX(KIND=8), DIMENSION(m_max_pad, 2, bloc_size)             :: cu
    REAL(KIND=8),    DIMENSION(2*m_max_pad-1, 2, bloc_size)         :: ru
    COMPLEX(KIND=8), DIMENSION(m_max_pad, bloc_size)                :: prod_cu
    REAL(KIND=8),    DIMENSION(2*m_max_pad-1,bloc_size)             :: prod_ru
    REAL(KIND=8),    DIMENSION(SIZE(c1_in,3),4, bloc_size*nb_procs) :: dist_field, combined_field
    COMPLEX(KIND=8), DIMENSION(bloc_size,SIZE(c1_in,3)*nb_procs)    :: dist_prod_cu, combined_prod_cu
    COMPLEX(KIND=8), DIMENSION(bloc_size)                           :: intermediate
    INTEGER   :: np, np_tot, m_max, m_max_c, MPID,  N_r_pad
    INTEGER(KIND=8) :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    INTEGER ::   nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code
    REAL(KIND=8) :: t
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator

    IF (PRESENT(temps)) temps = 0.d0

    np      = SIZE(c1_in,1)
    m_max_c = SIZE(c1_in,3) ! Number of complex (cosines + sines) coefficients per point
    m_max = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    np_tot = nb_procs*bloc_size
    N_r_pad=2*m_max_pad-1

    IF (m_max_c==0) THEN
       WRITE(*,*) ' BUG FFT_PAR_PROD_DCL '
       STOP
    END IF

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()
    DO i = 1, m_max_c
       dist_field(i,1:2,1:np) = TRANSPOSE(c1_in(:,:,i))
       dist_field(i,3:4,1:np) = TRANSPOSE(c2_in(:,:,i))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100

    longueur_tranche=bloc_size*m_max_c*4

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, 2
             ! Put real and imaginary parts in a complex
             ! nf=1 => c1_in
             ! nf=2 => c2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*2

    t = MPI_WTIME()
    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    !PRODUCT
    prod_ru(:,:) = ru(:,1,:)*ru(:,2,:)
    !PRODUCT

    howmany = bloc_size*1
    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
         inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE)
    CALL dfftw_execute(fftw_plan_multi_r2c)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    !JLG, FEB 4, 2011
!!$       prod_cu = prod_cu/N_r !Scaling
    prod_cu = prod_cu*(1.d0/N_r_pad) !Scaling
    !JLG, FEB 4, 2011
    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

    !Now we need to redistribute the Fourier coefficients on each processor
    t = MPI_WTIME()
    combined_prod_cu(:,1)=prod_cu(1,:)
    DO n=2, m_max
       combined_prod_cu(:,n)=2*CONJG(prod_cu(n,:))
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    t = MPI_WTIME()
    longueur_tranche=bloc_size*m_max_c*2
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu,longueur_tranche,MPID, dist_prod_cu,longueur_tranche, &
         MPID,communicator,code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          intermediate = dist_prod_cu(:,shiftl+i)
          DO n = 1, bloc_size
             IF (n+shiftc > np ) CYCLE
             c_out(n+shiftc, 1, i) = REAL (intermediate(n),KIND=8)
             c_out(n+shiftc, 2 , i)  = AIMAG(intermediate(n))
          END DO
       END DO
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

  END SUBROUTINE FFT_PAR_PROD_DCL

  SUBROUTINE FFT_PAR_DIV_DCL(communicator, c1_in, c2_in, c_out, nb_procs, bloc_size, m_max_pad, temps)
    !FFT (FFT(-1) V1 . FFT(-1) V2) = c_out
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    USE my_util
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: c_1in(1:np,1:2,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: c1_in, c2_in
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT) :: c_out
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    INTEGER, INTENT(IN)                                 :: nb_procs, bloc_size, m_max_pad
    COMPLEX(KIND=8), DIMENSION(m_max_pad, 2, bloc_size)             :: cu
    REAL(KIND=8),    DIMENSION(2*m_max_pad-1, 2, bloc_size)         :: ru
    COMPLEX(KIND=8), DIMENSION(m_max_pad, bloc_size)                :: prod_cu
    REAL(KIND=8),    DIMENSION(2*m_max_pad-1,bloc_size)             :: prod_ru
    REAL(KIND=8),    DIMENSION(SIZE(c1_in,3),4, bloc_size*nb_procs) :: dist_field, combined_field
    COMPLEX(KIND=8), DIMENSION(bloc_size,SIZE(c1_in,3)*nb_procs)    :: dist_prod_cu, combined_prod_cu
    COMPLEX(KIND=8), DIMENSION(bloc_size)                           :: intermediate
    INTEGER   :: np, np_tot, m_max, m_max_c, MPID,  N_r_pad
    INTEGER(KIND=8) :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    INTEGER ::   nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code
    REAL(KIND=8) :: t
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator

    IF (PRESENT(temps)) temps = 0.d0

    np      = SIZE(c1_in,1)
    m_max_c = SIZE(c1_in,3) ! Number of complex (cosines + sines) coefficients per point
    m_max = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    np_tot = nb_procs*bloc_size
    N_r_pad=2*m_max_pad-1

    IF (m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_PAR_DIV_DCL '
       STOP
    END IF

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()
    DO i = 1, m_max_c
       dist_field(i,1:2,1:np) = TRANSPOSE(c1_in(:,:,i))
       dist_field(i,3:4,1:np) = TRANSPOSE(c2_in(:,:,i))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100

    longueur_tranche=bloc_size*m_max_c*4

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, 2
             ! Put real and imaginary parts in a complex
             ! nf=1 => c1_in
             ! nf=2 => c2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*2

    t = MPI_WTIME()
    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    !DIVISION
    prod_ru(:,:) = ru(:,1,:)/ru(:,2,:)
    !DIVISION

    howmany = bloc_size*1
    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
         inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE)
    CALL dfftw_execute(fftw_plan_multi_r2c)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    !JLG, FEB 4, 2011
!!$       prod_cu = prod_cu/N_r !Scaling
    prod_cu = prod_cu*(1.d0/N_r_pad) !Scaling
    !JLG, FEB 4, 2011
    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

    !Now we need to redistribute the Fourier coefficients on each processor
    t = MPI_WTIME()
    combined_prod_cu(:,1)=prod_cu(1,:)
    DO n=2, m_max
       combined_prod_cu(:,n)=2*CONJG(prod_cu(n,:))
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    t = MPI_WTIME()
    longueur_tranche=bloc_size*m_max_c*2
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu,longueur_tranche,MPID, dist_prod_cu,longueur_tranche, &
         MPID,communicator,code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          intermediate = dist_prod_cu(:,shiftl+i)
          DO n = 1, bloc_size
             IF (n+shiftc > np ) CYCLE
             c_out(n+shiftc, 1, i) = REAL (intermediate(n),KIND=8)
             c_out(n+shiftc, 2 , i)  = AIMAG(intermediate(n))
          END DO
       END DO
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

  END SUBROUTINE FFT_PAR_DIV_DCL

  SUBROUTINE FFT_HEAVISIDE_DCL(communicator, V1_in, values, V_out, nb_procs, bloc_size, &
       m_max_pad, coeff_tanh, temps)
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:,:),INTENT(IN)  :: V1_in
    REAL(KIND=8), DIMENSION(:),      INTENT(IN)  :: values
    REAL(KIND=8),                    INTENT(IN)  :: coeff_tanh
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    INTEGER, INTENT(IN)                                 :: bloc_size, m_max_pad, nb_procs
    INTEGER          :: np, np_tot, nb_field, m_max, m_max_c, MPID, N_r_pad
    INTEGER(KIND=8)  :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    COMPLEX(KIND=8), DIMENSION(m_max_pad, bloc_size)                 :: cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1, bloc_size)                :: ru
    COMPLEX(KIND=8), DIMENSION(m_max_pad,bloc_size)                  :: prod_cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size)                 :: prod_ru
    COMPLEX(KIND=8), DIMENSION(bloc_size)                            :: intermediate
    REAL(KIND=8), DIMENSION(SIZE(V1_in,4),SIZE(V1_in,3),bloc_size*nb_procs)      :: dist_field, combined_field
    COMPLEX(KIND=8), DIMENSION(bloc_size,SIZE(V1_in,4)*nb_procs)                 :: combined_prod_cu
    COMPLEX(KIND=8), DIMENSION(bloc_size,SIZE(V1_in,4)*nb_procs)                 :: dist_prod_cu
    REAL(KIND=8), DIMENSION(SIZE(values),2*m_max_pad-1, bloc_size)               :: V1_real_loc
    INTEGER :: n1, n2
    INTEGER :: nb, shiftc, shiftl, jindex, longueur_tranche, i, n, code, interface_nb
    REAL(KIND=8) :: t
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator

    IF (PRESENT(temps)) temps = 0.d0

    np      = SIZE(V1_in,2)
    nb_field= SIZE(V1_in,3) ! Number of fields
    m_max_c = SIZE(V1_in,4) ! Number of complex (cosines + sines) coefficients per point
    m_max = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    N_r_pad=2*m_max_pad-1
    np_tot = nb_procs*bloc_size

    IF (MOD(nb_field,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_HEAVISIDE_DCL '
       STOP
    END IF

    ! Bloc_size is the number of points that are handled by one processor
    ! once the Fourier modes are all collected
    ! Computation of bloc_size and np_tot
    ! fin de la repartition des points

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()

    DO interface_nb = 1, SIZE(values)-1
       dist_field = 0.d0
       DO i = 1, m_max_c
          dist_field(i,1:nb_field,1:np) = TRANSPOSE(V1_in(interface_nb,:,:,i))
       END DO
       IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

       IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100

       longueur_tranche=bloc_size*m_max_c*nb_field

       t = MPI_WTIME()
       MPID=MPI_DOUBLE_PRECISION
       CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
            MPID, communicator, code)

       IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

       t = MPI_WTIME()
       !JLG, FEB 4, 2011
       cu = 0.d0
       !JLG, FEB 4, 2011
       DO n = 1, bloc_size
          DO nb = 1, nb_procs
             shiftc = (nb-1)*bloc_size
             shiftl = (nb-1)*m_max_c
             jindex = n + shiftc
             ! Put real and imaginary parts in a complex
             ! nf=1,2,3 => V1_in
             ! nf=4,5,6 => V2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,n) = 0.5d0*CMPLX(combined_field(:,1,jindex),&
                  -combined_field(:,2,jindex),KIND=8)
          END DO
       END DO
       cu(1,:) = 2*CMPLX(REAL(cu(1,:),KIND=8),0.d0,KIND=8)


       IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

       ! Set the parameters for dfftw
       fft_dim=1; istride=1; ostride=1;
       !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
       idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
       odist=m_max_pad; onembed(1)=m_max_pad
       !JLG, FEB 4, 2011

       howmany=bloc_size*nb_field/2

       t = MPI_WTIME()
       CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
            onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
       !write(*,*) ' FFT_PAR_CROSS_PROD: fftw_plan_multi_c2r', fftw_plan_multi_c2r
       CALL dfftw_execute(fftw_plan_multi_c2r)

       ! clean up
       CALL dfftw_destroy_plan(fftw_plan_multi_c2r)
       V1_real_loc(interface_nb,:,:)=ru
    END DO

    ! CROSS PRODUCT
    IF (nb_field==2) THEN
       DO n1 = 1, 2*m_max_pad-1
          DO n2 = 1, bloc_size
             prod_ru(n1,n2) = values(1)
             DO interface_nb = 1, SIZE(values)-1
                !reconstruction convex
                prod_ru(n1,n2) = prod_ru(n1,n2)*(1.d0-regul(V1_real_loc(interface_nb,n1,n2),coeff_tanh)) + &
                     values(interface_nb+1)*regul(V1_real_loc(interface_nb,n1,n2),coeff_tanh)
             END DO
          END DO
       END DO
    END IF
    ! CROSS PRODUCT

    howmany = howmany
    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
         inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_CROSS_PROD: fftw_plan_multi_r2c', fftw_plan_multi_r2c
    CALL dfftw_execute(fftw_plan_multi_r2c)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    !JLG, FEB 4, 2011
!!$       prod_cu = prod_cu/N_r !Scaling
    prod_cu = prod_cu*(1.d0/N_r_pad) !Scaling
    !JLG, FEB 4, 2011
    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t
    !Now we need to redistribute the Fourier coefficients on each processor

    t = MPI_WTIME()
    combined_prod_cu(:,1)=prod_cu(1,:)
    DO n=2, m_max
       !combined_prod_cu(:,:,n)=prod_cu(n,:,:)
       combined_prod_cu(:,n)=2*CONJG(prod_cu(n,:))
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    t = MPI_WTIME()
    longueur_tranche=bloc_size*m_max_c*nb_field
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu,longueur_tranche,MPID, dist_prod_cu,longueur_tranche, &
         MPID,communicator,code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t
    ! dimensions:
    t = MPI_WTIME()
    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          intermediate = dist_prod_cu(:,shiftl+i)
          DO n = 1, bloc_size
             IF (n+shiftc > np ) CYCLE
             v_out(n+shiftc, 1, i) = REAL (intermediate(n),KIND=8)
             v_out(n+shiftc, 2, i)  = AIMAG(intermediate(n))
          END DO
       END DO
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

  END SUBROUTINE FFT_HEAVISIDE_DCL

  SUBROUTINE FFT_SCALAR_VECT_DCL(communicator,V1_in, V2_in, V_out, pb, nb_procs, &
       bloc_size, m_max_pad, temps)
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    USE my_util
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V1_in ! VECTOR
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V2_in ! SCALAR
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    INTEGER, INTENT(IN)                                 :: bloc_size, m_max_pad, nb_procs
    INTEGER, INTENT(IN)                                 :: pb
    INTEGER          :: np, np_tot, nb_field1, nb_field2, m_max, m_max_c, MPID, N_r_pad
    INTEGER(KIND=8)  :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    COMPLEX(KIND=8), DIMENSION(m_max_pad, (SIZE(V1_in,2)+SIZE(V2_in,2))/2, bloc_size)  :: cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,(SIZE(V1_in,2)+SIZE(V2_in,2))/2,bloc_size)   :: ru
    COMPLEX(KIND=8), DIMENSION(m_max_pad,SIZE(V1_in,2)/2,bloc_size)  :: prod_cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,SIZE(V1_in,2)/2,bloc_size) :: prod_ru
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2, bloc_size)           :: intermediate
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),SIZE(V1_in,2)+SIZE(V2_in,2),bloc_size*nb_procs) :: dist_field
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),SIZE(V1_in,2)+SIZE(V2_in,2),bloc_size*nb_procs) :: combined_field
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: dist_prod_cu
    INTEGER :: i_field
    INTEGER :: nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code
    REAL(KIND=8) :: t
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator

    IF (PRESENT(temps)) temps = 0.d0

    np      = SIZE(V1_in,1)
    nb_field1= SIZE(V1_in,2) ! Number of fields
    nb_field2= SIZE(V2_in,2) ! Number of fields
    m_max_c = SIZE(V1_in,3) ! Number of complex (cosines + sines) coefficients per point
    m_max = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    N_r_pad=2*m_max_pad-1
    np_tot = nb_procs*bloc_size

    IF (MOD(nb_field1,2)/=0 .OR. MOD(nb_field2,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_SCALAR_VECT_DCL '
       STOP
    END IF

    ! Bloc_size is the number of points that are handled by one processor
    ! once the Fourier modes are all collected
    ! Computation of bloc_size and np_tot
    ! fin de la repartition des points

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()

    DO i = 1, m_max_c
       dist_field(i,1:nb_field1,1:np) = TRANSPOSE(V1_in(:,:,i))
       dist_field(i,nb_field1+1:nb_field1+nb_field2,1:np) = TRANSPOSE(V2_in(:,:,i))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100

    longueur_tranche=bloc_size*m_max_c*(nb_field1+nb_field2)

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, (nb_field1+nb_field2)/2
             ! Put real and imaginary parts in a complex
             ! nf=1,2,3 => V1_in
             ! nf=4 => V2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*(nb_field1+nb_field2)/2

    t = MPI_WTIME()
    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_CROSS_PROD: fftw_plan_multi_c2r', fftw_plan_multi_c2r
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    IF (nb_field1 == 6 .AND. nb_field2 == 2) THEN
       IF (pb==1) THEN
          prod_ru(:,1,:) = ru(:,4,:)*ru(:,1,:)
          prod_ru(:,2,:) = ru(:,4,:)*ru(:,2,:)
          prod_ru(:,3,:) = ru(:,4,:)*ru(:,3,:)
       ELSE IF (pb==2) THEN
          prod_ru(:,1,:) = ru(:,1,:)/ru(:,4,:)
          prod_ru(:,2,:) = ru(:,2,:)/ru(:,4,:)
          prod_ru(:,3,:) = ru(:,3,:)/ru(:,4,:)
       ELSE
          CALL error_petsc('error in problem type while calling FFT_SCALAR_VECT_DCL ')
       END IF
    ELSE
       CALL error_petsc('error in dimension of INPUT variables while calling FFT_SCALAR_VECT_DCL ')
    END IF

    howmany = bloc_size*nb_field1/2

    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
         inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_CROSS_PROD: fftw_plan_multi_r2c', fftw_plan_multi_r2c
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)



    !JLG, FEB 4, 2011
!!$       prod_cu = prod_cu/N_r !Scaling
    prod_cu = prod_cu*(1.d0/N_r_pad) !Scaling
    !JLG, FEB 4, 2011
    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

    !Now we need to redistribute the Fourier coefficients on each processor

    t = MPI_WTIME()
    combined_prod_cu(:,:,1)=prod_cu(1,:,:)
    DO n=2, m_max
       !combined_prod_cu(:,:,n)=prod_cu(n,:,:)
       combined_prod_cu(:,:,n)=2*CONJG(prod_cu(n,:,:))
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    longueur_tranche=bloc_size*m_max_c*nb_field1

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu,longueur_tranche,MPID, dist_prod_cu,longueur_tranche, &
         MPID,communicator,code)

    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t
    ! dimensions:
    t = MPI_WTIME()

    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          intermediate = dist_prod_cu(:,:,shiftl+i)
          DO n = 1, bloc_size
             IF (n+shiftc > np ) CYCLE
             DO i_field = 1, nb_field1/2
                v_out(n+shiftc, i_field*2-1, i) = REAL (intermediate(i_field,n),KIND=8)
                v_out(n+shiftc, i_field*2 , i)  = AIMAG(intermediate(i_field,n))
             END DO
          END DO
       END DO
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

  END SUBROUTINE FFT_SCALAR_VECT_DCL

  SUBROUTINE FFT_SCALAR_VECT_NO_OVERSHOOT(communicator, scalar_bounds, V1_in, V2_in, V_out, pb, nb_procs, &
       bloc_size, m_max_pad, temps)
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    USE my_util
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V1_in ! VECTOR
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V2_in ! SCALAR
    REAL(KIND=8), DIMENSION(:),      INTENT(IN)  :: scalar_bounds
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    INTEGER, INTENT(IN)                                 :: bloc_size, m_max_pad, nb_procs
    INTEGER, INTENT(IN)                                 :: pb
    INTEGER          :: np, np_tot, nb_field1, nb_field2, m_max, m_max_c, MPID, N_r_pad
    INTEGER(KIND=8)  :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    COMPLEX(KIND=8), DIMENSION(m_max_pad, (SIZE(V1_in,2)+SIZE(V2_in,2))/2, bloc_size)  :: cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,(SIZE(V1_in,2)+SIZE(V2_in,2))/2,bloc_size)   :: ru
    COMPLEX(KIND=8), DIMENSION(m_max_pad,SIZE(V1_in,2)/2,bloc_size)  :: prod_cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,SIZE(V1_in,2)/2,bloc_size) :: prod_ru
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2, bloc_size)           :: intermediate
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),SIZE(V1_in,2)+SIZE(V2_in,2),bloc_size*nb_procs) :: dist_field
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),SIZE(V1_in,2)+SIZE(V2_in,2),bloc_size*nb_procs) :: combined_field
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: dist_prod_cu
    INTEGER :: i_field
    INTEGER :: nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code
    REAL(KIND=8) :: t, scalar_min, scalar_max
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator

    IF (PRESENT(temps)) temps = 0.d0

    np      = SIZE(V1_in,1)
    nb_field1= SIZE(V1_in,2) ! Number of fields
    nb_field2= SIZE(V2_in,2) ! Number of fields
    m_max_c = SIZE(V1_in,3) ! Number of complex (cosines + sines) coefficients per point
    m_max = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    N_r_pad=2*m_max_pad-1
    np_tot = nb_procs*bloc_size

    IF (MOD(nb_field1,2)/=0 .OR. MOD(nb_field2,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_SCALAR_VECT_NO_OVERSHOOT '
       STOP
    END IF

    ! Bloc_size is the number of points that are handled by one processor
    ! once the Fourier modes are all collected
    ! Computation of bloc_size and np_tot
    ! fin de la repartition des points

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()

    DO i = 1, m_max_c
       dist_field(i,1:nb_field1,1:np) = TRANSPOSE(V1_in(:,:,i))
       dist_field(i,nb_field1+1:nb_field1+nb_field2,1:np) = TRANSPOSE(V2_in(:,:,i))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100

    longueur_tranche=bloc_size*m_max_c*(nb_field1+nb_field2)

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, (nb_field1+nb_field2)/2
             ! Put real and imaginary parts in a complex
             ! nf=1,2,3 => V1_in
             ! nf=4 => V2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*(nb_field1+nb_field2)/2

    t = MPI_WTIME()
    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_CROSS_PROD: fftw_plan_multi_c2r', fftw_plan_multi_c2r
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    IF (nb_field1 == 6 .AND. nb_field2 == 2) THEN
       scalar_min = MINVAL(scalar_bounds)
       scalar_max = MAXVAL(scalar_bounds)
       IF (pb==1) THEN
          prod_ru(:,1,:) =  MAX(scalar_min,MIN(ru(:,4,:),scalar_max))*ru(:,1,:)
          prod_ru(:,2,:) =  MAX(scalar_min,MIN(ru(:,4,:),scalar_max))*ru(:,2,:)
          prod_ru(:,3,:) =  MAX(scalar_min,MIN(ru(:,4,:),scalar_max))*ru(:,3,:)
       ELSE IF (pb==2) THEN
          prod_ru(:,1,:) = ru(:,1,:)/MAX(scalar_min,MIN(ru(:,4,:),scalar_max))
          prod_ru(:,2,:) = ru(:,2,:)/MAX(scalar_min,MIN(ru(:,4,:),scalar_max))
          prod_ru(:,3,:) = ru(:,3,:)/MAX(scalar_min,MIN(ru(:,4,:),scalar_max))
       ELSE
          CALL error_petsc('error in problem type while calling FFT_SCALAR_VECT_DCL ')
       END IF
    ELSE
       CALL error_petsc('error in dimension of INPUT variables while calling FFT_SCALAR_VECT_DCL ')
    END IF

    howmany = bloc_size*nb_field1/2

    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
         inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_CROSS_PROD: fftw_plan_multi_r2c', fftw_plan_multi_r2c
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)



    !JLG, FEB 4, 2011
!!$       prod_cu = prod_cu/N_r !Scaling
    prod_cu = prod_cu*(1.d0/N_r_pad) !Scaling
    !JLG, FEB 4, 2011
    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

    !Now we need to redistribute the Fourier coefficients on each processor

    t = MPI_WTIME()
    combined_prod_cu(:,:,1)=prod_cu(1,:,:)
    DO n=2, m_max
       !combined_prod_cu(:,:,n)=prod_cu(n,:,:)
       combined_prod_cu(:,:,n)=2*CONJG(prod_cu(n,:,:))
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    longueur_tranche=bloc_size*m_max_c*nb_field1

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu,longueur_tranche,MPID, dist_prod_cu,longueur_tranche, &
         MPID,communicator,code)

    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t
    ! dimensions:
    t = MPI_WTIME()

    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          intermediate = dist_prod_cu(:,:,shiftl+i)
          DO n = 1, bloc_size
             IF (n+shiftc > np ) CYCLE
             DO i_field = 1, nb_field1/2
                v_out(n+shiftc, i_field*2-1, i) = REAL (intermediate(i_field,n),KIND=8)
                v_out(n+shiftc, i_field*2 , i)  = AIMAG(intermediate(i_field,n))
             END DO
          END DO
       END DO
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

  END SUBROUTINE FFT_SCALAR_VECT_NO_OVERSHOOT

  SUBROUTINE FFT_MAX_VEL_DCL(communicator, V1_in, V_out, nb_procs, bloc_size, m_max_pad)
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    USE my_util
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V1_in ! VECTOR
    REAL(KIND=8),                    INTENT(OUT) :: V_out
    INTEGER, INTENT(IN)                          :: bloc_size, m_max_pad, nb_procs
    INTEGER          :: np, np_tot, nb_field1, m_max, m_max_c, MPID, N_r_pad
    INTEGER(KIND=8)  :: fftw_plan_multi_c2r
    COMPLEX(KIND=8), DIMENSION(m_max_pad, SIZE(V1_in,2)/2, bloc_size)  :: cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1, SIZE(V1_in,2)/2,bloc_size)  :: ru
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),SIZE(V1_in,2),bloc_size*nb_procs) :: dist_field, combined_field
    INTEGER :: nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code
    REAL(KIND=8), DIMENSION(2*m_max_pad-1, bloc_size) ::  norm_vel_loc
    REAL(KIND=8) :: max_velocity
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator

    np        = SIZE(V1_in,1)
    nb_field1 = SIZE(V1_in,2) ! Number of fields
    m_max_c   = SIZE(V1_in,3) ! Number of complex (cosines + sines) coefficients per point
    m_max     = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    N_r_pad   = 2*m_max_pad-1
    np_tot    = nb_procs*bloc_size

    IF (MOD(nb_field1,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_MAX_VEL_DCL '
       STOP
    END IF

    ! Bloc_size is the number of points that are handled by one processor
    ! once the Fourier modes are all collected
    ! Computation of bloc_size and np_tot
    ! fin de la repartition des points

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT

    DO i = 1, m_max_c
       dist_field(i,1:nb_field1,1:np) = TRANSPOSE(V1_in(:,:,i))
    END DO

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 0.d0

    longueur_tranche=bloc_size*m_max_c*nb_field1

    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)

    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, nb_field1/2
             ! Put real and imaginary parts in a complex
             ! nf=1,2,3 => V1_in
             ! nf=4 => V2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*nb_field1/2

    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_CROSS_PROD: fftw_plan_multi_c2r', fftw_plan_multi_c2r
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    IF (nb_field1==6) THEN
       norm_vel_loc(:,:) = SQRT(ru(:,1,:)**2 + ru(:,2,:)**2 + ru(:,3,:)**2)
    ELSE
       norm_vel_loc(:,:) = SQRT(ru(:,1,:)**2)
    END IF
    max_velocity = MAXVAL(norm_vel_loc)
    CALL MPI_ALLREDUCE(max_velocity, V_out, 1, MPI_DOUBLE_PRECISION, &
         MPI_MAX, communicator, code)

  END SUBROUTINE FFT_MAX_VEL_DCL

  SUBROUTINE FFT_TENSOR_DCL(communicator,V1_in, V2_in, V_out, nb_procs, bloc_size, m_max_pad, temps, opt_tension)
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    USE my_util
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),    INTENT(IN)  :: V1_in, V2_in
    REAL(KIND=8), DIMENSION(:,:,:,:),  INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    LOGICAL,                    OPTIONAL, INTENT(IN)    :: opt_tension
    INTEGER, INTENT(IN)                                 :: bloc_size, m_max_pad, nb_procs
    INTEGER          :: np, np_tot, nb_field, m_max, m_max_c, MPID, N_r_pad
    INTEGER(KIND=8)  :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    COMPLEX(KIND=8), DIMENSION(m_max_pad, SIZE(V1_in,2), bloc_size)   :: cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,SIZE(V1_in,2),bloc_size)    :: ru
    COMPLEX(KIND=8), DIMENSION(m_max_pad,SIZE(V1_in,2)/2,bloc_size)   :: prod_cu1, prod_cu2, prod_cu3
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,SIZE(V1_in,2)/2,bloc_size)  :: prod_ru1, prod_ru2, prod_ru3
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),2*SIZE(V1_in,2),bloc_size*nb_procs)    :: dist_field
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),2*SIZE(V1_in,2),bloc_size*nb_procs)    :: combined_field
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu1
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu2
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu3
    COMPLEX(KIND=8), DIMENSION(3,SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu_tot
    COMPLEX(KIND=8), DIMENSION(3,SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: dist_prod_cu_tot
    COMPLEX(KIND=8), DIMENSION(3,SIZE(V1_in,2)/2,bloc_size)                        :: intermediate_tot
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size)                  :: norm_grad_tension
    INTEGER :: i_field
    INTEGER :: nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code
    REAL(KIND=8) :: t
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator

    IF (PRESENT(temps)) temps = 0.d0

    np      = SIZE(V1_in,1)
    nb_field= SIZE(V1_in,2) ! Number of fields
    m_max_c = SIZE(V1_in,3) ! Number of complex (cosines + sines) coefficients per point
    m_max = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    N_r_pad=2*m_max_pad-1
    np_tot = nb_procs*bloc_size

    IF (MOD(nb_field,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_TENSOR_DCL '
       STOP
    END IF

    ! Bloc_size is the number of points that are handled by one processor
    ! once the Fourier modes are all collected
    ! Computation of bloc_size and np_tot
    ! fin de la repartition des points

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()

    DO i = 1, m_max_c
       dist_field(i,1:nb_field,1:np) = TRANSPOSE(V1_in(:,:,i))
       dist_field(i,nb_field+1:2*nb_field,1:np) = TRANSPOSE(V2_in(:,:,i))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100

    longueur_tranche=bloc_size*m_max_c*nb_field*2

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, nb_field
             ! Put real and imaginary parts in a complex
             ! nf=1,2,3 => V1_in
             ! nf=4,5,6 => V2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*nb_field


    t = MPI_WTIME()
    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_CROSS_PROD: fftw_plan_multi_c2r', fftw_plan_multi_c2r
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)


    ! TENSOR PRODUCT
    IF(.NOT.PRESENT(opt_tension)) THEN
       IF (nb_field==6) THEN
          prod_ru1(:,1,:) = ru(:,1,:)*ru(:,4,:)
          prod_ru1(:,2,:) = ru(:,1,:)*ru(:,5,:)
          prod_ru1(:,3,:) = ru(:,1,:)*ru(:,6,:)

          prod_ru2(:,1,:) = ru(:,2,:)*ru(:,4,:)
          prod_ru2(:,2,:) = ru(:,2,:)*ru(:,5,:)
          prod_ru2(:,3,:) = ru(:,2,:)*ru(:,6,:)

          prod_ru3(:,1,:) = ru(:,3,:)*ru(:,4,:)
          prod_ru3(:,2,:) = ru(:,3,:)*ru(:,5,:)
          prod_ru3(:,3,:) = ru(:,3,:)*ru(:,6,:)
       END IF
    ELSE IF (opt_tension) THEN
       IF (nb_field==6) THEN
          norm_grad_tension = SQRT(ru(:,4,:)**2 + ru(:,5,:)**2 + ru(:,6,:)**2)
          prod_ru1(:,1,:) = ru(:,1,:)*ru(:,4,:)/(norm_grad_tension + 1.d-14) &
               - norm_grad_tension
          prod_ru1(:,2,:) = ru(:,1,:)*ru(:,5,:)/(norm_grad_tension + 1.d-14)
          prod_ru1(:,3,:) = ru(:,1,:)*ru(:,6,:)/(norm_grad_tension + 1.d-14)

          prod_ru2(:,1,:) = ru(:,2,:)*ru(:,4,:)/(norm_grad_tension + 1.d-14)
          prod_ru2(:,2,:) = ru(:,2,:)*ru(:,5,:)/(norm_grad_tension + 1.d-14) &
               - norm_grad_tension
          prod_ru2(:,3,:) = ru(:,2,:)*ru(:,6,:)/(norm_grad_tension + 1.d-14)

          prod_ru3(:,1,:) = ru(:,3,:)*ru(:,4,:)/(norm_grad_tension + 1.d-14)
          prod_ru3(:,2,:) = ru(:,3,:)*ru(:,5,:)/(norm_grad_tension + 1.d-14)
          prod_ru3(:,3,:) = ru(:,3,:)*ru(:,6,:)/(norm_grad_tension + 1.d-14) &
               - norm_grad_tension
       END IF
    ELSE
       CALL error_petsc('BUG in FFT_TENSOR_DCL : opt_tension should be true if used')
    END IF
    ! TENSOR PRODUCT

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*nb_field/2

    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru1, &
         inembed, istride, idist, prod_cu1, onembed, ostride, odist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_CROSS_PROD: fftw_plan_multi_r2c', fftw_plan_multi_r2c
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru2, &
         inembed, istride, idist, prod_cu2, onembed, ostride, odist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_CROSS_PROD: fftw_plan_multi_r2c', fftw_plan_multi_r2c
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru3, &
         inembed, istride, idist, prod_cu3, onembed, ostride, odist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_CROSS_PROD: fftw_plan_multi_r2c', fftw_plan_multi_r2c
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    !JLG, FEB 4, 2011
!!$       prod_cu = prod_cu/N_r !Scaling
!!$    prod_cu = prod_cu/N_r_pad !Scaling

    prod_cu1 = prod_cu1*(1.d0/N_r_pad) !Scaling
    prod_cu2 = prod_cu2*(1.d0/N_r_pad) !Scaling
    prod_cu3 = prod_cu3*(1.d0/N_r_pad) !Scaling

    !JLG, FEB 4, 2011
    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

    !Now we need to redistribute the Fourier coefficients on each processor
    t = MPI_WTIME()
    combined_prod_cu1(:,:,1)=prod_cu1(1,:,:)
    combined_prod_cu2(:,:,1)=prod_cu2(1,:,:)
    combined_prod_cu3(:,:,1)=prod_cu3(1,:,:)
    DO n=2, m_max
       !combined_prod_cu(:,:,n)=prod_cu(n,:,:)
       combined_prod_cu1(:,:,n)=2*CONJG(prod_cu1(n,:,:))
       combined_prod_cu2(:,:,n)=2*CONJG(prod_cu2(n,:,:))
       combined_prod_cu3(:,:,n)=2*CONJG(prod_cu3(n,:,:))
    END DO

    combined_prod_cu_tot(1,:,:,:) = combined_prod_cu1(:,:,:)
    combined_prod_cu_tot(2,:,:,:) = combined_prod_cu2(:,:,:)
    combined_prod_cu_tot(3,:,:,:) = combined_prod_cu3(:,:,:)

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    t = MPI_WTIME()
    longueur_tranche=bloc_size*m_max_c*nb_field*3
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu_tot,longueur_tranche,MPID, dist_prod_cu_tot,longueur_tranche, &
         MPID,communicator,code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t
    ! dimensions:
    t = MPI_WTIME()
    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          intermediate_tot = dist_prod_cu_tot(:,:,:,shiftl+i)
          DO n = 1, bloc_size
             IF (n+shiftc > np ) CYCLE
             DO i_field = 1, nb_field/2
                v_out(:, n+shiftc, i_field*2-1, i) = REAL (intermediate_tot(:,i_field,n),KIND=8)
                v_out(:, n+shiftc, i_field*2 , i)  = AIMAG(intermediate_tot(:,i_field,n))
             END DO
          END DO
       END DO
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t
  END SUBROUTINE FFT_TENSOR_DCL

  SUBROUTINE FFT_PAR_VAR_ETA_PROD_T_DCL(communicator, eta, H_mesh, &
       c_in, c_out, nb_procs, bloc_size, m_max_pad, time, temps)
    !FFT (FFT(-1) V1 . FFT(-1) V2) = c_out
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    USE my_util
    USE def_type_mesh
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    TYPE(mesh_type)                                           :: H_mesh
    INTERFACE
       FUNCTION eta(H_mesh,angles,nb_angles,nb,ne,time) RESULT(vv)
         USE def_type_mesh
         USE input_data
         USE my_util
         IMPLICIT NONE
         TYPE(mesh_type), INTENT(IN)                :: H_mesh
         REAL(KIND=8), DIMENSION(:), INTENT(IN)     :: angles
         INTEGER, INTENT(IN)                        :: nb_angles
         INTEGER, INTENT(IN)                        :: nb, ne
         REAL(KIND=8), INTENT(IN)                   :: time
         REAL(KIND=8), DIMENSION(nb_angles,ne-nb+1) :: vv
       END FUNCTION eta
    END INTERFACE
    INTEGER, INTENT(IN)                                  :: nb_procs, bloc_size, m_max_pad
    REAL(KIND=8), DIMENSION(2*m_max_pad-1)                    :: angles
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size)          :: eta_n
    REAL(KIND=8)                                              :: pi
    INTEGER                  :: rank_F
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: c_in
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT) :: c_out
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT)  :: temps
    REAL(KIND=8)                               :: time
    COMPLEX(KIND=8), DIMENSION(m_max_pad, SIZE(c_in,2)/2, bloc_size)          :: cu
    REAL(KIND=8),    DIMENSION(2*m_max_pad-1, SIZE(c_in,2)/2, bloc_size)      :: ru
    COMPLEX(KIND=8), DIMENSION(m_max_pad, SIZE(c_in,2)/2, bloc_size)          :: prod_cu
    REAL(KIND=8),    DIMENSION(2*m_max_pad-1, SIZE(c_in,2)/2,bloc_size)       :: prod_ru
    REAL(KIND=8),    DIMENSION(SIZE(c_in,3),SIZE(c_in,2), bloc_size*nb_procs) :: dist_field, combined_field
    COMPLEX(KIND=8), DIMENSION(SIZE(c_in,2)/2,bloc_size,SIZE(c_in,3)*nb_procs):: dist_prod_cu
    COMPLEX(KIND=8), DIMENSION(SIZE(c_in,2)/2,bloc_size,SIZE(c_in,3)*nb_procs):: combined_prod_cu
    COMPLEX(KIND=8), DIMENSION(SIZE(c_in,2)/2,bloc_size)                      :: intermediate
    INTEGER   :: np, np_tot, m_max, m_max_c, MPID,  N_r_pad, nb_field, i_field
    INTEGER(KIND=8) :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    INTEGER ::   nb, ne, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code
    REAL(KIND=8) :: t
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    PetscErrorCode :: ierr
    MPI_Comm :: communicator

    IF (PRESENT(temps)) temps = 0.d0

    np       = SIZE(c_in,1)
    nb_field = SIZE(c_in,2)
    m_max_c  = SIZE(c_in,3) ! Number of complex (cosines + sines) coefficients per point
    m_max    = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    np_tot   = nb_procs*bloc_size
    N_r_pad  = 2*m_max_pad-1

    IF (m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_PAR_VAR_ETA_PROD_T_DCL '
       STOP
    END IF

    pi = ACOS(-1.d0)
    DO n = 1, N_r_pad
       angles(n) = 2*pi*(n-1)/N_r_pad
    END DO

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()
    DO i = 1, m_max_c
       dist_field(i,:,1:np) = TRANSPOSE(c_in(:,:,i))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100

    longueur_tranche=bloc_size*m_max_c*nb_field

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, nb_field/2
             ! Put real and imaginary parts in a complex
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*nb_field/2

    t = MPI_WTIME()
    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

!!!!!!!!!!!!!!!!!!!!!!!!!! PRODUCT with eta !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check if the proc is the last one to know if you need to add 1.d0 to eta
    CALL MPI_Comm_rank(communicator, rank_F, ierr)
    IF (rank_F==nb_procs-1) THEN
       IF (np-rank_F*bloc_size.GE.1) THEN
          nb = rank_F*bloc_size+1
          ne = np
          eta_n(:,1:np-rank_F*bloc_size) = ETA(H_mesh,angles,N_r_pad,nb,ne,time)
          eta_n(:,np-rank_F*bloc_size+1:bloc_size)  = 1.d0
       ELSE
          eta_n=1.d0
       END IF
    ELSE
       nb = rank_F*bloc_size+1
       ne = rank_F*bloc_size+bloc_size
       eta_n(:,:) = ETA(H_mesh,angles,N_r_pad,nb,ne,time)
    END IF

    IF (nb_field==2) THEN
       prod_ru(:,1,:) = eta_n*ru(:,1,:)
    ELSE IF (nb_field==6) THEN
       prod_ru(:,1,:) = eta_n*ru(:,1,:)
       prod_ru(:,2,:) = eta_n*ru(:,2,:)
       prod_ru(:,3,:) = eta_n*ru(:,3,:)
    ELSE
       CALL error_petsc('error in nb_field of c_in while calling FFT_PAR_VAR_ETA_PROD_T_DCL')
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!! PRODUCT with eta  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    howmany = bloc_size*nb_field/2

    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
         inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE)
    CALL dfftw_execute(fftw_plan_multi_r2c)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    !JLG, FEB 4, 2011
!!$       prod_cu = prod_cu/N_r !Scaling
    prod_cu = prod_cu*(1.d0/N_r_pad) !Scaling
    !JLG, FEB 4, 2011
    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

    !Now we need to redistribute the Fourier coefficients on each processor
    t = MPI_WTIME()
    combined_prod_cu(:,:,1)=prod_cu(1,:,:)
    DO n=2, m_max
       combined_prod_cu(:,:,n)=2*CONJG(prod_cu(n,:,:))
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    t = MPI_WTIME()
    longueur_tranche=bloc_size*m_max_c*nb_field
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu,longueur_tranche,MPID, dist_prod_cu,longueur_tranche, &
         MPID,communicator,code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          intermediate = dist_prod_cu(:,:,shiftl+i)
          DO n = 1, bloc_size
             IF (n+shiftc > np ) CYCLE
             DO i_field = 1, nb_field/2
                c_out(n+shiftc, 2*i_field-1, i) = REAL (intermediate(i_field,n),KIND=8)
                c_out(n+shiftc, 2*i_field, i)   = AIMAG(intermediate(i_field,n))
             END DO
          END DO
       END DO
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

  END SUBROUTINE FFT_PAR_VAR_ETA_PROD_T_DCL

  SUBROUTINE FFT_PAR_VAR_ETA_PROD_GAUSS_DCL(communicator, eta, H_mesh, &
       c_in, c_out, nb_procs, bloc_size, m_max_pad, rr_gauss, time, temps)
    !FFT (FFT(-1) V1 . FFT(-1) V2) = c_out
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    USE my_util
    USE def_type_mesh
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    TYPE(mesh_type)                                           :: H_mesh
    INTERFACE
       FUNCTION eta(H_mesh,rr_gauss,angles,nb_angles,nb,ne,time) RESULT(vv)
         USE def_type_mesh
         USE input_data
         USE my_util
         IMPLICIT NONE
         TYPE(mesh_type)                            :: H_mesh
         REAL(KIND=8), DIMENSION(:,:), INTENT(IN)   :: rr_gauss
         REAL(KIND=8), DIMENSION(:), INTENT(IN)     :: angles
         INTEGER,                    INTENT(IN)     :: nb_angles
         INTEGER,                    INTENT(IN)     :: nb, ne
         REAL(KIND=8),               INTENT(IN)     :: time
         REAL(KIND=8), DIMENSION(nb_angles,ne-nb+1) :: vv
       END FUNCTION eta
    END INTERFACE
    INTEGER, INTENT(IN)                                  :: nb_procs, bloc_size, m_max_pad
    REAL(KIND=8), DIMENSION(2*m_max_pad-1)                    :: angles
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size)          :: eta_n
    REAL(KIND=8)                                              :: pi
    INTEGER                                                   :: rank_F
    REAL(KIND=8)                                 :: time
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: c_in
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT) :: c_out
    REAL(KIND=8), DIMENSION(:,:),    INTENT(IN)  :: rr_gauss
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT)  :: temps
    COMPLEX(KIND=8), DIMENSION(m_max_pad, SIZE(c_in,2)/2, bloc_size)          :: cu
    REAL(KIND=8),    DIMENSION(2*m_max_pad-1, SIZE(c_in,2)/2, bloc_size)      :: ru
    COMPLEX(KIND=8), DIMENSION(m_max_pad, SIZE(c_in,2)/2, bloc_size)          :: prod_cu
    REAL(KIND=8),    DIMENSION(2*m_max_pad-1, SIZE(c_in,2)/2,bloc_size)       :: prod_ru
    REAL(KIND=8),    DIMENSION(SIZE(c_in,3),SIZE(c_in,2), bloc_size*nb_procs) :: dist_field, combined_field
    COMPLEX(KIND=8), DIMENSION(SIZE(c_in,2)/2,bloc_size,SIZE(c_in,3)*nb_procs):: dist_prod_cu
    COMPLEX(KIND=8), DIMENSION(SIZE(c_in,2)/2,bloc_size,SIZE(c_in,3)*nb_procs):: combined_prod_cu
    COMPLEX(KIND=8), DIMENSION(SIZE(c_in,2)/2,bloc_size)                      :: intermediate
    INTEGER   :: np, np_tot, m_max, m_max_c, MPID,  N_r_pad, nb_field, i_field
    INTEGER(KIND=8) :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    INTEGER ::   nb, ne, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code
    REAL(KIND=8) :: t
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    PetscErrorCode :: ierr
    MPI_Comm :: communicator

    IF (PRESENT(temps)) temps = 0.d0

    np       = SIZE(c_in,1)
    nb_field = SIZE(c_in,2)
    m_max_c  = SIZE(c_in,3) ! Number of complex (cosines + sines) coefficients per point
    m_max    = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    np_tot   = nb_procs*bloc_size
    N_r_pad  = 2*m_max_pad-1

    IF (m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_PAR_VAR_ETA_PROD_GAUSS_DCL '
       STOP
    END IF

    pi = ACOS(-1.d0)
    DO n = 1, N_r_pad
       angles(n) = 2*pi*(n-1)/N_r_pad
    END DO

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()
    DO i = 1, m_max_c
       dist_field(i,:,1:np) = TRANSPOSE(c_in(:,:,i))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100

    longueur_tranche=bloc_size*m_max_c*nb_field

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, nb_field/2
             ! Put real and imaginary parts in a complex
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*nb_field/2

    t = MPI_WTIME()
    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

!!!!!!!!!!!!!!!!!!!!!!!!!! PRODUCT with eta !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check if the proc is the last one to know if you need to add 1.d0 to eta
    CALL MPI_Comm_rank(communicator, rank_F, ierr)
    IF (rank_F==nb_procs-1) THEN
       IF (np-rank_F*bloc_size.GE.1) THEN
          nb = rank_F*bloc_size+1
          ne = np
          eta_n(:,1:np-rank_F*bloc_size) = ETA(H_mesh,rr_gauss(:,nb:ne),angles,N_r_pad,nb,ne,time)
          eta_n(:,np-rank_F*bloc_size+1:bloc_size)  = 1.d0
       ELSE
          eta_n=1.d0
       END IF
    ELSE
       nb = rank_F*bloc_size+1
       ne = rank_F*bloc_size+bloc_size
       eta_n(:,:) = ETA(H_mesh,rr_gauss(:,nb:ne),angles,N_r_pad,nb,ne,time)
    END IF

    IF (nb_field==2) THEN
       prod_ru(:,1,:) = eta_n*ru(:,1,:)
    ELSE IF (nb_field==6) THEN
       prod_ru(:,1,:) = eta_n*ru(:,1,:)
       prod_ru(:,2,:) = eta_n*ru(:,2,:)
       prod_ru(:,3,:) = eta_n*ru(:,3,:)
    ELSE
       CALL error_petsc('error in nb_field of c_in while calling FFT_PAR_VAR_ETA_PROD_GAUSS_DCL')
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!! PRODUCT with eta  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    howmany = bloc_size*nb_field/2

    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
         inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE)
    CALL dfftw_execute(fftw_plan_multi_r2c)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    !JLG, FEB 4, 2011
!!$       prod_cu = prod_cu/N_r !Scaling
    prod_cu = prod_cu*(1.d0/N_r_pad) !Scaling
    !JLG, FEB 4, 2011
    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

    !Now we need to redistribute the Fourier coefficients on each processor
    t = MPI_WTIME()
    combined_prod_cu(:,:,1)=prod_cu(1,:,:)
    DO n=2, m_max
       combined_prod_cu(:,:,n)=2*CONJG(prod_cu(n,:,:))
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    t = MPI_WTIME()
    longueur_tranche=bloc_size*m_max_c*nb_field
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu,longueur_tranche,MPID, dist_prod_cu,longueur_tranche, &
         MPID,communicator,code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          intermediate = dist_prod_cu(:,:,shiftl+i)
          DO n = 1, bloc_size
             IF (n+shiftc > np ) CYCLE
             DO i_field = 1, nb_field/2
                c_out(n+shiftc, 2*i_field-1, i) = REAL (intermediate(i_field,n),KIND=8)
                c_out(n+shiftc, 2*i_field, i)   = AIMAG(intermediate(i_field,n))
             END DO
          END DO
       END DO
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

  END SUBROUTINE FFT_PAR_VAR_ETA_PROD_GAUSS_DCL

  FUNCTION regul(phi,eps) RESULT(r)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: phi, eps
    REAL(KIND=8)             :: x, r
    x = phi - 0.5d0
    IF (x .LE. -eps) THEN
       r = 0.d0
    ELSE IF (x .LE. eps) THEN
       r = (1+ x*(x**2 - 3*eps**2)/(-2*eps**3))/2
    ELSE
       r = 1.d0
    END IF
  END FUNCTION regul

  FUNCTION regul_tab(phi,eps) RESULT(r)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phi
    REAL(KIND=8), INTENT(IN) :: eps
    REAL(KIND=8) :: x
    REAL(KIND=8), DIMENSION(SIZE(phi)) :: r
    INTEGER :: n
    DO n = 1, SIZE(phi)
       x = phi(n) - 0.5d0
       IF (x .LE. -eps) THEN
          r(n) = 0.d0
       ELSE IF (x .LE. eps) THEN
          r(n) = (1.d0 + x*(x**2 - 3*eps**2)/(-2*eps**3))/2
       ELSE
          r(n) = 1.d0
       END IF
    END DO
  END FUNCTION regul_tab

  SUBROUTINE FFT_CHECK_INTERFACE(communicator, V1_in, nb_fluids, interface_ok, &
       nb_procs, bloc_size, m_max_pad, temps)
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    USE my_util
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:,:),INTENT(IN)  :: V1_in
    INTEGER,                         INTENT(IN)  :: nb_fluids
    INTEGER,                         INTENT(OUT) :: interface_ok
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    INTEGER, INTENT(IN)                                 :: bloc_size, m_max_pad, nb_procs
    INTEGER          :: np, np_tot, nb_field, m_max, m_max_c, MPID, N_r_pad
    INTEGER(KIND=8)  :: fftw_plan_multi_c2r
    COMPLEX(KIND=8), DIMENSION(m_max_pad, bloc_size)                 :: cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1, bloc_size)                :: ru
    REAL(KIND=8), DIMENSION(SIZE(V1_in,4),SIZE(V1_in,3),bloc_size*nb_procs) :: dist_field, combined_field
    REAL(KIND=8), DIMENSION(nb_fluids-1,2*m_max_pad-1, bloc_size)          :: V1_real_loc
    INTEGER :: n1, n2
    INTEGER :: nb, shiftc, shiftl, jindex, longueur_tranche, i, n, code, interface_nb
    REAL(KIND=8) :: t
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator

    IF (PRESENT(temps)) temps = 0.d0

    np      = SIZE(V1_in,2)
    nb_field= SIZE(V1_in,3) ! Number of fields
    m_max_c = SIZE(V1_in,4) ! Number of complex (cosines + sines) coefficients per point
    m_max = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    N_r_pad=2*m_max_pad-1
    np_tot = nb_procs*bloc_size

    IF (MOD(nb_field,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_CHECK_INTERFACE '
       STOP
    END IF

    ! Bloc_size is the number of points that are handled by one processor
    ! once the Fourier modes are all collected
    ! Computation of bloc_size and np_tot
    ! fin de la repartition des points

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()

    DO interface_nb = 1, nb_fluids-1
       dist_field = 0.d0
       DO i = 1, m_max_c
          dist_field(i,1:nb_field,1:np) = TRANSPOSE(V1_in(interface_nb,:,:,i))
       END DO
       IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

       !IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100
       !TEST LC
       IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 0.d0

       longueur_tranche=bloc_size*m_max_c*nb_field

       t = MPI_WTIME()
       MPID=MPI_DOUBLE_PRECISION
       CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
            MPID, communicator, code)

       IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

       t = MPI_WTIME()
       !JLG, FEB 4, 2011
       cu = 0.d0
       !JLG, FEB 4, 2011
       DO n = 1, bloc_size
          DO nb = 1, nb_procs
             shiftc = (nb-1)*bloc_size
             shiftl = (nb-1)*m_max_c
             jindex = n + shiftc
             ! Put real and imaginary parts in a complex
             ! nf=1,2,3 => V1_in
             ! nf=4,5,6 => V2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,n) = 0.5d0*CMPLX(combined_field(:,1,jindex),&
                  -combined_field(:,2,jindex),KIND=8)
          END DO
       END DO
       cu(1,:) = 2*CMPLX(REAL(cu(1,:),KIND=8),0.d0,KIND=8)


       IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

       ! Set the parameters for dfftw
       fft_dim=1; istride=1; ostride=1;
       !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
       idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
       odist=m_max_pad; onembed(1)=m_max_pad
       !JLG, FEB 4, 2011

       howmany=bloc_size*nb_field/2

       t = MPI_WTIME()
       CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
            onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
       !write(*,*) ' FFT_PAR_CROSS_PROD: fftw_plan_multi_c2r', fftw_plan_multi_c2r
       CALL dfftw_execute(fftw_plan_multi_c2r)

       ! clean up
       CALL dfftw_destroy_plan(fftw_plan_multi_c2r)
       V1_real_loc(interface_nb,:,:)=ru
    END DO


    ! CHECK IF INTERFACES CROSS EACH OTHER
    ! no crossing > interface_ok=1 / crossing > interface_ok=0
    interface_ok = 1
    IF (nb_field==2) THEN
       IF (nb_fluids==2) THEN
          RETURN
       ELSE
          DO interface_nb = 1, nb_fluids-2
             DO n1 = 1, 2*m_max_pad-1
                DO n2 = 1, bloc_size
                   IF((0.1d0.LE.V1_real_loc(interface_nb,n1,n2)).AND. &
                        (V1_real_loc(interface_nb,n1,n2).LE.0.9d0)) THEN
                      IF(V1_real_loc(interface_nb,n1,n2).LT.V1_real_loc(interface_nb+1,n1,n2)) THEN
                         interface_ok = 0
                         RETURN
                      END IF
                   END IF
                END DO
             END DO
          END DO
       END IF
    ELSE
       CALL error_petsc('BUG in FFT_CHECK_INTERFACE: pb with V1_in dimensions')
    END IF
    ! CHECK IF INTERFACE CROSS EACH OTHER

  END SUBROUTINE FFT_CHECK_INTERFACE

  SUBROUTINE FFT_COMPUTE_ENTROPY_VISC(communicator, communicator_S, V1_in, V2_in, V3_in, V4_in, V5_in, &
       hloc_gauss, V1_out, V2_out, V3_out, &
       nb_procs, bloc_size, m_max_pad, residual_normalization, l_G, temps)
    !FFT (FFT(-1) V1 . FFT(-1) V2) = c_out
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    USE my_util
    USE input_data
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V1_in, V2_in, V3_in, V4_in, V5_in
    REAL(KIND=8), DIMENSION(:),      INTENT(IN)  :: hloc_gauss
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT) :: V1_out, V2_out, V3_out
    REAL(KIND=8),                    INTENT(IN)  :: residual_normalization
    INTEGER,                         INTENT(IN)  :: nb_procs, bloc_size, m_max_pad, l_G
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    COMPLEX(KIND=8), DIMENSION(m_max_pad, SIZE(V1_in,2)*5/2, bloc_size)  :: cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,SIZE(V1_in,2)*5/2,bloc_size)   :: ru
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size) :: visco_entro
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,SIZE(V1_in,2)/2,bloc_size) :: prod_ru_1, prod_ru_2, prod_ru_3
    COMPLEX(KIND=8), DIMENSION(m_max_pad,SIZE(V1_in,2)/2,bloc_size)  :: prod_cu_1, prod_cu_2, prod_cu_3
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size)            :: intermediate_1
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size)            :: intermediate_2
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size)            :: intermediate_3
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size) :: norm_vel
    REAL(KIND=8), DIMENSION(2*m_max_pad-1, bloc_size/l_G)                        :: norm_vel_int
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),5*SIZE(V1_in,2),bloc_size*nb_procs)    :: dist_field
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),5*SIZE(V1_in,2),bloc_size*nb_procs)    :: combined_field
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu_1
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu_2
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu_3
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: dist_prod_cu_1
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: dist_prod_cu_2
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: dist_prod_cu_3
    REAL(KIND=8), DIMENSION(bloc_size*nb_procs)                  :: hloc_gauss_tot
    INTEGER  :: np, np_tot, nb_field,  m_max, m_max_c, MPID,  N_r_pad
    INTEGER(KIND=8) :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    INTEGER ::   nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code, rank, l, i_field
    REAL(KIND=8) :: t, hh, x
    REAL(KIND=8) :: max_norm_vel_loc, max_norm_vel_loc_F, max_norm_vel_tot
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator
    MPI_Comm :: communicator_S
    CALL MPI_COMM_RANK(communicator, rank, code)

    IF (PRESENT(temps)) temps = 0.d0

    np      = SIZE(V1_in,1)
    nb_field= SIZE(V1_in,2) ! Number of fields
    m_max_c = SIZE(V1_in,3) ! Number of complex (cosines + sines) coefficients per point
    m_max = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    np_tot = nb_procs*bloc_size
    N_r_pad=2*m_max_pad-1

    IF (MOD(nb_field,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_COMPUTE_ENTROPY_VISC '
       STOP
    END IF

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()

    DO i = 1, m_max_c
       dist_field(i,           1:nb_field,  1:np) = TRANSPOSE(V1_in(:,:,i))
       dist_field(i,  nb_field+1:2*nb_field,1:np) = TRANSPOSE(V2_in(:,:,i))
       dist_field(i,2*nb_field+1:3*nb_field,1:np) = TRANSPOSE(V3_in(:,:,i))
       dist_field(i,3*nb_field+1:4*nb_field,1:np) = TRANSPOSE(V4_in(:,:,i))
       dist_field(i,4*nb_field+1:5*nb_field,1:np) = TRANSPOSE(V5_in(:,:,i))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    !IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100
    !LC 2016/02/16 (we do a max on vel later)
    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 0.d0
    !LC 2016/02/16

    hloc_gauss_tot(1:np) = hloc_gauss
    IF (np/=np_tot) hloc_gauss_tot(np+1:np_tot) = 1.d100

    longueur_tranche=bloc_size*m_max_c*nb_field*5

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, nb_field*5/2
             ! Put real and imaginary parts in a complex
             ! nf=1,2,3 => V1_in
             ! nf=4,5,6 => V2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*nb_field*5/2

    t = MPI_WTIME()
    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_DOT_PROD: fftw_plan_multi_c2r', fftw_plan_multi_c2r
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    !===ru(:,1:3,:)   is velocity
    !===ru(:,4:6,:)   is res_ns
    !===ru(:,7:9,:)   is first line of velocity gradient
    !===ru(:,10:12,:) is second line of velocity gradient
    !===ru(:,10:12,:) is third line of velocity gradient

    !===Compute local velocity
    norm_vel(:,:) = SQRT(ru(:,1,:)**2 + ru(:,2,:)**2 + ru(:,3,:)**2)
    DO i = 1, 2*m_max_pad - 1
       DO l = 1, bloc_size/l_G
          x = MAXVAL(norm_vel(i,(l-1)*l_G+1:l*l_G))
          norm_vel_int(i,l) = x
       END DO
    END DO
    IF (2*m_max_pad - 1 .GE. 3) THEN
       DO l = 1, bloc_size/l_G
          DO i = 2, 2*m_max_pad - 2
             norm_vel(i,(l-1)*l_G+1:l*l_G) = MAXVAL(norm_vel_int(i-1:i+1,l))
          END DO
          norm_vel(1,(l-1)*l_G+1:l*l_G) = MAX(norm_vel_int(1,l),norm_vel_int(2,l),norm_vel_int(2*m_max_pad - 1,l))
          norm_vel(2*m_max_pad - 1,(l-1)*l_G+1:l*l_G) = &
               MAX(norm_vel_int(2*m_max_pad - 2,l),norm_vel_int(2*m_max_pad - 1,l),norm_vel_int(1,l))
       END DO
    ELSE
       DO l = 1, bloc_size/l_G
          norm_vel(1,(l-1)*l_G+1:l*l_G) = norm_vel_int(1,l)
       END DO
    END IF
    !===End compute local velocity

    !==Compute maximum of velocity
    max_norm_vel_loc=MAXVAL(norm_vel)
    CALL MPI_ALLREDUCE(max_norm_vel_loc,max_norm_vel_loc_F,1,MPI_DOUBLE_PRECISION, MPI_MAX, communicator, code)
    CALL MPI_ALLREDUCE(max_norm_vel_loc_F,max_norm_vel_tot,1,MPI_DOUBLE_PRECISION, MPI_MAX, communicator_S, code)
    !==End compute maximum of velocity

    !===Compute entropy viscosity
    DO n = 1, bloc_size
       hh = hloc_gauss_tot(n+rank*bloc_size)
       visco_entro(:,n) = -inputs%LES_coeff1 + MIN(inputs%LES_coeff4*norm_vel(:,n), &
            inputs%LES_coeff2*hh*ABS(ru(:,1,n)*ru(:,4,n) + ru(:,2,n)*ru(:,5,n) + ru(:,3,n)*ru(:,6,n))&
            /(residual_normalization+1.d-14))
    END DO
    !===End compute entropy viscosity

    !===Compute visco_entro times gradient(velocity)
    prod_ru_1(:,1,:) = visco_entro*ru(:,7,:)
    prod_ru_1(:,2,:) = visco_entro*ru(:,8,:)
    prod_ru_1(:,3,:) = visco_entro*ru(:,9,:)

    prod_ru_2(:,1,:) = visco_entro*ru(:,10,:)
    prod_ru_2(:,2,:) = visco_entro*ru(:,11,:)
    prod_ru_2(:,3,:) = visco_entro*ru(:,12,:)

    prod_ru_3(:,1,:) = visco_entro*ru(:,13,:)
    prod_ru_3(:,2,:) = visco_entro*ru(:,14,:)
    prod_ru_3(:,3,:) = visco_entro*ru(:,15,:)
    !===End compute visco_entro times gradient(velocity)

    howmany = bloc_size*nb_field/2
    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru_1, &
         inembed, istride, idist, prod_cu_1, onembed, ostride, odist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_DOT_PROD: fftw_plan_multi_r2c', fftw_plan_multi_r2c
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru_2, &
         inembed, istride, idist, prod_cu_2, onembed, ostride, odist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_DOT_PROD: fftw_plan_multi_r2c', fftw_plan_multi_r2c
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru_3, &
         inembed, istride, idist, prod_cu_3, onembed, ostride, odist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_DOT_PROD: fftw_plan_multi_r2c', fftw_plan_multi_r2c
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    prod_cu_1 = prod_cu_1*(1.d0/N_r_pad) !Scaling
    prod_cu_2 = prod_cu_2*(1.d0/N_r_pad) !Scaling
    prod_cu_3 = prod_cu_3*(1.d0/N_r_pad) !Scaling

    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

    !Now we need to redistribute the Fourier coefficients on each processor
    t = MPI_WTIME()
    combined_prod_cu_1(:,:,1)=prod_cu_1(1,:,:)
    combined_prod_cu_2(:,:,1)=prod_cu_2(1,:,:)
    combined_prod_cu_3(:,:,1)=prod_cu_3(1,:,:)

    DO n=2, m_max
       combined_prod_cu_1(:,:,n)=2*CONJG(prod_cu_1(n,:,:))
       combined_prod_cu_2(:,:,n)=2*CONJG(prod_cu_2(n,:,:))
       combined_prod_cu_3(:,:,n)=2*CONJG(prod_cu_3(n,:,:))
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    t = MPI_WTIME()
    longueur_tranche=bloc_size*m_max_c*nb_field
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu_1,longueur_tranche,MPID, dist_prod_cu_1,longueur_tranche, &
         MPID,communicator,code)
    CALL MPI_ALLTOALL (combined_prod_cu_2,longueur_tranche,MPID, dist_prod_cu_2,longueur_tranche, &
         MPID,communicator,code)
    CALL MPI_ALLTOALL (combined_prod_cu_3,longueur_tranche,MPID, dist_prod_cu_3,longueur_tranche, &
         MPID,communicator,code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          intermediate_1 = dist_prod_cu_1(:,:,shiftl+i)
          intermediate_2 = dist_prod_cu_2(:,:,shiftl+i)
          intermediate_3 = dist_prod_cu_3(:,:,shiftl+i)
          DO n = 1, bloc_size
             IF (n+shiftc > np ) CYCLE
             DO i_field = 1, nb_field/2
                V1_out(n+shiftc, i_field*2-1, i) = REAL (intermediate_1(i_field,n),KIND=8)
                V1_out(n+shiftc, i_field*2 , i)  = AIMAG(intermediate_1(i_field,n))
                V2_out(n+shiftc, i_field*2-1, i) = REAL (intermediate_2(i_field,n),KIND=8)
                V2_out(n+shiftc, i_field*2 , i)  = AIMAG(intermediate_2(i_field,n))
                V3_out(n+shiftc, i_field*2-1, i) = REAL (intermediate_3(i_field,n),KIND=8)
                V3_out(n+shiftc, i_field*2 , i)  = AIMAG(intermediate_3(i_field,n))
             END DO
          END DO
       END DO
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

  END SUBROUTINE FFT_COMPUTE_ENTROPY_VISC

  SUBROUTINE FFT_COMPUTE_ENTROPY_VISC_MOM(communicator,communicator_S, V1_in, V2_in, V3_in, c1_in, hloc_gauss, &
       c1_real_out, nb_procs, bloc_size, m_max_pad, l_G, opt_c2_real_out, temps)
    !FFT (FFT(-1) V1 . FFT(-1) V2) = c_out
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    USE my_util
    USE input_data
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V1_in, V2_in, V3_in
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: c1_in
    REAL(KIND=8), DIMENSION(:),      INTENT(IN)  :: hloc_gauss
    REAL(KIND=8), DIMENSION(:,:),    INTENT(OUT) :: c1_real_out
    INTEGER,                         INTENT(IN)  :: nb_procs, bloc_size, m_max_pad, l_G
    REAL(KIND=8), DIMENSION(:,:), OPTIONAL,  INTENT(OUT)   :: opt_c2_real_out
    REAL(KIND=8), DIMENSION(:),   OPTIONAL,  INTENT(INOUT) :: temps
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),3*SIZE(V1_in,2)+SIZE(c1_in,2),bloc_size*nb_procs) :: dist_field
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),3*SIZE(V1_in,2)+SIZE(c1_in,2),bloc_size*nb_procs) :: combined_field
    COMPLEX(KIND=8), DIMENSION(m_max_pad, (3*SIZE(V1_in,2)+SIZE(c1_in,2))/2, bloc_size)     :: cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,(3*SIZE(V1_in,2)+SIZE(c1_in,2))/2,bloc_size)      :: ru
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size)             :: norm_vel, norm_mom
    REAL(KIND=8), DIMENSION(2*m_max_pad-1, bloc_size/l_G)        :: norm_vel_int, norm_mom_int
    REAL(KIND=8), DIMENSION(bloc_size*nb_procs)                  :: hloc_gauss_tot
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size)             :: visco_entro
    INTEGER  :: np, np_tot, nb_field, nb_field2, m_max, m_max_c, MPID,  N_r_pad
    INTEGER(KIND=8) :: fftw_plan_multi_c2r
    INTEGER ::   nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code, rank, l
    !   REAL(KIND=8) :: t, hh, x, max_vel_loc, max_vel_F, max_vel_tot
    ! CN 21/01/2018
    REAL(KIND=8) :: t, hh, x, max_vel_loc, max_vel_F, max_vel_tot, max_mom_loc, max_mom_F, max_mom_tot, max_v_m
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator
    MPI_Comm :: communicator_S

    CALL MPI_COMM_RANK(communicator, rank, code)

    IF (PRESENT(temps)) temps = 0.d0

    np        = SIZE(V1_in,1)
    nb_field  = SIZE(V1_in,2) ! Number of fields for vector
    nb_field2 = SIZE(c1_in,2) ! Number of fields for scalar
    m_max_c   = SIZE(V1_in,3) ! Number of complex (cosines + sines) coefficients per point
    m_max     = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    np_tot    = nb_procs*bloc_size
    N_r_pad   = 2*m_max_pad-1

    IF (MOD(nb_field,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_COMPUTE_ENTROPY_VISC_MOM '
       STOP
    END IF

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()

    DO i = 1, m_max_c
       dist_field(i,           1:nb_field,  1:np) = TRANSPOSE(V1_in(:,:,i))
       dist_field(i,  nb_field+1:2*nb_field,1:np) = TRANSPOSE(V2_in(:,:,i))
       dist_field(i,2*nb_field+1:3*nb_field,1:np) = TRANSPOSE(V3_in(:,:,i))
       dist_field(i,3*nb_field+1:3*nb_field+nb_field2,1:np) = TRANSPOSE(c1_in(:,:,i))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 0.d0

    hloc_gauss_tot(1:np) = hloc_gauss
    IF (np/=np_tot) hloc_gauss_tot(np+1:np_tot) = 1.d100

    longueur_tranche=bloc_size*m_max_c*(nb_field*3+nb_field2)

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, (3*nb_field+nb_field2)/2
             ! Put real and imaginary parts in a complex
             ! nf=1,2,3 => V1_in
             ! nf=4,5,6 => V2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*(nb_field*3+nb_field2)/2

    t = MPI_WTIME()
    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_DOT_PROD: fftw_plan_multi_c2r', fftw_plan_multi_c2r
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    !===ru(:,1:3,:)   is velocity
    !===ru(:,4:6,:)   is momentum
    !===ru(:,7:9,:)   is res_ns
    !===ru(:,10,:)    is res_mass

    !===Compute norm_vel and norm_mom
    norm_vel(:,:) = SQRT(ru(:,1,:)**2 + ru(:,2,:)**2 + ru(:,3,:)**2)
    norm_mom(:,:) = SQRT(ru(:,4,:)**2 + ru(:,5,:)**2 + ru(:,6,:)**2)
    DO i = 1, 2*m_max_pad - 1
       DO l = 1, bloc_size/l_G
          x = MAXVAL(norm_vel(i,(l-1)*l_G+1:l*l_G))
          norm_vel_int(i,l) = x
          x = MAXVAL(norm_mom(i,(l-1)*l_G+1:l*l_G))
          norm_mom_int(i,l) = x
       END DO
    END DO
    IF (2*m_max_pad - 1 .GE. 3) THEN
       DO l = 1, bloc_size/l_G
          DO i = 2, 2*m_max_pad - 2
             norm_vel(i,(l-1)*l_G+1:l*l_G) = MAXVAL(norm_vel_int(i-1:i+1,l))
             norm_mom(i,(l-1)*l_G+1:l*l_G) = MAXVAL(norm_mom_int(i-1:i+1,l))
          END DO
          norm_vel(1,(l-1)*l_G+1:l*l_G) = MAX(norm_vel_int(1,l),norm_vel_int(2,l),norm_vel_int(2*m_max_pad - 1,l))
          norm_vel(2*m_max_pad - 1,(l-1)*l_G+1:l*l_G) = &
               MAX(norm_vel_int(2*m_max_pad - 2,l),norm_vel_int(2*m_max_pad - 1,l),norm_vel_int(1,l))
          norm_mom(1,(l-1)*l_G+1:l*l_G) = MAX(norm_mom_int(1,l),norm_mom_int(2,l),norm_mom_int(2*m_max_pad - 1,l))
          norm_mom(2*m_max_pad - 1,(l-1)*l_G+1:l*l_G) = &
               MAX(norm_mom_int(2*m_max_pad - 2,l),norm_mom_int(2*m_max_pad - 1,l),norm_mom_int(1,l))
       END DO
    ELSE
       DO l = 1, bloc_size/l_G
          norm_vel(1,(l-1)*l_G+1:l*l_G) = norm_vel_int(1,l)
          norm_mom(1,(l-1)*l_G+1:l*l_G) = norm_mom_int(1,l)
       END DO
    END IF
    !===End compute norm_vel and norm_mom

    !===Compute MAX(norm_vel) on domain
    max_vel_loc=MAXVAL(norm_vel)
    CALL MPI_ALLREDUCE(max_vel_loc,max_vel_F, 1, MPI_DOUBLE_PRECISION, &
         MPI_MAX, communicator, code)
    CALL MPI_ALLREDUCE(max_vel_F, max_vel_tot, 1, MPI_DOUBLE_PRECISION, &
         MPI_MAX, communicator_S, code)
    !   !===End compute MAX(norm_vel) on domain
    ! CN 21/01/2018
    max_mom_loc=MAXVAL(norm_mom)
    CALL MPI_ALLREDUCE(max_mom_loc,max_mom_F, 1, MPI_DOUBLE_PRECISION, &
         MPI_MAX, communicator, code)
    CALL MPI_ALLREDUCE(max_mom_F, max_mom_tot, 1, MPI_DOUBLE_PRECISION, &
         MPI_MAX, communicator_S, code)
    max_v_m = max_vel_tot*max_mom_tot
    !===End compute MAX(norm_vel), MAX(norm_mom) on domain

    !===Compute entropy viscosities for momentum and level set equations
    DO n = 1, bloc_size
       hh = hloc_gauss_tot(n+rank*bloc_size)

       visco_entro(:,n) = inputs%LES_coeff2*hh* &
            MAX(ABS(ru(:,1,n)*ru(:,7,n) + ru(:,2,n)*ru(:,8,n) + ru(:,3,n)*ru(:,9,n)), &
            ABS(ru(:,10,n)*(ru(:,1,n)**2 + ru(:,2,n)**2 + ru(:,3,n)**2)))
    END DO

    !JLG Sept 21 2016
    !===Average entropy viscosity
    IF (2*m_max_pad - 1 .GE. 3) THEN
       DO l = 1, bloc_size/l_G
          DO i = 2, 2*m_max_pad - 2
             norm_vel_int(i,l) = MAXVAL(visco_entro(i-1:i+1,(l-1)*l_G+1:l*l_G))
          END DO
          norm_vel_int(1,l) =  MAX(MAXVAL(visco_entro(1,(l-1)*l_G+1:l*l_G)),MAXVAL(visco_entro(2,(l-1)*l_G+1:l*l_G)),&
               MAXVAL(visco_entro(2*m_max_pad - 1,(l-1)*l_G+1:l*l_G)))
          norm_vel_int(2*m_max_pad - 1,l) =   MAX(MAXVAL(visco_entro(2*m_max_pad - 2,(l-1)*l_G+1:l*l_G)),&
               MAXVAL(visco_entro(2*m_max_pad - 1,(l-1)*l_G+1:l*l_G)),MAXVAL(visco_entro(1,(l-1)*l_G+1:l*l_G)))
       END DO
    ELSE
       DO l = 1, bloc_size/l_G
          norm_vel_int(1,l) =  MAXVAL(visco_entro(1,(l-1)*l_G+1:l*l_G))
       END DO
    END IF
    DO i = 1, 2*m_max_pad - 1
       DO l = 1, bloc_size/l_G
          visco_entro(i,(l-1)*l_G+1:l*l_G)=norm_vel_int(i,l)
       END DO
    END DO
    !===End Average entropy viscosity
    !End JLG Sept 21 2016

    DO n = 1, bloc_size
       IF (l_G==3) THEN !===Level set viscosity
          c1_real_out(:,n) = MIN(4*inputs%LES_coeff4*max_vel_tot, &
               16*visco_entro(:,n)/(norm_vel(:,n)*norm_mom(:,n)+1.d-30))
       ELSE !===Momentum viscosity
          c1_real_out(:,n) = MIN(inputs%LES_coeff4*norm_vel(:,n), &
               visco_entro(:,n)/(norm_vel(:,n)*norm_mom(:,n)+1.d-30))
          IF (PRESENT(opt_c2_real_out)) THEN
             opt_c2_real_out(:,n) = MIN(inputs%LES_coeff4*max_vel_tot, &
                  visco_entro(:,n)/(norm_vel(:,n)*norm_mom(:,n)+1.d-30))
          END IF
       END IF
    END DO
    !===End compute entropy viscosities for momentum and level set equations

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

  END SUBROUTINE FFT_COMPUTE_ENTROPY_VISC_MOM

  SUBROUTINE FFT_COMPUTE_DIFFU_MOM(communicator, V1_in, V2_in, V3_in, &
       V1_out, V2_out, nb_procs, bloc_size, m_max_pad, temps)
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    USE my_util
    USE input_data
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),       INTENT(IN)    :: V1_in, V2_in !vectors
    REAL(KIND=8), DIMENSION(:,:,:),       INTENT(IN)    :: V3_in !scalar
    REAL(KIND=8), DIMENSION(:,:,:),       INTENT(OUT)   :: V1_out
    REAL(KIND=8), DIMENSION(:,:,:),       INTENT(OUT)   :: V2_out
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    INTEGER,                              INTENT(IN)    :: bloc_size, m_max_pad, nb_procs
    INTEGER          :: np, np_tot, nb_field1, nb_field2, nb_field3
    INTEGER          :: m_max, m_max_c, MPID, N_r_pad
    INTEGER(KIND=8)  :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),2*SIZE(V1_in,2)+SIZE(V3_in,2), bloc_size*nb_procs) :: dist_field
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),2*SIZE(V1_in,2)+SIZE(V3_in,2), bloc_size*nb_procs) :: combined_field
    COMPLEX(KIND=8), DIMENSION(m_max_pad, (2*SIZE(V1_in,2)+SIZE(V3_in,2))/2, bloc_size)  :: cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,(2*SIZE(V1_in,2)+SIZE(V3_in,2))/2,bloc_size)   :: ru
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,SIZE(V1_in,2)/2,bloc_size) :: prod_ru_1, prod_ru_2
    COMPLEX(KIND=8), DIMENSION(m_max_pad,SIZE(V1_in,2)/2,bloc_size)  :: prod_cu_1, prod_cu_2
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size)            :: intermediate_1, intermediate_2
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu_1
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu_2
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: dist_prod_cu_1
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: dist_prod_cu_2
    INTEGER :: i_field, rank
    INTEGER :: nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code
    REAL(KIND=8) :: t
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator
    PetscErrorCode         :: ierr

    IF (PRESENT(temps)) temps = 0.d0

    CALL MPI_Comm_rank(communicator, rank, ierr)

    np       = SIZE(V1_in,1)
    nb_field1= SIZE(V1_in,2)   ! Number of fields for vector
    nb_field2= SIZE(V2_in,2)   ! Number of fields
    nb_field3= SIZE(V3_in,2)   ! Number of fields
    m_max_c  = SIZE(V1_in,3)   ! Number of complex (cosines + sines) coefficients per point
    m_max    = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    N_r_pad  = 2*m_max_pad-1
    np_tot   = nb_procs*bloc_size

    IF (MOD(nb_field1,2)/=0 .OR. MOD(nb_field2,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_COMPUTE_DIFFU_MOM '
       STOP
    END IF

    ! Bloc_size is the number of points that are handled by one processor
    ! once the Fourier modes are all collected
    ! Computation of bloc_size and np_tot
    ! fin de la repartition des points

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()

    DO i = 1, m_max_c
       dist_field(i,1:nb_field1,1:np)                         = TRANSPOSE(V1_in(:,:,i))
       dist_field(i,nb_field1+1:2*nb_field1,1:np)             = TRANSPOSE(V2_in(:,:,i))
       dist_field(i,2*nb_field1+1:2*nb_field1+nb_field3,1:np) = TRANSPOSE(V3_in(:,:,i))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 0.d0

    longueur_tranche=bloc_size*m_max_c*(nb_field1+nb_field2+nb_field3)

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, (nb_field1+nb_field2+nb_field3)/2
             ! Put real and imaginary parts in a complex
             ! nf=1,2,3 => V1_in
             ! nf=4 => V2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
    !idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
    !odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*(nb_field1+nb_field2+nb_field3)/2

    t = MPI_WTIME()
    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_CROSS_PROD: fftw_plan_multi_c2r', fftw_plan_multi_c2r
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    !===ru(:,1:3,:)   is first line of strain rate tensor
    !===ru(:,4:6,:)   is 2nd line/column 2 and 3 of strain rate tensor + 3rd line/3rd column
    !===ru(:,7,:)    is dynamical viscosity

    IF  (nb_field1==6 .AND. nb_field2==6 .AND. nb_field3==2) THEN
       !===Compute visc_dyn times strain rate tensor
       prod_ru_1(:,1,:) = ru(:,7,:)*ru(:,1,:)
       prod_ru_1(:,2,:) = ru(:,7,:)*ru(:,2,:)
       prod_ru_1(:,3,:) = ru(:,7,:)*ru(:,3,:)

       prod_ru_2(:,1,:) = ru(:,7,:)*ru(:,4,:)
       prod_ru_2(:,2,:) = ru(:,7,:)*ru(:,5,:)
       prod_ru_2(:,3,:) = ru(:,7,:)*ru(:,6,:)
       !===End Compute visc_dyn times strain rate tensor
    ELSE
       CALL error_petsc('error on inputs while calling FFT_COMPUTE_DIFFU_MOM')
    END IF

    howmany = bloc_size*nb_field1/2 !vector
    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru_1, &
         inembed, istride, idist, prod_cu_1, onembed, ostride, odist, FFTW_ESTIMATE)
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    howmany = bloc_size*nb_field1/2 !vector
    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru_2, &
         inembed, istride, idist, prod_cu_2, onembed, ostride, odist, FFTW_ESTIMATE)
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    !JLG, FEB 4, 2011
    !prod_cu = prod_cu/N_r !Scaling
    prod_cu_1 = prod_cu_1*(1.d0/N_r_pad) !Scaling
    prod_cu_2 = prod_cu_2*(1.d0/N_r_pad) !Scaling
    !JLG, FEB 4, 2011
    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

    !Now we need to redistribute the Fourier coefficients on each processor

    t = MPI_WTIME()
    combined_prod_cu_1(:,:,1)=prod_cu_1(1,:,:)
    combined_prod_cu_2(:,:,1)=prod_cu_2(1,:,:)
    DO n=2, m_max
       combined_prod_cu_1(:,:,n)=2*CONJG(prod_cu_1(n,:,:))
       combined_prod_cu_2(:,:,n)=2*CONJG(prod_cu_2(n,:,:))
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    t = MPI_WTIME()
    longueur_tranche=bloc_size*m_max_c*nb_field1 !vectors
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu_1,longueur_tranche,MPID, dist_prod_cu_1,longueur_tranche, &
         MPID,communicator,code)
    CALL MPI_ALLTOALL (combined_prod_cu_2,longueur_tranche,MPID, dist_prod_cu_2,longueur_tranche, &
         MPID,communicator,code)

    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          intermediate_1 = dist_prod_cu_1(:,:,shiftl+i)
          intermediate_2 = dist_prod_cu_2(:,:,shiftl+i)
          DO n = 1, bloc_size
             IF (n+shiftc > np ) CYCLE
             DO i_field = 1, nb_field1/2
                V1_out(n+shiftc, i_field*2-1, i) = REAL (intermediate_1(i_field,n),KIND=8)
                V1_out(n+shiftc, i_field*2  , i) = AIMAG(intermediate_1(i_field,n))
                V2_out(n+shiftc, i_field*2-1, i) = REAL (intermediate_2(i_field,n),KIND=8)
                V2_out(n+shiftc, i_field*2  , i) = AIMAG(intermediate_2(i_field,n))
             END DO
          END DO
       END DO
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

  END SUBROUTINE FFT_COMPUTE_DIFFU_MOM

  SUBROUTINE FFT_PAR_COMPR_ENTRO_VISC_DCL(communicator, V1_in, V2_in, c1_in, c_in_real, hloc_gauss, &
       coeff1_in_level, V_out, c_out, nb_procs, bloc_size, m_max_pad, temps)
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    USE my_util
    USE input_data
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)         :: V1_in, V2_in !vector
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)         :: c1_in  !scalar
    REAL(KIND=8), DIMENSION(:,:),    INTENT(IN)         :: c_in_real
    REAL(KIND=8), DIMENSION(:),      INTENT(IN)         :: hloc_gauss
    REAL(KIND=8),                    INTENT(IN)         :: coeff1_in_level
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT)        :: V_out
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT)        :: c_out
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    INTEGER, INTENT(IN)                                 :: bloc_size, m_max_pad, nb_procs
    INTEGER          :: np, np_tot, nb_field1, nb_field2, nb_field3
    INTEGER          :: m_max, m_max_c, MPID, N_r_pad
    INTEGER(KIND=8)  :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),2*SIZE(V1_in,2)+SIZE(c1_in,2),bloc_size*nb_procs) :: dist_field
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),2*SIZE(V1_in,2)+SIZE(c1_in,2),bloc_size*nb_procs) :: combined_field
    COMPLEX(KIND=8), DIMENSION(m_max_pad, (2*SIZE(V1_in,2)+SIZE(c1_in,2))/2, bloc_size)  :: cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,(2*SIZE(V1_in,2)+SIZE(c1_in,2))/2,bloc_size)   :: ru
    REAL(KIND=8), DIMENSION(2*m_max_pad-1, bloc_size)                :: norm_grad_phi
    REAL(KIND=8), DIMENSION(bloc_size*nb_procs)                      :: hloc_gauss_tot
    REAL(KIND=8), DIMENSION(2*m_max_pad-1, bloc_size)                :: visc_comp, visc_entro
    COMPLEX(KIND=8), DIMENSION(m_max_pad,SIZE(V1_in,2)/2,bloc_size)  :: prod_cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,SIZE(V1_in,2)/2,bloc_size) :: prod_ru
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2, bloc_size)           :: intermediate
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: dist_prod_cu
    COMPLEX(KIND=8), DIMENSION(m_max_pad,bloc_size)              :: prod_cu2
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size)             :: prod_ru2
    COMPLEX(KIND=8), DIMENSION(bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu2
    COMPLEX(KIND=8), DIMENSION(bloc_size,SIZE(V1_in,3)*nb_procs) :: dist_prod_cu2
    COMPLEX(KIND=8), DIMENSION(bloc_size)                        :: intermediate2
    INTEGER :: i_field
    INTEGER :: nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code
    REAL(KIND=8) :: t
    INTEGER :: rank
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator

    CALL MPI_COMM_RANK(communicator, rank, code)

    IF (PRESENT(temps)) temps = 0.d0

    np       = SIZE(V1_in,1)
    nb_field1= SIZE(V1_in,2) ! Number of fields
    nb_field2= SIZE(V2_in,2) ! Number of fields
    nb_field3= SIZE(c1_in,2) ! Number of fields
    m_max_c  = SIZE(V1_in,3) ! Number of complex (cosines + sines) coefficients per point
    m_max    = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    N_r_pad  = 2*m_max_pad-1
    np_tot   = nb_procs*bloc_size

    IF (MOD(nb_field1,2)/=0 .OR. MOD(nb_field2,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_PAR_COMPR_ENTRO_VISC_DCL '
       STOP
    END IF

    ! Bloc_size is the number of points that are handled by one processor
    ! once the Fourier modes are all collected
    ! Computation of bloc_size and np_tot
    ! fin de la repartition des points

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()

    DO i = 1, m_max_c
       dist_field(i,1:nb_field1,1:np) = TRANSPOSE(V1_in(:,:,i))
       dist_field(i,nb_field1+1:2*nb_field1,1:np) = TRANSPOSE(V2_in(:,:,i))
       dist_field(i,2*nb_field1+1:2*nb_field1+nb_field3,1:np) = TRANSPOSE(c1_in(:,:,i))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100

    hloc_gauss_tot(1:np) = hloc_gauss
    IF (np/=np_tot) hloc_gauss_tot(np+1:np_tot) = 1.d100

    longueur_tranche=bloc_size*m_max_c*(nb_field1+nb_field2+nb_field3)

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, (nb_field1+nb_field2+nb_field3)/2
             ! Put real and imaginary parts in a complex
             ! nf=1,2,3 => V1_in
             ! nf=4 => V2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
    !idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
    !odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*(nb_field1+nb_field2+nb_field3)/2

    t = MPI_WTIME()
    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_CROSS_PROD: fftw_plan_multi_c2r', fftw_plan_multi_c2r
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    IF (nb_field1 == 6 .AND. nb_field2 == 6 .AND. nb_field3 == 2) THEN
       !===ru(:,1:3,:) = h*Grad(phi)
       !===ru(:,4:6,:) = h*Grad(phi_reg)
       !===ru(:,7,:)   = phi
       !===c_in_real   = visc_entro(:,:)

       !===Compute MAX(0,phi*(1-phi))
       prod_ru2(:,:) = MAX(0.d0, ru(:,7,:)*(1.d0 - ru(:,7,:)))

       !===Compute |Grad(phi_reg)|
       norm_grad_phi(:,:) = SQRT(ru(:,4,:)**2 + ru(:,5,:)**2 + ru(:,6,:)**2) + 1.d-14

       !===Compute visc_comp
       DO n = 1, bloc_size
          visc_entro(:,n) = -coeff1_in_level + c_in_real(:,n)
          visc_comp(:,n)= c_in_real(:,n)*inputs%level_set_comp_factor*prod_ru2(:,n)/norm_grad_phi(:,n)
       END DO

       !===Compute (coeff1_in_level-visc_entro)*Grad(phi) + visc_comp*Grad(phi_reg)
       prod_ru(:,1,:) = -visc_entro*ru(:,1,:) + visc_comp*ru(:,4,:)
       prod_ru(:,2,:) = -visc_entro*ru(:,2,:) + visc_comp*ru(:,5,:)
       prod_ru(:,3,:) = -visc_entro*ru(:,3,:) + visc_comp*ru(:,6,:)
       !TEST COUPEZ
       !       !===Compute MAX(0,phi*(1-phi))
       !       prod_ru2(:,:) = MAX(0.d0, ru(:,7,:)*(1.d0 - ru(:,7,:)))
       !
       !       !===Compute visc_comp
       !       DO n = 1, bloc_size
       !          visc_entro(:,n) = -coeff1_in_level + c_in_real(:,n)
       !       END DO
       !
       !       !===Compute (coeff1_in_level-visc_entro)*(h*Grad(phi))
       !       prod_ru(:,1,:) = -visc_entro*ru(:,1,:)
       !       prod_ru(:,2,:) = -visc_entro*ru(:,2,:)
       !       prod_ru(:,3,:) = -visc_entro*ru(:,3,:)
       !TEST COUPEZ
    ELSE
       CALL error_petsc('error in problem type while calling FFT_PAR_COMPR_ENTRO_VISC_DCL')
    END IF

    howmany = bloc_size*nb_field1/2
    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
         inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE)
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    howmany = bloc_size*1
    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru2, &
         inembed, istride, idist, prod_cu2, onembed, ostride, odist, FFTW_ESTIMATE)
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    !JLG, FEB 4, 2011
    !prod_cu = prod_cu/N_r !Scaling
    prod_cu = prod_cu*(1.d0/N_r_pad) !Scaling
    prod_cu2 = prod_cu2*(1.d0/N_r_pad) !Scaling
    !JLG, FEB 4, 2011
    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

    !Now we need to redistribute the Fourier coefficients on each processor

    t = MPI_WTIME()
    combined_prod_cu(:,:,1)=prod_cu(1,:,:)
    combined_prod_cu2(:,1)=prod_cu2(1,:)
    DO n=2, m_max
       combined_prod_cu(:,:,n)=2*CONJG(prod_cu(n,:,:))
       combined_prod_cu2(:,n)=2*CONJG(prod_cu2(n,:))
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    t = MPI_WTIME()
    longueur_tranche=bloc_size*m_max_c*nb_field1 !vector
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu,longueur_tranche,MPID, dist_prod_cu,longueur_tranche, &
         MPID,communicator,code)

    longueur_tranche=bloc_size*m_max_c*2 !scalar
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu2,longueur_tranche,MPID, dist_prod_cu2,longueur_tranche, &
         MPID,communicator,code)

    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t
    ! dimensions:
    t = MPI_WTIME()

    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          intermediate = dist_prod_cu(:,:,shiftl+i)
          intermediate2 = dist_prod_cu2(:,shiftl+i)
          DO n = 1, bloc_size
             IF (n+shiftc > np ) CYCLE
             DO i_field = 1, nb_field1/2
                V_out(n+shiftc, i_field*2-1, i) = REAL (intermediate(i_field,n),KIND=8)
                V_out(n+shiftc, i_field*2  , i) = AIMAG(intermediate(i_field,n))
             END DO
             c_out(n+shiftc, 1, i) = REAL (intermediate2(n),KIND=8)
             c_out(n+shiftc, 2, i) = AIMAG(intermediate2(n))
          END DO
       END DO
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

  END SUBROUTINE FFT_PAR_COMPR_ENTRO_VISC_DCL

  SUBROUTINE FFT_COMPUTE_ENTROPY_VISC_GRAD_MOM(communicator, V1_in, V2_in, V3_in, c_in_real, &
       V_out, nb_procs, bloc_size, m_max_pad, temps)
    !FFT (FFT(-1) V1 . FFT(-1) V2) = c_out
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    USE my_util
    USE input_data
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V1_in, V2_in, V3_in !Gradient of momentum
    REAL(KIND=8), DIMENSION(:,:),    INTENT(IN)  :: c_in_real !entropy viscosity
    REAL(KIND=8), DIMENSION(:,:,:,:),INTENT(OUT) :: V_out
    INTEGER,                         INTENT(IN)  :: nb_procs, bloc_size, m_max_pad
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    COMPLEX(KIND=8), DIMENSION(m_max_pad, SIZE(V1_in,2)*3/2, bloc_size)  :: cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,SIZE(V1_in,2)*3/2,bloc_size)   :: ru
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size)                 :: visco_entro
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,SIZE(V1_in,2)/2,bloc_size) :: prod_ru_1, prod_ru_2, prod_ru_3
    COMPLEX(KIND=8), DIMENSION(m_max_pad,SIZE(V1_in,2)/2,bloc_size)  :: prod_cu_1, prod_cu_2, prod_cu_3

    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size)            :: intermediate_1, intermediate_2, intermediate_3
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),3*SIZE(V1_in,2),bloc_size*nb_procs)    :: dist_field, combined_field
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu_1
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu_2
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu_3
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: dist_prod_cu_1
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: dist_prod_cu_2
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: dist_prod_cu_3
    INTEGER  :: np, np_tot, nb_field,  m_max, m_max_c, MPID,  N_r_pad
    INTEGER(KIND=8) :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    INTEGER ::   nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code, rank, i_field
    REAL(KIND=8) :: t
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator
    CALL MPI_COMM_RANK(communicator, rank, code)

    IF (PRESENT(temps)) temps = 0.d0

    np      = SIZE(V1_in,1)
    nb_field= SIZE(V1_in,2) ! Number of fields
    m_max_c = SIZE(V1_in,3) ! Number of complex (cosines + sines) coefficients per point
    m_max = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    np_tot = nb_procs*bloc_size
    N_r_pad=2*m_max_pad-1

    IF (MOD(nb_field,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_COMPUTE_ENTROPY_VISC_GRAD_MOM '
       STOP
    END IF

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()

    DO i = 1, m_max_c
       dist_field(i,           1:nb_field,  1:np) = TRANSPOSE(V1_in(:,:,i))
       dist_field(i,  nb_field+1:2*nb_field,1:np) = TRANSPOSE(V2_in(:,:,i))
       dist_field(i,2*nb_field+1:3*nb_field,1:np) = TRANSPOSE(V3_in(:,:,i))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100

    longueur_tranche=bloc_size*m_max_c*nb_field*3

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, nb_field*3/2
             ! Put real and imaginary parts in a complex
             ! nf=1,2,3 => V1_in
             ! nf=4,5,6 => V2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*nb_field*3/2

    t = MPI_WTIME()
    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_DOT_PROD: fftw_plan_multi_c2r', fftw_plan_multi_c2r
    CALL dfftw_execute(fftw_plan_multi_c2r)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    !===ru(:,1:3,:)   is first line of momentum gradient
    !===ru(:,4:6,:)   is second line of momentum gradient
    !===ru(:,7:9,:)   is third line of momentum gradient

    !===Compute -LES_coeff1_mom + entropy viscosity
    DO n = 1, bloc_size
       visco_entro(:,n) = -inputs%LES_coeff1_mom + c_in_real(:,n)
    END DO
    !===End compute -LES_coeff1_mom + entropy viscosity

    !===Compute (-LES_coeff1_mom + visco_entro) times gradient(momentum)
    prod_ru_1(:,1,:) = visco_entro*ru(:,1,:)
    prod_ru_1(:,2,:) = visco_entro*ru(:,2,:)
    prod_ru_1(:,3,:) = visco_entro*ru(:,3,:)

    prod_ru_2(:,1,:) = visco_entro*ru(:,4,:)
    prod_ru_2(:,2,:) = visco_entro*ru(:,5,:)
    prod_ru_2(:,3,:) = visco_entro*ru(:,6,:)

    prod_ru_3(:,1,:) = visco_entro*ru(:,7,:)
    prod_ru_3(:,2,:) = visco_entro*ru(:,8,:)
    prod_ru_3(:,3,:) = visco_entro*ru(:,9,:)
    !===End compute (-LES_coeff1_mom + visco_entro) times gradient(momentum)

    howmany = bloc_size*nb_field/2
    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru_1, &
         inembed, istride, idist, prod_cu_1, onembed, ostride, odist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_DOT_PROD: fftw_plan_multi_r2c', fftw_plan_multi_r2c
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru_2, &
         inembed, istride, idist, prod_cu_2, onembed, ostride, odist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_DOT_PROD: fftw_plan_multi_r2c', fftw_plan_multi_r2c
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru_3, &
         inembed, istride, idist, prod_cu_3, onembed, ostride, odist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_DOT_PROD: fftw_plan_multi_r2c', fftw_plan_multi_r2c
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    prod_cu_1 = prod_cu_1*(1.d0/N_r_pad) !Scaling
    prod_cu_2 = prod_cu_2*(1.d0/N_r_pad) !Scaling
    prod_cu_3 = prod_cu_3*(1.d0/N_r_pad) !Scaling

    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

    !Now we need to redistribute the Fourier coefficients on each processor
    t = MPI_WTIME()
    combined_prod_cu_1(:,:,1)=prod_cu_1(1,:,:)
    combined_prod_cu_2(:,:,1)=prod_cu_2(1,:,:)
    combined_prod_cu_3(:,:,1)=prod_cu_3(1,:,:)

    DO n=2, m_max
       combined_prod_cu_1(:,:,n)=2*CONJG(prod_cu_1(n,:,:))
       combined_prod_cu_2(:,:,n)=2*CONJG(prod_cu_2(n,:,:))
       combined_prod_cu_3(:,:,n)=2*CONJG(prod_cu_3(n,:,:))
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    t = MPI_WTIME()
    longueur_tranche=bloc_size*m_max_c*nb_field
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu_1,longueur_tranche,MPID, dist_prod_cu_1,longueur_tranche, &
         MPID,communicator,code)
    CALL MPI_ALLTOALL (combined_prod_cu_2,longueur_tranche,MPID, dist_prod_cu_2,longueur_tranche, &
         MPID,communicator,code)
    CALL MPI_ALLTOALL (combined_prod_cu_3,longueur_tranche,MPID, dist_prod_cu_3,longueur_tranche, &
         MPID,communicator,code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          intermediate_1 = dist_prod_cu_1(:,:,shiftl+i)
          intermediate_2 = dist_prod_cu_2(:,:,shiftl+i)
          intermediate_3 = dist_prod_cu_3(:,:,shiftl+i)
          DO n = 1, bloc_size
             IF (n+shiftc > np ) CYCLE
             DO i_field = 1, nb_field/2
                V_out(1,n+shiftc, i_field*2-1, i) = REAL (intermediate_1(i_field,n),KIND=8)
                V_out(1,n+shiftc, i_field*2 , i)  = AIMAG(intermediate_1(i_field,n))
                V_out(2,n+shiftc, i_field*2-1, i) = REAL (intermediate_2(i_field,n),KIND=8)
                V_out(2,n+shiftc, i_field*2 , i)  = AIMAG(intermediate_2(i_field,n))
                V_out(3,n+shiftc, i_field*2-1, i) = REAL (intermediate_3(i_field,n),KIND=8)
                V_out(3,n+shiftc, i_field*2 , i)  = AIMAG(intermediate_3(i_field,n))
             END DO
          END DO
       END DO
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

  END SUBROUTINE FFT_COMPUTE_ENTROPY_VISC_GRAD_MOM

  SUBROUTINE FFT_NO_OVERSHOOT_LEVEL_SET(communicator, c1_inout, nb_procs, bloc_size, m_max_pad, temps)
    !FFT (FFT(-1) V1 . FFT(-1) V2) = c_out
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    USE my_util
    USE input_data
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(INOUT)  :: c1_inout
    INTEGER,                         INTENT(IN)     :: nb_procs, bloc_size, m_max_pad
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    REAL(KIND=8), DIMENSION(SIZE(c1_inout,3),SIZE(c1_inout,2),bloc_size*nb_procs) :: dist_field
    REAL(KIND=8), DIMENSION(SIZE(c1_inout,3),SIZE(c1_inout,2),bloc_size*nb_procs) :: combined_field
    COMPLEX(KIND=8), DIMENSION(m_max_pad, SIZE(c1_inout,2)/2, bloc_size)     :: cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,SIZE(c1_inout,2)/2,bloc_size)      :: ru
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size)             :: prod_ru_1
    COMPLEX(KIND=8), DIMENSION(m_max_pad,bloc_size)              :: prod_cu_1
    COMPLEX(KIND=8), DIMENSION(bloc_size,SIZE(c1_inout,3)*nb_procs) :: combined_prod_cu_1
    COMPLEX(KIND=8), DIMENSION(bloc_size,SIZE(c1_inout,3)*nb_procs) :: dist_prod_cu_1
    COMPLEX(KIND=8), DIMENSION(bloc_size)                        :: intermediate_1
    INTEGER  :: np, np_tot, nb_field, m_max, m_max_c, MPID,  N_r_pad
    INTEGER(KIND=8) :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    INTEGER ::   nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code, rank
    REAL(KIND=8) :: t
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator

    CALL MPI_COMM_RANK(communicator, rank, code)

    IF (PRESENT(temps)) temps = 0.d0

    np        = SIZE(c1_inout,1)
    nb_field  = SIZE(c1_inout,2) ! Number of fields for scalar
    m_max_c   = SIZE(c1_inout,3) ! Number of complex (cosines + sines) coefficients per point
    m_max     = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    np_tot    = nb_procs*bloc_size
    N_r_pad   = 2*m_max_pad-1

    IF (MOD(nb_field,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_NO_OVERSHOOT_LEVEL_SET '
       STOP
    END IF

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()

    DO i = 1, m_max_c
       dist_field(i, 1:nb_field,  1:np) = TRANSPOSE(c1_inout(:,:,i))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 0.d0

    longueur_tranche=bloc_size*m_max_c*nb_field

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, nb_field/2
             ! Put real and imaginary parts in a complex
             ! nf=1,2,3 => V1_in
             ! nf=4,5,6 => V2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*nb_field/2

    t = MPI_WTIME()
    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_DOT_PROD: fftw_plan_multi_c2r', fftw_plan_multi_c2r
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    !===Compute MAX(0, MIN(1, c1_inout))
    DO n = 1, bloc_size
       prod_ru_1(:,n) = MAX(0.d0, MIN(1.d0, ru(:,1,n)))
    END DO
    !===End compute MAX(0, MIN(1, c1_inout))

    howmany = bloc_size*1
    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru_1, &
         inembed, istride, idist, prod_cu_1, onembed, ostride, odist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_DOT_PROD: fftw_plan_multi_r2c', fftw_plan_multi_r2c
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    prod_cu_1 = prod_cu_1*(1.d0/N_r_pad) !Scaling

    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

    !Now we need to redistribute the Fourier coefficients on each processor
    t = MPI_WTIME()
    combined_prod_cu_1(:,1)=prod_cu_1(1,:)

    DO n=2, m_max
       combined_prod_cu_1(:,n)=2*CONJG(prod_cu_1(n,:))
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    t = MPI_WTIME()
    longueur_tranche=bloc_size*m_max_c*2
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu_1,longueur_tranche,MPID, dist_prod_cu_1,longueur_tranche, &
         MPID,communicator,code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          intermediate_1 = dist_prod_cu_1(:,shiftl+i)
          DO n = 1, bloc_size
             IF (n+shiftc > np ) CYCLE
             c1_inout(n+shiftc, 1, i) = REAL (intermediate_1(n),KIND=8)
             c1_inout(n+shiftc, 2, i) = AIMAG(intermediate_1(n))
          END DO
       END DO
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

  END SUBROUTINE FFT_NO_OVERSHOOT_LEVEL_SET

  SUBROUTINE FFT_COMPRESSION_LEVEL_SET_DCL(communicator_F,communicator_S, V1_in, V2_in, c_in, c_out, &
       hloc_gauss, l_G, nb_procs, bloc_size, m_max_pad, temps)
    USE input_data
    !FFT (FFT(-1) V1 . FFT(-1) V2) = c_out
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V1_in, V2_in
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: c_in
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT) :: c_out
    REAL(KIND=8), DIMENSION(:),      INTENT(IN)  :: hloc_gauss
    INTEGER,                         INTENT(IN)  :: l_G
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    INTEGER, INTENT(IN)                                 :: nb_procs, bloc_size, m_max_pad
    COMPLEX(KIND=8), DIMENSION(m_max_pad, SIZE(V1_in,2)+SIZE(c_in,2)/2, bloc_size)  :: cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,SIZE(V1_in,2)+SIZE(c_in,2)/2,bloc_size)   :: ru
    COMPLEX(KIND=8), DIMENSION(m_max_pad,bloc_size)  :: prod_cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size) :: prod_ru
    COMPLEX(KIND=8), DIMENSION(bloc_size)            :: intermediate
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),2*SIZE(V1_in,2)+SIZE(c_in,2),bloc_size*nb_procs) :: dist_field
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),2*SIZE(V1_in,2)+SIZE(c_in,2),bloc_size*nb_procs) :: combined_field
    COMPLEX(KIND=8), DIMENSION(bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu
    COMPLEX(KIND=8), DIMENSION(bloc_size,SIZE(V1_in,3)*nb_procs) :: dist_prod_cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size)      :: norm_vel
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size/l_G)  :: norm_vel_int
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size)      :: norm_grad_phi
    REAL(KIND=8), DIMENSION(bloc_size*nb_procs)           :: hloc_gauss_tot
    REAL(KIND=8) :: hh, x, max_norm_vel_loc, max_norm_vel_loc_F, max_norm_vel_tot
    INTEGER :: rank, l
    INTEGER  :: np, np_tot, nb_field, nb_field2,  m_max, m_max_c, MPID,  N_r_pad
    INTEGER(KIND=8) :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    INTEGER ::   nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code
    REAL(KIND=8) :: t
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator_F
    MPI_Comm :: communicator_S

    IF (PRESENT(temps)) temps = 0.d0

    CALL MPI_COMM_RANK(communicator_F, rank, code)

    np       = SIZE(V1_in,1)
    nb_field = SIZE(V1_in,2) ! Number of fields
    nb_field2= SIZE(c_in,2) ! Number of fields
    m_max_c  = SIZE(V1_in,3) ! Number of complex (cosines + sines) coefficients per point
    m_max    = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    np_tot   = nb_procs*bloc_size
    N_r_pad  = 2*m_max_pad-1

    IF (MOD(nb_field,2)/=0 .OR. MOD(nb_field2,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_COMPRESSION_LEVEL_SET_DCL '
       STOP
    END IF

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()

    DO i = 1, m_max_c
       dist_field(i,1:nb_field,1:np) = TRANSPOSE(V1_in(:,:,i))
       dist_field(i,nb_field+1:2*nb_field,1:np) = TRANSPOSE(V2_in(:,:,i))
       dist_field(i,2*nb_field+1:2*nb_field+nb_field2,1:np) = TRANSPOSE(c_in(:,:,i))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 0.d0

    hloc_gauss_tot(1:np) = hloc_gauss
    IF (np/=np_tot) hloc_gauss_tot(np+1:np_tot) = 1.d100

    longueur_tranche=bloc_size*m_max_c*(nb_field*2+nb_field2)

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator_F, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, (nb_field+nb_field2/2)
             ! Put real and imaginary parts in a complex
             ! nf=1,2,3 => V1_in
             ! nf=4,5,6 => V2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*(nb_field+nb_field2/2)


    t = MPI_WTIME()
    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_DOT_PROD: fftw_plan_multi_c2r', fftw_plan_multi_c2r
    CALL dfftw_execute(fftw_plan_multi_c2r)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    !===ru(:,1:3,:)   is gradient(level_set)
    !===ru(:,4:6,:)   is velocity
    !===ru(:,7,:)     is level set

    !===Compute local velocity
    norm_vel(:,:) = SQRT(ru(:,4,:)**2 + ru(:,5,:)**2 + ru(:,6,:)**2)
    DO i = 1, 2*m_max_pad - 1
       DO l = 1, bloc_size/l_G
          x = MAXVAL(norm_vel(i,(l-1)*l_G+1:l*l_G))
          norm_vel_int(i,l) = x
       END DO
    END DO
    IF (2*m_max_pad - 1 .GE. 3) THEN
       DO l = 1, bloc_size/l_G
          DO i = 2, 2*m_max_pad - 2
             norm_vel(i,(l-1)*l_G+1:l*l_G) = MAXVAL(norm_vel_int(i-1:i+1,l))
          END DO
          norm_vel(1,(l-1)*l_G+1:l*l_G) = MAX(norm_vel_int(1,l),norm_vel_int(2,l),norm_vel_int(2*m_max_pad - 1,l))
          norm_vel(2*m_max_pad - 1,(l-1)*l_G+1:l*l_G) = &
               MAX(norm_vel_int(2*m_max_pad - 2,l),norm_vel_int(2*m_max_pad - 1,l),norm_vel_int(1,l))
       END DO
    ELSE
       DO l = 1, bloc_size/l_G
          norm_vel(1,(l-1)*l_G+1:l*l_G) = norm_vel_int(1,l)
       END DO
    END IF
    !===End compute local velocity

    !==Compute maximum of velocity
    max_norm_vel_loc=MAXVAL(norm_vel)
    CALL MPI_ALLREDUCE(max_norm_vel_loc,max_norm_vel_loc_F,1,MPI_DOUBLE_PRECISION, MPI_MAX, communicator_F, code)
    CALL MPI_ALLREDUCE(max_norm_vel_loc_F,max_norm_vel_tot,1,MPI_DOUBLE_PRECISION, MPI_MAX, communicator_S, code)
    !==End compute maximum of velocity

    !===Compute |Grad(level_set)|
    norm_grad_phi(:,:) = SQRT(ru(:,1,:)**2 + ru(:,2,:)**2 + ru(:,3,:)**2)

    !===Compute compression term
    !JLG Sept 30 2016
    hh = inputs%h_min_distance
    max_norm_vel_tot=min(1.d0,max_norm_vel_tot)
    DO n = 1, bloc_size
       !TEST JLG LC 2016
       prod_ru(:,n) = inputs%level_set_comp_factor*4*inputs%LES_coeff4* &
            max_norm_vel_tot*(2.d0*regul_tab(ru(:,7,n),0.01d0)-1.d0)* &
            (MAX((1.d0 - (2*ru(:,7,n)-1.d0)**2),0.d0)/(2*hh)-norm_grad_phi(:,n))
       !TEST JLG LC 2016
    END DO
    !=End compute compression term

    howmany = bloc_size*1
    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
         inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_DOT_PROD: fftw_plan_multi_r2c', fftw_plan_multi_r2c
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    !JLG, FEB 4, 2011
!!$       prod_cu = prod_cu/N_r !Scaling
    prod_cu = prod_cu*(1.d0/N_r_pad) !Scaling
    !JLG, FEB 4, 2011
    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

    !Now we need to redistribute the Fourier coefficients on each processor
    t = MPI_WTIME()
    combined_prod_cu(:,1)=prod_cu(1,:)
    DO n=2, m_max
       combined_prod_cu(:,n)=2*CONJG(prod_cu(n,:))
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    t = MPI_WTIME()
    longueur_tranche=bloc_size*m_max_c*2
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu,longueur_tranche,MPID, dist_prod_cu,longueur_tranche, &
         MPID,communicator_F,code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          intermediate = dist_prod_cu(:,shiftl+i)
          DO n = 1, bloc_size
             IF (n+shiftc > np ) CYCLE
             c_out(n+shiftc, 1, i) = REAL (intermediate(n),KIND=8)
             c_out(n+shiftc, 2 , i)  = AIMAG(intermediate(n))
          END DO
       END DO
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

  END SUBROUTINE FFT_COMPRESSION_LEVEL_SET_DCL

  SUBROUTINE FFT_PAR_ENTRO_VISC_DCL(communicator, V1_in, c1_in, c_in_real, hloc_gauss, &
       coeff1_in_level, V_out, c_out, nb_procs, bloc_size, m_max_pad, temps)
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    USE my_util
    USE input_data
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)         :: V1_in  !vector
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)         :: c1_in  !scalar
    REAL(KIND=8), DIMENSION(:,:),    INTENT(IN)         :: c_in_real
    REAL(KIND=8), DIMENSION(:),      INTENT(IN)         :: hloc_gauss
    REAL(KIND=8),                    INTENT(IN)         :: coeff1_in_level
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT)        :: V_out
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT)        :: c_out
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    INTEGER, INTENT(IN)                                 :: bloc_size, m_max_pad, nb_procs
    INTEGER          :: np, np_tot, nb_field1, nb_field2
    INTEGER          :: m_max, m_max_c, MPID, N_r_pad
    INTEGER(KIND=8)  :: fftw_plan_multi_c2r, fftw_plan_multi_r2c

    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),SIZE(V1_in,2)+SIZE(c1_in,2),bloc_size*nb_procs) :: dist_field
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),SIZE(V1_in,2)+SIZE(c1_in,2),bloc_size*nb_procs) :: combined_field
    COMPLEX(KIND=8), DIMENSION(m_max_pad, (SIZE(V1_in,2)+SIZE(c1_in,2))/2, bloc_size)  :: cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,(SIZE(V1_in,2)+SIZE(c1_in,2))/2,bloc_size)   :: ru
    REAL(KIND=8), DIMENSION(bloc_size*nb_procs)                      :: hloc_gauss_tot
    REAL(KIND=8), DIMENSION(2*m_max_pad-1, bloc_size)                :: visc_entro
    COMPLEX(KIND=8), DIMENSION(m_max_pad,SIZE(V1_in,2)/2,bloc_size)  :: prod_cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,SIZE(V1_in,2)/2,bloc_size) :: prod_ru
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2, bloc_size)           :: intermediate
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu
    COMPLEX(KIND=8), DIMENSION(SIZE(V1_in,2)/2,bloc_size,SIZE(V1_in,3)*nb_procs) :: dist_prod_cu
    COMPLEX(KIND=8), DIMENSION(m_max_pad,bloc_size)              :: prod_cu2
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size)             :: prod_ru2
    COMPLEX(KIND=8), DIMENSION(bloc_size,SIZE(V1_in,3)*nb_procs) :: combined_prod_cu2
    COMPLEX(KIND=8), DIMENSION(bloc_size,SIZE(V1_in,3)*nb_procs) :: dist_prod_cu2
    COMPLEX(KIND=8), DIMENSION(bloc_size)                        :: intermediate2
    INTEGER :: i_field
    INTEGER :: nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code
    REAL(KIND=8) :: t
    INTEGER :: rank

    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator

    CALL MPI_COMM_RANK(communicator, rank, code)

    IF (PRESENT(temps)) temps = 0.d0

    np       = SIZE(V1_in,1)
    nb_field1= SIZE(V1_in,2) ! Number of fields
    nb_field2= SIZE(c1_in,2) ! Number of fields
    m_max_c  = SIZE(V1_in,3) ! Number of complex (cosines + sines) coefficients per point
    m_max    = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    N_r_pad  = 2*m_max_pad-1
    np_tot   = nb_procs*bloc_size

    IF (MOD(nb_field1,2)/=0 .OR. MOD(nb_field2,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_PAR_ENTRO_VISC_DCL '
       STOP
    END IF

    ! Bloc_size is the number of points that are handled by one processor
    ! once the Fourier modes are all collected
    ! Computation of bloc_size and np_tot
    ! fin de la repartition des points

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()

    DO i = 1, m_max_c
       dist_field(i,1:nb_field1,1:np) = TRANSPOSE(V1_in(:,:,i))
       dist_field(i,nb_field1+1:nb_field1+nb_field2,1:np) = TRANSPOSE(c1_in(:,:,i))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100

    hloc_gauss_tot(1:np) = hloc_gauss
    IF (np/=np_tot) hloc_gauss_tot(np+1:np_tot) = 1.d100

    longueur_tranche=bloc_size*m_max_c*(nb_field1+nb_field2)

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, (nb_field1+nb_field2)/2
             ! Put real and imaginary parts in a complex
             ! nf=1,2,3 => V1_in
             ! nf=4 => V2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
    !idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
    !odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*(nb_field1+nb_field2)/2

    t = MPI_WTIME()
    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_CROSS_PROD: fftw_plan_multi_c2r', fftw_plan_multi_c2r
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    IF (nb_field1 == 6 .AND. nb_field2 == 2) THEN
       !===ru(:,1:3,:) = h*Grad(phi)
       !===ru(:,4,:)   = phi
       !===c_in_real   = visc_entro(:,:)

       !===Compute MAX(0,phi*(1-phi))
       prod_ru2(:,:) = MAX(0.d0, ru(:,4,:)*(1.d0 - ru(:,4,:)))

       !===Compute visc_entro
       DO n = 1, bloc_size
          visc_entro(:,n) = -coeff1_in_level + c_in_real(:,n)
       END DO

       !===Compute (coeff1_in_level-visc_entro)*(h*Grad(phi))
       prod_ru(:,1,:) = -visc_entro*ru(:,1,:)
       prod_ru(:,2,:) = -visc_entro*ru(:,2,:)
       prod_ru(:,3,:) = -visc_entro*ru(:,3,:)
    ELSE
       CALL error_petsc('error in problem type while calling FFT_PAR_ENTRO_VISC_DCL')
    END IF

    howmany = bloc_size*nb_field1/2
    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
         inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE)
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    howmany = bloc_size*1
    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru2, &
         inembed, istride, idist, prod_cu2, onembed, ostride, odist, FFTW_ESTIMATE)
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    !JLG, FEB 4, 2011
    !prod_cu = prod_cu/N_r !Scaling
    prod_cu = prod_cu*(1.d0/N_r_pad) !Scaling
    prod_cu2 = prod_cu2*(1.d0/N_r_pad) !Scaling
    !JLG, FEB 4, 2011
    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

    !Now we need to redistribute the Fourier coefficients on each processor

    t = MPI_WTIME()
    combined_prod_cu(:,:,1)=prod_cu(1,:,:)
    combined_prod_cu2(:,1)=prod_cu2(1,:)
    DO n=2, m_max
       combined_prod_cu(:,:,n)=2*CONJG(prod_cu(n,:,:))
       combined_prod_cu2(:,n)=2*CONJG(prod_cu2(n,:))
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    t = MPI_WTIME()
    longueur_tranche=bloc_size*m_max_c*nb_field1 !vector
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu,longueur_tranche,MPID, dist_prod_cu,longueur_tranche, &
         MPID,communicator,code)

    longueur_tranche=bloc_size*m_max_c*2 !scalar
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu2,longueur_tranche,MPID, dist_prod_cu2,longueur_tranche, &
         MPID,communicator,code)

    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t
    ! dimensions:
    t = MPI_WTIME()

    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          intermediate = dist_prod_cu(:,:,shiftl+i)
          intermediate2 = dist_prod_cu2(:,shiftl+i)
          DO n = 1, bloc_size
             IF (n+shiftc > np ) CYCLE
             DO i_field = 1, nb_field1/2
                V_out(n+shiftc, i_field*2-1, i) = REAL (intermediate(i_field,n),KIND=8)
                V_out(n+shiftc, i_field*2  , i) = AIMAG(intermediate(i_field,n))
             END DO
             c_out(n+shiftc, 1, i) = REAL (intermediate2(n),KIND=8)
             c_out(n+shiftc, 2, i) = AIMAG(intermediate2(n))
          END DO
       END DO
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

  END SUBROUTINE FFT_PAR_ENTRO_VISC_DCL

  SUBROUTINE FFT_PAR_SCAL_FUNCT(communicator, c1_inout, funct, nb_procs, bloc_size, m_max_pad, temps)
    ! MODIFICATION: based on FFT_NO_OVERSHOOT_LEVEL_SET
    ! used to compute a scalar field which is a function of another one
    ! uses an user defined function of the condlim
    USE my_util
    USE input_data
    USE boundary
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(INOUT)  :: c1_inout
    INTERFACE
       FUNCTION funct(x) RESULT(vv)
         IMPLICIT NONE
         REAL(KIND=8) :: x
         REAL(KIND=8) :: vv
       END FUNCTION funct
    END INTERFACE
    INTEGER,                         INTENT(IN)     :: nb_procs, bloc_size, m_max_pad
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    REAL(KIND=8), DIMENSION(SIZE(c1_inout,3),SIZE(c1_inout,2),bloc_size*nb_procs) :: dist_field
    REAL(KIND=8), DIMENSION(SIZE(c1_inout,3),SIZE(c1_inout,2),bloc_size*nb_procs) :: combined_field
    COMPLEX(KIND=8), DIMENSION(m_max_pad, SIZE(c1_inout,2)/2, bloc_size)     :: cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,SIZE(c1_inout,2)/2,bloc_size)      :: ru
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size)             :: prod_ru_1
    COMPLEX(KIND=8), DIMENSION(m_max_pad,bloc_size)              :: prod_cu_1
    COMPLEX(KIND=8), DIMENSION(bloc_size,SIZE(c1_inout,3)*nb_procs) :: combined_prod_cu_1
    COMPLEX(KIND=8), DIMENSION(bloc_size,SIZE(c1_inout,3)*nb_procs) :: dist_prod_cu_1
    COMPLEX(KIND=8), DIMENSION(bloc_size)                        :: intermediate_1
    INTEGER  :: np, np_tot, nb_field, m_max, m_max_c, MPID,  N_r_pad
    INTEGER(KIND=8) :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    INTEGER ::   nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code, rank
    REAL(KIND=8) :: t
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator

    CALL MPI_COMM_RANK(communicator, rank, code)

    IF (PRESENT(temps)) temps = 0.d0

    np        = SIZE(c1_inout,1)
    nb_field  = SIZE(c1_inout,2) ! Number of fields for scalar
    m_max_c   = SIZE(c1_inout,3) ! Number of complex (cosines + sines) coefficients per point
    m_max     = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    np_tot    = nb_procs*bloc_size
    N_r_pad   = 2*m_max_pad-1

    IF (MOD(nb_field,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_PAR_SCAL_FUNCT '
       STOP
    END IF

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()

    DO i = 1, m_max_c
       dist_field(i, 1:nb_field,  1:np) = TRANSPOSE(c1_inout(:,:,i))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 0.d0

    longueur_tranche=bloc_size*m_max_c*nb_field

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    cu = 0.d0
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, nb_field/2
             ! Put real and imaginary parts in a complex
             ! nf=1,2,3 => V1_in
             ! nf=4,5,6 => V2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad

    howmany=bloc_size*nb_field/2

    t = MPI_WTIME()
    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    !===Computation in the physical space
    DO i = 1, 2*m_max_pad-1
       DO n = 1, bloc_size
          prod_ru_1(i,n) = funct(ru(i,1,n))
       END DO
    END DO
    !===End computation in the physical space

    howmany = bloc_size*1
    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru_1, &
         inembed, istride, idist, prod_cu_1, onembed, ostride, odist, FFTW_ESTIMATE)
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    prod_cu_1 = prod_cu_1*(1.d0/N_r_pad) !Scaling

    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

    !Now we need to redistribute the Fourier coefficients on each processor
    t = MPI_WTIME()
    combined_prod_cu_1(:,1)=prod_cu_1(1,:)

    DO n=2, m_max
       combined_prod_cu_1(:,n)=2*CONJG(prod_cu_1(n,:))
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    t = MPI_WTIME()
    longueur_tranche=bloc_size*m_max_c*2
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu_1,longueur_tranche,MPID, dist_prod_cu_1,longueur_tranche, &
         MPID,communicator,code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          intermediate_1 = dist_prod_cu_1(:,shiftl+i)
          DO n = 1, bloc_size
             IF (n+shiftc > np ) CYCLE
             c1_inout(n+shiftc, 1, i) = REAL (intermediate_1(n),KIND=8)
             c1_inout(n+shiftc, 2, i) = AIMAG(intermediate_1(n))
          END DO
       END DO
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

  END SUBROUTINE FFT_PAR_SCAL_FUNCT

  SUBROUTINE FFT_MAX_MIN_VEL_DCL(communicator, V1_in, V_out, nb_procs, bloc_size, m_max_pad)
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    USE my_util
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V1_in ! VECTOR
    REAL(KIND=8), DIMENSION(2),      INTENT(OUT) :: V_out ! V_out(1)=max, V_out(2) = min
    INTEGER, INTENT(IN)                          :: bloc_size, m_max_pad, nb_procs
    INTEGER          :: np, np_tot, nb_field1, m_max, m_max_c, MPID, N_r_pad
    INTEGER(KIND=8)  :: fftw_plan_multi_c2r
    COMPLEX(KIND=8), DIMENSION(m_max_pad, SIZE(V1_in,2)/2, bloc_size)  :: cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1, SIZE(V1_in,2)/2,bloc_size)  :: ru
    REAL(KIND=8), DIMENSION(SIZE(V1_in,3),SIZE(V1_in,2),bloc_size*nb_procs) :: dist_field, combined_field
    INTEGER :: nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code
    REAL(KIND=8), DIMENSION(2*m_max_pad-1, bloc_size) ::  norm_vel_loc
    REAL(KIND=8) :: max_velocity, min_velocity
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator

    np        = SIZE(V1_in,1)
    nb_field1 = SIZE(V1_in,2) ! Number of fields
    m_max_c   = SIZE(V1_in,3) ! Number of complex (cosines + sines) coefficients per point
    m_max     = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    N_r_pad   = 2*m_max_pad-1
    np_tot    = nb_procs*bloc_size

    IF (MOD(nb_field1,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_MAX_MIN_VEL_DCL '
       STOP
    END IF

    ! Bloc_size is the number of points that are handled by one processor
    ! once the Fourier modes are all collected
    ! Computation of bloc_size and np_tot
    ! fin de la repartition des points

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT

    DO i = 1, m_max_c
       dist_field(i,1:nb_field1,1:np) = TRANSPOSE(V1_in(:,:,i))
    END DO

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 0.d0

    longueur_tranche=bloc_size*m_max_c*nb_field1

    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)

    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, nb_field1/2
             ! Put real and imaginary parts in a complex
             ! nf=1,2,3 => V1_in
             ! nf=4 => V2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*nb_field1/2

    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_CROSS_PROD: fftw_plan_multi_c2r', fftw_plan_multi_c2r
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    IF (nb_field1==6) THEN
       norm_vel_loc(:,:) = SQRT(ru(:,1,:)**2 + ru(:,2,:)**2 + ru(:,3,:)**2)
    ELSE
       norm_vel_loc(:,:) = SQRT(ru(:,1,:)**2)
    END IF
    max_velocity = MAXVAL(norm_vel_loc)
    !   min_velocity = MINVAL(norm_vel_loc)
    DO i=1, N_r_pad
       DO n=1, bloc_size
          min_velocity = norm_vel_loc(i,n)
          IF (min_velocity.GT.0.d0) THEN
             EXIT
          END IF
       END DO
    END DO
    DO i=1, N_r_pad
       DO n=1, bloc_size
          IF ((norm_vel_loc(i,n).LT.min_velocity).AND.(norm_vel_loc(i,n).GT.0.d0)) THEN
             min_velocity=norm_vel_loc(i,n)
          END IF
       END DO
    END DO
    CALL MPI_ALLREDUCE(max_velocity, V_out(1), 1, MPI_DOUBLE_PRECISION, &
         MPI_MAX, communicator, code)
    CALL MPI_ALLREDUCE(min_velocity, V_out(2), 1, MPI_DOUBLE_PRECISION, &
         MPI_MIN, communicator, code)
  END SUBROUTINE FFT_MAX_MIN_VEL_DCL

  SUBROUTINE FFT_SQUAREROOT_SCALAR(communicator, c1_inout, nb_procs, bloc_size, m_max_pad, temps)
    !FFT (FFT(-1) V1 . FFT(-1) V2) = c_out
    !This a de-aliased version of the code, FEB 4, 2011, JLG
    USE my_util
    USE input_data
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    ! Format: V_1in(1:np,1:6,1:m_max_c)
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(INOUT)  :: c1_inout
    INTEGER,                         INTENT(IN)     :: nb_procs, bloc_size, m_max_pad
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    REAL(KIND=8), DIMENSION(SIZE(c1_inout,3),SIZE(c1_inout,2),bloc_size*nb_procs) :: dist_field
    REAL(KIND=8), DIMENSION(SIZE(c1_inout,3),SIZE(c1_inout,2),bloc_size*nb_procs) :: combined_field
    COMPLEX(KIND=8), DIMENSION(m_max_pad, SIZE(c1_inout,2)/2, bloc_size)     :: cu
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,SIZE(c1_inout,2)/2,bloc_size)      :: ru
    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size)             :: prod_ru_1
    COMPLEX(KIND=8), DIMENSION(m_max_pad,bloc_size)              :: prod_cu_1
    COMPLEX(KIND=8), DIMENSION(bloc_size,SIZE(c1_inout,3)*nb_procs) :: combined_prod_cu_1
    COMPLEX(KIND=8), DIMENSION(bloc_size,SIZE(c1_inout,3)*nb_procs) :: dist_prod_cu_1
    COMPLEX(KIND=8), DIMENSION(bloc_size)                        :: intermediate_1
    INTEGER  :: np, np_tot, nb_field, m_max, m_max_c, MPID,  N_r_pad
    INTEGER(KIND=8) :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    INTEGER ::   nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code, rank
    REAL(KIND=8) :: t
    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters
    !#include "petsc/finclude/petsc.h"
    MPI_Comm :: communicator

    CALL MPI_COMM_RANK(communicator, rank, code)

    IF (PRESENT(temps)) temps = 0.d0

    np        = SIZE(c1_inout,1)
    nb_field  = SIZE(c1_inout,2) ! Number of fields for scalar
    m_max_c   = SIZE(c1_inout,3) ! Number of complex (cosines + sines) coefficients per point
    m_max     = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    np_tot    = nb_procs*bloc_size
    N_r_pad   = 2*m_max_pad-1

    IF (MOD(nb_field,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG in FFT_NO_OVERSHOOT_LEVEL_SET '
       STOP
    END IF

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part
    ! on nodal points
    ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    t = MPI_WTIME()

    DO i = 1, m_max_c
       dist_field(i, 1:nb_field,  1:np) = TRANSPOSE(c1_inout(:,:,i))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 0.d0

    longueur_tranche=bloc_size*m_max_c*nb_field

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, communicator, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    !JLG, FEB 4, 2011
    cu = 0.d0
    !JLG, FEB 4, 2011
    DO n = 1, bloc_size
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, nb_field/2
             ! Put real and imaginary parts in a complex
             ! nf=1,2,3 => V1_in
             ! nf=4,5,6 => V2_in
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),&
                  -combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    !JLG, FEB 4, 2011
    !Padding is done by initialization of cu: cu = 0
    !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
    !JLG, FEB 4, 2011

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1;
    !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
    idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
    odist=m_max_pad; onembed(1)=m_max_pad
    !JLG, FEB 4, 2011

    howmany=bloc_size*nb_field/2

    t = MPI_WTIME()
    CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_DOT_PROD: fftw_plan_multi_c2r', fftw_plan_multi_c2r
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    !===Compute SQRT( MAX(0, MIN(1, c1_inout)) )
    DO n = 1, bloc_size
       prod_ru_1(:,n) = SQRT(MAX(0.d0, MIN(1.d0, ru(:,1,n))))
    END DO
    !===End compute SQRT( MAX(0, MIN(1, c1_inout)) )

    howmany = bloc_size*1
    CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru_1, &
         inembed, istride, idist, prod_cu_1, onembed, ostride, odist, FFTW_ESTIMATE)
    !write(*,*) ' FFT_PAR_DOT_PROD: fftw_plan_multi_r2c', fftw_plan_multi_r2c
    CALL dfftw_execute(fftw_plan_multi_r2c)
    ! clean up
    CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    prod_cu_1 = prod_cu_1*(1.d0/N_r_pad) !Scaling

    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

    !Now we need to redistribute the Fourier coefficients on each processor
    t = MPI_WTIME()
    combined_prod_cu_1(:,1)=prod_cu_1(1,:)

    DO n=2, m_max
       combined_prod_cu_1(:,n)=2*CONJG(prod_cu_1(n,:))
    END DO

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    t = MPI_WTIME()
    longueur_tranche=bloc_size*m_max_c*2
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu_1,longueur_tranche,MPID, dist_prod_cu_1,longueur_tranche, &
         MPID,communicator,code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()
    DO i = 1, m_max_c
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size
          shiftl = (nb-1)*m_max_c
          intermediate_1 = dist_prod_cu_1(:,shiftl+i)
          DO n = 1, bloc_size
             IF (n+shiftc > np ) CYCLE
             c1_inout(n+shiftc, 1, i) = REAL (intermediate_1(n),KIND=8)
             c1_inout(n+shiftc, 2, i) = AIMAG(intermediate_1(n))
          END DO
       END DO
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

  END SUBROUTINE FFT_SQUAREROOT_SCALAR
  
END MODULE sft_parallele
