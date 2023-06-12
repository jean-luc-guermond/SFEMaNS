MODULE interpolation_tools

CONTAINS
  SUBROUTINE backup_suite(type_fic, filename,  n_it, rank, is_sequential)
    USE chaine_caractere
    IMPLICIT NONE
    CHARACTER(len=3),       INTENT(IN) :: type_fic
    CHARACTER(len=200),     INTENT(IN) :: filename
    INTEGER,                INTENT(IN) :: n_it, rank
    LOGICAL,                INTENT(IN) :: is_sequential

    CHARACTER(len=3)                   :: tit_s, tit
    CHARACTER(len=5)                   :: name_it, name_dom
    INTEGER                            :: l, lblank
    CHARACTER(len=200)                 :: cmd, dir_out

    WRITE(tit_s,'(i3)') rank
    lblank = eval_blank(3,tit_s)
    DO l = 1, lblank - 1
       tit_s(l:l) = '0'
    END DO
    IF (is_sequential) THEN
       name_dom=''
       dir_out='NON_PETSC_OUT/'
    ELSE
       name_dom='_S'//tit_s
       dir_out='PETSC_OUT'
    END IF
    WRITE(tit,'(i3)') n_it
    lblank = eval_blank(3,tit)
    DO l = 1, lblank - 1
       tit(l:l) = '0'
    END DO
    name_it='_I'//tit
    IF (type_fic=='mxw') THEN
       cmd = 'mv suite_maxwell'
    ELSE IF (type_fic=='nst') THEN
       cmd = 'mv suite_ns'
    ELSE
       WRITE(*,*) 'WARNING: exit without doing anything (backup_suite)'
       RETURN
    END IF

    cmd = trim(adjustl(cmd))//trim(adjustl(name_dom))//name_it//'.'//trim(adjustl(filename))
    WRITE(*,*) '('//trim(adjustl(cmd))//')'
    !JLG+LC Jan 6 2014
    IF (is_sequential) THEN
       IF (rank==0) CALL system(trim(adjustl(cmd))//' '//trim(adjustl(dir_out)))
    ELSE
       CALL system(trim(adjustl(cmd))//' '//trim(adjustl(dir_out)))
    END IF

  END SUBROUTINE backup_suite

  SUBROUTINE max_controle(cont, comm)
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    INTEGER, DIMENSION(:)                      :: cont
    INTEGER, DIMENSION(:), ALLOCATABLE         :: tmp
    INTEGER                                    :: np
    MPI_Comm           :: comm
    PetscErrorCode     :: ierr

    np =  SIZE(cont)
    ALLOCATE(tmp(np))
    tmp = 0
    CALL MPI_ALLREDUCE(cont, tmp, np, MPI_INTEGER, MPI_SUM, comm, ierr)
    cont = tmp
    DEALLOCATE(tmp)

  END SUBROUTINE max_controle

  SUBROUTINE somme_spatiale(field, comm, cont)
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:,:)               :: field
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE  :: tmp
    INTEGER                                      :: np, i, j, n2, n3
    INTEGER, DIMENSION(:), OPTIONAL              :: cont

    MPI_Comm        :: comm
    PetscErrorCode  :: ierr

    np = SIZE(field,1)
    n2 = SIZE(field,2)
    n3 = SIZE(field,3)

    ALLOCATE(tmp(np, n2, n3))
    tmp = 0.d0

    DO j = 1, n3
       DO i = 1, n2
          CALL MPI_ALLREDUCE(field(:,i,j), tmp(:,i,j), np, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
       END DO
    END DO
    IF (present(cont)) THEN
       DO j = 1, n3
          DO i = 1, n2
             field(:,i,j) = tmp(:,i,j)/cont(:)
          END DO
       END DO
    ELSE
       field = tmp
    END IF

    DEALLOCATE(tmp)

  END SUBROUTINE somme_spatiale

  SUBROUTINE inter_mesh_loc_to_glob(mesh_in, mesh_out, in_field, out_field, l_t_g, is_in, comm)
    !> Inject local field (in_field) into global field (out_field) if (is_in=true)
    !> Project global field (in_field) into local field (out_field) if (is_in=false)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                  :: mesh_in, mesh_out
    REAL(KIND=8), DIMENSION(:,:,:)   :: in_field, out_field
    INTEGER, DIMENSION(:)            :: l_t_g
    LOGICAL                          :: is_in
    INTEGER                          :: i2, i3, n1, n2, n3, m1, n
    MPI_Comm                         :: comm

    n  = mesh_out%np
    n1 = SIZE(in_field,1)
    m1 = SIZE(out_field, 1)
    n2 = SIZE(in_field, 2)
    n3 = SIZE(in_field, 3)

    IF (is_in) THEN
       DO i2 = 1, n2
          DO i3 = 1, n3
             out_field(l_t_g(1:mesh_in%dom_np),i2,i3) = in_field(1:mesh_in%dom_np,i2,i3)
          END DO
       END DO
       CALL somme_spatiale(out_field, comm)
    ELSE
       DO i2 = 1, n2
          DO i3 = 1, n3
             out_field(1:m1,i2,i3) = in_field(l_t_g(1:m1),i2,i3)
          END DO
       END DO
    END IF

  END SUBROUTINE inter_mesh_loc_to_glob

  SUBROUTINE loc_to_glob(mesh_loc, mesh_glob, l_t_g)
    USE def_type_mesh
    IMPLICIT NONE

    TYPE(mesh_type)                     :: mesh_loc, mesh_glob
    INTEGER, DIMENSION(mesh_loc%np)     :: l_t_g

    REAL(KIND=8)                        :: epsilon = 1.d-10
    INTEGER                             :: i, j
    REAL(KIND=8), DIMENSION(2)          :: r_loc, r_glob

    DO i=1, mesh_loc%np
       r_loc = mesh_loc%rr(:,i)
       DO j=1, mesh_glob%np
          r_glob = mesh_glob%rr(:,j)
          IF (SUM((r_loc-r_glob)**2) < epsilon) THEN
             l_t_g(i) = j
             EXIT
          END IF
       END DO
    END DO

    IF (MINVAL(l_t_g)==0) WRITE(*,*) 'BUG in loc_to_glob', mesh_loc%rr(:,MINLOC(l_t_g))


  END SUBROUTINE loc_to_glob


  SUBROUTINE interp_mesh(mesh_in, mesh_out, in_field, out_field, controle, type_fe)
    USE def_type_mesh
    IMPLICIT NONE

    TYPE(mesh_type)                    :: mesh_in, mesh_out
    REAL(KIND=8), DIMENSION(:,:,:)     :: in_field, out_field
    INTEGER, DIMENSION(mesh_out%np)    :: controle
    INTEGER                            :: type_fe

    INTEGER                             :: m, i, j, k, ni, l
    REAL(KIND=8), DIMENSION(mesh_in%gauss%n_w)  :: ff
    REAL(KIND=8), DIMENSION(mesh_in%gauss%n_ws)  :: ffe
    REAL(KIND=8), DIMENSION(3)          :: abc
    REAL(KIND=8), DIMENSION(2)          :: ab

    controle = 0
    out_field = 0.d0

    DO i = 1, mesh_out%np
       CALL find_elem(mesh_in, mesh_out%rr(:,i), abc, m)
       IF (m == 0) CYCLE
       CALL gauss_ff(abc, type_fe, ff)
       controle(i) = 1

       DO j = 1, SIZE(in_field,2)
          DO k = 1, SIZE(in_field,3)
             out_field(i,j,k)  = SUM(ff*in_field(mesh_in%jj(:,m),j,k))
          END DO
       END DO
       m = 0
    END DO

    DO j = 1, mesh_out%mes
       DO ni = 1, SIZE(mesh_out%jjs,1)
          i = mesh_out%jjs(ni,j)
          IF (controle(i)>0) CYCLE
          CALL find_edge(mesh_in, mesh_out%rr(:,i), m, ab)
          IF (m==0) CYCLE
          CALL gauss_ff_edge(ab, type_fe, ffe)
          controle(i) = 1
          DO l = 1, SIZE(in_field, 2)
             DO k = 1, SIZE(in_field, 3)
                out_field(i,l,k) = SUM(ffe*in_field(mesh_in%jjs(:,m),l,k))
             END DO
          END DO
       END DO
    END DO

    IF (MAXVAL(controle) > 1) WRITE(*,*) 'BUG in interp_mesh'

  END SUBROUTINE interp_mesh

  SUBROUTINE gauss_ff(abc, type_fe, ff)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(3)        :: abc
    INTEGER, INTENT(IN)               :: type_fe
    REAL(KIND=8), DIMENSION(3*type_fe):: ff

    IF (ABS(1.d0-SUM(abc)) > 1.d-12) THEN
       WRITE(*,*) 'bug in gauss_ff'
       STOP
    END IF

    IF (type_fe == 1) THEN
       ff = abc
    ELSE
       ff(1:3) = abc*(2*abc - 1)
       ff(4) = 4*abc(2)*abc(3)
       ff(5) = 4*abc(3)*abc(1)
       ff(6) = 4*abc(1)*abc(2)
    END IF
  END SUBROUTINE  gauss_ff

  SUBROUTINE find_elem(mesh, rr, abc, m)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)              :: mesh
    REAL(KIND=8), DIMENSION(2)   :: rr
    REAL(KIND=8), DIMENSION(3)   :: abc
    INTEGER                      :: m, n
    REAL(KIND=8), DIMENSION(2)   :: X1, X2, X3, Y12, Y23, Y31, R1, R2, R3

    m = 0

    DO n = 1, mesh%me
       X1 = mesh%rr(:,mesh%jj(1,n)) - rr
       X2 = mesh%rr(:,mesh%jj(2,n)) - rr
       X3 = mesh%rr(:,mesh%jj(3,n)) - rr
       Y23 = mesh%rr(:,mesh%jj(3,n))-mesh%rr(:,mesh%jj(2,n))
       Y31 = mesh%rr(:,mesh%jj(1,n))-mesh%rr(:,mesh%jj(3,n))
       Y12 = mesh%rr(:,mesh%jj(2,n))-mesh%rr(:,mesh%jj(1,n))
       R1 = pd_vect(Y23)
       R2 = pd_vect(Y31)
       R3 = pd_vect(Y12)
       abc(1) = pd_scal(X2,R1)
       abc(2) = pd_scal(X3,R2)
       abc(3) = pd_scal(X1,R3)

       IF (MINVAL(abc) < -1.d-12) CYCLE
       m=n
       abc = abc/SUM(abc)
       EXIT

    END DO

  END SUBROUTINE  find_elem

  SUBROUTINE gauss_ff_edge(ab, type_fe, ff)
    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(2)        :: ab
    INTEGER, INTENT(IN)               :: type_fe
    REAL(KIND=8), DIMENSION(1+type_fe):: ff

    IF (ABS(1.d0-SUM(ab)) > 1.d-12) THEN
       WRITE(*,*) 'bug in gauss_ff_edge'
       STOP
    END IF

    IF (type_fe == 1) THEN
       ff = ab
    ELSE
       ff(1) = ab(1)*(ab(1)-ab(2))
       ff(2) = ab(2)*(ab(2)-ab(1))
       ff(3) = 4*ab(1)*ab(2)
    END IF

  END SUBROUTINE gauss_ff_edge

  SUBROUTINE find_edge(mesh, rr, m, ab)
    USE def_type_mesh
    IMPLICIT NONE

    TYPE(mesh_type)                :: mesh
    REAL(KIND=8), DIMENSION(2)     :: rr, ab, abt
    INTEGER                        :: m
    REAL(KIND=8)                   :: x, y, h, hr
    INTEGER                        :: ms

    x = 1
    m = 1
    DO ms = 1, mesh%mes
       h = SUM((mesh%rr(:,mesh%jjs(1,ms))-mesh%rr(:,mesh%jjs(2,ms)))**2)
       h = SQRT(h)
       CALL dist(rr, mesh%rr(:,mesh%jjs(1,ms)), mesh%rr(:,mesh%jjs(2,ms)), y, abt)
       IF (y < x) THEN
          x  = y
          ab = abt
          m  = ms
          hr = h
       END IF
    END DO

    IF (x > hr) THEN
       m = 0
    END IF


  END SUBROUTINE find_edge

  SUBROUTINE dist(rr, rr1, rr2, y, abt)
    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(2)   :: rr, rr1, rr2, abt
    REAl(KIND=8)                 :: y
    REAL(KIND=8), DIMENSION(2)   :: Y12, X1, X2, R

    X1 = rr1-rr
    X2 = rr2-rr
    Y12 = X2-X1
    R = pd_vect(Y12)

    y = ABS(pd_scal(X1,R)/pd_scal(R,R))
    abt(1) = -pd_scal(X1,Y12)/pd_scal(Y12,Y12)
    abt(2) =  pd_scal(X2,Y12)/pd_scal(Y12,Y12)

    IF (abt(1)*abt(2) < -1.d-12) THEN
       y = SQRT(MIN(pd_scal(X1,X1), pd_scal(X2,X2)))
       IF (pd_scal(X1,X1) < pd_scal(X2,X2)) THEN
          abt(1) = 1.d0
          abt(2) = 0.d0
       ELSE
          abt(1) = 0.d0
          abt(2) = 1.d0
       END IF
    END IF

  END SUBROUTINE dist

  FUNCTION pd_vect(X) RESULT(Y)
    IMPLICIT NONE
    REAl(KIND=8), DIMENSION(2)      :: X, Y

    Y(1) = X(2)
    Y(2) = -X(1)
  END FUNCTION pd_vect

  FUNCTION pd_scal(X,Y) RESULT(pp)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:)      :: X,Y
    REAL(KIND=8)                    :: pp

    pp = SUM(X*Y)

  END FUNCTION pd_scal




END MODULE interpolation_tools
