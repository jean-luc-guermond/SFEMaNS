!
!Authors: Jean-Luc Guermond, Raphael Laguerre, Caroline Nore, Copyright 2005
!
MODULE bessel

  PUBLIC :: bessk, bessk0, bessk1
  PUBLIC :: bessi, bessi0, bessi1
  PUBLIC :: bessj0, bessj1
  PRIVATE

CONTAINS
  ! ----------------------------------------------------------------------

  FUNCTION BESSK(N,X) RESULT(vv)
    REAL *8 X,vv,TOX,BK,BKM,BKP
    INTEGER :: N, J
    ! ------------------------------------------------------------------------
    !     CE SOUS-PROGRAMME CALCULE LA FONCTION BESSEL MODIFIFIEE 3E ESPECE
    !     D'ORDRE N ENTIER POUR TOUT X REEL POSITIF > 0.  ON UTILISE ICI LA
    !     FORMULE DE RECURRENCE CLASSIQUE EN PARTANT DE BESSK0 ET BESSK1.
    !
    !     THIS ROUTINE CALCULATES THE MODIFIED BESSEL FUNCTION OF THE THIRD
    !     KIND OF INTEGER ORDER, N FOR ANY POSITIVE REAL ARGUMENT, X. THE
    !     CLASSICAL RECURSION FORMULA IS USED, STARTING FROM BESSK0 AND BESSK1.
    ! ------------------------------------------------------------------------
    !     REFERENCE:
    !     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
    !     MATHEMATICAL TABLES, VOL.5, 1962.
    ! ------------------------------------------------------------------------
    IF (N.EQ.0) THEN
       vv = BESSK0(X)
       RETURN
    ENDIF
    IF (N.EQ.1) THEN
       vv = BESSK1(X)
       RETURN
    ENDIF
    IF (X.EQ.0.D0) THEN
       vv = 1.D30
       RETURN
    ENDIF
    TOX = 2.D0/X
    BK  = BESSK1(X)
    BKM = BESSK0(X)
    DO 11 J=1,N-1
       BKP = BKM+DFLOAT(J)*TOX*BK
       BKM = BK
       BK  = BKP
11     CONTINUE
       vv = BK
       RETURN
     END FUNCTION BESSK

     ! ----------------------------------------------------------------------
     FUNCTION BESSK0(X) RESULT(vv)
       !     CALCUL DE LA FONCTION BESSEL MODIFIEE DU 3EME ESPECE D'ORDRE 0
       !     POUR TOUT X REEL NON NUL.
       !
       !     CALCULATES THE THE MODIFIED BESSEL FUNCTION OF THE THIRD KIND OF
       !     ORDER ZERO FOR ANY POSITIVE REAL ARGUMENT, X.
       ! ----------------------------------------------------------------------
       REAL*8 X,vv,Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7
       DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,0.23069756D0, &
            0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
       DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,0.2189568D-1, &
            -0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/
       IF(X.EQ.0.D0) THEN
          vv=1.D30
          RETURN
       ENDIF
       IF(X.LE.2.D0) THEN
          Y=X*X/4.D0
          AX=-LOG(X/2.D0)*BESSI0(X)
          vv=AX+(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
       ELSE
          Y=(2.D0/X)
          AX=EXP(-X)/SQRT(X)
          vv=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
       ENDIF
       RETURN
     END FUNCTION BESSK0
     ! ----------------------------------------------------------------------
     FUNCTION BESSK1(X) RESULT(vv)
       !     CALCUL DE LA FONCTION BESSEL MODIFIEE DE 3EME ESPECE D'ORDRE 1
       !     POUR TOUT X REEL POSITF NON NUL.
       !
       !     CALCULATES THE MODIFIED BESSEL FUNCTION OF THE THIRD KIND OF
       !     ORDER ONE FOR ANY POSITIVE REAL ARGUMENT, X.
       ! ----------------------------------------------------------------------
       REAL*8 X,vv,Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7
       DATA P1,P2,P3,P4,P5,P6,P7/1.D0,0.15443144D0,-0.67278579D0,  &
            -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
       DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,-0.3655620D-1, &
            0.1504268D-1,-0.780353D-2,0.325614D-2,-0.68245D-3/
       IF(X.EQ.0.D0) THEN
          vv=1.D32
          RETURN
       ENDIF
       IF(X.LE.2.D0) THEN
          Y=X*X/4.D0
          AX=LOG(X/2.D0)*BESSI1(X)
          vv=AX+(1.D0/X)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
       ELSE
          Y=(2.D0/X)
          AX=EXP(-X)/SQRT(X)
          vv=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
       ENDIF
       RETURN
     END FUNCTION BESSK1


     FUNCTION BESSI(N,X) RESULT(vv)
       !
       !     This subroutine calculates the first kind modified Bessel function
       !     of integer order N, for any REAL X. We use here the classical
       !     recursion formula, when X > N. For X < N, the Miller's algorithm
       !     is used to avoid overflows.
       !     REFERENCE:
       !     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
       !     MATHEMATICAL TABLES, VOL.5, 1962.
       !
       IMPLICIT NONE
       INTEGER, PARAMETER      :: IACC = 40
       REAL(KIND=8), PARAMETER :: BIGNO = 1.D10, BIGNI = 1.D-10
       REAL(KIND=8) ::  X,vv,TOX,BIM,BI,BIP
       INTEGER :: n, m, j
       IF (N.EQ.0) THEN
          vv = BESSI0(X)
          RETURN
       ENDIF
       IF (N.EQ.1) THEN
          vv = BESSI1(X)
          RETURN
       ENDIF
       IF(X.EQ.0.D0) THEN
          vv=0.D0
          RETURN
       ENDIF
       TOX = 2.D0/X
       BIP = 0.D0
       BI  = 1.D0
       vv = 0.D0
       M = 2*((N+INT(SQRT(FLOAT(IACC*N)))))
       DO 12 J = M,1,-1
          BIM = BIP+DFLOAT(J)*TOX*BI
          BIP = BI
          BI  = BIM
          IF (ABS(BI).GT.BIGNO) THEN
             BI  = BI*BIGNI
             BIP = BIP*BIGNI
             vv = vv*BIGNI
          ENDIF
          IF (J.EQ.N) vv = BIP
12        CONTINUE
          vv = vv*BESSI0(X)/BI
          RETURN
        END FUNCTION BESSI
        ! ----------------------------------------------------------------------
        ! ----------------------------------------------------------------------
        ! Auxiliary Bessel functions for N=0, N=1
        FUNCTION BESSI0(X) RESULT(vv)
          IMPLICIT NONE
          REAL(KIND=8) ::  X,vv,Y,P1,P2,P3,P4,P5,P6,P7,  &
               Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
          DATA P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,1.2067429D0,  &
               0.2659732D0,0.360768D-1,0.45813D-2/
          DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1, &
               0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,  &
               0.2635537D-1,-0.1647633D-1,0.392377D-2/
          IF(ABS(X).LT.3.75D0) THEN
             Y=(X/3.75D0)**2
             vv=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
          ELSE
             AX=ABS(X)
             Y=3.75D0/AX
             BX=EXP(AX)/SQRT(AX)
             AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
             vv=AX*BX
          ENDIF
          RETURN
        END FUNCTION BESSI0
        ! ----------------------------------------------------------------------
        FUNCTION BESSI1(X) RESULT(vv)
          IMPLICIT NONE
          REAL(KIND=8) :: X,vv,Y,P1,P2,P3,P4,P5,P6,P7,  &
               Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
          DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,  &
               0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
          DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1, &
               -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1, &
               -0.2895312D-1,0.1787654D-1,-0.420059D-2/
          IF(ABS(X).LT.3.75D0) THEN
             Y=(X/3.75D0)**2
             vv=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
          ELSE
             AX=ABS(X)
             Y=3.75D0/AX
             BX=EXP(AX)/SQRT(AX)
             AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
             vv=AX*BX
          ENDIF
          RETURN
        END FUNCTION BESSI1
        !-----------------------------------------------------------------------
        FUNCTION BESSJ0(x) RESULT(vv)
          IMPLICIT NONE
          REAL(KIND=8) :: x, y, ax, xx,vv , z
          REAL(KIND=8) ::  p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6
          DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4, &
               -.2073370639d-5,.2093887211d-6/
          DATA q1,q2,q3,q4,q5/-.1562499995d-1, &
               .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/

          DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0, &
               651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/
          DATA s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0,  &
               59272.64853d0,267.8532712d0,1.d0/
          !
          IF(ABS(x).LT.8.d0)THEN
             y=x**2
             vv=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*  &
                  (s4+y*(s5+y*s6)))))
          ELSE
             ax=ABS(x)
             z=8./ax
             y=z**2
             xx=ax-.785398164
             vv=SQRT(.636619772/ax)*(COS(xx)*(p1+y*(p2+y*(p3+y*(p4+y*  &
                  p5))))-z*SIN(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
          ENDIF
          RETURN
        END FUNCTION BESSJ0
        !-----------------------------------------------------------------------

        FUNCTION BESSJ1(x) RESULT(vv)
          IMPLICIT NONE
          REAL(KIND=8) :: x, y, ax, xx,vv , z
          REAL(KIND=8) ::  p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6
          !April 18 2008, JLG
          REAL(KIND=8) :: one
          DATA one/1.d0/
          !April 18 2008
          DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0, &
               242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/
          DATA s1,s2,s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0, &
               99447.43394d0,376.9991397d0,1.d0/
          DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4, &
               .2457520174d-5,-.240337019d-6/
          DATA  q1,q2,q3,q4,q5/.04687499995d0, &
               -.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
          !
          IF(ABS(x).LT.8.)THEN
             y=x**2
             vv=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+ &
                  y*(s4+y*(s5+y*s6)))))
          ELSE
             ax=ABS(x)
             z=8./ax
             y=z**2
             xx=ax-2.356194491
             !April 18 2008: 1.d0 -> one
             vv=SQRT(.636619772/ax)*(COS(xx)*(p1+y*(p2+y*(p3+y*(p4+y* &
                  p5))))-z*SIN(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*SIGN(one,x)
             !April 18 2008
          ENDIF
          RETURN
        END FUNCTION BESSJ1


      END MODULE bessel
