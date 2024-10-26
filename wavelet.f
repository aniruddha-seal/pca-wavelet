
C--- DECLARATION 

      IMPLICIT NONE
      INTEGER I,J,K,N,M,BS,BS1,FRR,ICOUNTER,JJ,IFLAG,LOC,DIM3
      PARAMETER(N=__FILL__,BS=10,BS1=2000,FRR=100) !BS1 = NO. OF POINTS TO BE COMPUTED
      DOUBLE PRECISION TI(N),CONS,A1,A2,A3,A(10000),B(N),RE(N),IM(N),
     *ABSWAV,PI,SIG,LAM,ABS1,SCL(N),ABSWA,AA,DIFF,A1P,A2P,A3P,FRE(BS1),
     *ABSWAV3(FRR),ABS11,ABS12,ABS13,FREQ,A0,SORT(FRR),TEMP

      DOUBLE COMPLEX Z(N)
      LOGICAL FIRST

      COMMON /TIME/ Z
      COMMON /PI/ PI,SIG,LAM
      COMMON /DIFF/ DIFF
      COMMON /TI/ TI

C--- PARAMETERS 

      DIFF =  1.5D-4 
      CONS =  1.0D0
      PI   =  4.0D0*DATAN(1.0D0)
      SIG  =  2.D0
      LAM  =  1.D0
C--      SIG  =  2.D0
C--      LAM  =  1.5D0

C--- RANGE OF FREQUENCIES (IN GENERAL, DEL-\OMEGA = + OR - 200CM^{-1} OF ZERO      ORDER OMEGA) 
C--- COMPUTE FREQUENCY VARIABLE FOR MORLET GROSSMAN WAVELET

      READ(5,*)DIM3  ! READ IF THE 3D DATA NEEDS TO PRINTED
      READ(5,*)A0    ! READ THE INITIAL FREQUENCY

      DO J = 1,FRR 
         A(J) = (1.0D0/(2.0D0*(A0+DBLE(10*J))*2.99792458D-02))*(LAM+
     *DSQRT(LAM*LAM+(1.0D0/(2.0D0*PI*PI*SIG*SIG))))
      END DO

C--- READING THE TIME SERIES
      DO J = 1,N
         READ(5,*) TI(J),RE(J)

C--- CONVERT FS TO PS
         TI(J)=TI(J)/1000.D0
         IM(J) = 0.D0
         Z(J) = DCMPLX(RE(J),IM(J))
         B(J) = TI(J) 
      END DO

C--- COMPUTING WAVELET TRANSFORM AT INITIAL FEW TIME POINTS

       DO J = 1,BS1
c      DO J = 10,10 
      
         DO K = 1,FRR
            CALL ABSWAV1(.TRUE.,A(K),B((J-1)*BS+1),ABSWAV3(K))
            FREQ=(1.0D0/(2.0D0*A(K)*2.99792458D-02))*(LAM+
     *DSQRT(LAM*LAM+(1.0D0/(2.0D0*PI*PI*SIG*SIG))))     
            IF (DIM3.EQ.1) WRITE(66,*)  B((J-1)*BS+1),FREQ,ABSWAV3(K)
         END DO 


C--- SORT THE ABSWAV3 IN DESCENDING ORDER AND STORE IT IN SORT(K)
         DO K=1,FRR
            SORT(K)=ABSWAV3(K)
         END DO

         DO JJ= 1,FRR-1
            DO K= JJ+1,FRR
               IF(SORT(JJ).LT.SORT(K)) THEN
               TEMP = SORT(JJ)
               SORT(JJ)  = SORT(K)
               SORT(K)  = TEMP
               ENDIF
            END DO
         END DO

C--- FIND THE TRUE MAXIMUM
         ICOUNTER=1
100      CONTINUE
         IFLAG=0
         DO WHILE (IFLAG.EQ.0.AND.ICOUNTER.LE.FRR)

C--- FIND LOCATION OF THE SORTED VALUES IN THE FREQ RANGE
            DO I = 1, FRR
               IF (ABSWAV3(I) .EQ. SORT(ICOUNTER)) LOC = I
            END DO

C--- CHECK WHICH LOCATION IS A MAXIMUM
            IF (LOC.EQ.1.OR.LOC.EQ.FRR) THEN
               ICOUNTER=ICOUNTER+1
            ELSE 
               ABS11 = ABSWAV3(LOC-1)
               ABS12 = ABSWAV3(LOC)
               ABS13 = ABSWAV3(LOC+1)
               IF (ABS12.GT.ABS11.AND.ABS12.GT.ABS13) THEN
                  A1=A(LOC-1)
                  A2=A(LOC+1)
                  IFLAG=1
                  CALL MNBRAK(B((J-1)*BS+1),A1,A2,A3)
                  IF (A3.GT.A(1).OR.A3.LT.A(FRR)) THEN
                     ICOUNTER=ICOUNTER+1
                     GOTO 100
                  END IF                              
               ELSE
                  ICOUNTER=ICOUNTER+1  
               END IF
            END IF
         END DO
      
         CALL FMAXIM(B((J-1)*BS+1),A1,A2,A3,AA,ABSWA)           
 
         WRITE(6,*) B((J-1)*BS+1), (1.0D0/(2.0D0*AA*2.99792458D-02))*
     *(LAM+DSQRT(LAM*LAM+(1.0D0/(2.0D0*PI*PI*SIG*SIG))))

C--- PRINT THE 3D NORMALIZED DATA
         DO K=1,FRR 
            FREQ=(1.0D0/(2.0D0*A(K)*2.99792458D-02))*(LAM+
     *DSQRT(LAM*LAM+(1.0D0/(2.0D0*PI*PI*SIG*SIG))))     
            IF (DIM3.EQ.1) WRITE(8,25) B((J-1)*BS+1),FREQ,
     *ABSWAV3(K)/SORT(1)
25          FORMAT(3(F15.5))
         END DO
         IF (DIM3.EQ.1) WRITE(8,*)

      END DO 

      END

C--- MAIN OVER 

C********************************************************************

C THIS SUBROUTINE RETURNS THE ABSOLUTE VALUE OF THE WAVELET TRANSFORM 

      SUBROUTINE ABSWAV1(FIRST,A,B,ABSWAV)
      IMPLICIT NONE
      INTEGER J,K,N,MI,MA,MIP,MAP,M
      PARAMETER(N=20001)
      DOUBLE PRECISION A,B,ABSWAV,TI(N),T1,T2,H,PI,SIG,LAM,
     *DIFF,A1,B1,TI1(N)
      DOUBLE COMPLEX Z(N),FUN,L1,L2,L3,L4,L5,Z1(N)
      LOGICAL FIRST

      COMMON /TIME/ Z
      COMMON /PI/ PI,SIG,LAM
      COMMON /DIFF/ DIFF
      COMMON /TI/ TI

C--- FINDING T1 AND  T2

      T1 = B - 3.0D0*DSQRT(2.0D0)*SIG*A
      T2 = B + 3.0D0*DSQRT(2.0D0)*SIG*A

C--- CHECKING T1 AND T2 ARE IN THE 'INPUT' RANGE

      IF (T1.LT.TI(1)) T1 = TI(1)
      IF (T2.GT.TI(N)) T2 = TI(N)

C--- SETTING THE INTERVAL TO BE EVEN 

      DO J=1,N 
         IF (DABS(TI(J)-T1).LT.DIFF) MIP = J
         IF (DABS(TI(J)-T2).LT.DIFF) MAP = J
      END DO  

      IF (MOD(MAP-MIP,2).EQ.0) THEN
         MI = MIP
         MA = MAP
      ELSE 
         MI = MIP+1
         MA = MAP
      END IF

      M  =  MA-MI

C--- DEFINING THE NEW ARRAY 

      DO J = 1,M+1
         TI1(J) = TI(MI+J-1)
         Z1(J)  = Z(MI+J-1)
      END DO 

      A1 =  TI1(1)
      B1 =  TI1(M+1)  
      H  =  (B1-A1)/M

C--- CALCULATING F(X0) AND F(XN)

      CALL FUNC(A,B,TI1(1),Z1(1),L1)
      CALL FUNC(A,B,TI1(M+1),Z1(M+1),L4)

C--- CALCULATING F(X(2J))

      L2 = (0.0D0,0.0D0)
      DO J = 1,M/2 - 1
         CALL FUNC(A,B,TI1(2*J),Z1(2*J),FUN)
         L2 = L2 + FUN
      END DO

C--- CALCULATING F(X(2J-1))

      L3 = (0.0D0,0.0D0)
      DO J = 1,M/2 
         CALL FUNC(A,B,TI1(2*J-1),Z1(2*J-1),FUN)
         L3 = L3 + FUN
      END DO

C--- SUMMING UP 

      L5 = (0.0D0,0.0D0)
      L5 = L5 + DCMPLX(H/(A*3.0D0),0.0D0)*(L1+(2.0D0,0.0D0)*L2
     *+(4.0D0,0.0D0)*L3+L4)
       
C--- CALCULATING 'ABSWAV' TO BE THE MODULUS 

      ABSWAV = 0.0D0
      IF (FIRST) THEN
         ABSWAV = ABSWAV + DSQRT(DBLE(L5*CONJG(L5)))
      ELSE
C--- NEGATIVE SIGN BECAUSE WE ARE LOOKING FOR A MINIMUM
         ABSWAV = ABSWAV - DSQRT(DBLE(L5*CONJG(L5)))
      END IF

C--- REFERENCE FOR INTEGATION: HTTP://EN.WIKIPEDIA.ORG/WIKI/SIMPSONS_RULE 
     
      RETURN
      END

C********************************************************************

C THIS SUBROUTINE COMPUTES THE COEFFICIENT OF EXPANSION FOR WAVELET

      SUBROUTINE FUNC(A,B,T,Z,FUN)

      INTEGER J,K,N
      PARAMETER(N=20001)
      DOUBLE PRECISION TI(N),T1,T,PI,SIG,LAM,A1,B1,A,B
      DOUBLE COMPLEX Z,PSI,FUN,FUN1
      COMMON /PI/ PI,SIG,LAM

      T1 = (T-B)/A
      
      FUN1= DCMPLX((DCOS(2.0D0*PI*LAM*T1)*EXP((-T1*T1)/
     *(2.0D0*SIG*SIG)))/(SIG*DSQRT(2.0D0*PI)),
     *-(DSIN(2.0D0*PI*LAM*T1)*EXP((-T1*T1)/
     *(2.0D0*SIG*SIG)))/(SIG*DSQRT(2.0D0*PI)))
      
      FUN = Z * FUN1 
      RETURN 

      END 

C********************************************************************

C THIS SUBROUTINE BRACKETS THE MAX OF A FUNCTION GIVEN TWO INITIALS
C CHECK NUMERICAL RECIPES

      SUBROUTINE MNBRAK(BB,AX,BX,CX)

      IMPLICIT NONE
      DOUBLE PRECISION AX,BX,CX,FA,FB,FC,GOLD,GLIMIT,TIN,BB
      PARAMETER (GOLD=1.618034D0,GLIMIT=100.0D0,TIN=1.0D-20)
      DOUBLE PRECISION DUM,FU,Q,R,U,ULIM
      DOUBLE PRECISION PI,SIG,LAM
      LOGICAL FIRST
      COMMON /PI/ PI,SIG,LAM
      PI   =  4.0D0*DATAN(1.0D0)
      SIG  =  2.D0
      LAM  =  1.D0

      CALL ABSWAV1(.FALSE.,AX,BB,FA)
      CALL ABSWAV1(.FALSE.,BX,BB,FB)

      IF (FB.GT.FA) THEN 
         DUM = AX
         AX = BX
         BX = DUM
         DUM = FB
         FB = FA
         FA = DUM
      END IF

      CX = BX + GOLD*(BX-AX)
      
      CALL ABSWAV1(.FALSE.,CX,BB,FC)

1     IF (FB.GE.FC) THEN 
          R = (BX-AX)*(FB-FC)
          Q = (BX-CX)*(FB-FA)
          U = BX-((BX-CX)*Q-(BX-AX)*R)/(2.0D0*SIGN(MAX(DABS(Q-R),TIN),
     *Q-R))
          ULIM = BX + GLIMIT*(CX-BX)

          IF ((BX-U)*(U-CX).GT.0.0D0) THEN 
             CALL ABSWAV1(.FALSE.,U,BB,FU)

             IF (FU.LT.FC) THEN 
                 AX = BX
                 FA = FB
                 BX = U
                 FB = FU
                 RETURN
             ELSE IF (FU.GT.FB) THEN 
                 CX = U
                 FC = FU
                 RETURN
             END IF

             U = CX + GOLD*(CX-BX)

             CALL ABSWAV1(.FALSE.,U,BB,FU)

          ELSE IF ((CX-U)*(U-ULIM).GT.0.0D0) THEN 

             CALL ABSWAV1(.FALSE.,U,BB,FU)

             IF (FU.LT.FC) THEN 
                 BX = CX
                 CX = U
                 U = CX + GOLD*(CX - BX)
                 FB = FC
                 FC = FU
                 CALL ABSWAV1(.FALSE.,U,BB,FU)
             END IF

          ELSE IF ((U-ULIM)*(ULIM-CX).GE.0.0D0) THEN 
              U = ULIM
              CALL ABSWAV1(.FALSE.,U,BB,FU)
 
          ELSE 
              U = CX + GOLD*(CX-BX)
              CALL ABSWAV1(.FALSE.,U,BB,FU)
          END IF
        
          AX = BX
          BX = CX
          CX = U
          FA = FB
          FB = FC
          FC = FU
        
          GOTO 1

      END IF

      RETURN
      END 

C****************************************************************** 

C THIS SUBROUTINE USES BRENT'S METHOD TO FIND MAX OF A FUNCTION
C CHECK NUMERICAL RECIPES

      SUBROUTINE FMAXIM(BB,AX,BX,CX,XMIN,BRENT)

      IMPLICIT NONE
      INTEGER ITMAX
      DOUBLE PRECISION BRENT,AX,BX,CX,TOL,XMIN,CGOLD,ZEPS
      PARAMETER(ITMAX=100,CGOLD=0.3819660D0,ZEPS=1.0D-10)
      INTEGER ITER
      DOUBLE PRECISION A,B,D,E,ETEMP,FU,FV,FW,FX,P,Q,R,
     *TOL1,TOL2,U,V,W,X,XM,BB
      LOGICAL FIRST

      TOL = 1.0D-08

      A = MIN(AX,CX)
      B = MAX(AX,CX)
      V = BX
      W = V
      X = V
      E = 0.0D0

      CALL ABSWAV1(.FALSE.,X,BB,FX)

      FV = FX
      FW = FX

      DO ITER = 1,ITMAX
         XM = 0.5D0*(A+B)
         TOL1 = TOL*DABS(X) + ZEPS
         TOL2 = 2.0D0*TOL1
          
         IF (DABS(X-XM).LE.(TOL2-0.5D0*(B-A))) GOTO 3
         IF (DABS(E).GT.TOL1) THEN 
             R = (X-W)*(FX-FV)
             Q = (X-V)*(FX-FW)
             P = (X-V)*Q-(X-W)*R
             Q = 2.0D0*(Q-R)

             IF (Q.GT.0.0D0) P = -P
             Q = DABS(Q)
             ETEMP = E
             E = D

             IF (DABS(P).GE.DABS(0.5D0*Q*ETEMP).OR.P.LE.Q*(A-X).OR.
     *P.GE.Q*(B-X)) GOTO 1

             D = P/Q
             U = X+D
             IF (U-A.LT.TOL2.OR.B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
             GOTO 2
         END IF

1        IF (X.GE.XM) THEN 
             E = A - X
         ELSE 
             E = B - X
         END IF

         D = CGOLD*E

2        IF (DABS(D).GE.TOL1) THEN 
             U = X + D
         ELSE 
             U = X + SIGN(TOL1,D)
         END IF 

         CALL ABSWAV1(.FALSE.,U,BB,FU)
 
         IF (FU.LE.FX) THEN 
        
            IF (U.GE.X) THEN 
               A = X
            ELSE 
               B = X
            END IF

            V = W
            FV = FW
            W = X
            FW = FX
            X = U
            FX = FU

         ELSE

            IF (U.LT.X) THEN 
               A = U
            ELSE
               B = U
            END IF

            IF (FU.LE.FW.OR.W.EQ.X) THEN 
               V = W
               FV = FW
               W = U
               FW = FU
            ELSE IF (FU.LE.FV.OR.V.EQ.X.OR.V.EQ.W) THEN 
               V = U
               FV = FU
            END IF

         END IF

      END DO

3     XMIN = X
      BRENT = FX

      RETURN 
      END 

C--- THE END ---
