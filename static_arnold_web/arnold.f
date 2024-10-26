      PROGRAM ARNOLD_WEB
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DO I= -7,7
         DO J= -7,7
            DO K= -7,7
               IF (ABS(I)+ABS(J)+ABS(K).LE.7) THEN
                  IF (J.NE.0) THEN
                     WRITE(*,25) I,J,K
                  ELSEIF (J.EQ.0) THEN
                     WRITE(70,*) I,J,K
                  ENDIF
               ENDIF
            END DO
         END DO
      END DO

25    FORMAT(3X,3(I2,3X),2(F7.4,3X))

      END PROGRAM
