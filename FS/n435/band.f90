PROGRAM BAND

        USE INTERFACE

        IMPLICIT NONE

        INTEGER :: KX,KY,IR1,IR2,I1,I2

!       READ HOPPINGS SENT BY WAN-SHENG

        OPEN(UNIT=1,FILE='input-ws.dat',STATUS='UNKNOWN')

        DO KX = -NNX, NNX
          DO KY = -NNY, NNY
            READ(1,*)IR1,IR2
            DO I1 = 1, NORBIT
              READ(1,*)(HOPIJ(IR1,IR2,I1,I2),I2=1, NORBIT) 
            END DO
          END DO
        END DO

        CLOSE(1)

        CALL ALPHABETA

END PROGRAM BAND

!#############################################################################################

