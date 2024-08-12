SUBROUTINE ALPHABETA

        USE INTERFACE

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NPOCT,DIMEN) :: AKHOLE        !       CENTER OF EACH POCKET (1 - KX, 2 - KY)
        REAL (KIND=8), DIMENSION(NPOCT,0:NPOINT) :: PHI              !       POLAR ANGLE FOR THE POCKETS
        REAL (KIND=8), DIMENSION(NPOCT,0:NPOINT) :: AKXPOC         !       X COMPONENT OF FERMI WAVE VECTOR (HOLE)
        REAL (KIND=8), DIMENSION(NPOCT,0:NPOINT) :: AKYPOC        !       Y COMPONENT OF FERMI WAVE VECTOR (HOLE)
        REAL (KIND=8), DIMENSION(DIMEN) :: AKONE           !       K POINT WHERE ENERGIES ARE CALCULATED
        REAL (KIND=8), DIMENSION(2,NORBIT) :: ENERPROD    !       ENERGIES FOR TWO SUCCESSIVE k'S
        REAL (KIND=8), DIMENSION(NPOCT) :: AKRADHOLE         !       RADIAL COORDINATE OF THE HOLE POCKET

        REAL (KIND=8) :: AK1X,AK1Y,PROD,AKMOD,AKX,AKY
        INTEGER :: I0,I1,I2,K1,K2,K3,L1,IBAND

!       COORDINATES OF THE CENTER OF EACH POCKET

        AKHOLE(1,1) = 0.D0
        AKHOLE(1,2) = 0.D0
        AKHOLE(2,1) = PI
        AKHOLE(2,2) = 0.D0
        AKHOLE(3,1) = 0.D0
        AKHOLE(3,2) = PI
        AKHOLE(4,1) = PI/2.D0
        AKHOLE(4,2) = -PI/2.D0
        AKHOLE(5,1) = -PI/2.D0
        AKHOLE(5,2) = -PI/2.D0
        AKHOLE(6,1) = -PI/2.D0
        AKHOLE(6,2) = PI/2.D0
        AKHOLE(7,1) = PI/2.D0
        AKHOLE(7,2) = PI/2.D0

!       ANGULAR COORDINATES FOR ALL POCKETS BUT THE PI/2,PI/2 ONES

        DO I2 = 1, NPOCT - 4
          DO I1 = 0, NP1
            PHI(I2,I1) = DBLE(I1)*PI/2.D0/DBLE(NP1)
          END DO

          DO I1 = NP1 + 1, NP2
            PHI(I2,I1) = PHI(I2,I1-NP1) + PI/2.D0
          END DO

          DO I1 = NP2 + 1, NP3
            PHI(I2,I1) = PHI(I2,I1-NP2) + PI
          END DO

          DO I1 = NP3 + 1, NPOINT
            PHI(I2,I1) = PHI(I2,I1-NP3) + 1.5D0*PI
          END DO
        END DO

!       ANGULAR COORDINATES FOR 4 PI/2,PI/2 POCKETS (ADJUSTED FOR SYMMETRY)
!       THE SEQUENCE (4 TO 7) IS FROM THE FOURTH TO THE FIRST QUADRANT

        DO I0 = 0, 3
          DO I1 = 0, NP1
            PHI(NPOCT-I0,I1) = DBLE(I1)*PI/2.D0/DBLE(NP1) + PI/4.D0 + DBLE(I0)*PI/2.D0
          END DO

          DO I1 = NP1 + 1, NP2
            PHI(NPOCT-I0,I1) = PHI(NPOCT-I0,I1-NP1) + PI/2.D0
          END DO

          DO I1 = NP2 + 1, NP3
            PHI(NPOCT-I0,I1) = PHI(NPOCT-I0,I1-NP2) + PI
          END DO

          DO I1 = NP3 + 1, NPOINT
            PHI(NPOCT-I0,I1) = PHI(NPOCT-I0,I1-NP3) + 1.5D0*PI
          END DO
        END DO

        CALL CHEMPOT
        
        DO K1 = 1, NPOCT                         !     'DO-LOOP' FOR HOLE POCKETS
          DO K2 = 0, NPOINT                     !     RUNNING THROUGH PHI
            AKONE(1) = AKHOLE(K1,1)
            AKONE(2) = AKHOLE(K1,2)

            CALL TIGHTB(AKONE(1),AKONE(2))

            DO I1 = 1, NORBIT
              ENERPROD(1,I1) = ENER(I1)
            END DO
            PROD = 1.D0
            AKMOD = 0.D0
            DO WHILE (PROD.GT.0.D0)
              AKMOD = AKMOD + STEP
              AK1X = AKMOD*DCOS(PHI(K1,K2))
              AK1Y = AKMOD*DSIN(PHI(K1,K2))
              AKONE(1) = AKHOLE(K1,1) + AK1X
              AKONE(2) = AKHOLE(K1,2) + AK1Y

              CALL TIGHTB(AKONE(1),AKONE(2))

              DO I1 = 1, NORBIT
                ENERPROD(2,I1) = ENER(I1)
              END DO

              DO I1 = 1, NORBIT
                PROD = ENERPROD(1,I1)*ENERPROD(2,I1)
                IF (PROD.LT.0.D0) THEN
                  IBAND = I1
                  DO K3 = 1, NORBIT
                    AWEIGHT(K1,K3,K2) = DCMPLX(VECR(K3,I1),VECI(K3,I1))
                    AORB(K1,K3,K2) = DREAL(AWEIGHT(K1,K3,K2)*DCONJG(AWEIGHT(K1,K3,K2)))
                  END DO

                  AKXPOC(K1,K2) = AK1X
                  AKYPOC(K1,K2) = AK1Y
                  AKX = AKONE(1)
                  AKY = AKONE(2)

                  CALL VELOCITY(AKX,AKY,AKMOD,IBAND)

                  VFERMI(K1,K2) = VELO
                  GO TO 100
                ENDIF
              END DO
100           CONTINUE
              DO I2 = 1, NORBIT
                ENERPROD(1,I2) = ENERPROD(2,I2)   ! SET FIRST POINT AS LAST
              END DO
            END DO      !       WHILE
          END DO        !       K2
        END DO          !       K1

        OPEN(UNIT=7,FILE='pocket00n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        OPEN(UNIT=8,FILE='pocketpi0n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        OPEN(UNIT=9,FILE='pocket0pin035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        OPEN(UNIT=10,FILE='pocketpi2-pi2n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        OPEN(UNIT=11,FILE='pocket-pi2-pi2n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        OPEN(UNIT=12,FILE='pocket-pi2pi2n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        OPEN(UNIT=13,FILE='pocketpi2pi2n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        OPEN(UNIT=17,FILE='pcp00n035.dat',STATUS='UNKNOWN')
        OPEN(UNIT=18,FILE='pcppi0n035.dat',STATUS='UNKNOWN')
        OPEN(UNIT=19,FILE='pcp0pin035.dat',STATUS='UNKNOWN')
        OPEN(UNIT=20,FILE='pcppi2-pi2n035.dat',STATUS='UNKNOWN')
        OPEN(UNIT=21,FILE='pcp-pi2-pi2n035.dat',STATUS='UNKNOWN')
        OPEN(UNIT=22,FILE='pcp-pi2pi2n035.dat',STATUS='UNKNOWN')
        OPEN(UNIT=23,FILE='pcppi2pi2n035.dat',STATUS='UNKNOWN')
        OPEN(UNIT=24,FILE='coeff00n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        OPEN(UNIT=25,FILE='coeffpi0n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        OPEN(UNIT=26,FILE='coeff0pin035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        OPEN(UNIT=27,FILE='coeffpi2-pi2n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        OPEN(UNIT=28,FILE='coeff-pi2-pi2n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        OPEN(UNIT=29,FILE='coeff-pi2pi2n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        OPEN(UNIT=30,FILE='coeffpi2pi2n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')

        !OPEN(UNIT=31, FILE='Big.dat',STATUS='UNKNOWN')
        !OPEN(UNIT=32,FILE='coeff-real-n035.txt',STATUS='UNKNOWN')
        !OPEN(UNIT=33,FILE='coeff-imag-n035.txt',STATUS='UNKNOWN')

        OPEN(UNIT=41,FILE='fermi-veln035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')

        !OPEN(UNIT=42,FILE='fermi-vel-n035.dat',STATUS='UNKNOWN')

        DO K2 = 0, NPOINT
          DO K1 = 1, NPOCT
            AKRADHOLE(K1) = DSQRT(AKXPOC(K1,K2)**2.D0 + AKYPOC(K1,K2)**2.D0)
            AK1X = AKHOLE(K1,1)+AKRADHOLE(K1)*DCOS(PHI(K1,K2))
            AK1Y = AKHOLE(K1,2)+AKRADHOLE(K1)*DSIN(PHI(K1,K2))
!       CHECK IF THE POCKET k-POINT IS OUT OF THE FIRST BRILLOUIN ZONE

            IF(AK1X.GT.PI) AK1X = AK1X - 2.D0 * PI
            IF(AK1X.LT.-PI) AK1X = AK1X + 2.D0 * PI
            IF(AK1Y.GT.PI) AK1Y = AK1Y - 2.D0 * PI
            IF(AK1Y.LT.-PI) AK1Y = AK1Y + 2.D0 * PI

            WRITE(K1+6)AK1X,AK1Y
            !write(*,*)AK1X,AK1Y
            WRITE(16+K1,500)AK1X,AK1Y
            WRITE(23+K1)PHI(K1,K2),((AWEIGHT(K1,L1,K2)),L1=1,NORBIT)

          END DO
          WRITE(41)(VFERMI(L1,K2),L1=1,NPOCT)
        END DO

        DO K1 = 1, NPOCT
          DO K2 = 0, NPOINT
            AKRADHOLE(K1) = DSQRT(AKXPOC(K1,K2)**2.D0 + AKYPOC(K1,K2)**2.D0)
            AK1X = AKHOLE(K1,1)+AKRADHOLE(K1)*DCOS(PHI(K1,K2))
            AK1Y = AKHOLE(K1,2)+AKRADHOLE(K1)*DSIN(PHI(K1,K2))
            WRITE(*,*) AK1X,AK1Y
          END DO
        END DO

500     FORMAT(2(F11.6,1X))

        CLOSE(7)
        CLOSE(8)
        CLOSE(9)
        CLOSE(10)
        CLOSE(11)
        CLOSE(12)
        CLOSE(13)
        CLOSE(17)
        CLOSE(18)
        CLOSE(19)
        CLOSE(20)
        CLOSE(21)
        CLOSE(22)
        CLOSE(23)
        CLOSE(24)
        CLOSE(25)
        CLOSE(26)
        CLOSE(27)
        CLOSE(28)
        CLOSE(29)
        CLOSE(30)
        CLOSE(41)

END SUBROUTINE ALPHABETA

SUBROUTINE VELOCITY(AKX,AKY,AKMOD,IBAND)

        USE INTERFACE

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(2,4) :: AK
        REAL (KIND=8), INTENT(IN) :: AKX,AKY,AKMOD
        INTEGER, INTENT(IN) :: IBAND

        REAL (KIND=8) :: STEPHF,AK1X,AK1Y,VX,VY
        INTEGER :: I1

        STEPHF = STEP/2.D0
        AK(1,1) = AKX + STEPHF
        AK(2,1) = AKY
        AK(1,2) = AKX - STEPHF
        AK(2,2) = AKY
        AK(1,3) = AKX
        AK(2,3) = AKY + STEPHF
        AK(1,4) = AKX
        AK(2,4) = AKY - STEPHF
!
        DO I1 = 1, 4
          AK1X = AK(1,I1)
          AK1Y = AK(2,I1)

          CALL TIGHTB(AK1X,AK1Y)
              
          EGRAD(I1) = ENER(IBAND)
        END DO                  !       I1

        VX = EGRAD(1) - EGRAD(2)
        VY = EGRAD(3) - EGRAD(4)

        VELO = DSQRT(VX*VX + VY*VY)/STEP

END SUBROUTINE VELOCITY
