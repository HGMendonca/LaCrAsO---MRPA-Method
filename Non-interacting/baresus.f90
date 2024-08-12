PROGRAM SUSPOL

        USE INTERFACE

        IMPLICIT NONE

        INTEGER, DIMENSION(NORBIT,NORBIT) :: INDEX

        REAL (KIND=8), DIMENSION(-NGRID2:NGRID2,-NGRID2:NGRID2) :: BZWEIGHT      !       WEIGHTS FOR THE EDGES OF THE BZ
        REAL (KIND=8), DIMENSION(NORBIT) :: EMUK  !       EIGENVALUES AT WAVE VECTOR K AND BAND MU
        REAL (KIND=8), DIMENSION(NORBIT) :: ENUKPQ  !       EIGENVALUES AT WAVE VECTOR K+q AND BAND NU
        REAL (KIND=8), DIMENSION(NORBIT,-NGRID2:NGRID2,-NGRID2:NGRID2) :: EMUK0  !       EIGENVALUES AT WAVE VECTOR K AND BAND MU
        REAL (KIND=8), DIMENSION(NORBIT,-NGRID2:NGRID2,-NGRID2:NGRID2,-NGRID2:NGRID2,-NGRID2:NGRID2) :: ENUKPQ0  !       EIGENVALUES AT WAVE VECTOR K+q AND BAND NU

        COMPLEX (KIND=8), DIMENSION(NORBIT,NORBIT,NORBIT,NORBIT,-NGRID2:NGRID2,-NGRID2:NGRID2) :: CHI      !       BARE SUSCEPTIBILITY
        COMPLEX (KIND=8), DIMENSION(NORBIT*NORBIT,NORBIT*NORBIT,-NGRID2:NGRID2,-NGRID2:NGRID2) :: CHINEW      !       BARE SUSCEPTIBILITY
        COMPLEX (KIND=8), DIMENSION(NORBIT,NORBIT,-NGRID2:NGRID2,-NGRID2:NGRID2) :: AMUSPK0 !       PROJECTION OF EIGENVECTOR IN K ON ORBITAL s
        COMPLEX (KIND=8), DIMENSION(NORBIT,NORBIT,-NGRID2:NGRID2,-NGRID2:NGRID2,-NGRID2:NGRID2,-NGRID2:NGRID2) :: ANUQTKPQ0 !       PROJECTION OF EIGENVECTOR IN K+q ON ORBITAL t
        COMPLEX (KIND=8), DIMENSION(NORBIT) :: AMUSK !       PROJECTION OF EIGENVECTOR IN K ON ORBITAL s
        COMPLEX (KIND=8), DIMENSION(NORBIT) :: AMUPK !       PROJECTION OF EIGENVECTOR IN K ON ORBITAL p
        COMPLEX (KIND=8), DIMENSION(NORBIT) :: ANUQKPQ !       PROJECTION OF EIGENVECTOR IN K+q ON ORBITAL q
        COMPLEX (KIND=8), DIMENSION(NORBIT) :: ANUTKPQ !       PROJECTION OF EIGENVECTOR IN K+q ON ORBITAL t
        COMPLEX (KIND=8), DIMENSION(-NGRID2:NGRID2,-NGRID2:NGRID2) :: ASOMA        !       FINAL RESULT

        REAL (KIND=8), DIMENSION(-NGRID2:NGRID2) :: AKX       !       WAVE VECTOR IN THE X DIRECTION
        REAL (KIND=8), DIMENSION(-NGRID2:NGRID2) :: AKY       !       WAVE VECTOR IN THE Y DIRECTION
        REAL (KIND=8), DIMENSION(-NGRID2:NGRID2) :: BKX       !       WAVE VECTOR IN THE X DIRECTION
        REAL (KIND=8), DIMENSION(-NGRID2:NGRID2) :: BKY       !       WAVE VECTOR IN THE Y DIRECTION
        REAL (KIND=8), DIMENSION(-NGRID2:NGRID2) :: A1024_1024    ! AUXILIARY

        COMPLEX (KIND=8) APROD        !       TAKES THE PRODUCT OF THE EIGENVECTOR PROJECTIONS
        COMPLEX (KIND=8) ENERDEM      !       TAKES THE DENOMINATOR FOR THE SUMMATION

        REAL (KIND=8) :: D                      !     MEASURES MOMENTUM IN x AXIS
        REAL (KIND=8) :: DKX,DKY
        REAL (KIND=8) :: QX,QY
        REAL (KIND=8) :: FENUKPQ,FEMUK
        REAL (KIND=8) :: DIFF

        INTEGER :: I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I,itn1,itn2             !       DO LOOP COUNTERS
        INTEGER :: K1,K2,K3,K4,K5,K6,L1,L2             !       DO LOOP COUNTERS
        INTEGER :: KX,KY,IR1,IR2

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

        A1024_1024  = (/ (DBLE(I), I = -NGRID2, NGRID2)  /)

        CALL CHEMPOT    !       INITIALIZES THE CHEMICAL POTENTIAL AMU

!       DEFINING WEIGHTS FOR THE BZ SUMMATION

        BZWEIGHT(-NGRID2:NGRID2,-NGRID2:NGRID2) = 1.D0

        DO I1 = -NGRID2, NGRID2
          DO I2 = -NGRID2, NGRID2
            IF(IABS(I1).EQ.NGRID2.AND.IABS(I2).EQ.NGRID2) BZWEIGHT(I2,I1) = 1.D0/4.D0
            IF(IABS(I1).EQ.NGRID2.AND.IABS(I2).NE.NGRID2) BZWEIGHT(I2,I1) = 1.D0/2.D0
            IF(IABS(I2).EQ.NGRID2.AND.IABS(I1).NE.NGRID2) BZWEIGHT(I2,I1) = 1.D0/2.D0
          END DO
        END DO

!       WAVE VECTORS FOR THE GRID

        AKX(-NGRID2:NGRID2) = A1024_1024(-NGRID2:NGRID2) * PI/DBLE(NGRID2)

        AKY(-NGRID2:NGRID2) = AKX(-NGRID2:NGRID2)

!       WAVE VECTORS FOR Q

        BKX(-NGRID2:NGRID2) = A1024_1024(-NGRID2:NGRID2) * PI/DBLE(NGRID2)

        BKY(-NGRID2:NGRID2) = BKX(-NGRID2:NGRID2)

!       ZEROING THE CHI() ARRAY

        CHI(1:NORBIT,1:NORBIT,1:NORBIT,1:NORBIT,-NGRID2:NGRID2,-NGRID2:NGRID2) = DCMPLX(0.D0,0.D0)

!       TRANSLATION FROM ONE BASIS TO ANOTHER

        DO K1 = 1, NORBIT
          DO K2 = 1, NORBIT
            INDEX(K1,K2) = (K1-1)*NORBIT + K2
          END DO
        END DO

!	HERE WE CALCULATE THE BAND STRUCTURE IN THE GRID AND STORE 
!	IT IN THE APPROPRIATE VECTORS

        DO I7 = -NGRID2, NGRID2         !       KY IN THE GRID
          DO I8 = -NGRID2, NGRID2       !       KX IN THE GRID
            DKX = AKX(I8)
            DKY = AKY(I7)

            CALL TIGHTB(DKX,DKY)

            DO K1 = 1, NORBIT
              EMUK0(K1,I8,I7) = ENER(K1)
              DO I1 = 1, NORBIT
                AMUSPK0(I1,K1,I8,I7) = DCMPLX(VECR(I1,K1),VECI(I1,K1))
              END DO
            END DO

!   NOW WE DO THE SAME FOR K + Q

            DO I5 = -NGRID2,NGRID2    !       THIS INDEXES Qx
              DO I6 = -NGRID2,NGRID2 !       THIS INDEXES Qy
                QX = BKX(I5)
                QY = BKY(I6)
                DKX = AKX(I8) + QX
                IF(DKX.GT.PI) DKX = DKX - 2.D0 * PI
                IF(DKX.LT.-PI) DKX = DKX + 2.D0 * PI
                DKY = AKY(I7) + QY
                IF(DKY.GT.PI) DKY = DKY - 2.D0 * PI
                IF(DKY.LT.-PI) DKY = DKY + 2.D0 * PI
 
                CALL TIGHTB(DKX,DKY)

                DO K2 = 1, NORBIT
                  ENUKPQ0(K2,I6,I5,I8,I7) = ENER(K2)
                  DO I2 = 1, NORBIT
                    ANUQTKPQ0(I2,K2,I6,I5,I8,I7) = DCMPLX(VECR(I2,K2),VECI(I2,K2))
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

!	NOW WE CALCULATE THE SUMMATION

        DO I7 = -NGRID2, NGRID2         !       KY IN THE GRID
          DO I8 = -NGRID2, NGRID2       !       KX IN THE GRID
            DO I5 = -NGRID2, NGRID2    !       Qx
              DO I6 = -NGRID2, NGRID2 !       Qy
                DO I1 = 1, NORBIT            !       THIS INDEXES THE ORB s
                  DO I2 = 1, NORBIT          !       THIS INDEXES THE ORB p
                    DO I3 = 1, NORBIT        !       THIS INDEXES THE ORB q
                      DO I4 = 1, NORBIT      !       THIS INDEXES THE ORB t
                        DO K1 = 1, NORBIT    !       THIS INDEXES THE BAND
                          EMUK(K1) = EMUK0(K1,I8,I7)
                          AMUSK(K1) = AMUSPK0(I1,K1,I8,I7)
                          AMUPK(K1) = AMUSPK0(I2,K1,I8,I7)
                        END DO
!
!   NOW WE DO THE SAME FOR K + Q
!
                        DO K2 = 1, NORBIT
                          ENUKPQ(K2) = ENUKPQ0(K2,I6,I5,I8,I7)
                          ANUQKPQ(K2) = ANUQTKPQ0(I3,K2,I6,I5,I8,I7)
                          ANUTKPQ(K2) = ANUQTKPQ0(I4,K2,I6,I5,I8,I7)
                        END DO
!
!       NOW WE HAVE EVERYTHING WE NEED TO START THE SUMMATIONS
!
                        DO I9 = 1, NORBIT    !       MU BANDS
                          DO I10 = 1, NORBIT !       NU BANDS
                              APROD = AMUSK(I9)*DCONJG(AMUPK(I9))*ANUQKPQ(I10)*DCONJG(ANUTKPQ(I10))
                              FENUKPQ = 1.D0/(DEXP((ENUKPQ(I10))/TEMP)+1.D0) !  CHEMICAL POTENTIAL WAS USED ALREADY IN THE TIGHT BINDING STEP
                              FEMUK = 1.D0/(DEXP((EMUK(I9))/TEMP)+1.D0) !  CHEMICAL POTENTIAL WAS USED ALREADY IN THE TIGHT BINDING STEP
                              DIFF = FENUKPQ - FEMUK
                              ENERDEM = DCMPLX(ENUKPQ(I10)-EMUK(I9),ETA)
                              CHI(I4,I3,I2,I1,I6,I5) = CHI(I4,I3,I2,I1,I6,I5) - DIFF*APROD/ENERDEM/DBLE(NPOINTS)*BZWEIGHT(I8,I7)
                          END DO!       I10
                        END DO  !       I9
                      END DO    !       I4
                    END DO      !       I3
                  END DO        !       I2
                END DO          !       I1
              END DO            !       I6
            END DO              !       I5
          END DO                !       I8
        END DO                  !       I7

        OPEN(UNIT=1,FILE='baresu158x158n435.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')

        DO I5 = -NGRID2, NGRID2          !       THIS INDEXES Qx
          DO I6 = -NGRID2, NGRID2        !       THIS INDEXES Qy
            DO I1 = 1, NORBIT            !       THIS INDEXES THE ORB s
              DO I2 = 1, NORBIT          !       THIS INDEXES THE ORB p
                DO I3 = 1, NORBIT        !       THIS INDEXES THE ORB q
                  DO I4 = 1, NORBIT      !       THIS INDEXES THE ORB t
                    WRITE(1)CHI(I4,I3,I2,I1,I6,I5)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

        CLOSE(1)

        DO K1 = -NGRID2, NGRID2
          DO K2 = -NGRID2, NGRID2
            ASOMA(K1,K2) = DCMPLX(0.D0,0.D0)
            DO K3 = 1, NORBIT
              DO K4 = 1, NORBIT
                ASOMA(K1,K2) = ASOMA(K1,K2) + CHI(K3,K4,K4,K3,K2,K1)/2.D0
              END DO
            END DO
          END DO
        END DO

        OPEN(UNIT=2,FILE='baresus158n435.dat',STATUS='UNKNOWN')

        DO K1 = -NGRID2, NGRID2
          WRITE(2,100)(DREAL(ASOMA(K1,K2)),K2 = -NGRID2, NGRID2)
        END DO

100     FORMAT(300(F10.6,1X))

        CLOSE(2)

END PROGRAM SUSPOL

