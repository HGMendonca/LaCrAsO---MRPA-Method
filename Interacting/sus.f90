SUBROUTINE LINDHART(QX,QY,CHI)

        USE INTERFACE

        IMPLICIT NONE

        INTEGER, DIMENSION(NORBIT,NORBIT) :: INDEX

        REAL (KIND=8), INTENT(IN) :: QX,QY 
        COMPLEX (KIND=8), DIMENSION(NORBIT*NORBIT,NORBIT*NORBIT), INTENT(OUT) :: CHI      !       BARE SUSCEPTIBILITY

        REAL (KIND=8), DIMENSION(-NGRID2:NGRID2,-NGRID2:NGRID2) :: BZWEIGHT      !       WEIGHTS FOR THE EDGES OF THE BZ
        REAL (KIND=8), DIMENSION(NORBIT) :: EMUK  !       EIGENVALUES AT WAVE VECTOR K AND BAND MU
        REAL (KIND=8), DIMENSION(NORBIT) :: ENUKPQ  !       EIGENVALUES AT WAVE VECTOR K+q AND BAND NU
        REAL (KIND=8), DIMENSION(NORBIT,-NGRID2:NGRID2,-NGRID2:NGRID2) :: EMUK0  !       EIGENVALUES AT WAVE VECTOR K AND BAND MU
        REAL (KIND=8), DIMENSION(NORBIT,-NGRID2:NGRID2,-NGRID2:NGRID2) :: ENUKPQ0  !       EIGENVALUES AT WAVE VECTOR K+q AND BAND NU

        COMPLEX (KIND=8), DIMENSION(NORBIT,NORBIT,-NGRID2:NGRID2,-NGRID2:NGRID2) :: AMUSPK0 !       PROJECTION OF EIGENVECTOR IN K ON ORBITAL s
        COMPLEX (KIND=8), DIMENSION(NORBIT,NORBIT,-NGRID2:NGRID2,-NGRID2:NGRID2) :: ANUQTKPQ0 !       PROJECTION OF EIGENVECTOR IN K+q ON ORBITAL t
        COMPLEX (KIND=8), DIMENSION(NORBIT) :: AMUSK !       PROJECTION OF EIGENVECTOR IN K ON ORBITAL s
        COMPLEX (KIND=8), DIMENSION(NORBIT) :: AMUPK !       PROJECTION OF EIGENVECTOR IN K ON ORBITAL p
        COMPLEX (KIND=8), DIMENSION(NORBIT) :: ANUQKPQ !       PROJECTION OF EIGENVECTOR IN K+q ON ORBITAL q
        COMPLEX (KIND=8), DIMENSION(NORBIT) :: ANUTKPQ !       PROJECTION OF EIGENVECTOR IN K+q ON ORBITAL t

        REAL (KIND=8), DIMENSION(-NGRID2:NGRID2) :: AKX       !       WAVE VECTOR IN THE X DIRECTION
        REAL (KIND=8), DIMENSION(-NGRID2:NGRID2) :: AKY       !       WAVE VECTOR IN THE Y DIRECTION
        REAL (KIND=8), DIMENSION(-NGRID2:NGRID2) :: A1024_1024    ! AUXILIARY

        COMPLEX (KIND=8) APROD        !       TAKES THE PRODUCT OF THE EIGENVECTOR PROJECTIONS
        COMPLEX (KIND=8) ENERDEM      !       TAKES THE DENOMINATOR FOR THE SUMMATION

        REAL (KIND=8) :: D                      !     MEASURES MOMENTUM IN x AXIS
        REAL (KIND=8) :: DKX,DKY
        REAL (KIND=8) :: FENUKPQ,FEMUK
        REAL (KIND=8) :: DIFF

        INTEGER :: I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I,itn1,itn2             !       DO LOOP COUNTERS
        INTEGER :: K1,K2,K3,K4,K5,K6,L1,L2             !       DO LOOP COUNTERS

        A1024_1024  = (/ (DBLE(I), I = -NGRID2, NGRID2)  /)

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

!       TRANSLATION FROM ONE BASIS TO ANOTHER

        DO K1 = 1, NORBIT
          DO K2 = 1, NORBIT
            INDEX(K1,K2) = (K1-1)*NORBIT + K2
          END DO
        END DO

!       ZEROING THE CHI() ARRAY

        CHI = DCMPLX(0.D0,0.D0)

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

            DKX = AKX(I8) + QX
            IF(DKX.GT.PI) DKX = DKX - 2.D0 * PI
            IF(DKX.LT.-PI) DKX = DKX + 2.D0 * PI
            DKY = AKY(I7) + QY
            IF(DKY.GT.PI) DKY = DKY - 2.D0 * PI
            IF(DKY.LT.-PI) DKY = DKY + 2.D0 * PI
 
            CALL TIGHTB(DKX,DKY)
            DO K2 = 1, NORBIT
              ENUKPQ0(K2,I8,I7) = ENER(K2)
              DO I2 = 1, NORBIT
                ANUQTKPQ0(I2,K2,I8,I7) = DCMPLX(VECR(I2,K2),VECI(I2,K2))
              END DO
            END DO
          END DO
        END DO

!	NOW WE CALCULATE THE SUMMATION

        DO I1 = 1, NORBIT            !       THIS INDEXES THE ORB s
          DO I2 = 1, NORBIT          !       THIS INDEXES THE ORB p
            DO I3 = 1, NORBIT        !       THIS INDEXES THE ORB q
              DO I4 = 1, NORBIT      !       THIS INDEXES THE ORB t
                L1 = INDEX(I1,I4)
                L2 = INDEX(I2,I3)
                DO I7 = -NGRID2, NGRID2         !       KY IN THE GRID
                  DO I8 = -NGRID2, NGRID2       !       KX IN THE GRID
                    DO K1 = 1, NORBIT    !       THIS INDEXES THE BAND
                      EMUK(K1) = EMUK0(K1,I8,I7)
                      AMUSK(K1) = AMUSPK0(I1,K1,I8,I7)
                      AMUPK(K1) = AMUSPK0(I2,K1,I8,I7)
                    END DO
!
!   NOW WE DO THE SAME FOR K + Q
!
                    DO K2 = 1, NORBIT
                      ENUKPQ(K2) = ENUKPQ0(K2,I8,I7)
                      ANUQKPQ(K2) = ANUQTKPQ0(I3,K2,I8,I7)
                      ANUTKPQ(K2) = ANUQTKPQ0(I4,K2,I8,I7)
                    END DO
!
!       NOW WE HAVE EVERYTHING WE NEED TO START THE SUMMATIONS
!
                    DO I9 = 1, NORBIT    !       MU BANDS
                      DO I10 = 1, NORBIT !       NU BANDS
                          APROD = AMUSK(I9)*DCONJG(AMUPK(I9))*ANUQKPQ(I10)*DCONJG(ANUTKPQ(I10))
                          FENUKPQ = 1.D0/(DEXP((ENUKPQ(I10))/TEMP)+1.D0) ! CHEMIC POTENT WAS USED ALREADY IN THE TIGHT BINDING STEP
                          FEMUK = 1.D0/(DEXP((EMUK(I9))/TEMP)+1.D0) !  CHEMICAL POTENTIAL WAS USED ALREADY IN THE TIGHT BINDING STEP
                          DIFF = FENUKPQ - FEMUK
                          ENERDEM = DCMPLX(ENUKPQ(I10)-EMUK(I9),ETA)
                          CHI(L1,L2) = CHI(L1,L2) - DIFF*APROD/ENERDEM/DBLE(NPOINTS)*BZWEIGHT(I8,I7)
                      END DO!       I10
                    END DO  !       I9
                  END DO    !       I4
                END DO      !       I3
              END DO        !       I2
            END DO          !       I1
          END DO            !       I8
        END DO              !       I7

END SUBROUTINE LINDHART
 
SUBROUTINE QVALUE(K1,K2,BKX,BKY,QX,QY,NKS)
 
        IMPLICIT NONE
 
        INTEGER, INTENT(IN) :: NKS

        REAL (KIND=8), DIMENSION(0:NKS), INTENT(IN) :: BKX       !       WAVE VECTOR IN THE X DIRECTION (USED FOR qx)
        REAL (KIND=8), DIMENSION(0:NKS), INTENT(IN) :: BKY       !       WAVE VECTOR IN THE X DIRECTION (USED FOR qx)
 
        REAL (KIND=8), INTENT(OUT) :: QX
        REAL (KIND=8), INTENT(OUT) :: QY

        REAL (KIND=8) :: PI
        INTEGER, INTENT(IN) :: K1
        INTEGER, INTENT(IN) :: K2
        INTEGER :: K3

        PI = DACOS(-1.D0)

        GO TO (10,20,30) K1
 
10      CONTINUE        !       LINE (0,0) --> (PI,0)
        QY = 0.D0
        QX = BKX(K2)
        GO TO 40
 
20      CONTINUE        !       LINE (PI,0) --> (PI,PI)
        QX = PI
        QY = BKY(K2)
        GO TO 40
 
30      CONTINUE        !       LINE (PI,PI) --> (0,0)
        K3 = NKS - K2
        QX = BKX(K3)
        QY = QX
 
40      CONTINUE

END SUBROUTINE QVALUE
 


