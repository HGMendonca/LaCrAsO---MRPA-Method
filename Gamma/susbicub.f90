SUBROUTINE SUSBICUB(QX,QY,CHIOUT,CHIINT)
      
        USE INTERFACE

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NORBIT2,NORBIT2,-NGRID2:NGRID2,-NGRID2:NGRID2), INTENT(IN) :: CHIOUT! ROTATED INPUT FROM SUSPOL.F
        REAL (KIND=8), DIMENSION(NORBIT2,NORBIT2), INTENT(OUT) :: CHIINT     !       THE INTERPOLATED MATRIX
        REAL (KIND=8), INTENT(IN) :: QX,QY

        INTEGER, PARAMETER :: NP = 4    !       NUMBER OF POINTS FOR INTERPOLATION (4)

        REAL (KIND=8) :: AKXL,AKXU,AKYL,AKYU,ANSY,ANSY1,ANSY2

        REAL (KIND=8), DIMENSION(NP) :: CHIP        !      INPUT FROM SUSPOL.F
        REAL (KIND=8), DIMENSION(NP) :: CHIX        !      DERIVATIVE WITH x
        REAL (KIND=8), DIMENSION(NP) :: CHIY        !      DERIVATIVE WITH y
        REAL (KIND=8), DIMENSION(NP) :: CHIXY       !      SECOND DERIVATIVE WITH x AND y

        INTEGER, DIMENSION(NP) :: KDXP1,KDXM1,KDYP1,KDYM1,KPX,KPY

        INTEGER :: II1,KX,II2,KY,K1,I1,I2,K,J,JM1,JP1,KM1,KP1

!       WAVE VECTORS FOR THE GRID

        DO I1 = -NGRID2, NGRID2
          AKX(I1) = DBLE(I1)*PI/DBLE(NGRID2)
          AKY(I1) = AKX(I1)
        END DO

!	FIND THE CELL, CALCULATE THE DERIVATIVES, PASS 
!	THE PARAMETERS, GET THE INTERPOLATION

        DO II1 = -NGRID2, NGRID2
          IF((AKX(II1) - QX).GE.0.D0)THEN
            KX = II1
            GO TO 100
          END IF
        END DO

100     CONTINUE
        DO II2 = -NGRID2, NGRID2
          IF((AKY(II2) - QY).GE.0.D0)THEN
            KY = II2
            GO TO 200
          END IF
        END DO
200     CONTINUE

!	INDEXES FOR ALL POINTS

        KPX(1) = KX - 1
        KPY(1) = KY - 1

        KPX(2) = KX
        KPY(2) = KY - 1

        KPX(3) = KX
        KPY(3) = KY

        KPX(4) = KX - 1
        KPY(4) = KY

!       NOW WE TEST TO SEE IF THE INDEXES ARE OUT OF THE GRID

        DO K1 = 1, NP
          IF(KPX(K1).GT.NGRID2) KPX(K1) = KPX(K1) - NGRID
          IF(KPY(K1).GT.NGRID2) KPY(K1) = KPY(K1) - NGRID
          IF(KPX(K1).LT.-NGRID2) KPX(K1) = KPX(K1) + NGRID
          IF(KPY(K1).LT.-NGRID2) KPY(K1) = KPY(K1) + NGRID
        END DO

!       write(*,*)i1,i2,(kpx(l1),l1=1,NP),(kpy(l1),l1=1,NP)

!	INDEXES FOR THE DERIVATIVES

        DO K1 = 1, NP
          KDXP1(K1) = KPX(K1) + 1
          KDXM1(K1) = KPX(K1) - 1
          KDYP1(K1) = KPY(K1) + 1
          KDYM1(K1) = KPY(K1) - 1

!	NOW WE TEST TO SEE IF THE INDEXES ARE OUT OF THE GRID

          IF(KDXP1(K1).GT.NGRID2) KDXP1(K1) = KDXP1(K1) - NGRID
          IF(KDYP1(K1).GT.NGRID2) KDYP1(K1) = KDYP1(K1) - NGRID
          IF(KDXM1(K1).LT.-NGRID2) KDXM1(K1) = KDXM1(K1) + NGRID
          IF(KDYM1(K1).LT.-NGRID2) KDYM1(K1) = KDYM1(K1) + NGRID
        END DO

!       HERE WE HAVE TO INTRODUCE THE INDEXES FOR
!       THE PAIRS OF ORBITALS

        DO I1 = 1, NORBIT2         !      FIRST INDEX FOR THE LINDHARD FUNCTION
          DO I2 = 1, NORBIT2       !      SECOND INDEX FOR THE LINDHARD FUNCTION
            DO K1 = 1, NP
              J = KPX(K1)
              K = KPY(K1)
              JP1 = KDXP1(K1)
              JM1 = KDXM1(K1)
              KP1 = KDYP1(K1)
              KM1 = KDYM1(K1)
              CHIP(K1) = CHIOUT(I1,I2,J,K)
              CHIX(K1) = (CHIOUT(I1,I2,JP1,K) - CHIOUT(I1,I2,JM1,K))/AKXCELL/2.D0
              CHIY(K1) = (CHIOUT(I1,I2,J,KP1) - CHIOUT(I1,I2,J,KM1))/AKYCELL/2.D0
              CHIXY(K1) = (CHIOUT(I1,I2,JP1,KP1)-CHIOUT(I1,I2,JP1,KM1)-CHIOUT(I1,I2,JM1,KP1)&
                           +CHIOUT(I1,I2,JM1,KM1))/AKXCELL/AKYCELL/4.D0
            END DO

            AKXL = AKX(KPX(1))
            AKXU = AKX(KPX(2))
            AKYL = AKY(KPY(1))
            AKYU = AKY(KPY(3))

!           write(*,*)'call bcuint'

            CALL BCUINT(CHIP,CHIX,CHIY,CHIXY,AKXL,AKXU,AKYL,AKYU,QX,QY,ANSY,ANSY1,ANSY2)

!           write(*,*)'back from bcuint'

            CHIINT(I1,I2) = ANSY

          END DO
        END DO

END SUBROUTINE SUSBICUB

