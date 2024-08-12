SUBROUTINE TIGHTB(BKX,BKY)

        USE INTERFACE

        IMPLICIT NONE

        REAL (KIND=8), INTENT(IN) :: BKX        !       X-COMPONENT OF WAVE VECTOR
        REAL (KIND=8), INTENT(IN) :: BKY        !       Y-COMPONENT OF WAVE VECTOR

        COMPLEX (KIND=8), DIMENSION(NORBIT,NORBIT) :: HAMILT

        INTEGER :: IX,IY,K1,K2
        REAL (KIND=8) :: PHASE

        CHARACTER (KIND=1), PARAMETER :: JOBZ = 'V'
        CHARACTER (KIND=1), PARAMETER :: UPLO = 'U'
        INTEGER, PARAMETER :: LDAB = NORBIT
        INTEGER, PARAMETER :: LWORK = 2 * NORBIT - 1
        COMPLEX (KIND=8), DIMENSION(NORBIT,NORBIT) :: AB
        COMPLEX (KIND=8), DIMENSION(LWORK) :: WORK
        REAL (KIND=8), DIMENSION(NORBIT) :: W
        REAL (KIND=8), DIMENSION(3*NORBIT-2) :: RWORK
        INTEGER :: INFO

        HAMILT = DCMPLX(0.D0,0.D0)
        AB = DCMPLX(0.D0,0.D0)

        DO IX = -NNX, NNX
          DO IY = -NNY, NNY
            PHASE = DBLE(IX)*BKX+DBLE(IY)*BKY
            DO K1 = 1, NORBIT
              DO K2 = 1, NORBIT
                HAMILT(K1,K2) = HAMILT(K1,K2) + HOPIJ(IX,IY,K1,K2)*DCMPLX(DCOS(PHASE),-DSIN(PHASE))
              END DO
            END DO
          END DO
        END DO

        DO K1 = 1, NORBIT
          DO K2 = K1, NORBIT
            AB(K1,K2) = HAMILT(K1,K2)
          END DO
        END DO

        CALL ZHEEV (JOBZ,UPLO,NORBIT,AB,LDAB,W,WORK,LWORK,RWORK,INFO)

        DO K1 = 1, NORBIT
          ENER(K1) = W(K1) - AMU
        END DO

        DO K1 = 1, NORBIT
          DO K2 = 1, NORBIT       
            VECR(K1,K2) = DREAL(AB(K1,K2))
            VECI(K1,K2) = DIMAG(AB(K1,K2))
          END DO
        END DO

END SUBROUTINE tightb
