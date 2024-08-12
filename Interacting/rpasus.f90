PROGRAM RPASUS

        USE INTERFACE

        IMPLICIT NONE

        INTEGER, DIMENSION(NORBIT2) :: IPIV     ! USED BY LAPACK
        COMPLEX (KIND=8), DIMENSION(NORBIT2) :: WORK     ! USED BY LAPACK
        INTEGER :: INFOC1,INFOC2,INFOS1,INFOS2     ! USED BY LAPACK

        REAL (KIND=8) :: QX,QY 
        REAL (KIND=8) :: AJK,AJPK,AUK,AVK           !       KUBO PARAMETERS (PRB 75, 224509)
        REAL (KIND=8) :: AJ,AJP,AU,AV           !       GRASER PARAMETERS
        REAL (KIND=8) :: D
        REAL (KIND=8), DIMENSION(NFILES,NORBIT2,NORBIT2) :: UC
        REAL (KIND=8), DIMENSION(NFILES,NORBIT2,NORBIT2) :: US
        REAL (KIND=8), DIMENSION(0:NLINE) :: BKX
        REAL (KIND=8), DIMENSION(0:NLINE) :: BKY
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: CHI
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: RPACHIC
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: RPACHIS
        COMPLEX (KIND=8), DIMENSION(NORBIT,NORBIT) :: RPACPART          !       DIAGONAL PART OF CHARGE SUSCEPTIBILIYY
        COMPLEX (KIND=8), DIMENSION(NORBIT,NORBIT) :: RPASPART          !       DIAGONAL PART OF SPIN SUSCEPTIBILIYY
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: ATOINVC       !       MATRIX TO BE INVERTED           
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: ATOINVS       !       MATRIX TO BE INVERTED
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: AINVEDC       !       INVERTED MATRIX 
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: AINVEDS       !       INVERTED MATRIX
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: CPROD       !       CHECKING INVERSION
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: SPROD       !       CHECKING INVERSION
        COMPLEX (KIND=8) :: ASOMAC               !       FINAL RESULT FOR CHARGE
        COMPLEX (KIND=8) :: ASOMAS               !       FINAL RESULT FOR SPIN
        CHARACTER(LEN=30), DIMENSION(NFILES) :: FILEC          !       NAMES FOR THE FILES

        INTEGER :: I1,I2,I3,I4,K1,K2,L1,L2,M1             !       DO LOOP COUNTERS
        INTEGER :: KX,KY,IR1,IR2

        INTEGER :: IUNITD

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

        CALL CHEMPOT            !       INITIALIZES THE CHEMICAL POTENTIAL AMU

        FILEC(1) =  'susn435g512a02-c1.dat'
        IF(LUVAL.EQ.0) FILEC(1) =  'susn435g512u0.dat'
        FILEC(2) =  'susn435g512a02-c2.dat'
        FILEC(3) =  'susn435g512a02-c3.dat'
        FILEC(4) =  'susn435g512a02-c4.dat'
        FILEC(5) =  'susn435g512a02-c5.dat'
        FILEC(6) =  'susn435g512a02-c6.dat'
        FILEC(7) =  'susn435g512a02-c7.dat'
        FILEC(8) =  'susn435g512a02-c8.dat'
        FILEC(9) =  'susn435g512a02-c9.dat'
        FILEC(10) =  'susn435g512a02-c10.dat'
        FILEC(11) =  'susn435g512a02-c11.dat'

!       ZEROING BOTH MATRICES UC AND US

          UC = 0.D0
          US = 0.D0

!       INTRODUCING THE INTERACTION (UNITS OF eV)

        OPEN(UNIT=42,FILE='uvaln435g512-c.dat',STATUS='UNKNOWN')

        IF(IPART.EQ.1) THEN
          OPEN(UNIT=43,FILE='susn435g512-part1.dat',STATUS='UNKNOWN')
          OPEN(UNIT=44,FILE='susn435g512-part2.dat',STATUS='UNKNOWN')
          OPEN(UNIT=45,FILE='susn435g512-part3.dat',STATUS='UNKNOWN')
          OPEN(UNIT=46,FILE='susn435g512-part4.dat',STATUS='UNKNOWN')
          OPEN(UNIT=47,FILE='susn435g512-part5.dat',STATUS='UNKNOWN')
        ENDIF

        DO M1 = 1, LUVAL + 1

          OPEN(UNIT=M1,FILE=FILEC(M1),STATUS='UNKNOWN')

!       WE FIRST INITIALIZE KUBO'S PARAMETERS

          AUK = AUINIT + DBLE(M1-1)/DBLE(LUVAL)*AUINCR
          IF(LUVAL.EQ.0) AUK = AUINIT
          AJK = AUK * ALFA
          AJPK = AJK
          AVK = AUK - 2.0D0 * AJK

          WRITE(42,300)M1,AUK,AJK,AJPK,AVK
300       FORMAT(I2,1X,4(F10.5,1X))

!       NOW WE TRANSFORM FROM KUBO'S TO GRASER'S NOTATION

          AU = AUK
          AV = AVK - AJK / 2.D0
          AJ = 2.D0 * AJK
          AJP = AJPK

!       CREATING MATRIX UC

          DO K1 = 1, NORBIT2
            UC(M1,K1,K1) = 3.D0 * AJ / 4.D0 - AV
          END DO
!
          DO K1 = 1, NORBIT
            L1 = (NORBIT+1) * K1 - NORBIT
            UC(M1,L1,L1) = AU
          END DO
!
          DO K1 = 1, NORBIT
            L1 = (NORBIT+1) * K1 - NORBIT
            DO K2 = 1, 5
              L2 = (NORBIT+1) * K2 - NORBIT
              IF(L1.NE.L2) UC(M1,L1,L2) = 2.D0 * AV
            END DO
          END DO
!
          UC(M1,2,6) = AJP
          UC(M1,3,11) = AJP
          UC(M1,4,16) = AJP
          UC(M1,5,21) = AJP
          UC(M1,8,12) = AJP
          UC(M1,9,17) = AJP
          UC(M1,10,22) = AJP
          UC(M1,14,18) = AJP
          UC(M1,15,23) = AJP
          UC(M1,20,24) = AJP
!
          UC(M1,6,2) = AJP
          UC(M1,11,3) = AJP
          UC(M1,16,4) = AJP
          UC(M1,21,5) = AJP
          UC(M1,12,8) = AJP
          UC(M1,17,9) = AJP
          UC(M1,22,10) = AJP
          UC(M1,18,14) = AJP
          UC(M1,23,15) = AJP
          UC(M1,24,20) = AJP
!       
!       CREATING MATRIX US
!
          DO K1 = 1, NORBIT2
            US(M1,K1,K1) = AJ / 4.D0 + AV
          END DO
!
          DO K1 = 1, NORBIT
            L1 = (NORBIT+1) * K1 - NORBIT
            US(M1,L1,L1) = AU
          END DO
!
          DO K1 = 1, NORBIT
            L1 = (NORBIT+1) * K1 - NORBIT
            DO K2 = 1, NORBIT
              L2 = (NORBIT+1) * K2 - NORBIT
              IF(L1.NE.L2) US(M1,L1,L2) = AJ / 2.d0
            END DO
          END DO
!
          US(M1,2,6) = AJP
          US(M1,3,11) = AJP
          US(M1,4,16) = AJP
          US(M1,5,21) = AJP
          US(M1,8,12) = AJP
          US(M1,9,17) = AJP
          US(M1,10,22) = AJP
          US(M1,14,18) = AJP
          US(M1,15,23) = AJP
          US(M1,20,24) = AJP
!
          US(M1,6,2) = AJP
          US(M1,11,3) = AJP
          US(M1,16,4) = AJP
          US(M1,21,5) = AJP
          US(M1,12,8) = AJP
          US(M1,17,9) = AJP
          US(M1,22,10) = AJP
          US(M1,18,14) = AJP
          US(M1,23,15) = AJP
          US(M1,24,20) = AJP
        END DO          !       M1 (FOR U VALUES)

!       WAVE VECTORS FOR THE SYMMETRY LINES

        DO I1 = 0, NLINE
          BKX(I1) = DBLE(I1)*PI/DBLE(NLINE)
          BKY(I1) = BKX(I1)
        END DO

        D = 0.D0
        DO I1 = 1, NSYM    !       THIS INDEXES THE SYMMETRY LINE
          DO I2 = 0, NLINE !       THIS INDEXES A POINT IN THE SYMMETRY LINE
            CALL QVALUE(I1,I2,BKX,BKY,QX,QY,NLINE)      !  RETURNS QX AND QY
            CALL LINDHART(QX,QY,CHI)    !       CALCULATES THE BARE SUSCEPTIBILITY

            DO M1= 1, LUVAL + 1
              ASOMAC = DCMPLX(0.D0,0.D0)
              ASOMAS = DCMPLX(0.D0,0.D0)
              ATOINVC = DCMPLX(0.D0,0.D0)
              ATOINVS = DCMPLX(0.D0,0.D0)
              RPACHIC = DCMPLX(0.D0,0.D0)
              RPACHIS = DCMPLX(0.D0,0.D0)
              CPROD = DCMPLX(0.D0,0.D0)
              SPROD = DCMPLX(0.D0,0.D0)

              DO I3 = 1, NORBIT2
                DO I4 = 1, NORBIT2
                  DO K1 = 1, NORBIT2
                    ATOINVC(I3,I4)=ATOINVC(I3,I4)+UC(M1,I3,K1)*CHI(K1,I4)
                    ATOINVS(I3,I4)=ATOINVS(I3,I4)-US(M1,I3,K1)*CHI(K1,I4)
                  END DO
                END DO
              END DO

              DO I3 = 1, NORBIT2
                ATOINVC(I3,I3) = 1.D0 + ATOINVC(I3,I3)
                ATOINVS(I3,I3) = 1.D0 + ATOINVS(I3,I3)
              END DO

              CALL ZGETRF(NORBIT2,NORBIT2,ATOINVC,NORBIT2,IPIV,INFOC1)
              IF(INFOC1.NE.0)THEN
                WRITE(*,*)'SOMETHING WRONG HERE C1',INFOC1
                STOP
              ENDIF

              CALL ZGETRI(NORBIT2,ATOINVC,NORBIT2,IPIV,WORK,NORBIT2,INFOC2)
              IF(INFOC2.NE.0)THEN
                WRITE(*,*)'SOMETHING WRONG HERE C2',INFOC2
                STOP
              ENDIF

              CALL ZGETRF(NORBIT2,NORBIT2,ATOINVS,NORBIT2,IPIV,INFOS1)
              IF(INFOS1.NE.0)THEN
                WRITE(*,*)'SOMETHING WRONG HERE S1',INFOS1
                STOP
              ENDIF

              CALL ZGETRI(NORBIT2,ATOINVS,NORBIT2,IPIV,WORK,NORBIT2,INFOS2)
              IF(INFOS2.NE.0)THEN
                WRITE(*,*)'SOMETHING WRONG HERE S2',INFOS2
                STOP
              ENDIF

              DO I3 = 1, NORBIT2
                DO I4 = 1, NORBIT2
                  DO K1 = 1, NORBIT2
                    RPACHIC(I3,I4)=RPACHIC(I3,I4)+CHI(I3,K1)*ATOINVC(K1,I4)
                    RPACHIS(I3,I4)=RPACHIS(I3,I4)+CHI(I3,K1)*ATOINVS(K1,I4)
                  END DO
                END DO
              END DO

!           CALCULATING THE SUM OF CHI'S COMPONENTS

              DO K1 = 1, NORBIT
                L1 = (NORBIT+1) * K1 - NORBIT
                DO K2 = 1, NORBIT
                  L2 = (NORBIT+1) * K2 - NORBIT
                  ASOMAC = ASOMAC + RPACHIC(L1,L2)/2.d0
                  ASOMAS = ASOMAS + RPACHIS(L1,L2)/2.d0
                  RPASPART(K1,K2) = RPACHIS(L1,L2)
                  RPACPART(K1,K2) = RPACHIC(L1,L2)
                END DO
              END DO
!
              D = PI*DBLE(I1-1) + BKX(I2)
              WRITE(M1,100)D,ASOMAC,ASOMAS

              IF(IPART.EQ.1)THEN
                DO K1 = 1, NORBIT
                  IUNITD = 42 + K1
                  WRITE(IUNITD,200)D,(RPASPART(K1,K2),K2=1,NORBIT)
                END DO
              ENDIF

              END DO              !       (END DO FOR DIFFERNT U VALUES)
            END DO                !       THIS DO IS FOR A POINT IN THE SYMMETRY LINE
          END DO                  !       THIS DO IS FOR THE SYMMETRY LINE

100     FORMAT(5(E13.6,1X))
200     FORMAT(E13.6,1X,5(E13.6,1X,E13.6))
        CLOSE(42)
        IF(IPART.EQ.1)THEN
          CLOSE(43)
          CLOSE(44)
          CLOSE(45)
          CLOSE(46)
          CLOSE(47)
        ENDIF

END PROGRAM RPASUS

