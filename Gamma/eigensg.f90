PROGRAM EIGENSG

        USE INTERFACE
!
        IMPLICIT NONE
!
        COMPLEX (KIND=8), DIMENSION(NORBIT,NORBIT,NORBIT,NORBIT,-NGRID2:NGRID2,-NGRID2:NGRID2) :: CHIINPUT!  INPUT FROM SUSPOL.F
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2,-NGRID2:NGRID2,-NGRID2:NGRID2) :: CHIOUT    !    ROTATED INPUT FROM SUSPOL.F
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: CHI
        COMPLEX (KIND=8), DIMENSION(1:NTOT2,1:NTOT2,NORBIT2,NORBIT2) :: CHISAVE
        REAL (KIND=8), DIMENSION(NORBIT2,NORBIT2,-NGRID2:NGRID2,-NGRID2:NGRID2) :: CHIRE
        REAL (KIND=8), DIMENSION(NORBIT2,NORBIT2,-NGRID2:NGRID2,-NGRID2:NGRID2) :: CHIIMG
        REAL (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: CHIREINT
        REAL (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: CHIIMGINT
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: RPACHIC
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: RPACHIS
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: ATOINVC       !       MATRIX TO BE INVERTED           
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: ATOINVS       !       MATRIX TO BE INVERTED
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: AINVEDC       !       INVERTED MATRIX 
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: AINVEDS       !       INVERTED MATRIX
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: CPROD       !       CHECKING INVERSION
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: SPROD       !       CHECKING INVERSION
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: AUX1       !      AUXILIARY 1 
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: AUX2       !      AUXILIARY 2
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: AUX3       !      AUXILIARY 3
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: AUX4       !      AUXILIARY 4
        COMPLEX (KIND=8), DIMENSION(NORBIT2,NORBIT2) :: AUX5       !      AUXILIARY 5

        COMPLEX (KIND=8), DIMENSION(NPOC,0:NPSURF,NPOC,0:NPSURF,NORBIT2,NORBIT2) :: GAMMA         !  INTERACTION VERTEX
        COMPLEX (KIND=8), DIMENSION(NPOC,0:NPSURF,NPOC,0:NPSURF) :: GAMMAIJ         !  INTERACTION KERNEL
        COMPLEX (KIND=8), DIMENSION(NPOC,0:NPSURF,NPOC,0:NPSURF) :: GAMMAIJSYM         !  SYMMETRIC INTERACTION KERNEL
        COMPLEX (KIND=8), DIMENSION(1:NTOT2,1:NTOT2) :: GAMMA2SAVE
        COMPLEX (KIND=8), DIMENSION(NPOC,NORBIT,0:NPSURF) :: AWEIGHT !       COEFFICIENT OF EACH ORBITAL AT THE POCKET SURFACE

        INTEGER, DIMENSION(2,NORBIT2) :: MAT1           !       MAKES THE INDEX INTERCHANGE AT THE END
        INTEGER, DIMENSION(NORBIT,NORBIT) :: MAT2            !       MAKES THE INDEX INTERCHANGE AT THE END
        INTEGER, DIMENSION(NORBIT,NORBIT) :: INDEX      !       CONTAINS THE 'ROTATION' FROM ONE BASIS TO ANOTHER

        INTEGER, DIMENSION(NORBIT2) :: IPIV     ! USED BY LAPACK
        COMPLEX (KIND=8), DIMENSION(NORBIT2) :: WORK     ! USED BY LAPACK
        INTEGER :: INFOC1,INFOC2,INFOS1,INFOS2     ! USED BY LAPACK

        REAL (KIND=8), DIMENSION(NTOT2) :: EIGENR    !       EIGENVALUES (REAL PART)
        REAL (KIND=8), DIMENSION(NTOT2) :: EIGENI    !       EIGENVALUES (IMAGINARY PART)
        REAL (KIND=8), DIMENSION(NTOT2,NTOT2) :: VECR    !       EIGENVECTORS (REAL PART)
        REAL (KIND=8), DIMENSION(NTOT2,NTOT2) :: VECI    !       EIGENVECTORS (IMAGINARY PART)
        COMPLEX (KIND=8), DIMENSION(NPOC,NPOC) :: LIJ     !       PAIRING STRENGTH \lambda_ij
        COMPLEX (KIND=8) :: NZERO                       !        \SUM_j N_j(0)

        COMPLEX (KIND=8) :: AAA1,BBB1

        INTEGER :: I1,I2,I3,I4,I5,I6,I7,I8
        INTEGER :: L1,L2,L3,L4,LL1,LL2,LLL1
        INTEGER :: K1,K2,K3,K4,K5,K6,KK1,KK2,KK3,KK4,K
        INTEGER :: J2,J4
        INTEGER :: M1,M3,M4,N2
        INTEGER :: NQ,NP,NT,NS

        REAL (KIND=8) :: QX,QY
        REAL (KIND=8), DIMENSION(-NTOT2/2:NTOT2/2, -NTOT2/2:NTOT2/2) :: SAVEQX
        REAL (KIND=8), DIMENSION(-NTOT2/2:NTOT2/2, -NTOT2/2:NTOT2/2) :: SAVEQY

        integer :: indice1, indice2

!       INTRODUCING THE INTERACTION (UNITS OF eV)

        DO I1 = 1, NITER
          IF(NITER.EQ.1) AUK(I1) = AUINIT + AUINTER*DBLE(I1-1)/DBLE(NITER)
          IF(NITER.GT.1) AUK(I1) = AUINIT + AUINTER*DBLE(I1-1)/DBLE(NITER-1)
          AJK(I1) = AUK(I1) * ALFA
          AJPK(I1) = AJK(I1)
          AVK(I1) = AUK(I1) - 2.0D0 * AJK(I1)
        END DO

!       NOW WE TRANSFORM FROM KUBO'S TO GRASER'S NOTATION

        DO I1 = 1, NITER
          AU(I1) = AUK(I1)
          AV(I1) = AVK(I1) - AJK(I1) / 2.D0
          AJ(I1) = 2.D0 * AJK(I1)
          AJP(I1) = AJPK(I1)
        END DO

!       START INPUT

        IF(IDOP.EQ.-5) OPEN(UNIT=99,FILE='baresu32x32n-05.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        IF(IDOP.EQ.24) OPEN(UNIT=99,FILE='baresu104x104n024.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        IF(IDOP.EQ.35) OPEN(UNIT=99,FILE='baresu158x158n435.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        IF(IDOP.EQ.45) OPEN(UNIT=99,FILE='baresu104x104n045.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
!
        DO K1 = -NGRID2, NGRID2
          DO K2 = -NGRID2, NGRID2
            DO K3 = 1, NORBIT
              DO K4 = 1, NORBIT
                DO K5 = 1, NORBIT
                  DO K6 = 1, NORBIT
                    READ(99)CHIINPUT(K6,K5,K4,K3,K2,K1)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

        CLOSE(99)

!       ZEROING THE CHIOUT() ARRAY

        DO K1 = -NGRID2 , NGRID2
          DO K2 = -NGRID2 , NGRID2
            DO K3 = 1, NORBIT2
              DO K4 = 1, NORBIT2
                CHIOUT(K4,K3,K2,K1) = DCMPLX(0.D0,0.D0)
              END DO
            END DO
          END DO
        END DO

!       ROTATE THE INPUT

        DO K1 = 1, NORBIT
          DO K2 = 1, NORBIT
            INDEX(K1,K2) = (K1-1)*NORBIT + K2
          END DO
        END DO

       DO I1 = -NGRID2, NGRID2         !       KX IN THE GRID
         DO I2 = -NGRID2, NGRID2       !       KY IN THE GRID
            DO I3 = 1, NORBIT            !       THIS INDEXES THE ROW s
              DO I4 = 1, NORBIT          !       THIS INDEXES THE ROW p
                DO I5 = 1, NORBIT        !       THIS INDEXES THE ROW q
                  DO I6 = 1, NORBIT      !       THIS INDEXES THE ROW t
                    L1 = INDEX(I6,I3)
                    L2 = INDEX(I5,I4)
                    CHIOUT(L1,L2,I2,I1) = CHIINPUT(I6,I5,I4,I3,I2,I1)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

        DO I1 = -NGRID2, NGRID2
          DO I2 = -NGRID2, NGRID2
            DO i3 = 1, NORBIT2
              DO i4 = 1, NORBIT2
              IF(DABS(DIMAG(CHIOUT(I4,I3,I2,I1))).GT.TOL) CHIOUT(I4,I3,I2,I1) = DCMPLX(DREAL(CHIOUT(I4,I3,I2,I1)),0.D0)
              end do
            end do
          end do
        end do

!       READING THE K VECTORS FOR THE ELECTRON AND HOLE POCKETS


        IF(IDOP.EQ.-5)THEN
          OPEN(UNIT=1,FILE='pocket00n-05.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=2,FILE='pocketpipin-05.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=3,FILE='pocketpi2pi2n-05.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        ENDIF

        IF(IDOP.EQ.24)THEN
          OPEN(UNIT=1,FILE='pocket00n024.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=2,FILE='pocketpi0n024.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=3,FILE='pocket0pin024.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=4,FILE='pocketpi2pi2n024.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        ENDIF

        IF(IDOP.EQ.35)THEN
          OPEN(UNIT=7,FILE='../pocket00n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=8,FILE='../pocketpi0n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=9,FILE='../pocket0pin035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=10,FILE='../pocketpi2-pi2n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=11,FILE='../pocket-pi2-pi2n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=12,FILE='../pocket-pi2pi2n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=13,FILE='../pocketpi2pi2n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        ENDIF

        IF(IDOP.EQ.45)THEN
          OPEN(UNIT=1,FILE='pocketpipin045.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=2,FILE='pocket00n045-1.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=3,FILE='pocket00n045-2.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        ENDIF

        DO I3 = 1, NPOC
          DO I1 = 0, NPSURF
            READ(I3+6)(AK(I3,I2,I1), I2=1,2)
          END DO
        END DO

        CLOSE(7)
        CLOSE(8)
        CLOSE(9)
        CLOSE(10)
        CLOSE(11)
        CLOSE(12)
        CLOSE(13)
        IF(IDOP.EQ.24.OR.IDOP.EQ.35) CLOSE(4)

!       ZEROING THE COORDINATES ALONG THE AXES

        DO I1 = 1, NPOC
          DO I2 = 0, NPSURF
            IF(DABS(AK(I1,1,I2)).LT.1.D-3) AK(I1,1,I2) = 0.D0
            IF(DABS(AK(I1,2,I2)).LT.1.D-3) AK(I1,2,I2) = 0.D0
          END DO
        END DO

        IF(IDOP.EQ.-5)THEN
          OPEN(UNIT=7,FILE='coeff00n-05.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=8,FILE='coeffpipin-05.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=9,FILE='coeffpi2pi2n-05.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        ENDIF

        IF(IDOP.EQ.24)THEN
          OPEN(UNIT=7,FILE='coeff00n024.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=8,FILE='coeffpi0n024.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=9,FILE='coeff0pin024.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=10,FILE='coeffpi2pi2n024.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        ENDIF

        IF(IDOP.EQ.35)THEN
          OPEN(UNIT=8,FILE='../coeff00n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=9,FILE='../coeffpi0n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=10,FILE='../coeff0pin035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=11,FILE='../coeffpi2-pi2n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=12,FILE='../coeff-pi2-pi2n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=13,FILE='../coeff-pi2pi2n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=14,FILE='../coeffpi2pi2n035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        ENDIF

        IF(IDOP.EQ.45)THEN
          OPEN(UNIT=7,FILE='coeffpipin045.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=8,FILE='coeff00n045-1.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
          OPEN(UNIT=9,FILE='coeff00n045-2.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        ENDIF

        DO K1 = 1, NPOC
          DO K2 = 0, NPSURF
            READ(7+K1)PHIREAD(K1,K2),(AWEIGHT(K1,L1,K2),L1=1,NORBIT)
          END DO
        END DO

        CLOSE(8)
        CLOSE(9)
        CLOSE(10)
        CLOSE(11)
        CLOSE(12)
        CLOSE(13)
        CLOSE(14)
        IF(IDOP.EQ.24.OR.IDOP.EQ.35) CLOSE(10)

        DO K2 = 0, NPSURF
          PHI(K2) = PHIREAD(1,K2)
        END DO

!       INPUT IS OVER

!       THESE MATRICES WILL CONTAIN THE INDEX INTERCHANGE

        K = 0
        DO I1 = 1, NORBIT
          DO I2 = 1, NORBIT
            K = K + 1
            MAT1(1,K) = I1
            MAT1(2,K) = I2
            MAT2(I1,I2) = K
          END DO
        END DO

        IF(IFLAG.EQ.1)THEN
          IF(IDOP.EQ.-5)THEN
            OPEN(UNIT=27,FILE='vecar64sgn-05a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=28,FILE='vecbr64sgn-05a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=29,FILE='veccr64sgn-05a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=30,FILE='vecdr64sgn-05a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=37,FILE='vecar64sgphin-05a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=38,FILE='vecbr64sgphin-05a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=39,FILE='veccr64sgphin-05a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=40,FILE='vecdr64sgphin-05a.dat',STATUS='UNKNOWN')

            OPEN(UNIT=41,FILE='strength64sn-05a.dat',STATUS='UNKNOWN')
          ENDIF

          IF(IDOP.EQ.24)THEN
            OPEN(UNIT=27,FILE='vecar64sgn024a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=28,FILE='vecbr64sgn024a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=29,FILE='veccr64sgn024a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=30,FILE='vecdr64sgn024a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=37,FILE='vecar64sgphin024a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=38,FILE='vecbr64sgphin024a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=39,FILE='veccr64sgphin024a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=40,FILE='vecdr64sgphin024a.dat',STATUS='UNKNOWN')

            OPEN(UNIT=41,FILE='strength64sn024a.dat',STATUS='UNKNOWN')
          ENDIF

          IF(IDOP.EQ.35)THEN
            OPEN(UNIT=27,FILE='vecar32-test.dat',STATUS='UNKNOWN')
            OPEN(UNIT=28,FILE='vecbr32-test.dat',STATUS='UNKNOWN')
            OPEN(UNIT=29,FILE='veccr32-test.dat',STATUS='UNKNOWN')
            OPEN(UNIT=30,FILE='vecdr32-test.dat',STATUS='UNKNOWN')
            OPEN(UNIT=37,FILE='vecar32sgphi-test.dat',STATUS='UNKNOWN')
            OPEN(UNIT=38,FILE='vecbr32sgphi-test.dat',STATUS='UNKNOWN')
            OPEN(UNIT=39,FILE='veccr32sgphi-test.dat',STATUS='UNKNOWN')
            OPEN(UNIT=40,FILE='vecdr32sgphi-test.dat',STATUS='UNKNOWN')

            OPEN(UNIT=41,FILE='strength64sn035a.dat',STATUS='UNKNOWN')
          ENDIF

          IF(IDOP.EQ.45)THEN
            OPEN(UNIT=27,FILE='vecar64sgn045a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=28,FILE='vecbr64sgn045a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=29,FILE='veccr64sgn045a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=30,FILE='vecdr64sgn045a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=37,FILE='vecar64sgphin045a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=38,FILE='vecbr64sgphin045a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=39,FILE='veccr64sgphin045a.dat',STATUS='UNKNOWN')
            OPEN(UNIT=40,FILE='vecdr64sgphin045a.dat',STATUS='UNKNOWN')

            OPEN(UNIT=41,FILE='strength64sn045a.dat',STATUS='UNKNOWN')
          ENDIF
        ENDIF          

        IF(IDOP.EQ.-5)THEN
          OPEN(UNIT=25,FILE='eigr32sgn-05b.dat',STATUS='UNKNOWN')
          OPEN(UNIT=26,FILE='eigi32sgn-05b.dat',STATUS='UNKNOWN')
        ENDIF

        IF(IDOP.EQ.24)THEN
          OPEN(UNIT=25,FILE='eigr104sgn024b.dat',STATUS='UNKNOWN')
          OPEN(UNIT=26,FILE='eigi104sgn024b.dat',STATUS='UNKNOWN')
        ENDIF

        IF(IDOP.EQ.35)THEN
          OPEN(UNIT=25,FILE='eigr32-test.dat',STATUS='UNKNOWN')
          OPEN(UNIT=26,FILE='eigi32-test.dat',STATUS='UNKNOWN')
        ENDIF

        IF(IDOP.EQ.45)THEN
          OPEN(UNIT=25,FILE='eigr104sgn045c.dat',STATUS='UNKNOWN')
          OPEN(UNIT=26,FILE='eigi104sgn045c.dat',STATUS='UNKNOWN')
        ENDIF


        DO M1 = 1 , NITER

!       ZEROING BOTH MATRICES UC AND US

        DO I1 = 1, NORBIT2
          DO I2 = 1, NORBIT2
            UC(I1,I2) = 0.D0
            US(I1,I2) = 0.D0
          END DO
        END DO

!       ZEROING \GAMMA

        DO I1 = 1, NORBIT2
          DO I2 = 1, NORBIT2
            DO I3 = 0, NPSURF
              DO I4 = 1, NPOC
                DO I5 = 0, NPSURF
                  DO I6 = 1, NPOC
                    GAMMA(I6,I5,I4,I3,I2,I1) = DCMPLX(0.D0,0.d0)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

!       CREATING THE MATRIX UC

        DO K1 = 1, NORBIT2
          UC(K1,K1) = 3.D0 * AJ(M1) / 4.D0 - AV(M1)
        END DO

        DO K1 = 1, NORBIT
          L1 = (NORBIT+1) * K1 - NORBIT
          UC(L1,L1) = AU(M1)
        END DO

        DO K1 = 1, NORBIT
          L1 = (NORBIT+1) * K1 - NORBIT
          DO K2 = 1, NORBIT
            L2 = (NORBIT+1) * K2 - NORBIT
            IF(L1.NE.L2) UC(L1,L2) = 2.D0 * AV(M1)
          END DO
        END DO

        UC(2,6) = AJP(M1)
        UC(3,11) = AJP(M1)
        UC(4,16) = AJP(M1)
        UC(5,21) = AJP(M1)
        UC(8,12) = AJP(M1)
        UC(9,17) = AJP(M1)
        UC(10,22) = AJP(M1)
        UC(14,18) = AJP(M1)
        UC(15,23) = AJP(M1)
        UC(20,24) = AJP(M1)

        UC(6,2) = AJP(M1)
        UC(11,3) = AJP(M1)
        UC(16,4) = AJP(M1)
        UC(21,5) = AJP(M1)
        UC(12,8) = AJP(M1)
        UC(17,9) = AJP(M1)
        UC(22,10) = AJP(M1)
        UC(18,14) = AJP(M1)
        UC(23,15) = AJP(M1)
        UC(24,20) = AJP(M1)

!       CREATING THE MATRIX US

        DO K1 = 1, NORBIT2
          US(K1,K1) = AJ(M1) / 4.D0 + AV(M1)
        END DO

        DO K1 = 1, NORBIT
          L1 = (NORBIT+1) * K1 - NORBIT
          US(L1,L1) =  AU(M1)
        END DO

        DO K1 = 1, NORBIT
          L1 = (NORBIT+1) * K1 - NORBIT
          DO K2 = 1, NORBIT
            L2 = (NORBIT+1) * K2 - NORBIT
            IF(L1.NE.L2) US(L1,L2) = AJ(M1) / 2.d0
          END DO
        END DO

        US(2,6) = AJP(M1)
        US(3,11) = AJP(M1)
        US(4,16) = AJP(M1)
        US(5,21) = AJP(M1)
        US(8,12) = AJP(M1)
        US(9,17) = AJP(M1)
        US(10,22) = AJP(M1)
        US(14,18) = AJP(M1)
        US(15,23) = AJP(M1)
        US(20,24) = AJP(M1)

        US(6,2) = AJP(M1)
        US(11,3) = AJP(M1)
        US(16,4) = AJP(M1)
        US(21,5) = AJP(M1)
        US(12,8) = AJP(M1)
        US(17,9) = AJP(M1)
        US(22,10) = AJP(M1)
        US(18,14) = AJP(M1)
        US(23,15) = AJP(M1)
        US(24,20) = AJP(M1)

        !OPEN(UNIT=1,FILE='gbaresuNxNnD-plus.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')

        DO I1 = 1, NPOC    !       THIS INDEXES THE POCKETS
          DO I2 = 1, NPOC  !       THIS INDEXES THE POCKETS
            DO I3 = 0, NPSURF       !       A POINT IN FS 1
              DO I4 = 0, NPSURF     !       A POINT IN FS 2
        !DO I1 = 1, NPOC
          !DO I3 = 0, NPSURF
            ! Ajuste para incluir NTOT2/2 nos índices
            !indice1 = (I1 - 1) * (NPSURF + 1) + I3 + 1
            !DO I2 = 1, NPOC
              !DO I4 = 0, NPSURF
                ! Ajuste para incluir NTOT2/2 nos índices
                !indice2 = (I2 - 1) * (NPSURF + 1) + I4 + 1

                !PRINT*, indice1, indice2

!       CALCULATE THE VALUE OF QX AND QY

                QX = AK(I1,1,I3) - AK(I2,1,I4)
                QY = AK(I1,2,I3) - AK(I2,2,I4)

!       CHECK THAT THE VECTORS ARE INSIDE THE GRID

                IF(QX.GT.PI) QX = QX - 2.D0 * PI
                IF(QX.LT.-PI) QX = QX + 2.D0 * PI
                IF(QY.GT.PI) QY = QY - 2.D0 * PI
                IF(QY.LT.-PI) QY = QY + 2.D0 * PI

!       SAVE Q POINTS
                !SAVEQX(INDICE1,INDICE2) = QX
                !SAVEQY(INDICE1,INDICE2) = QY

                DO KK1 = -NGRID2, NGRID2
                  DO KK2 = -NGRID2, NGRID2
                    DO KK3 = 1, NORBIT2
                      DO KK4 = 1, NORBIT2
                        CHIRE(KK4,KK3,KK2,KK1) = DREAL(CHIOUT(KK4,KK3,KK2,KK1))
                      END DO
                    END DO
                  END DO
                END DO

                CALL SUSBICUB(QX,QY,CHIRE,CHIREINT)    !       INTERPOLATES THE BARE SUSCETIBILITY

                DO KK1 = -NGRID2, NGRID2
                  DO KK2 = -NGRID2, NGRID2
                    DO KK3 = 1, NORBIT2
                      DO KK4 = 1, NORBIT2
                        CHIIMG(KK4,KK3,KK2,KK1) = DIMAG(CHIOUT(KK4,KK3,KK2,KK1))
                      END DO
                    END DO 
                  END DO
                END DO

                CALL SUSBICUB(QX,QY,CHIIMG,CHIIMGINT)    !  INTERPOLATES THE BARE SUSCETIBILITY

                DO KK1 = 1, NORBIT2
                  DO KK2 = 1, NORBIT2
                    CHI(KK2,KK1) = DCMPLX(CHIREINT(KK2,KK1),CHIIMGINT(KK2,KK1))
                    CHISAVE(indice1, indice2, KK1, KK2) = CHI(KK2,KK1)
                  END DO
                END DO

!       SAVE \CHI
                !DO KK1 = 1, NORBIT2
                !  DO KK2 = 1, NORBIT2
                !    WRITE(1)CHISAVE(indice1, indice2, KK1, KK2)
                !  END DO
                !END DO

                DO L3 = 1, NORBIT2
                  DO L4 = 1, NORBIT2
                    ATOINVC(L4,L3) = DCMPLX(0.D0,0.D0)
                    ATOINVS(L4,L3) = DCMPLX(0.D0,0.D0)
                    RPACHIC(L4,L3) = DCMPLX(0.D0,0.D0)
                    RPACHIS(L4,L3) = DCMPLX(0.D0,0.D0)
                    CPROD(L4,L3) = DCMPLX(0.D0,0.D0)
                    SPROD(L4,L3) = DCMPLX(0.D0,0.D0)
                    AUX1(L4,L3) = DCMPLX(0.D0,0.D0)
                    AUX2(L4,L3) = DCMPLX(0.D0,0.D0)
                    AUX3(L4,L3) = DCMPLX(0.D0,0.D0)
                    AUX4(L4,L3) = DCMPLX(0.D0,0.D0)
                    AUX5(L4,L3) = DCMPLX(0.D0,0.D0)
                  END DO
                END DO

!               THIS CALCULATES THE DRESSED SUSCEPTIBILITY 
!               FOR THE CHARGE AND FOR THE SPIN
!               THE FINAL RESULT IS STORED IN 
!               1) RPACHIC(4,4) FOR CHARGE
!               2) RPACHIS(4,4) FOR SPIN

                DO M3 = 1, NORBIT2
                  DO M4 = 1, NORBIT2
                    DO K1 = 1, NORBIT2
                      ATOINVC(M3,M4)=ATOINVC(M3,M4)+UC(M3,K1)*CHI(K1,M4)
                      ATOINVS(M3,M4)=ATOINVS(M3,M4)-US(M3,K1)*CHI(K1,M4)
                    END DO
                  END DO
                END DO

                !DO M3 = 1, NORBIT2
                !  WRITE(*,*)(real(CHI(M4,M3)), M4=1,NORBIT2)
                !END DO
                !STOP

                DO M3 = 1, NORBIT2
                  ATOINVC(M3,M3) = 1.D0 + ATOINVC(M3,M3)
                  ATOINVS(M3,M3) = 1.D0 + ATOINVS(M3,M3)
                END DO

                !DO M3 = 1, NORBIT2
                !    WRITE(*,*)(ATOINVS(M3,M4),M4=1,NORBIT2)
                !END DO
                !STOP

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

                DO M3 = 1, NORBIT2
                  DO M4 = 1, NORBIT2
                    DO K1 = 1, NORBIT2
                      RPACHIC(M3,M4)=RPACHIC(M3,M4)+CHI(M3,K1)*ATOINVC(K1,M4)
                      RPACHIS(M3,M4)=RPACHIS(M3,M4)+CHI(M3,K1)*ATOINVS(K1,M4)
                    END DO
!                   write(*,*)RPACHIC(M3,M4)
!                   write(*,*)RPACHIS(M3,M4)
                  END DO
                END DO
!               stop

                do ll1 = 1, NORBIT2
                  do ll2 = 1, NORBIT2
                    if(dabs(dimag(RPACHIC(ll2,ll1))).gt.tol) RPACHIC(ll2,ll1) = dcmplx(dreal(RPACHIC(ll2,ll1)),0.d0)
                  end do
                end do

                do ll1 = 1, NORBIT2
                  do ll2 = 1, NORBIT2
                    if(dabs(dimag(RPACHIS(ll2,ll1))).gt.tol) RPACHIS(ll2,ll1) = dcmplx(dreal(RPACHIS(ll2,ll1)),0.d0)
                  end do
                end do

!       STARTING THE MATRIX PRODUCTS THAT DEFINE \GAMMA

                DO M3 = 1, NORBIT2
                  DO M4 = 1, NORBIT2
                    DO K1 = 1, NORBIT2
                      AUX1(M4,M3)=AUX1(M4,M3)+US(M4,K1)*RPACHIS(K1,M3)
                      AUX2(M4,M3)=AUX2(M4,M3)+UC(M4,K1)*RPACHIC(K1,M3)
                    END DO
                  END DO
                END DO

                DO M3 = 1, NORBIT2
                  DO M4 = 1, NORBIT2
                    DO K1 = 1, NORBIT2
                      AUX3(M4,M3)=AUX3(M4,M3)+3.D0*AUX1(M4,K1)*US(K1,M3)
                      AUX4(M4,M3)=AUX4(M4,M3)-AUX2(M4,K1)*UC(K1,M3)
                    END DO
                  END DO
                END DO

                DO M3 = 1, NORBIT2
                  DO M4 = 1, NORBIT2
                    AUX5(M4,M3)=(AUX3(M4,M3)+US(M4,M3)+AUX4(M4,M3)+UC(M4,M3))/2.D0
                    NP = MAT1(1,M4)
                    NS = MAT1(2,M4)
                    NT = MAT1(1,M3)
                    NQ = MAT1(2,M3)
                    GAMMA(I1,I3,I2,I4,MAT2(NS,NT),MAT2(NP,NQ))=AUX5(M4,M3)
                  END DO      !       M4
                END DO        !       M3
              END DO            !       I4      A POINT IN FS 2
            END DO              !       I3      A POINT IN FS 1
          END DO                !       I2      THIS INDEXES THE POCKETS
        END DO                  !       I1      THIS INDEXES THE POCKETS

        !CLOSE(1)
!stop
        !OPEN(UNIT=2,FILE='qx_k-kprime.txt',STATUS='UNKNOWN')
        !OPEN(UNIT=3,FILE='qy_k-kprime.txt',STATUS='UNKNOWN')
        !  DO LL1 = -NTOT2/2, NTOT2/2
        !    WRITE(2,*)(SAVEQX(LL1,LL2) ,LL2= -NTOT2/2, NTOT2/2)
        !    WRITE(3,*)(SAVEQY(LL1,LL2) ,LL2= -NTOT2/2, NTOT2/2)
        !  END DO
        !CLOSE(2)
        !CLOSE(3)

!        open(unit=1,file='gamma1-real.dat',status='unknown')
!       do ll1 = 1, NORBIT2
!          write(1,*)(dreal(GAMMA(1,0,1,0,ll1,ll2)),ll2=1,norbit2)
!       end do
!               write(1,*)'========================================'
!111     format(25(f7.3,1x))
!        close(1)
!       stop
!       NOW THAT WE HAVE \GAMMA, AS CALCULATED ABOVE, LET'S DEFINE A \GAMMA 
!       THAT DEPENDS ON THE FERMI SURFACE POCKET AND THE POINT IN IT

        DO I1 = 1, NPOC
          DO I2 = 0, NPSURF
            indice1 = (I1 - 1) * (NPSURF + 1) + I2 + 1
            DO I3 = 1, NPOC
              DO I4 = 0, NPSURF
                indice2 = (I3 - 1) * (NPSURF + 1) + I4 + 1
                GAMMAIJ(I1,I2,I3,I4) = DCMPLX(0.D0,0.D0)
!           GIVEN k (I2), FIND -k (J2) (POCKETS 1 TO 3)
                IF(I1.LT.4) THEN 
                  J2 = I2 + (NPSURF+1)/2
                  IF(J2.GT.NPSURF) J2 = J2 - (NPSURF + 1)
                END IF
!           GIVEN k' (I4), FIND -k' (J4) (POCKETS 1 TO 3)
                IF(I3.LT.4) THEN
                  J4 = I4 + (NPSURF+1)/2
                  IF(J4.GT.NPSURF) J4 = J4 - (NPSURF + 1)
                END IF
!           GIVEN k (I2), FIND -k (J2) (POCKETS 4 TO 7)
                IF(I1.GT.3)THEN
                  IF(I1.EQ.4.OR.I1.EQ.5) K1 = I1 + 2
                  IF(I1.EQ.6.OR.I1.EQ.7) K1 = I1 - 2
                  J2 = I2
                END IF
!           GIVEN k' (I4), FIND -k' (J4) (POCKETS 4 TO 7)
                IF(I3.GT.3)THEN
                  IF(I3.EQ.4.OR.I3.EQ.5) K3 = I3 + 2
                  IF(I3.EQ.6.OR.I3.EQ.7) K3 = I3 - 2
                  J4 = I4
                END IF
                DO I5 = 1, NORBIT
                  DO I6 = 1, NORBIT
                    DO I7 = 1, NORBIT
                      DO I8 = 1, NORBIT
                        IF(I1.LT.4)THEN ! (POCKETS 1 TO 3)
                          AAA1 = DCONJG(AWEIGHT(I1,I5,J2))*DCONJG(AWEIGHT(I1,I6,I2))
                        ELSE ! (POCKETS 4 TO 7)
                          AAA1 = DCONJG(AWEIGHT(K1,I5,J2))*DCONJG(AWEIGHT(I1,I6,I2))
                        END IF
                        IF(I3.LT.4)THEN ! (POCKETS 1 TO 3)
                          BBB1 = AWEIGHT(I3,I8,I4)*AWEIGHT(I3,I7,J4)
                        ELSE ! (POCKETS 4 TO 7)
                          BBB1 = AWEIGHT(I3,I8,I4)*AWEIGHT(K3,I7,J4)
                        ENDIF
                        GAMMAIJ(I1,I2,I3,I4) = GAMMAIJ(I1,I2,I3,I4)+AAA1*BBB1*DREAL(GAMMA(I1,I2,I3,I4,MAT2(I6,I5),MAT2(I8,I7)))
                        GAMMA2SAVE(indice1,indice2) = GAMMAIJ(I1,I2,I3,I4)
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
        !open(unit=1,file='gamma2.0-real.dat',status='unknown')
        !      do ll1 = 1,NTOT2
        !        write(1,*)(dreal(GAMMA2SAVE(ll1,ll2)),ll2=1,NTOT2)
        !      end do
        !close(1)
        !stop
        DO I1 = 1, NPOC
          DO I2 = 0, NPSURF
            indice1 = (I1 - 1) * (NPSURF + 1) + I2 + 1
            DO I3 = 1, NPOC
              DO I4 = 0, NPSURF
                indice2 = (I3 - 1) * (NPSURF + 1) + I4 + 1
!           GIVEN k' (I4), FIND -k' (J4) (POCKETS 1 TO 3)
                IF(I3.LT.4) THEN
                  J4 = I4 + (NPSURF+1)/2
                  IF(J4.GT.NPSURF) J4 = J4 - (NPSURF + 1)
                END IF
!           GIVEN k' (I4), FIND -k' (J4) (POCKETS 4 TO 7)
                IF(I3.GT.3)THEN
                  IF(I3.EQ.4.OR.I3.EQ.5) K3 = I3 + 2
                  IF(I3.EQ.6.OR.I3.EQ.7) K3 = I3 - 2
                  J4 = I4
                END IF
                IF(I3.LT.4)THEN ! (POCKETS 1 TO 3)
                  GAMMAIJSYM(I1,I2,I3,I4) = (GAMMAIJ(I1,I2,I3,I4) + GAMMAIJ(I1,I2,I3,J4))/2.D0
                ELSE ! (POCKETS 4 TO 7)
                  GAMMAIJSYM(I1,I2,I3,I4) = (GAMMAIJ(I1,I2,I3,I4) + GAMMAIJ(I1,I2,K3,J4))/2.D0
                END IF

                !GAMMA2SAVE(indice1,indice2) = GAMMAIJSYM(I1,I2,I3,I4)
              END DO
            END DO
          END DO
        END DO

!       open(unit=1,file='gamma2-real.dat',status='unknown')
!       open(unit=2,file='gamma2-imag.dat',status='unknown')
!        do ll1 = -NTOT2/2,NTOT2/2-1
!          write(1,*)(dreal(GAMMA2SAVE(ll1,ll2)),ll2=-NTOT2/2,NTOT2/2-1)
!         write(2,*)(dimag(GAMMAIJSYM(ll1,ll2)),ll2=1,NTOT2)
!        end do
!              write(1,*)'========================================'
!              write(2,*)'========================================'
!111     format(25(f7.3,1x))
!       stop

        CALL AUTOVALOR(GAMMAIJSYM,EIGENR,EIGENI,VECR,VECI)

        WRITE(25,100)AU(M1),(EIGENR(LLL1),LLL1=1, 10)
        WRITE(26,100)AU(M1),(EIGENI(LLL1),LLL1=1, 10)

        IF(IFLAG.EQ.1)THEN     
          N2 = 0
          DO I1 = 1, NPOC
            DO I2 = 0, NPSURF
              N2 = N2 + 1
              WRITE(27,*)AK(I1,1,I2),AK(I1,2,I2),VECR(N2,1)
              WRITE(28,*)AK(I1,1,I2),AK(I1,2,I2),VECR(N2,2)
              WRITE(29,*)AK(I1,1,I2),AK(I1,2,I2),VECR(N2,3)
              WRITE(30,*)AK(I1,1,I2),AK(I1,2,I2),VECR(N2,4)
              WRITE(37,*)DBLE(I1-1)*2.D0*PI+PHI(I2),VECR(N2,1),VECI(N2,1)
              WRITE(38,*)DBLE(I1-1)*2.D0*PI+PHI(I2),VECR(N2,2),VECI(N2,2)
              WRITE(39,*)DBLE(I1-1)*2.D0*PI+PHI(I2),VECR(N2,3),VECI(N2,3)
              WRITE(40,*)DBLE(I1-1)*2.D0*PI+PHI(I2),VECR(N2,4),VECI(N2,4)
            END DO
          END DO
        END IF


        END DO          !       M1(EXTERNAL DO FOR THE U VALUES)

        CLOSE(25)
        CLOSE(26)
        IF(IFLAG.EQ.1)THEN     
          CLOSE(27)
          CLOSE(28)
          CLOSE(29)
          CLOSE(30)
          CLOSE(37)
          CLOSE(38)
          CLOSE(39)
          CLOSE(40)
          CLOSE(41)
        END IF

100     FORMAT(11(E13.7,1X))
200     FORMAT(3(E13.7,1X))
300     FORMAT(2(I2,1X),2(E13.7,1X))

END PROGRAM EIGENSG

SUBROUTINE AUTOVALOR(GAMMA,EIGENR,EIGENI,VECR,VECI)

        USE INTERFACE 

        IMPLICIT NONE

        COMPLEX (KIND=8), DIMENSION(NPOC,0:NPSURF,NPOC,0:NPSURF), INTENT(IN) :: GAMMA     !       SYMMETRIC VERTEX FUNCTION
        REAL (KIND=8), DIMENSION(NTOT2), INTENT(OUT) :: EIGENR    !       EIGENVALUES (REAL PART)
        REAL (KIND=8), DIMENSION(NTOT2), INTENT(OUT) :: EIGENI    !       EIGENVALUES (IMAGINARY PART)
        REAL (KIND=8), DIMENSION(NTOT2,NTOT2), INTENT(OUT) :: VECR    !       EIGENVECTORS (REAL PART)
        REAL (KIND=8), DIMENSION(NTOT2,NTOT2), INTENT(OUT) :: VECI    !       EIGENVECTORS (IMAGINARY PART)

        COMPLEX (KIND=8), DIMENSION(NTOT2,NTOT2) :: TDIAGR  !       MATRIX TO BE DIAGONALIZED
        COMPLEX (KIND=8), DIMENSION(-NTOT2/2:NTOT2/2-1,-NTOT2/2:NTOT2/2-1) :: GAMMA3SAVE
        REAL (KIND=8), DIMENSION(NPOC,2) :: AKCENTER

        INTEGER :: K2,L1,KC,I1,I2,I3,I4,KR, ll1,ll2, indice1, indice2
        REAL (KIND=8) :: AKMOD
        REAL (KIND=8) :: DK
        REAL (KIND=8) :: TEMPDK

        IF(IDOP.EQ.-5) OPEN(UNIT=12,FILE='fermi-veln-05.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        IF(IDOP.EQ.24) OPEN(UNIT=12,FILE='fermi-veln024.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        IF(IDOP.EQ.35) OPEN(UNIT=12,FILE='../fermi-veln035.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
        IF(IDOP.EQ.45) OPEN(UNIT=12,FILE='fermi-veln045.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')

        DO K2 = 0, NPSURF
          READ(12)(VFERMI(L1,K2),L1=1,NPOC)
        END DO

        CLOSE(12)

        IF(IDOP.EQ.-5) THEN
          AKCENTER(1,1) = PI
          AKCENTER(1,2) = PI
          AKCENTER(2,1) = 0.D0
          AKCENTER(2,2) = 0.D0
          AKCENTER(3,1) = PI/2.D0
          AKCENTER(3,2) = PI/2.D0
        ENDIF

        IF(IDOP.EQ.24.OR.IDOP.EQ.35) THEN
          AKCENTER(1,1) = 0.D0
          AKCENTER(1,2) = 0.D0
          AKCENTER(2,1) = PI
          AKCENTER(2,2) = 0.D0
          AKCENTER(3,1) = 0.D0
          AKCENTER(3,2) = PI
          AKCENTER(4,1) = PI/2.D0
          AKCENTER(4,2) = -PI/2.D0
          AKCENTER(5,1) = -PI/2.D0
          AKCENTER(5,2) = -PI/2.D0
          AKCENTER(6,1) = -PI/2.D0
          AKCENTER(6,2) = PI/2.D0
          AKCENTER(7,1) = PI/2.D0
          AKCENTER(7,2) = PI/2.D0
        ENDIF

        IF(IDOP.EQ.45) THEN
          AKCENTER(1,1) = PI
          AKCENTER(1,2) = PI
          AKCENTER(2,1) = 0.D0
          AKCENTER(2,2) = 0.D0
          AKCENTER(3,1) = 0.D0
          AKCENTER(3,2) = 0.D0
        ENDIF

        KC = 0
        DO I3 = 1, NPOC
          DO I4 = 0, NPSURF
            indice2 = (I3 - 1) * (NPSURF + 1) + I4 + 1
            KC = KC + 1
            AKMOD=DSQRT((AK(I3,1,I4)-AKCENTER(I3,1))**2.D0+(AK(I3,2,I4)-AKCENTER(I3,2))**2.D0)
            ! REGULARIZAÇÃO DO AKMOD DEVIDO AO LIMITE DA ZB
            IF (I3 == 2) THEN
              IF (AKMOD > 2) THEN
                AKMOD=DSQRT((AK(I3,1,I4)-(-PI))**2.D0+(AK(I3,2,I4)-(0.D0))**2.D0)
              END IF
            END IF
            IF (I3 == 3) THEN
              IF (AKMOD > 2) THEN
                AKMOD=DSQRT((AK(I3,1,I4)-(0.d0))**2.D0+(AK(I3,2,I4)-(-PI))**2.D0)
              END IF
            END IF
            KR = 0
            DO I1 = 1, NPOC
              DO I2 = 0, NPSURF
                indice1 = (I1 - 1) * (NPSURF + 1) + I2 + 1
                KR = KR + 1
                TDIAGR(KR,KC)=-GAMMA(I1,I2,I3,I4)/4.D0/PI/PI/VFERMI(I3,I4)!*DELTAFI*AKMOD
                
                !GAMMA3SAVE(indice1,indice2) = TDIAGR(KR, KC)
              END DO
            END DO
          END DO
        END DO

!        open(unit=99,file='gamma3-real.dat',status='unknown')
!       open(unit=2,file='gamma3-imag.dat',status='unknown')

!        do ll1 = -NTOT2/2,NTOT2/2-1
!          write(99,*)(dreal(GAMMA3SAVE(ll1,ll2)),ll2=-NTOT2/2,NTOT2/2-1)
!                write(2,111)(dimag(TDIAGR(ll2,ll1)),ll2=-NPP2, NPP2-1)
!          end do
!              write(99,*)'========================================'
!              write(2,*)'========================================'
!111     format(25(f6.3,1x))
!       stop
!       close(99)

        CALL DIAGONAL(TDIAGR,EIGENR,EIGENI,VECR,VECI)

END SUBROUTINE AUTOVALOR

SUBROUTINE DIAGONAL(TBR,EIGENRMAX,EIGENIMAX,VECRMAX,VECIMAX)

        USE INTERFACE

        IMPLICIT NONE

        COMPLEX (KIND=8), DIMENSION(NTOT2,NTOT2), INTENT(IN) :: TBR  !       MATRIX TO BE DIAGONALIZED

        REAL (KIND=8) :: AMINIMO

        COMPLEX (KIND=8), DIMENSION(NTOT2) :: EIGEN    !       EIGENVALUES FOR GAP FUNCTION 
        COMPLEX (KIND=8), DIMENSION(1,NTOT2) :: VL    !       NOT USED
        COMPLEX (KIND=8), DIMENSION(NTOT2,NTOT2) :: VR    !       RIGHT EIGENVECTORS
        COMPLEX (KIND=8), DIMENSION(3*NTOT2) :: WORK    !       EIGENVECTORS
        REAL (KIND=8), DIMENSION(2*NTOT2) :: RWORK    !       EIGENVALUES

        REAL (KIND=8), DIMENSION(NTOT2), INTENT(OUT) :: EIGENRMAX    !       EIGENVALUES (REAL PART)
        REAL (KIND=8), DIMENSION(NTOT2), INTENT(OUT) :: EIGENIMAX    !       EIGENVALUES (IMAGINARY PART)
        REAL (KIND=8), DIMENSION(NTOT2,NTOT2), INTENT(OUT) :: VECRMAX    !       EIGENVECTORS (REAL PART)
        REAL (KIND=8), DIMENSION(NTOT2,NTOT2), INTENT(OUT) :: VECIMAX    !       EIGENVEVTORS (IMAGINARY PART)

        CHARACTER(LEN=1) :: JOBVL,JOBVR

        INTEGER :: INFO
        INTEGER :: I1,I2,LL
        INTEGER :: IMAX

        JOBVL = 'N'
        JOBVR = 'N'
        IF(IFLAG.EQ.1) JOBVR = 'V'

        CALL ZGEEV(JOBVL,JOBVR,NTOT2,TBR,NTOT2,EIGEN,VL,1,VR,NTOT2,WORK,3*NTOT2,RWORK,INFO)

        IF(INFO.NE.0)THEN
          WRITE(*,*)'SOMETHING WRONG HERE: ZGEEV  ',INFO
          STOP
        ENDIF

        DO I1 = 1, NTOT2
          AMINIMO = -1.D6
          DO I2 = 1, NTOT2
            IF(DREAL(EIGEN(I2)).GT.AMINIMO) THEN
              AMINIMO = DREAL(EIGEN(I2))
              IMAX = I2
            ENDIF
          END DO
          EIGENRMAX(I1) = DREAL(EIGEN(IMAX))
          EIGENIMAX(I1) = DIMAG(EIGEN(IMAX))
          IF(IFLAG.EQ.1) THEN
            DO LL = 1, NTOT2
              VECRMAX(LL,I1) = DREAL(VR(LL,IMAX))
              VECIMAX(LL,I1) = DIMAG(VR(LL,IMAX))
            END DO
          ENDIF
          EIGEN(IMAX) = DCMPLX(-1.D7,0.D0)
        END DO

END SUBROUTINE DIAGONAL
