SUBROUTINE CHEMPOT
 
        USE INTERFACE

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(20) :: CHEM

        CHEM(11) = 0.1139217937D2  ! -0.5
        CHEM(10) = -0.8294944743D0  ! -0.4
        CHEM(9) = -0.7398260883D0  ! -0.2
        CHEM(8) = -0.6762833072D0  ! 0.0
        CHEM(7) =  11.539008998897742  ! 0.25
        CHEM(6) =  0.1154033177D2  ! 0.35
        CHEM(5) =  0.1159516784D2  ! 0.65
        CHEM(4) =  0.1156385487D2  ! 0.5
        CHEM(3) =  0.1157356590D2  ! 0.55
        CHEM(2) =  0.1155521928D2  ! 0.45
        CHEM(1) =  0.1154824675D2  ! 0.4

        IF(AN.EQ.-0.5D0) AMU = CHEM(11)
        IF(AN.EQ.-0.4D0) AMU = CHEM(10)
        IF(AN.EQ.-0.2D0) AMU = CHEM(9)
        IF(AN.EQ.0.D0) AMU = CHEM(8)
        IF(AN.EQ.0.25D0) AMU = CHEM(7)
        IF(AN.EQ.0.35D0) AMU = CHEM(6)
        IF(AN.EQ.0.65D0) AMU = CHEM(5)
        IF(AN.EQ.0.5D0) AMU = CHEM(4)
        IF(AN.EQ.0.55D0) AMU = CHEM(3)
        IF(AN.EQ.0.45D0) AMU = CHEM(2)
        IF(AN.EQ.0.4D0) AMU = CHEM(1)
        
END SUBROUTINE CHEMPOT
