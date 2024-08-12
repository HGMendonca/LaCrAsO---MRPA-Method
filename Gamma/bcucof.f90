SUBROUTINE bcucof(y,y1,y2,y12,d1,d2,c)

      IMPLICIT NONE 

        REAL (KIND=8), INTENT(IN) :: d1,d2
        REAL (KIND=8), DIMENSION(4), INTENT(IN) :: y,y1,y12,y2
        REAL (KIND=8), DIMENSION(4,4), INTENT(OUT) :: c

        INTEGER :: i,j,k,l
        REAL*8 d1d2,xx,cl(16),wt(16,16),x(16)
        SAVE wt
        DATA wt/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4,10*&
        0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4,&
        1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,&
        -6,4,2*0,3,-2,0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,&
        10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4,&
        -2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,&
        2,-2,2*0,-1,1/
        
        d1d2=d1*d2
!     write(*,*)(y(i),i=1,4)
!     write(*,*)(y1(i),i=1,4)
!     write(*,*)(y2(i),i=1,4)
!     write(*,*)(y12(i),i=1,4)
      do 11 i=1,4
        x(i)=y(i)
        x(i+4)=y1(i)*d1
        x(i+8)=y2(i)*d2
        x(i+12)=y12(i)*d1d2
11    continue
      do 13 i=1,16
        xx=0.
        do 12 k=1,16
          xx=xx+wt(i,k)*x(k)
12      continue
        cl(i)=xx
13    continue
      l=0
      do 15 i=1,4
        do 14 j=1,4
          l=l+1
          c(i,j)=cl(l)
14      continue
15    continue

END SUBROUTINE bcucof

