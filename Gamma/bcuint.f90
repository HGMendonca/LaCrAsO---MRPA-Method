SUBROUTINE bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,ansy1,ansy2)

        IMPLICIT NONE

        REAL (KIND=8), INTENT(OUT) :: ansy,ansy1,ansy2
        REAL (KIND=8), INTENT(IN) :: x1,x1l,x1u,x2,x2l,x2u
        REAL (KIND=8), DIMENSION(4), INTENT(IN) :: y,y1,y12,y2

!    USES bcucof

      INTEGER :: i
      REAL (KIND=8) :: t,u
      REAL (KIND=8), DIMENSION(4,4) :: c

!     write(*,*)(y(i),i=1,4)
!     write(*,*)(y1(i),i=1,4)
!     write(*,*)(y2(i),i=1,4)
!     write(*,*)(y12(i),i=1,4)
!     write(*,*)x1l,x1u,x2l,x2u

      call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)

      if(x1u.eq.x1l.or.x2u.eq.x2l)pause 'bad input in bcuint'
      t=(x1-x1l)/(x1u-x1l)
      u=(x2-x2l)/(x2u-x2l)
      ansy=0.D0
      ansy2=0.D0
      ansy1=0.D0
      do 11 i=4,1,-1
        ansy=t*ansy+((c(i,4)*u+c(i,3))*u+c(i,2))*u+c(i,1)
        ansy2=t*ansy2+(3.*c(i,4)*u+2.*c(i,3))*u+c(i,2)
        ansy1=u*ansy1+(3.*c(4,i)*t+2.*c(3,i))*t+c(2,i)
11    continue
      ansy1=ansy1/(x1u-x1l)
      ansy2=ansy2/(x2u-x2l)

END SUBROUTINE bcuint
