!=======================================================================
      FUNCTION FUNCBODY(X,Y,Z,T)
!=======================================================================
!
!     Function for an immersed body
!
!     Required condition to use secant method :
!        1. FUNCBODY(X,Y,Z)<0  inside the body
!           FUNCBODY(X,Y,Z)=0  at the body
!           FUNCBODY(X,Y,Z)>0  outside the body
!     Ex.
!        FUNCBODY=X**2+Y**2+Z**2-0.5**2       ! for sphere
!        FUNCBODY=X**2+Y**2-0.5**2            ! for infinite cylinder
!
!     Parameter T corresponds to non-dimensional time.
!     For static-body problems, T is set zero.
!-----------------------------------------------------------------------
      IMPLICIT NONE

      REAL*8     :: X,Y,Z
      REAL*8     :: T
      REAL*8     :: FUNCBODY

      REAL*8     :: R,P,PI,TT,RR
      REAL*8     :: YTMP, ZTMP

      PI = ACOS(-1.)

      ! THICK WALL CHANNEL
      ! IF (ABS(Y) .GE. 0.5) THEN
      !   FUNCBODY = -1.
      ! ELSE
      !   FUNCBODY = 1.
      ! ENDIF
 
      ! ! UNIT CIRCULAR SECTION PIPE
      ! IF (-Y**2.-Z**2.+0.5D0**2. .GT. 0.D0) THEN
      !   FUNCBODY = 1.
      ! ELSE
      !   FUNCBODY = -1.
      ! ENDIF

      ! ! RIBS
      ! IF (FUNCBODY .GT. 0.) THEN
      !   R = (Y**2.+Z**2.)**0.5D0
      !   P = ATAN2(Y, Z) / PI * 180. ! DEGREES
      !   IF (P .LT. 0) P = P + 360.
      !   P = MOD(P,6.)
      !   IF (P .GT. 3.) P = 6. - P

      !   TT = ATAN2(.5*SIN(3./180.*PI),.5*COS(3./180.*PI)-.5 + 6./127.) / PI * 180.
      !   RR = (.5 - 6./127.)*SIN(TT*PI/180.)/SIN((TT-P)*PI/180.)
      !   IF (R .GT. RR) FUNCBODY = -1.
      ! ENDIF

      P  = 1.02746*2 / 16.
      RR = P / 4. * 3.**0.5D0
      FUNCBODY = Y - 4.*RR/P*ABS(MOD(MOD(Z-P/4.,P)+P,P) - P/2.) + 1.


      RETURN
      END FUNCTION FUNCBODY

! COORD. ROTATING MODULES FOR ROTATIONERY BODY, MADE BY SANGJOON LEE
!=======================================================================
       SUBROUTINE ROTATE_ZAXIS(XR,YR,ZR,THETA)
!=======================================================================
!     ROTATE X AND Y COORDINATE TO COUNTERCLOCKWISE DIRECTION 
!     THETA IN DEGREE UNITS
      IMPLICIT NONE
      REAL*8   , intent(in)     :: THETA
      REAL*8   , intent(inout)  :: XR,YR,ZR
      REAL*8                    :: PI
      REAL*8                    :: XTEMPO,YTEMPO

      PI = ACOS(-1.)

      XTEMPO=XR
      YTEMPO=YR

      XR=COS(THETA*PI/180.)*XTEMPO-SIN(THETA*PI/180.)&
         *YTEMPO
      YR=SIN(THETA*PI/180.)*XTEMPO+COS(THETA*PI/180.)&
         *YTEMPO
      ZR=ZR

      RETURN
      END

!=======================================================================
       SUBROUTINE ROTATE_YAXIS(XR,YR,ZR,THETA)
!=======================================================================
!     ROTATE X AND Z COORDINATE TO COUNTERCLOCKWISE DIRECTION
!     THETA IN DEGREE UNITS
      IMPLICIT NONE
      REAL*8   , intent(in)     :: THETA
      REAL*8   , intent(inout)  :: XR,YR,ZR
      REAL*8                    :: PI
      REAL*8                    :: XTEMPO,ZTEMPO

      PI = ACOS(-1.)

      XTEMPO=XR
      ZTEMPO=ZR

      XR=COS(THETA*PI/180.)*XTEMPO-SIN(THETA*PI/180.)&
         *ZTEMPO
      ZR=-SIN(THETA*PI/180.)*XTEMPO-COS(THETA*PI/180.)&
         *ZTEMPO
      YR=YR

      RETURN
      END

!=======================================================================
       SUBROUTINE ROTATE_XAXIS(XR,YR,ZR,THETA)
!=======================================================================
!     ROTATE Y AND Z COORDINATE TO COUNTERCLOCKWISE DIRECTION
!     THETA IN DEGREE UNITS
      IMPLICIT NONE
      REAL*8   , intent(in)     :: THETA
      REAL*8   , intent(inout)  :: XR,YR,ZR
      REAL*8                    :: PI
      REAL*8                    :: YTEMPO,ZTEMPO

      PI = ACOS(-1.)
      
      YTEMPO=YR
      ZTEMPO=ZR

      ZR=COS(THETA*PI/180.)*ZTEMPO+SIN(THETA*PI/180.)&
         *YTEMPO
      YR=-SIN(THETA*PI/180.)*ZTEMPO+COS(THETA*PI/180.)&
         *YTEMPO
      XR=XR

      RETURN
      END
