!---------------------------------------------------!
!                SHARIF UNIVERSITY OF THECHNOLOGY              !
!        Hamid Abbaszadeh 83701556                        !
!        email: abbaszadeh.h@gmail.com                    !
!                                                   !
!        June 18, 2023                                                          !
!        Bouncing Balls (Problem Number 5) Version 2.0.0  !
!---------------------------------------------------!
PROGRAM BouncingBall
   IMPLICIT NONE
   INTEGER, PARAMETER :: N = 1000000
   REAL, PARAMETER :: g = 9.8, Mass = 1.0
   REAL, PARAMETER :: fixdt = 0.001, rt = 10.0, endtime = 100.0
   TYPE :: Object
      REAL, DIMENSION(N) :: x, v
      REAL :: m
   END TYPE Object
   REAL :: dt
   INTEGER :: i, j
   REAL :: m1, m2, ma, temp
   REAL, DIMENSION(N) :: t
   TYPE(Object) :: p1, p2
   LOGICAL, DIMENSION(N) :: b = .FALSE., c = .FALSE.

   OPEN (UNIT=11, FILE="alltxv.dat")
   OPEN (UNIT=12, FILE="tx2v2_1.dat")
   OPEN (UNIT=13, FILE="tx2v2_2.dat")

   p1%x(1) = 1.0
   p1%v(1) = 0.0
   p1%m = 1.0*Mass

   p2%x(1) = 3.0
   p2%v(1) = 0.0
   p2%m = 1.0*p1%m

   m1 = 2.0*p1%m/(p1%m + p2%m)
   m2 = 2.0*p2%m/(p1%m + p2%m)
   ma = (p1%m - p2%m)/(p1%m + p2%m)

   dt = fixdt
   t(1) = 0.0
   i = 1

   DO WHILE (t(i) <= endtime)
      p1%v(i + 1) = p1%v(i) - g*dt
      p1%x(i + 1) = p1%x(i) + (p1%v(i + 1) + p1%v(i))*dt/2.0

      IF (p1%x(i + 1) < 0.0) THEN
         p1%v(i) = -p1%v(i)
         b(i) = .TRUE.
         dt = dt/rt
      ELSE
         p2%v(i + 1) = p2%v(i) - g*dt
         p2%x(i + 1) = p2%x(i) + (p2%v(i + 1) + p2%v(i))*dt/2.0

         IF (p1%x(i + 1) > p2%x(i + 1)) THEN
            temp = ma*p1%v(i) + m2*p2%v(i)
            p2%v(i) = m1*p1%v(i) - ma*p2%v(i)
            p1%v(i) = temp
            c(i) = .TRUE.
            dt = dt/rt
         ELSE
            t(i + 1) = t(i) + dt
            dt = fixdt
            i = i + 1
         END IF
      END IF
   END DO

   DO j = 1, i
      WRITE (11, '(5(F9.5, 1X))') t(j), p1%x(j), p2%x(j)
      IF (b(j)) WRITE (12, '(2(F9.5, 1X))') p2%x(j), p2%v(j)
      IF (c(j)) WRITE (13, '(2(F9.5, 1X))') p2%x(j), p2%v(j)
   END DO
END PROGRAM BouncingBall
