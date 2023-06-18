!---------------------------------------------------------------!
!		SHARIF UNIVERSITY OF THECHNOLOGY		!
!	Hamid Abbaszadeh 83701556				!
!	email: abbaszadeh.h@gmail.com			!
!								!
!	Novwmber 8, 2005 					!
!	Bouncing Balls (Problem Number 5) Version 1.0.9		!
!---------------------------------------------------------------!
      PROGRAM BOUNCINGBALL
         IMPLICIT NONE
         INTEGER, PARAMETER :: N = 1000000
         REAL, PARAMETER :: g = 9.8, Mass = 1.0
         REAL, PARAMETER :: fixdt = 0.001, rt = 10, endtime = 100.0
         TYPE OBJECT
            REAL, DIMENSION(N) :: x, v
            REAL :: m
         END TYPE OBJECT
         REAL :: dt
         INTEGER :: i, j
         REAL :: m1, m2, ma, temp
         REAL, DIMENSION(N) :: t
         TYPE (OBJECT) :: p1, p2
         LOGICAL, DIMENSION(N) :: b = .FALSE., c = .FALSE.
!------------INIT-------------------------------
         OPEN(UNIT = 11, FILE = "alltxv.dat")
         OPEN(UNIT = 12, FILE = "tx2v2_1.dat")
         OPEN(UNIT = 13, FILE = "tx2v2_2.dat")
   90    FORMAT(F9.5, F9.5, F9.5, F9.5, F9.5)
   91    FORMAT(F9.5, F9.5)
         p1%x(1) = 1.0
         p1%v(1) = 0.0
         p1%m = 1.0 * Mass
         p2%x(1) = 3.0
         p2%v(1) = 0.0
         p2%m = 1.0 * p1%m
         m1 = 2 * p1%m / (p1%m + p2%m)
         m2 = 2 * p2%m / (p1%m + p2%m)
         ma = (p1%m - p2%m) / (p1%m + p2%m)

         dt = fixdt
         t(1) = 0.0
         i = 1
!----------------main----------------------------
         DO
            p1%v(i + 1) = p1%v(i) - g * dt
            p1%x(i + 1) = p1%x(i) + (p1%v(i + 1) + p1%v(i)) * dt / 2
            IF (p1%x(i + 1) .LT. 0) THEN
               p1%v(i)= -p1%v(i)
               b(i) = .TRUE.
               dt = dt / rt
            ELSE
               p2%v(i + 1) = p2%v(i) - g*dt
               p2%x(i + 1) = p2%x(i) + (p2%v(i + 1) + p2%v(i)) * dt / 2
               IF (p1%x(i + 1) .GT. p2%x(i + 1)) THEN
                  temp = ma * p1%v(i) + m2 * p2%v(i)
                  p2%v(i) = m1 * p1%v(i) - ma * p2%v(i)
                  p1%v(i) = temp
                  c(i) = .TRUE.
                  dt = dt / rt
               ELSE
                  t(i + 1) = t(i) + dt
                  IF (t(i + 1) .GT. endtime) EXIT
                  dt = fixdt
                  i = i + 1
               END IF
            END IF
         END DO
!----------------exporting data---------------------
         DO j = 1, i
            WRITE(11, 90) t(j), p1%x(j), p2%x(j), p1%v(j), p2%v(j)
            IF (b(j)) WRITE(12, 91) p2%x(j), p2%v(j)
            IF (c(j)) WRITE(13, 91) p2%x(j), p2%v(j)
         END DO
      END PROGRAM BOUNCINGBALL
