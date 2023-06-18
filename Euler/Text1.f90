!---------------------------------------------!
!		SHARIF UNIVERSITY OF THECHNOLOGY        !
!                                             !
!	Hamid Abbaszadeh 83701556                  !
!	email: abbaszadeh.h@gmail.com              !
!                                             !
!	Dr. Langari                                !
!                                             !
!	October 17, 2005                           !
!                                             !
!	Euler (Problem Number 4)                   !
!---------------------------------------------!
PROGRAM Euler
   IMPLICIT NONE

   INTEGER, PARAMETER :: M = 200
   INTEGER :: i
   REAL :: ys, r, dt
   REAL, DIMENSION(M) :: x, y
   CHARACTER*14 :: File_name


   File_name = "Euler_data.txt"
   y(1) = 400.0
   x(1) = 0.0
   ys = 300.0
   dt = 0.05
   r = 1.0

   DO i = 2, M
      x(i) = x(i - 1) + dt
      y(i) = y(i - 1) - r * (y(i - 1) - ys) * dt
   END DO

   OPEN(15, FILE = File_name)

   DO i = 1, M
      WRITE (15, *) x(i) , "," , y(i)
   END DO

   Close(15)

END PROGRAM Euler