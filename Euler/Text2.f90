!-------------------------------------------!
!		SHARIF UNIVERSITY OF THECHNOLOGY	     !
!											           !
!	Hamid Abbaszadeh 83701556                !
!	email: abbaszadeh.h@gmail.com            !
!                                           !
!	Dr. Langari                              !
!                                           !
!	October 17, 2005                         !
!                                           !
!	Euler (Problem Number 5)                 !
!-------------------------------------------!
PROGRAM Euler
   IMPLICIT NONE

   INTEGER, PARAMETER :: M=300
   INTEGER :: i
   REAL :: ys, r, dt
   REAL, DIMENSION(M) :: x, y
   CHARACTER*14 :: File_name
   CHARACTER :: answer


   PRINT *, "***************** Euler Method ******************"

   DO
      PRINT *, "Please give Initial Tempereture T (350.0) :"
      READ *, y(1)
      PRINT *, "Give the Room Tempereture RT (300.0) :"
      READ *, ys
      PRINT *, "Give the rate of colding r (0.1) :"
      READ *, r
      PRINT *, "Give Delta t (0.1 to 0.01) :"
      READ *, dt
      PRINT *, "Give a File name for exorting Data :"
      READ *, File_name

      x(1) = 0.0

      DO i = 2, M
         x(i) = x(i - 1) + dt
         y(i) = y(i - 1) - r * (y(i - 1) + ys) * dt
      END DO

      OPEN (15, FILE = File_name)

      DO i= 1, M
         WRITE (15, *) x(i), y(i)
      END DO

      Close(15)

      PRINT *, "Do you want to repeat program (Y or N)?"
      READ *, answer
      IF (answer == "N" .OR. answer == "n") EXIT

      File_name = ""

   END DO

   CLOSE (15)

END PROGRAM Euler
