!-------------------------------------------------!
!   SHARIF UNIVERSITY OF THECHNOLOGY              !
!                                                 !
!   Hamid Abbaszadeh 83701556                     !
!   email: abbaszadeh.h@gmail.com                 !
!                                                 !
!   Dr. Langari                                   !
!                                                 !
!   October 29, 2005                              !
!                                                 !
!   Euler (Problem Number 5) Version 1.0.1        !
!-------------------------------------------------!
PROGRAM cyclic
   IMPLICIT NONE
   TYPE OBJ
      REAL :: X
      REAL :: Y
      REAL :: Vx
      REAL :: Vy
   END TYPE OBJ
   INTEGER, PARAMETER :: N = 2500
   REAL, PARAMETER :: dt = 0.02, G = 1.0, M = 1.0
   INTEGER :: i, j
   REAL :: Radius, Radius_mid, T
   TYPE(OBJ) :: P1
   REAL, DIMENSION(N, 2):: Velocity, P_eu, P_rk
   REAL, DIMENSION(1, 2):: P_mid, V_mid

!--------------Initial Conditions-----------------------------!
   P_eu(1, :) = (/1.0, 0.0/)
   P_rk(1, :) = P_eu(1, :)
   Velocity(1, :) = (/0.0, 1.0/)

   OPEN (UNIT=1, FILE="c_eurk.txt")
!-------------------Euler Method------------------------------!
   DO i = 2, N

      Radius = (P_eu(i - 1, 1)**2 + P_eu(i - 1, 2)**2)**1.5

      Velocity(i, :) = Velocity(i - 1, :) - G*M*(P_eu(i - 1, :)/Radius)*dt
      P_eu(i, :) = P_eu(i - 1, :) + Velocity(i - 1, :)*dt

   END DO
!-------------------RK Method---------------------------------!
   DO i = 2, N

      Radius = (P_rk(i - 1, 1)**2 + P_rk(i - 1, 2)**2)**1.5

      V_mid(1, :) = Velocity(i - 1, :) - G*M*(P_rk(i - 1, :)/Radius)*dt/2
      P_mid(1, :) = P_rk(i - 1, :) + V_mid(1, :)*dt

      Radius_mid = (P_mid(1, 1)**2 + P_mid(1, 2)**2)**1.5

      Velocity(i, :) = Velocity(i - 1, :) - G*M*(P_mid(1, :)/Radius_mid)*dt
      P_rk(i, :) = P_rk(i - 1, :) + Velocity(i - 1, :)*dt

   END DO
!------------------------Exporting DATA-----------------------!
   DO i = 1, N
      WRITE (1, *) P_eu(i, :), P_rk(i, :)
   END DO
!------------------------Cacuating T^2/a^3--------------------!
   T = 0.0
   DO i = 2, N
      IF (1.0 == P_rk(i, 1)) EXIT
      T = T + dt
   END DO
   PRINT *, T**2/Radius**3
!-------------------------------------------------------------!
   ClOSE (1)

END PROGRAM cyclic
