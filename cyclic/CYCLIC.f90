!-----------------------------------------------!
!		SHARIF UNIVERSITY OF THECHNOLOGY		!
!												!
!	Hamid Abbaszadeh 83701556					!
!	email: abbaszadeh.h@gmail.com				!
!												!
!	Dr. Langari									!
!												!
!	October 24, 2005 							!
!												!
!	Euler (Problem Number 5)					! 
!-----------------------------------------------!
	PROGRAM cyclic
	IMPLICIT NONE

	INTEGER, PARAMETER	:: N=2000
	REAL, PARAMETER		:: dt=0.005,G=1.0,M=1.0
	INTEGER				:: i
	REAL				:: r,r_mid,x_mid,y_mid,V_x_mid,V_y_mid
	REAL, DIMENSION(N)	:: V_x,V_y,x,y
 
!--------------Initial Conditions-----------------------------!		
		x(1) = 1.0
		y(1) = 0.0
		V_x(1) = 0.0
		V_y(1) = 1.0

		OPEN (1,FILE="c_euler.dat")
		OPEN (2,FILE="c_RK.dat")
		
!-------------------Euler Method------------------------------!
		DO i=2,N

			r = (x(i-1)**2 + y(i-1)**2)**1.5
			V_x(i) = V_x(i-1) - G*M*(x(i-1)/r)*dt
			V_y(i) = V_y(i-1) - G*M*(y(i-1)/r)*dt
			x(i) = x(i-1) + V_x(i-1)*dt
			y(i) = y(i-1) + V_y(i-1)*dt

		END DO

		DO i=1,N
			WRITE (1,*) x(i),",",y(i)
		END DO
		
!-------------------RK Method---------------------------------!
		DO i=2,N
					
			r= (x(i-1)**2 + y(i-1)**2)**1.5
			
			V_x_mid = V_x(i-1) - G*M*(x(i-1)/r)*dt/2
			V_y_mid = V_y(i-1) - G*M*(y(i-1)/r)*dt/2
			x_mid = x(i-1) + V_x_mid*dt
			y_mid = y(i-1) + V_y_mid*dt

			r_mid = (x_mid**2 + y_mid**2)**1.5
			
			V_x(i) = V_x(i-1) - G*M*(x_mid/r_mid)*dt
			V_y(i) = V_y(i-1) - G*M*(y_mid/r_mid)*dt
			x(i) = x(i-1) + V_x(i-1)*dt
			y(i) = y(i-1) + V_y(i-1)*dt

		END DO

		DO i=1,N
			WRITE (2,*) x(i),",",y(i)
		END DO

		ClOSE(1)
		CLOSE(2)

	END PROGRAM cyclic