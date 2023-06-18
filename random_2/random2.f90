!---------------------------------------------!
!		SHARIF UNIVERSITY OF THECHNOLOGY	!
!											!
!	Hamid Abbaszadeh 83701556				!
!	email: abbaszadeh.h@gmail.com			!
!											!
!	Dr. Langari								!
!											!
!	October 15, 2005						!
!											!
!	RANDOM NUMBER GENRATOR (SECOND PROBLEM)	!
!---------------------------------------------!
	PROGRAM random2
	
	IMPLICIT NONE

	CHARACTER*14			:: FILEN1
	INTEGER, PARAMETER		:: M=128, C=5, N=256
	INTEGER					:: i
	INTEGER, DIMENSION(N)	:: r, k
	INTEGER, DIMENSION(N/2) :: x, y
	
	FILEN1='RAN2_DATA.TXT'
	OPEN(UNIT=15,FILE=FILEN1)

C	Genrating Random Number	
	k(1)=1
	DO i=2,N
		k(i)=MOD(k(i-1),M)*C
		r(i)=k(i)/M
	END DO
	
C	Spliting Random Number to x and y Array
	DO i=1,N/2
		x(i)=r(2*i-1)
		y(i)=2*r(i)
	END DO

	DO i=1,N/2
		PRINT *,x(i),y(i)		!Print Random Number in screen
		WRITE (15,*) x(i),y(i)	!Exprot Random Nuber to RAN2_DATA.TXT File
	END DO

	CLOSE (15)

	END PROGRAM random2