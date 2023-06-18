!---------------------------------------------!
!	 SHARIF UNIVERSITY OF THECHNOLOGY	    !
!											!
!	Hamid Abbaszadeh 83701556				!
!	email: abbaszadeh.h@gmail.com			!
!											!
!	Dr. Langari								!
!											!
!	October 15, 2005						!
!											!
!	RANDOM NUMBER GENRATOR (First Problem)	!
!---------------------------------------------!

      PROGRAM random1
         IMPLICIT NONE

         CHARACTER*15		:: File_Name
         INTEGER, PARAMETER	:: M=256
         INTEGER				:: n, i, IBM, iseed
         REAL, DIMENSION(M)	:: k
         REAL, PARAMETER		:: 	nor=4.65661287E-10

         File_Name="RAN1_DATA.TXT"
         OPEN (17,FILE=File_Name)

         PRINT *, 'Give number of Random number that you want ?'
         READ *, n
         PRINT *,'Give SEED number for genrating Random numbers ?'
         READ *, iseed


         IBM=2*iseed-1

         DO i=1,n
            IBM=IBM*16807
            IF (IBM.LT.0) THEN
               IBM=-IBM
            END IF
            k(i)=nor * REAL(IBM)
         END DO

         DO i=1,n
            PRINT *,k(i)
            WRITE (17,*) k(i)
         END DO

         CLOSE(17)

      END PROGRAM random1
