!---------------------------------------------!
!  SHARIF UNIVERSITY OF THECHNOLOGY           !
!                                             !
!        Hamid Abbaszadeh 83701556            !
!        email: abbaszadeh.h@gmail.com        !
!                                             !
!        Dr. Langari                          !
!                                             !
!        October 15, 2005                     !
!                                             !
!  SPIN GAME (Problem Number 3)               !
!---------------------------------------------!

PROGRAM SPIN_GAME
   IMPLICIT NONE

   CHARACTER*14 :: FILEN1
   INTEGER, PARAMETER :: M = 1000
   INTEGER :: rand1, rand2, Energy
   INTEGER :: OLD_SPIN, BeforeJ, NextJ
   INTEGER :: IBM, i, j, E1, E2, iseed, n, k
   INTEGER, DIMENSION(M) :: SPIN, NS
   INTEGER, DIMENSION(5, 20) :: NT

   FILEN1 = 'SPIN_DATA.TXT'
   OPEN (UNIT=15, FILE=FILEN1)

   PRINT *, "Please Give a number for SEED ?"
   READ *, iseed

   IBM = 2*iseed - 1

   k = 1
   DO k = 1, 5
      DO i = 1, M
         NS(i) = 0
         SPIN(i) = rand1(IBM)
      END DO
      DO n = 1, 20
         DO i = 1, M

            j = rand2(IBM)

            IF (j == 1) THEN
               BeforeJ = 1000
               NextJ = 2
            ELSE IF (j == 1000) THEN
               BeforeJ = 999
               NextJ = 1
            ELSE
               BeforeJ = j - 1
               NextJ = j + 1
            END IF
            IF (SPIN(BeforeJ) == SPIN(NextJ)) THEN
               IF (SPIN(j) /= SPIN(NextJ)) THEN
                  SPIN(j) = -SPIN(j)
                  NS(j) = NS(j) + 1
               END IF
            ELSE
               OLD_SPIN = SPIN(j)
               E1 = Energy(SPIN)
               SPIN(j) = -SPIN(j)
               E2 = Energy(SPIN)
               SPIN(j) = OLD_SPIN

               IF (E1 == E2) THEN
                  SPIN(j) = SPIN(j)*RAND1(IBM)
                  IF (SPIN(j) /= OLD_SPIN) THEN
                     NS(j) = NS(j) + 1
                  END IF
               ELSE IF (E2 < E1) THEN
                  SPIN(j) = -SPIN(j)
                  NS(j) = NS(j) + 1
               END IF
            END IF
         END DO
         DO i = 1, 1000
            IF (NS(i) == 0) THEN
               NT(k, n) = NT(k, n) + 1
            END IF
         END DO
      END DO
   END DO

   DO i = 1, 20
      WRITE (15, *) NT(1, i), NT(2, i), NT(3, i), NT(4, i), NT(5, i)
   END DO

   CLOSE (15)
END PROGRAM SPIN_GAME

!-------------- FUNCTIONS OF THIS PROGRAM --------------------
INTEGER FUNCTION rand1(BM)
   IMPLICIT NONE
   INTEGER BM

   BM = BM*16807
   IF (BM >= 0) THEN
      rand1 = +1
   ELSE
      rand1 = -1
   END IF

   RETURN
END FUNCTION rand1

INTEGER FUNCTION rand2(BM)
   IMPLICIT NONE
   INTEGER BM

   BM = BM*16807
   IF (BM < 0) THEN
      BM = -BM
   END IF
   rand2 = BM/2149633 + 1

   RETURN
END FUNCTION rand2

INTEGER FUNCTION Energy(SP)
   IMPLICIT NONE
   INTEGER SP(1000)
   INTEGER i

   Energy = 0
   DO i = 1, 999
      Energy = Energy - SP(i)*SP(i + 1)
   END DO
   Energy = Energy - SP(1000)*SP(1)

   RETURN
END FUNCTION Energy
