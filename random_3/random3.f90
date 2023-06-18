! Hamid Abbaszadeh 83701556
! email:abbaszadeh.h@gmail.com

program random1
   implicit none
   integer IBM, k, ICI, I, iseed
   real N(100)

   print *, 'Give iseed number'
   read *, iseed
   IBM = 2*iseed - 1
   do k = 0, 99
      ICI = 0
      do I = 1, 32
         ICI = ISHFT(ICI, 1)
         IBM = IBM*16807
         if (IBM .LT. 0) then
            ICI = ICI + 1
         end if
      end do
      N(k) = .5 + 2.328306e-10*ICI
   end do

   do k = 0, 99
      print *, N(K)
   end do

END program random1