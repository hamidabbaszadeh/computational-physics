program ising ! 2D Monte Carlo Simulation of Ising Model
use DFPORT
implicit none

! Variable declarations:
integer :: i, j, m, n, m2, n2 ! dummy integers
integer, allocatable :: A(:,:) ! matrix containing spins
integer :: nrows, ncols ! number of rows and cols of A
real :: temp, beta ! temperature, inverse temperature
integer :: ConfigType ! starting configuration type
integer :: npass ! number of passes for MC algorithm
integer :: ipass ! the current pass number
integer :: nequil ! number of equilibration steps
integer :: trial_spin ! values of changed spin
real :: high_temp ! starting temp for scan
real :: low_temp ! final temp for scan
real :: temp_interval ! interval between scan points
integer :: nscans ! number of scans (each at diff T)
integer :: iscan ! current scan number
logical :: MovieOn ! set to .true. to make movie of 1 temp
real :: deltaU ! change in energy between 2 configs
real :: deltaU1 ! energy changes for lattice gas
real :: log_eta ! log of random number to compare to
real :: magnetization ! magnetization of all spins in lattice
real :: magnetization_ave ! cumulative average magnetization
real :: magnetization2_ave ! cumulative average of mag. squared
real :: energy ! energy of all spins in lattice
real :: energy_ave ! cumulative average of energy
real :: energy2_ave ! cumulative average of energy squared
integer :: output_count ! # times things have been added to averages
integer,parameter :: L=100
 
print*, "________________MONTE CARLO 2D ISING MODEL________________"
print*, "Monte Carlo Statistics for 2D Ising Model with"
print*, " periodic boundary conditions."
print*, "The critical temperature is approximately 2.3, as seen on"
print*, " Chandler p. 123."

 
 
nrows=L
ncols=L
npass=20000000
nequil=1000000
high_temp=5.00
low_temp=0.50
temp_interval=0.1
ConfigType=3
MovieOn=.false.
 

! Set the dimensions of the matrix of spin arrays. This program uses
! periodic boundary conditions, so the first two rows and columns are
! the same as the last two.
allocate(A(nrows,ncols))

! Open output files:
open(unit=32,file='spin-array.dat',status='replace',action='write')
write(32,*) nrows
write(32,*) ncols
nscans = int((high_temp - low_temp)/temp_interval) + 1
if (MovieOn) then
	write(32,*) 51
	write(32,*) 1
else
	write(32,*) nscans
	write(32,*) 2
end if

open(unit=33,file='magnetization.dat',status='replace',action='write')
write(33,*) "temp ave_magnetization ave_magnetization^2 susceptibility"

open(unit=34,file='energy.dat',status='replace',action='write')
write(34,*) "temp ave_energy ave_energy^2 C_v"
 
open(unit=36,file='t-e-m40.dat',status='replace',action='write')
write(36,*) "time temp energy magnetization magnetization_ave"

open(unit=36,file='t-e-m30.dat',status='replace',action='write')
write(37,*) "time temp energy magnetization magnetization_ave"

open(unit=36,file='t-e-m25.dat',status='replace',action='write')
write(38,*) "time temp energy magnetization magnetization_ave"

open(unit=36,file='t-e-m23.dat',status='replace',action='write')
write(39,*) "time temp energy magnetization magnetization_ave"

open(unit=36,file='t-e-m21.dat',status='replace',action='write')
write(40,*) "time temp energy magnetization magnetization_ave"

open(unit=36,file='t-e-m15.dat',status='replace',action='write')
write(41,*) "time temp energy magnetization magnetization_ave"

open(unit=36,file='t-e-m10.dat',status='replace',action='write')
write(42,*) "time temp energy magnetization magnetization_ave"
 
  
scan_loop: do iscan = 1, nscans		!main loop
	temp = high_temp - temp_interval*(iscan-1)
	print*, "Running program for T =", temp

	! Initialize variables
	beta = 1.0/temp
	output_count = 0
	energy_ave = 0.0
	energy2_ave = 0.0
	magnetization_ave = 0.0
	magnetization2_ave = 0.0

	! Set up the initial spin configuration.
	select case(ConfigType)
	case(1) ! checkerboard setup
		A(1,1) = 1
		do i = 1, nrows-1
			A(i+1,1) = -A(i,1)
		enddo
		do j = 1, ncols-1
			A(:,j+1) = -A(:,j)
		enddo
		! (note: the requirement that nrows and ncols are even is to
		! ensure that the first two rows/cols start out the same as the
		! last two)
	case(2) !interface
		do i = 1, nrows
			do j = 1, (ncols+2)/2
				A(i,j) = 1
			enddo
			do j = (ncols)/2 + 1, ncols
				A(i,j) = -1
			enddo
		enddo
	case(3) ! unequal interface
		do i = 1, nrows
			do j = 1, (ncols)/4
				A(i,j) = 1
			enddo
			do j = (ncols)/4 + 1, ncols
				A(i,j) = -1
			enddo
		enddo
	case default
		print*, "Error! Check ConfigType parameter in ising.in"
		stop
	end select

	! Main loop containing Monte Carlo algorithm:
	MC_passes: do ipass = 0, npass

		!If MovieOn is .true., write the spin array to an output every
		!npass/50 steps.
		!if ((MovieOn) .and. (mod(ipass,npass/50) == 0)) then
		!	do i = 2, nrows+1
		!		do j = 2, ncols+1
		!			write(32,*) A(i,j)
		!		enddo
		!	enddo
		!endif

		! If ipass is greater than nequil (the number of equilibration steps),
		! calculate the magnetization and energy:
		if (ipass > nequil) then
			output_count = output_count + 1
			magnetization = sum(A(2:nrows+1,2:nrows+1))/(ncols*nrows*1.0)
			magnetization_ave = magnetization_ave + magnetization
			magnetization2_ave = magnetization2_ave + magnetization**2
			energy = 0.0
			do i = 1, nrows
				energy
				do j = 2, ncols + 1
					energy = energy - A(i,j)*(A(i-1,j)+A(i+1,j)+A(i,j-1)+A(i,j+1))
				enddo
			enddo
			! Divide the energy by the total number of spins to get the ave
			! energy per spin, and divide by 2 to account for double counting.
			energy = energy/(ncols*nrows*2.0)
			energy_ave = energy_ave + energy
			energy2_ave = energy2_ave + energy**2
			if(temp==4.0) write(36,*) output_count ,temp ,energy ,magnetization ,magnetization_ave 
			if(temp==3.0) write(37,*) output_count ,temp ,energy ,magnetization ,magnetization_ave
			if(temp==2.5) write(38,*) output_count ,temp ,energy ,magnetization ,magnetization_ave
			if(temp==2.3) write(39,*) output_count ,temp ,energy ,magnetization ,magnetization_ave
			if(temp==2.1) write(40,*) output_count ,temp ,energy ,magnetization ,magnetization_ave
			if(temp==1.5) write(41,*) output_count ,temp ,energy ,magnetization ,magnetization_ave
			if(temp==1.0) write(42,*) output_count ,temp ,energy ,magnetization ,magnetization_ave
 		endif

		! Randomly choose a spin to change:
		m = nint((nrows-1)*rand() + 2) ! choose a random row
		n = nint((ncols-1)*rand() + 2) ! choose a random column
		trial_spin = -A(m,n) ! trial spin value

		! Find change in energy (deltaU) due to trial move.
		! If exp(-beta*deltaU) > eta, where eta is random, accept move:
		deltaU = -trial_spin*(A(m-1,n)+A(m+1,n)+A(m,n-1)+A(m,n+1))*2
		log_eta = dlog(rand(5) + 1.0d-10) ! random number 0-1 (+ tiny offset)
		if (-beta*deltaU > log_eta) then
			A(m,n) = trial_spin
			if (m == 2) A(nrows+2,n) = trial_spin
			if (m == nrows+1) A(1,n) = trial_spin
			if (n == 2) A(m,ncols+2) = trial_spin
			if (n == ncols+1) A(m,1) = trial_spin
		endif
 
	enddo MC_passes
	close(36)

	!Write final spin array to output file
	!if (.not. MovieOn) then
	!	do i = 2, nrows + 1
	!		do j = 2, ncols + 1
	!			write(32,*) A(i,j)
	!		enddo
	!	enddo
	!endif
	write(33,*) temp, abs(magnetization_ave/output_count), &
	magnetization2_ave/output_count, &
	beta*(magnetization2_ave/output_count - (magnetization_ave/output_count)**2)
	write(34,*) temp, energy_ave/output_count, energy2_ave/output_count, &
	(beta**2)*(energy2_ave/output_count - (energy_ave/output_count)**2)
enddo scan_loop

close(32)
close(33)
close(34)

print*, "Program ising.f90 complete!"
print*, "Look at 'spin-array' with IDL program see_spins.pro"

end program ising