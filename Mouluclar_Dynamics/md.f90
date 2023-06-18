!---------------------------------------------------------------!
!		SHARIF UNIVERSITY OF THECHNOLOGY		!
!	Hamid Abbaszadeh 83701556				!
!	email: abbaszadeh.h@gmail.com			!
!								!
!	Novwmber 8, 2005 					!
!	Mouluclar Dynamics (Problem Number 6) Version 1.0.0	!
!---------------------------------------------------------------!
program md
use struc
use molecular
use export
implicit none
type(obj), dimension(:), allocatable	:: p_new,p_current,p_old, p_first
type(eng), dimension(:), allocatable	:: temp, temp_mean
type(eng)				:: temp0
type(conf)				:: config
real, dimension(1000)			:: gcum=0.0
real, dimension(1000)			:: pvx=0.0, pvy=0.0
real					:: clock_start, clock_finish
integer					:: i, j

 config%timebond 	= 10000
 config%n_particles	= 100
 config%step		= 1
 config%lennar		= 1			!1 for Lennard-Jones and 2 for 1/R Potential
 config%VV		= 1			!1 for Vedrlet Method and 2 for Velocity Verlet
 config%dt 		= 0.02
 config%dt2		= config%dt*config%dt
 config%dr 		= 0.1
 config%Cutoff		= 35.0
 config%L1		= 16.0
 config%L		= 0.0
 config%dL		= 2.0
 config%v_max		= 1.0
 config%temp0		= 1.0
 config%temp1		= 1.0
 config%R			= config%temp1/config%temp0

allocate(p_new(config%n_particles),p_current(config%n_particles))
allocate(p_first(config%n_particles),p_old(config%n_particles))
allocate(temp(config%timebond))
allocate(temp_mean(config%step))

print *,"What Potential Lennard-Jones[1], 1/r[2]   :"
read *,config%lennar
print *,"What Method Verlet[1], Velocity Verlet[2] :"
read *,config%VV

call print1(config)
do i=1,config%step
	call cpu_time(clock_start)
	temp=temp0
	gcum=0.0
	config%L = config%L1 + config%dL*i
!--------------------------------------------------------------------
	call init(p_old,p_current,p_first,config)
	if (config%VV==1) then
		do j=1,config%timebond
			!if (j==1000.or.j==2000.or.j==3000.or.j==4000.or.j==5000.or.j==6000.or.j==7000.or.j==8000.or.j==9000) config%R=1.2
			call verlet(p_old,p_current,p_new,p_first,temp(j),gcum,config)
			!config%R=1.0
		end do
	else
		do j=1,config%timebond
			call VVerlet(p_current,p_first,temp(j),gcum,config)
		end do
	end if
	call mean(temp,temp_mean(i),config)
!------------------ Exporting Data ----------------------------------
	call cpu_time(clock_finish)
	call printresult(temp_mean(i),clock_finish - clock_start,config%L)
	call exporting(i,p_first,p_current,temp,gcum,config)
end do
call exporting1(config,temp_mean)
end program md
