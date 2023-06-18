module export
use struc
implicit none

contains
subroutine exporting(j,p_first,p_old,temp,gcum,config)
	type(eng), intent(in), dimension(:)	:: temp
	real, intent(in), dimension(:)		:: gcum
	type(obj), intent(in), dimension(:)	:: p_old, p_first
	type(conf), intent(in)				:: config
	integer, intent(in)					:: j
	integer		:: i
	92 format(2(f14.5))
	94 format(4(f14.5))
	99 format(9(f14.5))
	open(10,file="data/" // config%LR2(config%lennar) //config%V1(config%VV) // "firend" // achar(j+64) // ".dat")
	open(11,file="data/" // config%LR2(config%lennar) //config%V1(config%VV) // "energy" // achar(j+64) // ".dat")
	open(12,file="data/" // config%LR2(config%lennar) //config%V1(config%VV) // "densig" // achar(j+64) // ".dat")
	do i=1,config%n_particles
		write(10,94),p_first(i)%x,p_first(i)%y,p_old(i)%x,p_old(i)%y
	end do
	do i=1,config%timebond
		write (11,99), i*config%dt,temp(i)
	end do
	i=1
	do while (config%dr*i<config%L*0.5)
		write (12,92), config%dr*i,gcum(i)/(3.1415*config%dr*(2*i+1)*config%timebond)
		i=i+1
	end do
	close(10)
	close(11)
	close(12)
end subroutine exporting
!----------------------------------------------------------
subroutine exporting1(config,tm)
	type(conf), intent(in)				:: config
	type(eng), intent(in), dimension(:)	:: tm
	integer		:: i
	90 format(8(f14.5))
	open(9,file="data/" // config%LR2(config%lennar)// config%V1(config%VV) // "config.dat")
	write(9,*),"Molecular Dynamics with The" // config%LR(config%Lennar)
	write(9,*),"with ",config%V2(config%VV)
	write(9,*),config%V2(config%VV)
	write(9,*),"Timebond  :",config%timebond,   " Cutoff     :",config%Cutoff
	write(9,*),"Particles :",config%n_particles," temp1      :",config%temp1
	write(9,*),"Temprature:",config%temp1
	write(9,*),"v max     :",config%v_max
	write(9,*),"dt        :",config%dt
	write(9,*),"dr        :",config%dr
	write(9,*),"-------------------------------------------------------------------------------"
	write(9,*),"L        MPE       MKE        MTE      MTemp     MVirial        Mr2        Mdr2"
	do i=1,config%step
		write(9,90),config%L1 + config%dL*i,tm(i)
	end do
	close(9)
end subroutine exporting1
!----------------------------------------------------------
subroutine print1(config)
	type(conf)	:: config
	print *,"Molecular Dynamics with The " // config%LR(config%Lennar)
	print *,"with ",config%V2(config%VV)
	print *,"Hamid Abbaszadeh, hamidabp@yahoo.com 2006/01/12"
	print *,"Number of Particles :",config%n_particles
	print *,"Total Time          :",config%timebond*config%dt
	print *,"Temprature          :",config%temp1
	print *,"dt                  :",config%dt
	print *,"dr                  :",config%dr
	print *,"----------------------------------------------------------------------------"
	print *,"         L        MTE      MTemp        MPR       Mdr2        Mr2   CPU time"
end subroutine print1
!----------------------------------------------------------
subroutine printresult(tm,t,L)
	type(eng)	:: tm
	real		:: t,L
	11 format(7(f11.5))
	write(*,11),L,tm%e,tm%t,tm%pr, tm%dr2,tm%r2,t
end subroutine printresult
end module export