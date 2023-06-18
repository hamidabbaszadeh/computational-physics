module struc
implicit none

type obj
	real :: x=0.0, y=0.0
	real :: vx=0.0, vy=0.0
	real :: ax=0.0, ay=0.0
end type obj
type eng
	real :: k=0.0
	real :: p=0.0
	real :: e=0.0
	real :: t=0.0
	real :: pr=0.0
	real :: virial=0.0
	real :: r2=0.0
	real :: dr2=0.0
end type eng
type conf
	character(23), dimension(2)	:: LR = (/"Lennard-Jones Potential","1/r Potential"/)
	character(23), dimension(2)	:: V2 = (/"Verlet Method","Velocity Verlet Method"/)
	character(1), dimension(2)	:: LR2 = (/ "L","R"/)
	character(1), dimension(2)  :: V1 = (/"V","S"/)
	integer	:: timebond=0
	integer	:: n_particles=0
	integer	:: step=0
	integer :: Lennar=1
	integer :: VV=1
	real	:: dt=0.0
	real	:: dt2=0.0
	real	:: dr=0.0
	real	:: Cutoff=0.0
	real	:: L1=0.0
	real	:: L=0.0
	real	:: dL=0.0
	real	:: v_max=0.0
	real	:: temp0=0.0
	real	:: temp1=0.0
	real	:: R=0.0
end type conf
end module struc