program ee
implicit none

integer, parameter			:: N=12, L=2**N

integer, dimension(0:L-1)	:: fi=0
real, dimension(0,L-1)		:: f0=0, f1=0
integer	:: t

do t=0,L-1
	fi(t)=t
end do

f0(100)=1

H0(f0,f1,N,L)
print *,f0
print *,f1


contains
function P(f,i,j)
	integer				:: P
	integer, intent(in)	:: f
	integer, intent(in)	::i,j
	integer :: b1,b2
	b1=0
	b2=0
	if (ibclr(f,i).ne.f) b1=1
	if (ibclr(f,j).ne.f) b2=1
	P=f+(b2-b1)*(2**i-2**j)
end function P
subroutine H0(f0,f1,N,L)
	real, dimension(:), intent(in)	:: f0
	real, dimension(:), intent(out)	:: f1
	integer, intent(in)				:: N, L
	integer		:: i, j, c
	f1=0
	do i=0,L-1
		do j=0,N-1
			if (j/=L-1) then
				c=P(f0(i),j,j+1)
			else
				c=P(f0(i),L-1,0)
			end if
			f1(c) = f0(i) + f1(c)
		end do
	end do
end subroutine H0

function Braket(f0,f1,N,L)
	real								:: Braket
	integer, dimension(:), intent(in)	:: f0, f1
	integer, intent(in)					:: N, L
	Braket = sum(f0*f1)
end function Braket

end program ee
