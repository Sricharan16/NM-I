program main
implicit none
real ,dimension(:),allocatable::x,f
integer :: i,n,c,j
n=2
allocate(x(n))
allocate(f(n))
read(*,*) x

do i=1,1000
	c=0
	call func(x,n,f)
	!write(*,20) x(1) ,x(2)
	!print *,""
	do j=1,n
	if (abs((f(j)-x(j))/x(j)) .le. 0.0000001) then
		c=c+1
	end if
	end do
	if (c .eq. n) EXIT
	x=f
end do
20 format(f10.4,f10.4)
print *, i
write(*,*) "Final =",f


end program main

subroutine func(x,n,f)
real,dimension(n),intent(inout)::x,f
integer,intent(in) ::n
f(1)= x(1) +  (-2*x(1)**3 -8*x(1) +4*x(2) +4)/20.0
f(2)= x(2) - (-4*x(1) + x(2)**2 + 3*x(2) +1)/24.0
end subroutine func