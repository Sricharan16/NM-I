program main
implicit none
real ::x,f,xprev
integer :: i
read(*,*) x
call func(x,f)
do i=1,1000
	x=f
	call func(x,f)
	write(*,*) i,")"," X=", x,"f =",f
	if (abs((f-x)/x) .le. 0.001) then
		EXIT
	end if
end do
write(*,*) "Final =",f


end program main

subroutine func(x,f)
real,intent(inout)::x
real,intent(out)::f
f= exp(x)/3
end subroutine func