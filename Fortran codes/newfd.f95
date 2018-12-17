program main
implicit none
real,dimension(:),allocatable ::x,y,dely,del2y,del3y
real ::alpha,givenval,yval
integer ::i,n

open(5,file="inp.txt")
read(5,*) n
allocate(x(n))
allocate(y(n))
allocate(dely(n-1))
allocate(del2y(n-2))
allocate(del3y(n-3))
do i=1,n
	read(5,*) x(i),y(i)
end do

do i=1,n-1
	dely(i) = y(i+1)-y(i)
	print *, dely(i)
end do


do i=1,n-2
	del2y(i)=dely(i+1)-dely(i)
	print *, del2y(i)
end do

do i=1,n-3
	del3y(i) = del2y(i+1)-del2y(i)
	print *, del3y(i)
end do
print *, "enter x"
read(5,*) givenval
print *,givenval
call alphfun(alpha,x,n,givenval)
call yfunc(y(1),alpha,dely(1),del2y(1),del3y(1),yval)
print *,"Y=",yval

close(5)
end program main


subroutine alphfun(alpha,x,n,val)
real,intent(in)::val
integer,intent(in) :: n
real,dimension(n),intent(in) :: x
real,intent(out) :: alpha
alpha = (val - x(1))/(x(2)-x(1))
print *,"alpha =", alpha
end subroutine alphfun
	

subroutine yfunc(y0,alpha,dely,del2y,del3y,yofx)
real,intent(in)::dely,del2y,del3y,y0,alpha
real,intent(out) :: yofx
yofx = (y0) + alpha*dely + alpha*(alpha-1)*del2y/2 + alpha*(alpha-1)*(alpha-2)*del3y/(6) 
end subroutine yfunc
