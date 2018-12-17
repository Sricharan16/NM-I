program main
implicit none

!input matrix
real*8,dimension(5,6) :: A
real :: tol=0.0001,c,d
integer :: i,j,k,l
real*8,dimension(5)::x
A(2,1)=0.0
A(2,2)=-1.0
A(2,3)=2.0
A(2,4)=-3.0
A(2,5)=4.0
A(1,1)=2.0
A(1,2)=3.0
A(1,3)=-1.0
A(1,4)=5.0
A(1,5)=-2.0
A(3,1)=-1.0
A(3,2)=3.0
A(3,3)=2.0
A(3,4)=-5.0
A(3,5)=1.0
A(4,1)=1.0
A(4,2)=2.0
A(4,3)=1.0
A(4,4)=2.0
A(4,5)=3.0
A(5,1)=-4.0
A(5,2)=-6.0
A(5,3)=-2.0
A(5,4)=8.0
A(5,5)=-1.0
x(1)=0.0
x(2)=0.0
x(3)=0.0
x(4)=0.0
x(5)=0.0
A(2,6)=-38.5
A(1,6)=32.4
A(3,6)=-17.9
A(4,6)=-13.9
A(5,6)=4.9

do i=1,5
	d=A(i,i)
	do k=1,6
		A(i,k)=A(i,k)/d
	end do
	do j=1,5
 		if(i.ne.j) then
 			c=A(j,i)
 			do k=1,6
 				A(j,k)=A(j,k)-(c*A(i,k))
 			end do
 		end if
 	end do
end do

 write (*,201) (A(i,6),i=1,5)
    201 format (6f12.6)
end program main