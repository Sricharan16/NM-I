program main
implicit none
real,dimension(:,:),allocatable :: a
real,dimension(:),allocatable :: x
integer :: n
n=5
allocate(a(n,n+1))
allocate(x(n))

a(1,1)=15
a(1,2)=-1
a(1,3)=2
a(1,4)=-3
a(1,5)=4
a(1,6)=8
a(2,1)=2
a(2,2)=23
a(2,3)=-1
a(2,4)=5
a(2,5)=-2
a(2,6)=82.4
a(3,1)=-1
a(3,2)=3
a(3,3)=92
a(3,4)=-5
a(3,5)=1
a(3,6)=-764.9
a(4,1)=1
a(4,2)=2
a(4,3)=1
a(4,4)=27
a(4,5)=3
a(4,6)=-8.9
a(5,1)=-4
a(5,2)=-6
a(5,3)=-2
a(5,4)=8
a(5,5)=41
a(5,6)=-201.9
x(1) =1
x(2) =1
x(3) =1
x(4) =1
x(5) =1
print *, "Augmented matrix"
call printmatrix(a,n,n+1)
call gausssiedel(a,n,x)
print *, "X:"
call printmatrix(x,n,1)
end program main

subroutine gausssiedel(a,n,x)
implicit none
integer,intent(in) ::n
real,dimension(n,n+1),intent(inout)::a
real,dimension(n),intent(inout)::x
real,dimension(n) :: x1
real:: s,epsilon,error
integer :: p,j,c,k
!call printmatrix(a,n,n+1) 
epsilon =0.0001
do k =1, 1000
	c=0
	x1=x
	do p = 1,n,1			!One iteration of all x
		s=0
		do j =1,n
			s = s+ a(p,j)*x(j)	!Using the latest values
		end do
		s = s - a(p,p)*x(p)
		x(p) = (a(p,n+1) - s)/a(p,p)
	end do				!end of iteration
	
	do p =1,n,1
		error =abs((x(p)-x1(p))/x1(p))
		if( error .lt. epsilon) then
			c =c+1
		end if
	end do
	
	if (c .eq. n) then
		EXIT
	end if
end do


!call printmatrix(a,3,4)
end subroutine gausssiedel

subroutine printmatrix(a,rows,col)
implicit none
integer,intent(in) :: rows,col
real,dimension(rows,col),intent(in) :: a
integer :: rowcount,colcount
10 format(f7.2)

do rowcount =1,rows
	do colcount =1,col
	write(*,10,advance ='no')a(rowcount,colcount)
	end do
	write(*,*)
end do
end subroutine printmatrix

