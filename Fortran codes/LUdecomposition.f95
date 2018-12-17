subroutine ludecomp(a,n,x)
implicit none
integer,intent(in) ::n
real,dimension(n,n+1),intent(inout)::a
real,dimension(n),intent(out)::x
real,dimension(n,n) :: l,u
real,dimension(n) :: b,y
real::pivot,factor,temp1,temp2
integer :: pivotrow,pivotcol,stepcount,rowcount,colcount,i,j
do i=1,n
	l(i,i) =1
	do j=i+1,n
		l(i,j)=0
	end do
end do
do i=1,n
	do j=1,n
		u(i,j)=a(i,j)
	end do
end do
do i=1,n
	b(i)=a(i,n+1)
end do
!call printmatrix(a,3,4)
	do stepcount = 1,n-1,1
		pivotrow =stepcount
		pivotcol = stepcount
		pivot = u(pivotrow,pivotcol)
		do while(pivot .eq. 0) 
			do i =1,n
				temp1 = u(pivotrow,i)
				u(pivotrow,i)=u(pivotrow+1,i)
				u(pivotrow+1,i)=temp1
			end do
			temp2=b(pivotrow)
			b(pivotrow)=b(pivotrow+1)
			b(pivotrow+1)=temp2
		pivot = u(pivotrow,pivotcol)
		end do
		
		do rowcount = pivotrow +1 , n,1
		l(rowcount,pivotcol) = u(rowcount,pivotcol)/pivot
			do colcount = pivotcol, n
				u(rowcount,colcount) = u(rowcount,colcount) - (l(rowcount,pivotcol) * u(pivotrow,colcount))
			end do
		end do
	end do
print *,"U"
call printmatrix(u,n,n)
print *,"L"
call printmatrix(l,n,n)
call forwardsubstitution(l,b,n,y)
call backsubstitution(u,y,n,x)
end subroutine ludecomp


subroutine forwardsubstitution(a,b,n,y)
implicit none
integer,intent(in) :: n
real,dimension(n,n),intent(in) :: a
real,dimension(n),intent(in) ::b
real,dimension(n),intent(out)::y
real :: s
integer :: rowcount , colcount
s=0
y(1) =b(1)/a(1,1)
do rowcount = 2 , n,1
	do colcount = 1, rowcount-1,1
		s = s + a(rowcount,colcount)*y(colcount)
	end do
	y(rowcount) = (b(rowcount) -s)/a(rowcount,rowcount)
	s=0
end do 
print *,"y"

call printmatrix(y,n,1)
end subroutine forwardsubstitution


subroutine backsubstitution(a,b,n,x)
implicit none
integer,intent(in) :: n
real,dimension(n,n),intent(in) :: a
real,dimension(n),intent(in) ::b
real,dimension(n),intent(out)::x
real :: s
integer :: rowcount , colcount
s=0
x(n) =b(n)/a(n,n)
do rowcount = n-1 , 1,-1
	do colcount = n, rowcount+1,-1
		s = s + a(rowcount,colcount)*x(colcount)
	end do
	x(rowcount) = (b(rowcount) -s)/a(rowcount,rowcount)
	s=0
end do
print *,"x" 
call printmatrix(x,n,1)
end subroutine backsubstitution

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

program main
implicit none
real,dimension(:,:),allocatable :: a
real,dimension(:),allocatable :: x
integer :: n
n=5
allocate(a(n,n+1))
allocate(x(n))

a(1,1)=0
a(1,2)=-1
a(1,3)=2
a(1,4)=-3
a(1,5)=4
a(1,6)=-38.5
a(2,1)=2
a(2,2)=3
a(2,3)=-1
a(2,4)=5
a(2,5)=-2
a(2,6)=32.4
a(3,1)=-1
a(3,2)=3
a(3,3)=2
a(3,4)=-5
a(3,5)=1
a(3,6)=-17.9
a(4,1)=1
a(4,2)=2
a(4,3)=1
a(4,4)=2
a(4,5)=3
a(4,6)=-13.9
a(5,1)=-4
a(5,2)=-6
a(5,3)=-2
a(5,4)=8
a(5,5)=-1
a(5,6)=4.9
print *,"Aug A"
call printmatrix(a,n,n+1)
call ludecomp(a,n,x)
!call printmatrix(x,3,1)
end program main

