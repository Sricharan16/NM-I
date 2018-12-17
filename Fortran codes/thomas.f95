subroutine thomas(a,n,x)
implicit none
integer,intent(in) ::n
real,dimension(n,n+1),intent(inout)::a
real,dimension(n),intent(out)::x
real::pivot,factor
integer :: k,i,row,col
!call printmatrix(a,3,4)
	do i=1,n-1
		pivot =a(i,i)
		do k= i,i+1
			a(i,k)=a(i,k)/pivot
		end do
		a(i,n+1) =a(i,n+1)/pivot
		factor = a(i+1,i)
		do k=i ,i+2
			a(i+1,k) = a(i+1,k) - a(i,k)*factor
		end do
		if (i .le. n-2) then
			a(i+1,n+1)=a(i+1,n+1) - a(i,n+1)*factor	
		end if		
	end do

call printmatrix(a,n,n+1)
call backsubstitution(a,n,x)
end subroutine thomas

subroutine backsubstitution(a,n,x)
implicit none
integer,intent(in) :: n
real,dimension(n,n+1),intent(in) :: a
real,dimension(n),intent(out)::x
real :: s
integer :: rowcount , colcount
s=0
x(n) =a(n,n+1)/a(n,n)
do rowcount = n-1 , 1,-1
	do colcount = n, rowcount+1,-1
		s = s + a(rowcount,colcount)*x(colcount)
	end do
	x(rowcount) = (a(rowcount,n+1) -s)/a(rowcount,rowcount)
	s=0
end do 
!call printmatrix(a,n,n+1)
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

a(1,1)=3
a(1,2)=-1
a(1,3)=0
a(1,4)=0
a(1,5)=0
a(1,6)=2
a(2,1)=-1
a(2,2)=3
a(2,3)=-1
a(2,4)=0
a(2,5)=0
a(2,6)=1
a(3,1)=0
a(3,2)=-1
a(3,3)=3
a(3,4)=-1
a(3,5)=0
a(3,6)=1
a(4,1)=0
a(4,2)=0
a(4,3)=-1
a(4,4)=3
a(4,5)=-1
a(4,6)=1
a(5,1)=0
a(5,2)=0
a(5,3)=0
a(5,4)=-1
a(5,5)=3
a(5,6)=2
print *,"Initial"
call printmatrix(a,n,n+1)
print *,"A"
call thomas(a,n,x)
print *,"X"
call printmatrix(x,n,1)
end program main

