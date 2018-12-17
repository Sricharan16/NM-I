subroutine gaussel(a,n,x)
implicit none
integer,intent(in) ::n
real,dimension(n,n+1),intent(inout)::a
real,dimension(n),intent(out)::x
real ::temp
real::pivot,factor
integer :: pivotrow,pivotcol,stepcount,rowcount,colcount,i
!call printmatrix(a,n,n+1)
	do stepcount = 1,n-1,1
		pivotrow =stepcount
		pivotcol = stepcount
		pivot = a(pivotrow,pivotcol)
		do while(pivot .eq. 0) 
			do i =1,n+1
				temp = a(pivotrow,i)
				a(pivotrow,i)=a(pivotrow+1,i)
				a(pivotrow+1,i)=temp
			end do
			pivot = a(pivotrow,pivotcol)
		end do		
		do rowcount = pivotrow +1 , n,1
				factor = a(rowcount,pivotcol)/pivot
				do colcount = pivotcol, n+1
					a(rowcount,colcount) = a(rowcount,colcount) - (factor * a(pivotrow,colcount))
				end do
		end do
	end do
print *, "Matrix after Gauss Elimiation :"
call printmatrix(a,n,n+1)
call backsubstitution(a,n,x)
end subroutine gaussel

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
print *, "Augmented Matrix"
call printmatrix(a,n,n+1)
call gaussel(a,n,x)
print *, "X"
call printmatrix(x,n,1)
end program main

