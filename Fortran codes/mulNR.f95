program main
implicit none
real,dimension(:,:),allocatable :: x,f,xprev
real,dimension(:,:),allocatable :: jac,jinv
real,dimension(:,:),allocatable ::pro
integer :: i,n,c,j
n=2
allocate(x(n,1))
allocate(xprev(n,1))
allocate(pro(n,1))
allocate(f(n,1))
allocate(jac(n,n))
allocate(jinv(n,n))

read(*,*) x
c=0
do i=1,1000
xprev =x
call func(x,n,f)
call jacobian(f,jac,n,x)
call inverse(jac,n,jinv)
call mul(jinv,f,n,1,pro) 
x =x-pro

do j=1,2
	if (abs((x(j,1)-xprev(j,1))/xprev(j,1)) .le. 0.0001 ) then
		c =c+1
	end if
end do
if (c .eq. n) EXIT
c=0
end do
print *, "No. of iterations =",i
print *, "x: "
call printmatrix(x,n,1)

end program main

subroutine func(x,n,f)
real,dimension(n),intent(inout)::x,f
integer,intent(in) ::n
f(1)= -2*x(1)**3 -8*x(1) +4*x(2) +4
f(2)= -4*x(1) + x(2)**2 + 3*x(2) +1
end subroutine func

subroutine jacobian(f,jac,n,x)
real,dimension(n),intent(inout)::f,x
real,dimension(n,n),intent(inout)::jac
jac(1,1) = -6*x(1)**2 -8
jac(1,2) = 4
jac(2,1) = -4
jac(2,2) = 2*x(2) +3

end subroutine jacobian


subroutine inverse(a,n,id)
implicit none
integer,intent(in) ::n
real,dimension(n,(n)),intent(inout)::a
real,dimension(n,n),intent(out)::id
real::pivot,factor,f,piv,temp,temp1
integer :: pivotrow,pivotcol,stepcount,rowcount,colcount,i,j,k
!call printmatrix(a,3,4)
do i=1,n
	id(i,i)=1
end do
!call printmatrix(id,n,n)
	do stepcount = 1,n-1,1
		pivotrow =stepcount
		pivotcol = stepcount
		pivot = a(pivotrow,pivotcol)
		do while(pivot .eq. 0) 
			do i =1,n
				temp = a(pivotrow,i)
				a(pivotrow,i)=a(pivotrow+1,i)
				a(pivotrow+1,i)=temp
			end do
			do i=1,n
				temp1 = id(pivotrow,i)
				id(pivotrow,i)=id(pivotrow+1,i)
				id(pivotrow+1,i)=temp1
			end do

			pivot = a(pivotrow,pivotcol)
		end do	
		do rowcount = pivotrow +1 , n,1
		factor = a(rowcount,pivotcol)/pivot
			do colcount = pivotcol, n
				a(rowcount,colcount) = a(rowcount,colcount) - (factor * a(pivotrow,colcount))
			end do
			do colcount = 1, n
				id(rowcount,colcount) = id(rowcount,colcount) - (factor * id(pivotrow,colcount))
			end do
		end do	
	end do
	do i=1 ,n
		f = a(i,i)
		do j=i,n
			a(i,j) =a(i,j)/f
			
		end do
		do j=1,n
			id(i,j)=id(i,j)/f
		end do
	end do
	
	do i = n,2,-1
		piv =a(i,i)
		do j=i-1,1,-1
			f =a(j,i)
			do k=1,n
				a(j,k) =a(j,k) - f*a(i,k)
			end do
			do k=1,n
			id(j,k) = id(j,k) - f*id(i,k)
			end do

		end do
	end do	
	
!print *, "Matrix after Gauss-Jordon Elimination :"
!call printmatrix(a,n,n+1)

!print *, "Inverse"
!c!all printmatrix(id,n,n)
!call mul(a1,id,n,pro)
!call printmatrix(pro,n,n)
end subroutine inverse

subroutine mul(a,b,n,m,pro)
implicit none
integer,intent(in) ::n,m
real,dimension(n,n),intent(inout)::a
real,dimension(n,m),intent(inout)::b
real,dimension(n,m),intent(out)::pro
real::pivot,factor,f,piv,temp,s
integer :: i,j,k
!call printmatrix(a,3,4)
	do i=1,n
		do j=1,m
			s=0
			do k=1,n
				s = s+ a(i,k)*b(k,j)
			end do
			pro(i,j) =s
		end do
	end do


end subroutine mul

subroutine printmatrix(a,rows,col)
implicit none
integer,intent(in) :: rows,col
real,dimension(rows,col),intent(in) :: a
integer :: rowcount,colcount
10 format(f7.5)

do rowcount =1,rows
	do colcount =1,col
	write(*,10,advance ='no')a(rowcount,colcount)
	end do
	write(*,*)
end do
end subroutine printmatrix