program main
implicit none
real,dimension(:),allocatable ::x,y,a
real,dimension(:,:),allocatable :: xm,xt,b,p,p1,a1
integer :: n,i,j
real :: sx,sy,sxy,sx2
open(5,file="input.txt")
read (5,*) n
allocate(x(n))
allocate(y(n))
allocate(b(2,1))
allocate(a1(2,1))
allocate(p(2,2))
allocate(p1(2,3))
allocate(xm(n,2))
allocate(xt(2,n))
allocate(a(2))
do i=1,n
	read (5,*) x(i),y(i)
end do
call sigma(x,y,n,sx,sy,sxy,sx2)
print *, sx,sy,sxy,sx2
a(2) = ((n*sxy) - (sx*sy))/((n*sx2) - sx**2)
a(1) = -1*(a(2)*sx - sy)/n

do i=1,n
	xm(i,1) = 1
	xm(i,2) = x(i)
	xt(1,i) = 1
	xt(2,i) = x(i)
end do

call matmul(xt,xm,2,n,2,p)
call matmul(xt,y,2,n,1,b)
print *, "X transpose"
call printmatrix(xt,2,n)
print *,"X"
call printmatrix(xm,n,2)
print *, "XtX"
call printmatrix(p,2,2)
print *, "XtY"
call printmatrix(b,2,1)
do i=1,2
	do j=1,2
		p1(i,j)=p(i,j)
	end do
	p1(i,3)= b(i,1)
end do

call gaussel(p1,2,a1)
print *, "A"
call printmatrix(a1,2,1)


write (*,10) "y=",a(1),"+",a(2),"x"
10 format(A,f7.4,A,f7.4,A)

close(5)
end program main


subroutine printmatrix(a,rows,col)
implicit none
integer,intent(in) :: rows,col
real,dimension(rows,col),intent(in) :: a
integer :: rowcount,colcount
10 format(f9.4)

do rowcount =1,rows
	do colcount =1,col
	write(*,10,advance ='no')a(rowcount,colcount)
	end do
	write(*,*)
end do
end subroutine printmatrix

subroutine matmul(a,b,m,n,k,c)
real,dimension(m,n),intent(in)::a
real,dimension(n,k),intent(in)::b
real,dimension(m,k),intent(out)::c
integer,intent(in)::m,n,k
integer :: i,j,l
do i=1,m
	do j=1,k
		do l=1,n
			c(i,j) = c(i,j) + a(i,l)*b(l,j)
		end do
	end do
end do
end subroutine matmul

subroutine sigma(x,y,n,sx,sy,sxy,sx2)
real,dimension(n),intent(inout)::x,y
integer,intent(in) :: n
real,intent(out) :: sx,sy,sxy,sx2
sx = 0
sy =0
sxy =0
sx2 =0
do i =1,n
	sx= sx + x(i)
	sx2 = sx2 + x(i)**2
	sy = sy + y(i)
	sxy = sxy + x(i)*y(i)
end do

end subroutine sigma


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