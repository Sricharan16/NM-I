!**********************************************************************************************************************************
!  Main program
!**********************************************************************************************************************************

program linreg

   implicit none                                                                    
   
   integer ::I
     integer:: n
    integer, parameter  :: dbl = kind (0.0d0)                                                                                    
   real(dbl), dimension(:,:), allocatable :: x          
    real(dbl), dimension(:,:), allocatable :: a    
     real(dbl), dimension(:,:), allocatable :: b                                     
   open (unit=99, file='data.txt', status='old', action='read')
   read(99, *), n
   allocate(a(n,2))
   allocate(b(n,1))
   allocate(x(2,1))

! printing matrix A and vector b after the elimination 
!   write (*,202)
!   do i = 1,n
!      write (*,201)  (a(i,j),j=1,n), b(i)
!   end do
! ! printing solutions
!   write (*,203)
!   write (*,201) (x(i),i=1,n)
! 200 format (' Gauss elimination ' &
!     ,/,/,' Matrix A and vector b')
! 201 format (6f12.6)
! 202 format (/,' Matrix A and vector b after elimination')
! 203 format (/,' Solutions x(n)')

!allocate(x(n,2))

   do I=1,n
    a(I,1)=1
      read(99,*) a(I,2),b(I,1)
      write(*,*) a(I,2),b(I,1)
   enddo
   call gauss_elimination(a,b,x,n)
    write (*,203)
  write (*,201) (x(i,1),i=1,n)
200 format (' Gauss elimination ' &
    ,/,/,' Matrix A and vector b')
201 format (6f12.6)
202 format (/,' Matrix A and vector b after elimination')
203 format (/,' Solutions x(n)')

end program linreg

subroutine gauss_elimination(a,b,x,n)
implicit none 
integer n
double precision a(n,2), b(n,1), x(2,1)
double precision s(n)
double precision c, pivot, store
integer i, j, k, l

!forward elimination
do k=1, n-1
  do i=k,2                  
    s(i) = 0.0
    do j=k,n               
      s(i) = max(s(i),abs(a(i,j)))
    end do
  end do
  pivot = abs(a(k,k)/s(k))
  l = k
  do j=k+1,n
    if(abs(a(j,k)/s(j)) > pivot) then
      pivot = abs(a(j,k)/s(j))
      l = j
    end if
  end do
! Checking if the given linear system has a sigular matrix
  if(pivot == 0.0) then
    write(*,*) ' The matrix is sigular '
    return
  end if
if (l /= k) then
  do j=k,n
     store = a(k,j)
     a(k,j) = a(l,j)
     a(l,j) = store
  end do
  store = b(k)
  b(k) = b(l)
  b(l) = store
end if
!elimination
   do i=k+1,n
      c=a(i,k)/a(k,k)
      a(i,k) = 0.0
      b(i)=b(i)- c*b(k)
      do j=k+1,n
         a(i,j) = a(i,j)-c*a(k,j)
      end do
   end do
end do
!back substiturion 
x(n) = b(n)/a(n,n)
do i=n-1,1,-1
   c=0.0
   do j=i+1,n
     c= c + a(i,j)*x(j)
   end do 
   x(i) = (b(i)- c)/a(i,i)
end do

end subroutine gauss_elimination

