program main
implicit none
integer, parameter :: n=5
double precision a(n,n), b(n), x(n)
integer i,j
! matrix a
  data (a(1,i), i=1,5) /  0,  -1,  2, -3, 4/
  data (a(2,i), i=1,5) /  2,  3,  -1, 5 , -2 /
  data (a(3,i), i=1,5) /  -1, 3, 2,-5,1 /
  data (a(4,i), i=1,5) /  1,2,1,2,3/
  data (a(5,i), i=1,5) /  -4,-6,-2 ,8,-1 /
! matrix b
  data (b(i),   i=1,5) /  -38.5,32.4,-17.9,-13.9,4.9 /

! printing  header and original equations
  write (*,200)
  do i=1,n
     write (*,201) (a(i,j),j=1,n), b(i)
  end do

  call gauss_elimination(a,b,x,n)

! printing matrix A and vector b after the elimination 
  write (*,202)
  do i = 1,n
     write (*,201)  (a(i,j),j=1,n), b(i)
  end do
! printing solutions
  write (*,203)
  write (*,201) (x(i),i=1,n)
200 format (' Gauss elimination ' &
    ,/,/,' Matrix A and vector b')
201 format (6f12.6)
202 format (/,' Matrix A and vector b after elimination')
203 format (/,' Solutions x(n)')
end program main

subroutine gauss_elimination(a,b,x,n)
implicit none 
integer n
double precision a(n,n), b(n), x(n)
double precision s(n)
double precision c, pivot, store
integer i, j, k, l

!forward elimination
do k=1, n-1
  do i=k,n                  
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

