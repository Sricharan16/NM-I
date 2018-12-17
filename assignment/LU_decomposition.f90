module LU
contains
   subroutine ludcmp(A,n,indx,D,code)
     implicit none
     integer, parameter :: nmax = 100
     real, parameter :: tiny = 1.5D-16

     real*8, intent(inout), dimension(N,N) :: A
     integer, intent(in) :: N
     integer, intent(out) :: D, code
     integer, intent(out), dimension(N) :: indx
     !f2py depend(N) A, indx

     real*8  :: AMAX, DUM, SUMM, vv(NMAX)
     integer :: i, j, k, imax

     D=1; code=0

     do I=1,N
       amax=0.d0
       do J=1,N
         if (DABS(A(I,J)).gt.amax) amax=DABS(A(I,J))
       end do 
       if(amax.LT.tiny) then
         CODE = 1
         return 
       end if
       vv(i) = 1.d0 / amax
     end do 

     do J=1,N
       do I=1,J-1
         summ = A(I,J)
         do K=1,I-1
           summ = summ - A(I,K)*A(K,J) 
         end do
         A(i,j) = summ
       end do 
       amax = 0.d0
       do i=j,n
         summ = A(i,j)
         do k=1,j-1
           summ = summ - A(i,k)*A(k,j) 
         end do
         A(i,j) = summ
         dum = vv(i)*DABS(summ)
         if(dum.ge.amax) then
           imax = i
           amax = dum
         end if 
       end do  
       
       if(j.ne.imax) then
         do k=1,n
           dum = A(imax,k)
           A(imax,k) = A(j,k)
           A(j,k) = dum
         end do
         d = -d
         vv(imax) = vv(j)
       end if

       indx(j) = imax
       if(dabs(A(j,j)) < tiny) A(j,j) = tiny

       if(j.ne.n) then
         dum = 1.d0 / A(j,j)
         do i=j+1,n
           A(i,j) = A(i,j)*dum
         end do
       end if
     end do

     return 
end subroutine ludcmp


 subroutine lubksb(A, N, indx, B)
 integer, intent(in) :: N 
 real*8, intent(in), dimension(N,N) :: A
 integer, intent(in), dimension(N) :: indx
 real*8, intent(inout), dimension(N) :: B
 real*8  summ
 ii = 0

 do i=1,N
   ll = indx(i)
   summ = B(ll)
   B(ll) = B(i)
   if(ii.ne.0) then
     do j=ii,i-1
       summ = summ - A(i,j)*B(j)
     end do
   else if(summ.ne.0.d0) then
     ii = i
   end if
   B(i) = summ
 end do

 do i=N,1,-1
   summ = B(i)
   if(i < N) then
     do J=i+1,N
       summ = summ - A(i,j)*B(j)
     end do
   end if
   B(i) = summ / A(i,i)
 end do
return 
end subroutine lubksb
end module LU
program lu_solver
  use LU
  implicit none
!input matrix A
!solution vector
!temporary vector (n+1)
  real*8, pointer ::  A(:,:)   
  real*8, pointer ::  B(:)     
  real*8, pointer ::  temp(:)  
  integer,pointer ::  indx(:)  !integer vector (n)
  integer :: i, d, rc, n = 5
  allocate(A(n,n))
  allocate(B(n))
  allocate(temp(n+1))
  allocate(indx(n))
A(1,1) =0
A(1,2) = -1
A(1,3) = 2
A(1,4) = -3
A(1,5) = 4
A(2,1) = 2
A(2,2) = 3
A(2,3) = -1
A(2,4) = 5
A(2,5) = -2
A(3,1) = -1
A(3,2) = 3
A(3,3) = 2
A(3,4) = -5
A(3,5) = 1
A(4,1) = 1
A(4,2) = 2
A(4,3) = 1
A(4,4) = 2
A(4,5) = 3
A(5,1) = -4
A(5,2) = -6
A(5,3) = -2
A(5,4) = 8
A(5,5) = -1
! matrix b
B(1)=-38.5
B(2)=32.4
B(3)=-17.9
B(4)=-13.9
B(5)=4.9
!call LU decompositions
  call LUDCMP(A,n,indx,D,rc)

!call solver
  if (rc.eq.0) then
    call LUBKSB(A,n,indx,B)
  endif

!print results or error message
  if (rc.eq.1) then
    write(*,*) ' The matrix is singular and hence no solution'
  else
    write(*,*) 'By solving using LU decomposition, we get:'
      write (*,201) (B(i),i=1,n)
    201 format (6f12.6)
  end if

end program lu_solver