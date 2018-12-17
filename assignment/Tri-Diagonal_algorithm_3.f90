      !n is the number of equations
      !b is the main diagonal
      !c is diagonal above main diagonal
      !a is diagonal below main diagonal 
      !d is right part of the linear equation system
      !x is the output Solution
subroutine solve_tridiag(a,b,c,d,x,n)
      implicit none
        integer,intent(in) :: n
        real(8),dimension(n),intent(in) :: a,b,c,d
        real(8),dimension(n),intent(out) :: x
        real(8),dimension(n) :: cp,dp
        real(8) :: m
        integer i
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
         do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           cp(i) = c(i)/m
           dp(i) = (d(i)-dp(i-1)*a(i))/m
         end do
         x(n) = dp(n)

        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)
        end do

    end subroutine solve_tridiag
 program main
implicit none
integer, parameter :: n=5
integer i
double precision a(n), b(n),c(n),d(n),x(n)
data (a(i), i=1,5) /  0,  -1,  -1, -1, -1/
data (b(i), i=1,5) /  3,  3,  3, 3, 3/
data (c(i), i=1,5) /  -1,  -1,  -1, -1, 0/
data (d(i), i=1,5) /  2,  1,  1, 1, 2/
call solve_tridiag(a,b,c,d,x,n)
write (*,203)
write (*,201) (x(i),i=1,n)
201 format (6f12.6)
203 format (/,' Solutions x(n)')
end program main