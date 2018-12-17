subroutine quadratic(a,b,c)
 
 IMPLICIT NONE
 INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(4)

 REAL(dp) ,intent(in):: a, b, c
 REAL(dp) :: e, discriminant, rroot1, rroot2
 COMPLEX(dp) :: croot1, croot2
 
 WRITE(*,"(3(A,E23.15))") "Coefficients are: a = ", a, "   b = ", b, "   c = ", c
 e = 1.0e-9_dp
 discriminant = b*b - 4.0_dp*a*c
 
 IF (ABS(discriminant) < e) THEN
    rroot1 = -b / (2.0_dp * a)
    WRITE(*,*) "The roots are real and equal:"
    WRITE(*,"(A,E23.15)") "Root = ", rroot1
 ELSE IF (discriminant > 0) THEN
    rroot1 = -(b + SIGN(SQRT(discriminant), b)) / (2.0_dp * a)
    rroot2 = c / (a * rroot1)
    WRITE(*,*) "The roots are real:"
    WRITE(*,"(2(A,E23.15))") "Root1 = ", rroot1, "  Root2 = ", rroot2
 ELSE
    croot1 = (-b + SQRT(CMPLX(discriminant))) / (2.0_dp*a) 
    croot2 = CONJG(croot1)
    WRITE(*,*) "The roots are complex:" 
    WRITE(*,"(2(A,2E23.15,A))") "Root1 = ", croot1, "j ", "  Root2 = ", croot2, "j"
 END IF
end subroutine quadratic
program mainfunc
    implicit none
integer :: n=4,i,iteration
real, dimension(:), allocatable :: a
real, dimension(:), allocatable :: b
real, dimension(:), allocatable :: c
real::r=0.1,s=0.1,newr,news,epsilonr,epsilons
real:: tolerance = 1.0e-4
logical :: solutionFound = .false.
allocate(a(n+1))
allocate(b(n+1))
allocate(c(n+1))
a(5)=1
a(4)=0
a(3)=0
a(2)=-1
a(1)=-10
iteration=0
do while(solutionFound .eqv. .false.)
        b(n+1)=a(n+1)
        b(n)=a(n)+r*b(n+1)
        i=n-1
         do while(i.gt.0)
            b(i)=a(i)+r*b(i+1)+s*b(i+2)
            i=i-1
        end do
        c(n+1)=b(n+1)
        c(n)=b(n)+r*c(n+1)
        i=n-1
         do while(i.gt.1)
            c(i)=b(i)+r*c(i+1)+s*c(i+2)
            i=i-1
        end do
        newr=(b(1)*c(4)-b(2)*c(3))/(c(3)*c(3)-c(2)*c(4))
        news=(c(2)*b(2)-b(1)*c(3))/(c(3)*c(3)-c(2)*c(4))
        write(*,*)"Iteration number=",iteration," r(k)=",newr," s(k)=",news
        r=r+newr
        s=s+news
        epsilons=(news)/s
        epsilonr=(newr)/r
        write(*,*)"Iteration number=",iteration," epsilonr=",epsilonr," epsilons=",epsilons
        if((abs(epsilonr).LT.tolerance).AND.(abs(epsilons).LT.tolerance)) then
        solutionFound=.true.
        end if
        iteration=iteration+1
    end do
    write(*,*)"THe roots of polynomial are"

    call quadratic(1.0,-1.0*r,-1.0*s)
    call quadratic(b(5),b(4),b(3))
end program mainfunc