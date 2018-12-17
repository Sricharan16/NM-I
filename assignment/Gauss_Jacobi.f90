logical function diagonallyDominant(matrix, rows, columns)
implicit none
integer, intent(in) :: rows, columns
real, dimension(rows, columns), intent(in) :: matrix
integer :: rowCount, columnCount
real :: diagonalElement
real :: sumOfNondiagonalElements
diagonallyDominant = .true.
do rowCount = 1, rows
diagonalElement = abs(matrix(rowCount, rowCount))
sumOfNondiagonalElements = 0
do columnCount = 1, columns
if(columnCount /= rowCount) then
sumOfNondiagonalElements = sumOfNondiagonalElements +abs(matrix(rowCount, columnCount))
end if
end do
if(sumOfNondiagonalElements > diagonalElement) then
diagonallyDominant = .false.
return
end if
end do
end function diagonallyDominant

subroutine printMatrix2D(matrix, rows, columns)
implicit none
integer, intent(in) :: rows, columns
real, dimension(rows, columns), intent(in) :: matrix
integer :: rowCounter, columnCounter
do rowCounter = 1, rows
do columnCounter = 1, columns
write(*,10, advance='no') matrix(rowCounter, columnCounter)
end do
write(*,*)
end do
write(*,*)
10 format(f7.2)
end subroutine printMatrix2D

subroutine jacobi(aIn, n, xGuess, tolerance, x)
    implicit none
    logical, external :: diagonallyDominant

    integer, intent(in) :: n
    real, dimension(n, (n+1)), intent(in) :: aIn
    real, dimension(n), intent(in) :: xGuess
    real, intent(in) :: tolerance

    real, dimension(n), intent(out) :: x

    real, dimension(n, (n+1)) :: a
    real, dimension(n, n) :: coefficientMatrix
    real, dimension(n) :: xPrevious
    real, dimension(n) :: error

    real :: factor
    logical :: solutionFound = .false.
    integer :: iteration = 0
    integer :: rowCount, columnCount

    a = aIn
    write(*,*) "The Augmented Matrix"
    call printMatrix2D(a, n, (n+1))

    do rowCount = 1, n
        do columnCount = 1, n
        coefficientMatrix(rowCount, columnCount) = a(rowCount, columnCount)
        end do
    end do
    if(diagonallyDominant(coefficientMatrix, n, n) .eqv. .false.) then
        stop "Error (Jacobi Method): The coefficient matrix is not diagonallydominant."
    end if

    x = xGuess
    xPrevious = x
    do while(solutionFound .eqv. .false.)
         do rowCount = 1, n
            factor = 0
            do columnCount = 1, n
                if(columnCount /= rowCount) then
                    factor = factor + a(rowCount, columnCount) *xPrevious(columnCount)
                end if
            end do
            x(rowCount) = (1 / a(rowCount, rowCount)) * (a(rowCount, (n+1)) - factor)
        end do
        do rowCount = 1, n
            error(rowCount) = abs((x(rowCount) - xPrevious(rowCount)) /(x(rowCount)))
        end do
        write(*,10) "Iteration #", iteration
        write(*,*) "Estimated Solution: "
        call printMatrix2D(x, n, 1)
        write(*,*) "Error: "
        call printMatrix2D(error, n, 1)
            if(maxval(error) < tolerance) then
                 solutionFound = .true.
            end if
        xPrevious = x
    end do
10 format(a11, i2)
end subroutine jacobi



program mainIterative
    implicit none
integer :: n
real, dimension(:,:), allocatable :: a
real, dimension(:), allocatable :: x
real, dimension(:), allocatable :: xGuess
real, parameter :: tolerance = 1.0e-6
n = 5
allocate(a(n, n+1))
allocate(x(n))
allocate(xGuess(n))
a(1,1) = 15
a(1,2) = -1
a(1,3) = 2
a(1,4) = -3
a(1,5) = 4
a(1,6) = 8
a(2,1) = 2
a(2,2) = 23
a(2,3) = -1
a(2,4) = 5
a(2,5) = -2
a(2,6) = 82.4
a(3,1) = -1
a(3,2) = 3
a(3,3) = 92
a(3,4) = -5
a(3,5) = 1
a(3,6) = -764.9
a(4,1) = 1
a(4,2) = 2
a(4,3) = 1
a(4,4) = 27
a(4,5) = 3
a(4,6) = -8.9
a(5,1) = -4
a(5,2) = -6
a(5,3) = -2
a(5,4) = 8
a(5,5) = 41

a(5,6) = -201.9
xGuess = 0
call jacobi(a, n, xGuess, tolerance, x)
write(*,*) "Solution by Jacobi Iterative Method:"
call printMatrix2D(x, n, 1)

end program mainIterative
