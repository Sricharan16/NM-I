subroutine newtonRaphson 
	implicit none
	real,external :: functionValue,firstDerivatuve
	real::x,xPrevious,slope,error,tolerance=1.0e-6
	integer ::iterationNumber=0
	print *,"Enter an intial value for x"
	read *,x
program mainProgram 
	implicit none 
	call newtonRaphson
end program mainProgram