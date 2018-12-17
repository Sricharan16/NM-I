program picard 
implicit none
	real::x,f,tol=0.00001
	integer::i
	print *,"Enter the initial guess for x:"
	read *,x

	do i=1,1000
		print *,"Iteration No.",i
		call my_fun(x,f)
		print *,"X=",x,"Functional value =",f
		if(abs(x-f).le.tol)exit
		x=f
	end do 
	print *,"The solution is x =",x
end program picard

subroutine my_fun(x,f)
	implicit none 
	real::x,f
	f=(1.0/3.0)*exp(x)
end subroutine my_fun