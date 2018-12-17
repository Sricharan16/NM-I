program multivariable_ss 
implicit none
	real::x1,x2,f1,f2,tol=1e-06
	integer::i
	print *,"Enter the initial guess for x1:"
	read *,x1
	print *,"Enter the initial guess for x2:"
	read *,x2


	do i=1,1000
		print *,"Iteration No.",i
		call my_fun2(x1,x2,f1,f2)
		write(*,20) x1,x2
		print *," "
		if((abs(x1-f1).le.tol).and.(abs(x2-f2).le.tol))exit
		x1=f1
		x2=f2
	end do 
	20 format(f10.4,f10.4)
	print *,"The solution is x1 =",x1,"x2=",x2
end program multivariable_ss

subroutine my_fun2(x1,x2,f1,f2)
	implicit none 
	real::x1,x2,f1,f2
	f1=x1+(1.0/20.0)*(4-8*x1+4*x2-2*x1**3)
	f2=x2-(1.0/24.0)*(1-4*x1+3*x2+x2**2)
end subroutine my_fun2
