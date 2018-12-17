subroutine my_fun2(x1,x2,f1,f2)
	implicit none 
	real::x1,x2,f1,f2
	f1=x1+(1.0/20.0)*(4-8*x1+4*x2-2*x1**3)
	f2=x2-(1.0/24.0)
end subroutine my_fun2
