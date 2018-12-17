!inverse and multiplication of matrices submodules to be written
subroutine mv_newtonraphson(n,x)
	implicit none 
	real::x1,f1,x2,f2
	real::det
	integer,intent(in)::n
	integer::i
	real,dimension(n,n)::a
	real,dimension(n,n),intent(inout)::x
	print *,"Enter the initial guess for x1:"
	read *,x1
	print *,"Enter the initial guess for x2:"
	read *,x2
	!call my_fun(x1,x2,f1,f2)
	a(1,1)=-8-6*x1**2
	a(1,2)=4
	a(2,1)=-4
	a(2,2)=3+2*x2
	det=a(1,1)*a(2,2)-a(2,1)*a(1,2)
	det=1/det
	x(1,1)=a(2,2)*det
	x(1,2)=-1*a(1,2)*det
	x(2,1)=-1*a(2,1)*det
	x(2,2)=a(1,1)*det
	do i=1,1000
		print *,"Iteration No.",i,"\n"
	end do 
end subroutine mv_newtonraphson

program main
implicit none 
	real,dimension(2,2)::x
	call mv_newtonraphson(2,x)
end program main


subroutine my_fun_f1(x1,x2,f1)
implicit none 
	real::x1,f1,x2
	f1=4-8*x1+4*x2-2*x1**3
	!f2=1-4*x1+3*x2+x2**2
end subroutine my_fun_f1

subroutine my_fun_f2(x1,x2,f2)
implicit none 
	real::x1,x2,f2
	!f1=4-8*x1+4*x2-2*x1**3
	f2=1-4*x1+3*x2+x2**2
end subroutine my_fun_f2