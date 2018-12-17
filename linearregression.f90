!**********************************************************************************************************************************
!  Main program
!**********************************************************************************************************************************

program linreg

   implicit none                                                                    
   
   integer :: n,i
    integer, parameter  :: dbl = kind (0.0d0)                                        

   real(dbl)           ::  b                                                        
   real(dbl)           ::  m                                                        
  ! real(dbl)           ::  n = 0.0d0   
   real(dbl), dimension(:,:), allocatable :: x                                        
   real(dbl)           ::  r                                                        
   character (len=80)  ::  str                                                      
   real(dbl)           ::  sumx  = 0.0d0                                            
   real(dbl)           ::  sumx2 = 0.0d0                                            
   real(dbl)           ::  sumxy = 0.0d0                                            
   real(dbl)           ::  sumy  = 0.0d0                                            
   real(dbl)           ::  sumy2 = 0.0d0                    
   open (unit=99, file='data.txt', status='old', action='read')
   read(99, *), n
   allocate(x(n,2))

   do I=1,n
      read(99,*) x(I,:)
      write(*,*) x(I,:)
   enddo
   write(*,*) x(1,1),x(1,2)

   do  i=1,n                                                                     
     ! n = n + 1.0d0                                                                 
      sumx  = sumx + x(i,1)                                                              
      sumx2 = sumx2 + x(i,1) * x(i,1)                                                         
      sumxy = sumxy + x(i,1) * x(i,2)                                                         
      sumy  = sumy + x(i,2)                                                            
      sumy2 = sumy2 + x(i,2) * x(i,2)                                                         
   end do
    write(*,*)sumx
    write(*,*)sumx2
     write(*,*)sumxy
      write(*,*)sumy
       write(*,*)sumy2
        write(*,*)n
   m = (n * sumxy  -  sumx * sumy) / (n * sumx2 - sumx**2)                          
   b = (sumy * sumx2  -  sumx * sumxy) / (n * sumx2  -  sumx**2)                   
   r = (sumxy - sumx * sumy / n) /                                     &            
                     sqrt((sumx2 - sumx**2/n) * (sumy2 - sumy**2/n))

   write (unit=*, fmt="(/a,es15.6)") " Slope        m = ", m                        
   write (unit=*, fmt="(a, es15.6)") " y-intercept  b = ", b
   write (unit=*, fmt="(a, es15.6)") " Correlation  r = ", r

end program linreg