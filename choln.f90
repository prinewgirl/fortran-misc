subroutine show_matrix(B,m,k)
    implicit none
    integer m,k,j
    real(kind=16), dimension(m,k) :: B
	do j=1,size(B,1) 
      write(*,'(20G12.4)') B(j,:) 
    end do      
    write(*,*) 
end subroutine show_matrix

subroutine choldc(a,n,np,p)
implicit none
integer n,np
real*16 a(np,np),p(1:n)

integer i,j,k
real*16 summ

do i=1,n
  do j=i,n
    summ=a(i,j)
    do k=i-1,1,-1
      summ = summ -a(i,k)*a(j,k)
    enddo
    
    if(i==j) then
      if(summ <= 0) write(*,*) 'A is not symmetric positive defined'
      p(i) = sqrt(summ)
      else
        a(j,i)=summ/p(i)
      endif
    enddo
  enddo

return 
endsubroutine choldc  

subroutine cholsl(a,n,np,p,b,x)
integer n,np
real*16 a(np,np),b(n),p(n),x(n)
integer i,k
real summ
do i=1,n
summ = b(i)
do k=i-1,1,-1
summ = summ-a(i,k)*x(k)
enddo
x(i) = summ/p(i)
enddo
do i=n,1,-1
summ = x(i)
do k=i+1,n
summ = summ-a(k,i)*x(k)
enddo
x(i) = summ/p(i)
enddo
return
end subroutine cholsl
  
program cholesky

      implicit none
      integer i, it,n
      real*16,allocatable ::A(:,:),r0(:,:),p(:),x(:)

      !define data 
        !A*x = b
        n = 3
        
       allocate(A(n,n))
       allocate(r0(n,1))
       allocate(p(n))
       allocate(x(n))


        A(1,1) = 10
        A(1,2) = 1
        A(1,3) = 0
        
        A(2,1) = 1
        A(2,2) = 10
        A(2,3) = 1
        
        A(3,1) = 0
        A(3,2) = 1
        A(3,3) = 10
        
        r0(1,1) = 11
        r0(2,1) = 11
        r0(3,1) = 1
        
        call choldc(A,n,n,p)
        call cholsl(a,n,n,p,r0,x)
        write(*,*) x
endprogram cholesky
