subroutine show_matrix(B,m,k)
    implicit none
    integer :: m,k,j
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
      if (a(j,i) /= 0) then 
        a(j,i)=summ/p(i)
      endif  
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
  

program pcg
      implicit none
      real(kind=16) :: tol,alphak,ro0, ro, betak
      integer :: n,i,k
      real*16,allocatable ::A(:,:),Mi(:,:), b(:,:),x0(:,:)
      real*16,allocatable :: r0(:,:),p0(:,:),p(:)
      real*16,allocatable :: zk(:,:),gk(:,:)
    
        n = 3
      
       allocate(A(n,n))
       allocate(b(n,1))
       allocate(Mi(n,n))
       allocate(x0(n,n))
       allocate(r0(n,1))
       allocate(p0(n,1))
       allocate(p(n))
       allocate(zk(n,1))
       allocate(gk(n,1))


       
  write(*,*) ' type  the maximum number of iterations'
        read (*,*) k
        
        !define data 
        !A*x = b
        n = 3
        A(1,1) = 10
        A(1,2) = 1
        A(1,3) = 0
        
        A(2,1) = 1
        A(2,2) = 10
        A(2,3) = 1
        
        A(3,1) = 0
        A(3,2) = 1
        A(3,3) = 10
        
        b(1,1) = 11
        b(2,1) = 11
        b(3,1) = 1
        
        Mi = A
        
      !algorithm
       
      call choldc(Mi,n,n,p)
     
      r0 = b - matmul(A,x0)
      p0 = matmul(Mi,r0)
      ro0 = sum(r0 * matmul(Mi,r0))
      
      
      do i=1,k
      zk = matmul(A,p0)
      alphak = sum(p0*r0)/sum(zk*p0)
      x0 = x0+alphak*p0
      r0 = r0 - alphak*zk
      
      
      call cholsl(Mi,n,n,p,r0,gk)
        
      ro = ro0
      ro0 = sum(r0*gk)
      betak = ro0/ro
      p0 = gk + betak*p0
      
      enddo
  call show_matrix(x0,n,1)
end program pcg
