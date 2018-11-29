subroutine show_matrix(B,m,k)
    implicit none
    integer :: m,k,j
    real(kind=16), dimension(m,k) :: B
	do j=1,size(B,1) 
      write(*,'(20G12.4)') B(j,:) 
    end do      
    write(*,*) 
end subroutine show_matrix

program pcg
      implicit none
      real(kind=16) :: tol,alphak,ro0, ro, betak
      integer :: n,i,k
      real*16,allocatable ::A(:,:),Mi(:,:), b(:,:),x0(:,:)
      real*16,allocatable :: r0(:,:),p0(:,:)
      real*16,allocatable :: zk(:,:),gk(:,:)
    
        n = 3
      
       allocate(A(n,n))
       allocate(b(n,1))
       allocate(Mi(n,n))
       allocate(x0(n,n))
       allocate(r0(n,1))
       allocate(p0(n,1))
       
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
        
      !algorithm
     
     do i=1,n
       Mi(i,i) = 1/A(i,i)
     enddo
     
      r0 = b - matmul(A,x0)
      p0 = matmul(Mi,r0)
      ro0 = sum(r0 * matmul(Mi,r0))
      
      
      do i=1,k
      zk = matmul(A,p0)
      alphak = sum(p0*r0)/sum(zk*p0)
      x0 = x0+alphak*p0
      r0 = r0 - alphak*zk
      gk = matmul(Mi,r0)
      ro = ro0
      ro0 = sum(r0*gk)
      betak = ro0/ro
      p0 = gk + betak*p0
      
      enddo
  call show_matrix(x0,n,1)
end program pcg
