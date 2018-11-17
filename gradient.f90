subroutine show_matrix(B,m,k)
    implicit none
    integer m,k,j
    real(kind=16), dimension(m,k) :: B
	do j=1,size(B,1) 
      write(*,'(20G12.4)') B(j,:) 
    end do      
    write(*,*) 
end subroutine show_matrix


program gradient
      implicit none
      real(kind=16) :: error, r , lambda
      integer :: n,i, it
      real*16,allocatable ::A(:,:),b(:,:)
      real*16,allocatable :: v(:,:), v_bef(:,:)
       write(*,*) 'type it'
       read (*,*) it
       n = 3
       allocate(A(n,n))
       allocate(b(n,1))
       allocate(v(n,1))
       allocate(v_bef(n,1))

        
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
        
        !first guess: x = 0
        b = -b
        v(1:n,1) = 0
        do i=1,it
          lambda = sum(b*b)/sum(matmul(A,b)*b)
          v_bef = v
          v =v -lambda*b
          b = b-lambda*matmul(A,b)
          
        enddo
        error = norm2(v - v_bef)/norm2(v)
        
        write (*,*) 'solution'
        call show_matrix(v,n,1)
        
        write (*,*) 'error' 
        write (*,'(20G12.4)') error
    
        deallocate(A,b,v)
end program gradient

