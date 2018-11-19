subroutine show_matrix(B,m,k)
    implicit none
    integer m,k,j
    real(kind=16), dimension(m,k) :: B
	do j=1,size(B,1) 
      write(*,'(20G12.4)') B(j,:) 
    end do      
    write(*,*) 
end subroutine show_matrix

subroutine choldc(a,n,np)
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
do i=1,n
  do j=1,n
    if(i<j) then
      a(i,j) = 0
    endif  
    if(i ==j) then
      a(i,j) = sqrt(a(i,j))
    endif
    enddo
enddo 

return 
endsubroutine choldc  


program cholesky

      implicit none
      integer i, it,n
      real*16,allocatable ::A(:,:),r0(:,:)    

      !define data 
        !A*x = b
        n = 3
        
       allocate(A(n,n))
       allocate(r0(n,1))
     

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
        
        call choldc(A,n,n)
        call show_matrix(A,n,n)
endprogram cholesky
