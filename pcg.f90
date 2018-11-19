subroutine show_matrix(B,m,k)
    implicit none
    integer m,k,j
    real(kind=16), dimension(m,k) :: B
	do j=1,size(B,1) 
      write(*,'(20G12.4)') B(j,:) 
    end do      
    write(*,*) 
end subroutine show_matrix

program pcg

endprogram pcg
