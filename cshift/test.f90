program test_cshift
    integer, dimension(4,10) :: a,b
    integer :: i, j,k
    a=0
    do i=1,4
      do j=1,10
        a(i,j) = j
      enddo
    enddo

    do j=1,10
      write(*,*) (a(i,j), i=1,4)
    enddo

    a=cshift(a,2,dim=2)
    do j=1,10
      write(*,*) (a(i,j), i=1,4)
    enddo

end program test_cshift
