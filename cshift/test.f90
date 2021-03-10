program test_cshift
    integer, dimension(4,10) :: a,b
    integer :: i, j,k
    a=0
    do i=1,4
      do j=1,10
        a(j,i) = j
      enddo
    enddo

    write(*,*) a(1,1),a(1,2), a(1,3), a(1,4)
    write(*,*) a(2,1),a(2,2), a(2,3), a(2,4)
    write(*,*) a(3,1),a(3,2), a(3,3), a(3,4)
    a=cshift(a,2,dim=1)
    write(*,*) a(1,1),a(1,2), a(1,3), a(1,4)
    write(*,*) a(2,1),a(2,2), a(2,3), a(2,4)
    write(*,*) a(3,1),a(3,2), a(3,3), a(3,4)
    open(12, file="test.txt")
    call check(a)
    close(12)

end program test_cshift

subroutine check(a)
  integer,dimension(4,10) :: a
  integer :: i, j ,k

  a=cshift(a,2,dim=2)
  do j=1,10
    write(12,*) (a(i,j), i=1,4)
  enddo

end subroutine check
