function constvec_r(n, val) result(v)
! Returns constant vector of val
! === Parameters ===
! n - vector dimenson
! val - (REAL) 2x1 array containing inteval
! === Output ===
! c - (array) nx1 array with all entries set to val
integer, intent(in) :: n
real(kind=dp), intent(in) :: val
real(kind=dp) :: v(n)
integer :: j
do j=1,n
    v(j) = val
enddo
end function constvec_r