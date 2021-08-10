function DEYE(n) result(I)
! Returns identity matrix
! === Parameters ===
! n - matrix dimension
! === Output =======
! I - (array) nxn matrix
integer, intent(in) :: n
real(kind=dp), allocatable :: I(:,:)
integer :: j
allocate(I(n,n))
I = 0.0_dp
do j=1,n
    I(j,j) = 1
enddo
end function DEYE