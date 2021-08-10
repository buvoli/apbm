function linspace(a,b,n) result(x)
! Returns row vector x of n points linearly spaced between [a,b]
! === Parameters ===
! a - (real) left endpoin
! b - (real) right endpoint
! === Output ===
! x - (array) array of equispaced points
real(kind=8), intent(in)  :: a,b
integer,      intent(in)  :: n
real(kind=8), allocatable :: x(:)
integer :: i
real(kind=8) :: s
allocate(x(n))

if(n == 1) then
    s = 0
else
    s = (b-a)/(n-1.d0) ! scaling factor
endif

do i=1,n
    x(i) = a + s*(i-1.d0)
enddo
end function linspace