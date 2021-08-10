function radaupts(n, I_in) result(x)
! Returns legendre points on interval I by applying the recurrence relation for the Legendre polynomials and their derivatives to perform Newton iteration on the WKB approximation to the roots.
!This code was adapted from Chebfun (https://www.chebfun.org) functions radaupts.m and jacpts.m
! === Parameters ===========================================
! n - number of radau points
! I_in - (array) 2x1 array containing inteval
! === Output ==============================================
! x - (array) nx1 array of radau points

integer, intent(in) :: n
real(dp), intent(in), optional :: I_in(2)
real(dp) :: I(2)
real(dp), allocatable :: x(:)

if(present(I_in)) then
    I = I_in
else
    I = (/ -1.d0, 1.d0 /)
endif

allocate(x(n))
x(1) = -1

if ( n > 1 ) then
    x(2:n) = jacpts(n - 1, 0, 1);
end if

! scale according to interval
x = (I(2) - I(1)) / 2.d0 * x + (I(2) + I(1)) / 2.d0

end function radaupts

function jacpts(n, a, b) result(x)

    integer, intent(in) :: n, a, b
    real(kind=8) :: I(2)
    real(kind=8), allocatable :: x(:), x1(:), x2(:)

    allocate(x(n))

    x1 = rec_main(n, a, b, .True.)
    x2 = rec_main(n, b, a, .False.)
    
    if( size(x1) + size(x2) .ne. n ) then
        call exit(1)
    end if

    x(1:size(x2)) = -1 * x2(size(x2):1:-1)
    x(size(x2)+1:n) = x1

end function jacpts

function rec_main(n, a, b, flag) result(x)

! parameters
integer, intent(in) ::  n, a, b
LOGICAL, intent(in)  :: flag

! local variables
integer :: num_x, counter, i
real(dp), allocatable :: x(:), P(:), Pp(:), dx(:), eps

! Asymptotic formula (WKB) - only positive x.
if ( flag ) then
    num_x = ceiling(real(n, dp) / 2.0_dp)
else
    num_x = floor(real(n, dp) / 2.0_dp) 
end if

allocate(x(num_x), P(num_x), Pp(num_x), dx(num_x))
do i = 1, num_x
    x(i) = num_x - i + 1
end do

x = (2.0_dp * x + a - 0.5_dp) * PI / (2.0_dp * n + a + b + 1.0_dp)
x = x + 1 / (2.0_dp * n + a + b + 1.0_dp)**2 * ((.25_dp - a**2_dp) * cotan(.5_dp * x) - (.25_dp - b**2_dp) * tan(.5_dp * x))
x = cos(x);

! Newton Interations
eps = .00000000000000023_dp
dx  = sqrt(eps) / 500
counter = 0

do while ( maxval(abs(dx)) > sqrt(eps) / 1000 .and. counter < 10 )
    counter = counter + 1
    call eval_Jac(x, n, a, b, P, PP)
    dx = -P / PP
    x = x + dx
enddo

end function rec_main

! EVALJAC   Evaluate Jacobi polynomial and derivative via recurrence relation.

subroutine eval_Jac(x, n, a, b, P, Pp)

! parameters
real(dp), intent(in) :: x(:)
integer, intent(in) :: n, a, b
real(dp), intent(out) :: P(:), Pp(:)

! local variables
integer :: ab
real(dp), allocatable :: P_m1(:), P_p1(:), PP_m1(:), PP_p1(:)
real(dp) :: C1, C2, C3, C4
integer k, num_x

num_x = size(x)
allocate(P_m1(num_x), P_p1(num_x), PP_m1(num_x), PP_p1(num_x))

! Initialise:
ab   = a + b
P    = .5 * (a - b + (ab + 2) * x)  
P_m1  = 1 
Pp   = .5 * (ab + 2)         
Pp_m1 = 0 

if ( n == 0 ) then
    P = P_m1 
    Pp = Pp_m1 
end if

do k = 1, n-1
    ! Constants
    C1 = 2*(k + 1)*(k + ab + 1)*(2*k + ab)
    C2 = (2*k + ab + 1)*(a**2 - b**2)
    C3 = (2*k + ab) * (2*k + ab + 1) * (2*k + ab + 2)
    C4 = 2*(k + a)*(k + b)*(2*k + ab + 2)

    ! Recurrence:
    P_p1 = ( (C2+C3*x)*P - C4*P_m1 ) / C1;
    Pp_p1 = ( (C2+C3*x)*Pp + C3*P - C4*Pp_m1 ) / C1;

    ! Update:
    P_m1 = P; 
    P = P_p1;  
    Pp_m1 = Pp; 
    Pp = Pp_p1;
end do
    
end subroutine eval_Jac