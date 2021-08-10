
! QUADW(x,a,b) returns quadrature weights for the nodes x such that
!
!   \int^b(i)_a(i) f(x) dx \approx \sum_{j=1}^n w(j,i) f(x(i))
!
! == Parameters ===============================================================
!   x (vector) - nodes
!   a (vector) - left endpoints
!   b (vector) - right endpoints
! =============================================================================

function quadW(x, a, b) result(w)

    real(kind=dp) :: x(:)
    real(kind=dp) :: a(:)
    real(kind=dp) :: b(:)

    real(kind=dp), allocatable :: V(:,:), rhs(:), pwrs(:), W(:,:)
    integer, allocatable :: P(:)
    integer :: dim, i, j, num_ab
    
    dim = size(x)
    num_ab = size(a)
    allocate(V(dim, dim), P(dim), rhs(dim), pwrs(dim), W(num_ab,dim))

    ! construct Vandermond matrix
    do i = 1, dim
        do j = 1, dim
            V(i,j) = x(j) ** (i - 1)
        end do
    end do
    call DLU(V, P) ! compute LU factorization of V
    
    pwrs = linspace(1.0_dp, real(dim, dp), dim)    
    do i = 1, num_ab
        rhs = (b(i)**pwrs - a(i)**pwrs) / pwrs
        w(i, :) = DLUS(V, P, rhs)
    end do

end function quadW