
! INTERPW (x, xp, e) computes polynomial interpolation weights w(i,j) for
!
! f(e(i)) \approx \sum_{j=1}^n w(i,j) f(x(j)) + \sum_{j=1}^m w(i,j+n) f'(xp(j))
!
! == Parameters ===============================================================
!   x (vector)  - nodes where function value is known
!   xp (vector) - nodes where function derivative is known
!   e (vector)  - evaluation points
! =============================================================================

function interpW(x, xp, e) result(w)

    real(kind=dp), intent(in)  :: x(:)
    real(kind=dp), intent(in)  :: xp(:)
    real(kind=dp), intent(in)  :: e(:)
    real(kind=dp), allocatable :: w(:,:)
    
    real(kind=dp), allocatable :: V(:,:), rhs(:)
    integer, allocatable :: P(:)
    integer :: n, m, dim, i, j, ne

    ne = size(e)

    n = size(x)
    m = size(xp)
    dim = n + m

    allocate(V(dim, dim), P(dim), rhs(dim), w(ne, dim))

    do i = 1, dim
        do j = 1, n
            V(i, j) = x(j)**(i - 1);
        enddo
        do j = 1, m
            V(i, j + n) = (i - 1) * xp(j) ** max(0,(i-2));
        enddo
    enddo

    call DLU(V, P) ! compute LU factorization of V
    rhs = linspace(0.0_dp, dim - 1.0_dp, dim)
    do i = 1, ne
        w(i, :) = DLUS(V, P, e(i) ** rhs)
    enddo

end function interpW