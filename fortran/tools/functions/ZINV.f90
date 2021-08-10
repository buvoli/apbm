! Simplified calling sequence for lapack solver using existing LU decomposing (DGETRF)
function ZINV(A) result(Ai)
    implicit none
    complex(kind=8), intent(in)  :: A(:,:)
    integer,         allocatable :: P(:)
    complex(kind=8), allocatable :: Ai(:,:), work(:)
    integer :: info, n

    n = size(A, 1)

    allocate(P(n))
    allocate(Ai(n,n), work(n))
    
    ! LU factorize
    Ai = A
    call ZGETRF(n, n, Ai, n, P, info) !http://www.netlib.no/netlib/lapack/complex16/zgetrf.f

    ! Invert
    call ZGETRI(n, Ai, n, P, work, n, info) !http://www.netlib.no/netlib/lapack/complex16/zgetri.f

end function ZINV