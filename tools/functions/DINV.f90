! Simplified calling sequence for lapack solver using existing LU decomposing (DGETRF)
function DINV(A) result(Ai)
    implicit none
    real(kind=dp), intent(in)  :: A(:,:)
    integer,         allocatable :: P(:)
    real(kind=dp), allocatable :: Ai(:,:), work(:)
    integer :: info, n

    n = size(A, 1)

    allocate(P(n))
    allocate(Ai(n,n), work(n))
    
    ! LU factorize
    Ai = A
    call DGETRF(n, n, Ai, n, P, info) !http://www.netlib.org/lapack/explore-3.1.1-html/dgetrf.f.html

    ! Invert
    call DGETRI(n, Ai, n, P, work, n, info) !http://www.netlib.org/lapack/explore-3.1.1-html/dgetri.f.html

end function DINV