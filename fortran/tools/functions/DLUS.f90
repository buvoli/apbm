! Simplified calling sequence for lapack solver using existing LU decomposing (DGETRF)
function DLUS(A,P,f) result(x)
    implicit none
    real(dp), intent(in)    :: A(:,:)
    real(dp), intent(in)    :: f(size(A,1))
    real(dp), allocatable   :: x(:)
    integer, intent(in)     :: P(size(A,1))
    integer                 :: info,n

    n = size(A,1)
    allocate(x(n))
    x = f;
    call DGETRS('N',n,1,A,n,P,x,n,INFO)  !http://www.netlib.no/netlib/lapack/complex16/zgetrs.f
    
end function DLUS