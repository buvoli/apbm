! ============================ Module Description =============================
! imexbdf_mod Implementation of IMEX-BDF method. Contains the subroutine:
!   1. imexbdf - IMEX-BDF timestepping code
!   1. initY   - computes initial conditions for imex integrator using 
!                fully-implicit polynomial iterator.
! =============================================================================

module imexbdf_mod

    ! Module Parameters
    use tools_mod, only: dp, tic, toc, linspace, constvec_r, deye, quadw, interpw, Zinv
    implicit none
    type imexbdf_settings
        integer :: order
    end type imexbdf_settings

    contains

    ! =========================================================================
    ! IMEXBDF   Implicit-Explicit BDF integrator
    !
    ! Arguments
    !
    !   L       (input) COMPLEX*16 array, dimension(n)
    !           vector cooresponding to PDE linear operator
    !
    !   N       (input) SUBROUTINE:
    !           cooresponds to PDE nonlinear operator. Must be of the form:
    !               subroutine N(t,y_in,N_out)
    !                   real(kind=8),    intent(in)  :: t
    !                   complex(kind=8), intent(in)  :: y_in(:)
    !                   complex(kind=8), intent(out) :: N_out(:)
    !               end subroutine N
    !
    !   tspan   (input) DOUBLE array, dimension(2)
    !           contains left and right integration bounds
    !
    !   y_in    (input) COMPLEX*16 array, dimension (N)
    !           Initial condition.
    !
    !   settings (input) IMEXBDF_SETTINGS
    !           derived data type that contains imex-bdf order
    !
    !   y_out   (output) COMPLEX*16 array, dimension (N)
    !           Initial condition.
    !
    !   times   (output) DOUBLE array, dimension(4)
    !           times(1:2) = [cputime,clocktime] for timestepping loop
    ! =========================================================================

    subroutine imexbdf(L, N, tspan, y_in, Nt, settings, y_out, times)

        ! Arguments
        complex(dp), intent(in)         :: L(:)
        real(dp),    intent(in)         :: tspan(2)
        complex(dp), intent(inout)      :: y_in(:)
        type(imexbdf_settings), intent(in)  :: settings
        integer,     intent(in)         :: Nt
        real(dp),    intent(out)        :: times(4)
        complex(dp), intent(out)        :: y_out(size(y_in))
        ! define interface subroutine function N
        interface
            subroutine N(t, yh_in, N_out, thread_id)
            import :: dp
            complex(dp), intent(in)  :: t
            complex(dp), intent(in)  :: yh_in(:)
            complex(dp), intent(out) :: N_out(:)
            integer, intent(in), optional :: thread_id
            end subroutine N
        end interface

        ! Local Variables
        complex(dp), allocatable :: Y(:,:), YN(:,:), RHS(:)
        real(dp), allocatable :: z(:), BDFW(:,:), EDW(:,:), empty(:)
        real(dp) :: h
        integer :: ode_dim, i, j, k, order
        complex(dp) :: t0
        
        order = settings%order
        ode_dim = size(L)
        allocate(Y(ode_dim, order), YN(ode_dim, order), RHS(ode_dim), empty(0))
        
        t0 = tspan(1)
        h  = (tspan(2)-tspan(1))/(real(Nt,8))

        call tic();
        
        call initY(L, N, h, settings, y_in, t0, Y, YN)

        ! initialize IMEX-BDF coefficients
        z = linspace(0.0_dp, order - 1.0_dp, order)
        BDFW = transpose(interpW(z, [ z(order) + 1 ], [ z(order) + 1 ])) ! weights
        EDW  = BDFW(order+1,1) * transpose(interpW(z, empty, [ z(order) + 1 ])); ! extrapolated derivative

        ! timestepping loop
        do i = 1, Nt
            RHS = matmul(Y, BDFW(1:order, 1)) + matmul(YN, h * EDW(:,1));
           
            Y(:,  1:order-1) = Y(:,  2:order)
            YN(:, 1:order-1) = YN(:, 2:order)
            
            Y(:, order) = RHS / (1 - h *  BDFW(order+1,1) * L)
            call N(t0 + h * z(order), Y(:,order), YN(:,order))
            
            t0 = tspan(1) + h * i;
        enddo

        call toc(times(1:2))

        y_out = Y(:, 1)
        
    end subroutine imexbdf

    ! =========================================================================
    ! INITY   Initializes starting values for IMEX-BDF integrator using a 
    !         fully-implicit polynomial iterator.
    !
    ! Arguments
    !
    !   L       (input) COMPLEX*16 array, dimension(n)
    !           vector cooresponding to PDE linear operator
    !
    !   N       (input) SUBROUTINE:
    !           cooresponds to PDE nonlinear operator. Must be of the form:
    !               subroutine N(t,y_in,N_out)
    !                   real(kind=8),    intent(in)  :: t
    !                   complex(kind=8), intent(in)  :: y_in(:)
    !                   complex(kind=8), intent(out) :: N_out(:)
    !               end subroutine N
    !
    !   h   (input) DOUBLE array, dimension(2)
    !           contains left and right integration bounds
    !
    !   y_in    (input) COMPLEX*16 array, dimension (N)
    !           Initial condition.
    !
    !   settings (input) IMEXDIRK_TABLEAU
    !           derived data type that contains the IMEX RK tableau
    !
    !   y_out   (output) COMPLEX*16 array, dimension (N)
    !           Initial condition.
    !
    !   times   (output) DOUBLE array, dimension(4)
    !           times(1:2) = [cputime,clocktime] for timestepping loop
    ! =========================================================================

    subroutine initY(L, N, h, settings, y0, t0, Y, YN)

        ! Arguments
        complex(dp), intent(in) :: L(:), t0
        real(dp),    intent(in) :: h
        type(imexbdf_settings), intent(in) :: settings
        complex(dp), intent(inout) :: y0(:)
        complex(dp), intent(inout) :: Y(:,:), YN(:,:)
        interface ! interface for subroutine N
            subroutine N(t, yh_in, N_out, thread_id)
            import :: dp
            complex(dp), intent(in)  :: t
            complex(dp), intent(in)  :: yh_in(:)
            complex(dp), intent(out) :: N_out(:)
            integer, intent(in), optional :: thread_id
            end subroutine N
        end interface
    
        ! Variables
        real(dp), allocatable :: B(:,:), z(:)
        complex(dp), allocatable :: ImhB(:,:,:)
        integer :: nL, q, i, j, kappa
        real(dp) :: r
        
        q  = settings%order;
        nL = size(L)
        
        z = linspace(0.0_dp, 2.0_dp, q)
        B = transpose(quadW(z, constvec_r(q, 0.0_dp), z)); ! iterator matrix
        r = h * (q - 1) / 2 ! convert from r to h, (node interval has width h (q-1), r = width/2 )

        allocate(ImhB(q, q, nL))
        do i = 1, nL
            ImhB(:, :, i) = Zinv( deye(q) - r * L(i) * B );
        enddo

        ! zeroth-order guess
        Y(:,1) = y0
        call N(t0, y0, YN(:,1))

        do i = 2, q
            Y(:, i)  = Y(:,1);
            YN(:, i) = YN(:,1)
        enddo

        ! iterator
        kappa = q
        do i = 1, kappa
            ! ---> form explicit rhs
            YN(:,2:q) = matmul(YN, r * B(:, 2:q))
            do j = 2, q
                Y(:, j) = Y(:,1) + YN(:, j)
            enddo
            ! ---> invert implicit component
            do j = 1, Nl
                Y(j, :) = matmul(Y(j, :), ImhB(:, :, j))
            enddo
            ! ---> eval nonlinear term
            do j = 2, q
                call N(t0 + r * z(j), Y(:,j), YN(:,j))
            enddo
        enddo

    end subroutine initY

end module imexbdf_mod