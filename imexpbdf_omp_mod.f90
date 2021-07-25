! ============================ Module Description =============================
! imexbdf_mod Implementation of IMEX-BDF method. Contains the subroutine:
!   1. imexbdf - IMEX-BDF timestepping code
!   1. initY   - computes initial conditions for imex integrator using 
!                fully-implicit polynomial iterator.
! =============================================================================

module imexpbdf_omp_mod

    ! Module Parameters
    use omp_lib
    use tools_mod, only: dp, tic, toc, linspace, constvec_r, deye, quadw, interpw, Zinv
    use imexpbdf_mod, only: imexpbdf_settings
    implicit none

    contains

    ! =========================================================================
    ! IMEXBDF   Implicit-Explicit PBDF integrator
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
    !           derived data type that contains nodes and alpha for imex-pbdf
    !
    !   y_out   (output) COMPLEX*16 array, dimension (N)
    !           Initial condition.
    !
    !   times   (output) DOUBLE array, dimension(4)
    !           times(1:2) = [cputime,clocktime] for timestepping loop
    ! =========================================================================

    subroutine imexpbdf_omp(L, N, tspan, y_in, Nt, settings, y_out, times)

        ! Arguments
        complex(dp), intent(in)         :: L(:)
        real(dp),    intent(in)         :: tspan(2)
        complex(dp), intent(inout)      :: y_in(:)
        type(imexpbdf_settings), intent(in) :: settings
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
        complex(dp), allocatable :: Y(:,:), YN(:,:), RHS(:,:)
        real(dp), allocatable :: PBDFW(:,:), EDW(:,:), empty(:)
        real(dp) :: h, r
        integer :: ode_dim, i, j, k, q, thread_id
        complex(dp) :: t0
        
        thread_id = -1;
        q = size(settings%z)
        ode_dim = size(L)
        allocate(Y(ode_dim, q), YN(ode_dim, q), RHS(ode_dim, q), empty(0))
        
        t0 = tspan(1)
        h  = (tspan(2)-tspan(1))/(real(Nt,8))
        r  = h / settings%alpha

        call tic();
        
        call initY(L, N, h, settings, y_in, t0, Y, YN)

        ! PBDF coefficients
        PBDFW = transpose(& 
            interpW( &
                settings%z, &
                [ settings%z(q) + settings%alpha ], &
                settings%z + settings%alpha &
            ) &
        ) 
        ! extrapolated derivative coefficients
        EDW = matmul( &
            transpose(interpW(settings%z, empty, [ settings%z(q) + settings%alpha ])), &
            PBDFW(q+1:q+1,:) &
        )

        ! timestepping loop
        !$ call omp_set_num_threads(q)
        !$omp parallel private(thread_id, t0)
        !$ thread_id = omp_get_thread_num() + 1
        do i = 1, Nt
            
            ! jth thread computes jth explicit method rhs
            RHS(:,thread_id) = matmul(Y, PBDFW(1:q, thread_id)) + matmul(YN, r * EDW(:,thread_id));
            !$omp barrier
            
            if(thread_id == q) then ! thread leader computes y_q^{[n+1]} (the only implicit output)
                Y(:, q) = RHS(:, q) / (1 - r * PBDFW(q+1,q) * L)
            end if
            !$omp barrier
            
            if(thread_id < q) then ! worker threads compute remaining outputs using y_q^{[n+1]}
                Y(:, thread_id) = RHS(:, thread_id) + r * L * Y(:,q) * PBDFW(q+1, thread_id)        
            end if
            
            ! jth thread computes jth nonlinear term
            call N(t0 + h * settings%z(thread_id), Y(:,thread_id), YN(:,thread_id), thread_id)      
            t0 = tspan(1) + h * i;
            !$omp barrier
        enddo
        !$omp end parallel

        call toc(times(1:2))

        y_out = Y(:, 1)
        
    end subroutine imexpbdf_omp

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
        type(imexpbdf_settings), intent(in) :: settings
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
        real(dp), allocatable :: B(:,:)
        complex(dp), allocatable :: ImhB(:,:,:)
        integer :: nL, q, i, j, kappa
        real(dp) :: r
        
        ! verify valid nodes
        if(settings%z(1) .ne. -0.0_dp) then
            write (*,*) "Invalid nodes provided to IMEX-PBDF; translate nodes so that z(1) = 0."
            call exit(1)
        endif

        q  = size(settings%z)
        nL = size(L)

        B = transpose(quadW(settings%z, constvec_r(q, settings%z(1)), settings%z)); ! iterator matrix
        r = h / settings%alpha; ! compute r for propagator

        allocate(ImhB(q, q, nL))
        do i = 1, nL
            ImhB(:, :, i) = Zinv( deye(q) - r * L(i) * B )
        enddo

        ! zeroth-order guess
        Y(:,1) = y0
        call N(t0, y0, YN(:,1))

        do i = 2, q
            Y(:, i)  = Y(:,1);
            YN(:, i) = YN(:,1)
        enddo

        ! iterator
        kappa = q + 1
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
                call N(t0 + r * settings%z(j), Y(:,j), YN(:,j))
            enddo
        enddo

    end subroutine initY

end module imexpbdf_omp_mod