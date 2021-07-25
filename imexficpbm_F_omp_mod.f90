! ============================ Module Description ==============================
! imexficpbm_F_mod -  Serial implementation of a fully-implicit, composite, 
!               IMEX, Adams PBMs with:
!                   Nodes: Left Sweeping real
!                   Explicit Active Index Set (AIS): PMFO
!                   Implicit Active Index Set (AIS): FI (outputs only)
!                   Endpoints: Fixed
!
! Contains subroutines:
!   1. iemx_pbm_omp - ETDPBM timestepping code
!   2. initY   - computes initial conditions for imex integrator using an
!                the iterator.
! ==============================================================================

module imexficpbm_F_omp_mod

    ! Module Parameters
    use omp_lib
    use tools_mod, only: dp, tic, toc, linspace, constvec_r, zeye, quadw, interpw, Zinv
    use imexficpbm_F_mod, only: imexfipbm_F_settings, imexficpbm_F_settings, imexpbm_F_coefficients, settings2Coefficient
    implicit none

    contains

    ! =========================================================================
    ! IMEXBDF   Implicit-Explicit IMEX PBM integrator
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

    subroutine imexficpbm_F_omp(L, N, tspan, y_in, Nt, settings, y_out, times)

        ! Arguments
        complex(dp), intent(in)         :: L(:)
        real(dp),    intent(in)         :: tspan(2)
        complex(dp), intent(inout)      :: y_in(:)
        type(imexficpbm_F_settings), intent(in) :: settings
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
        complex(dp), allocatable :: Y(:,:), YN(:,:)
        real(dp) :: h, r
        complex(dp) :: t0
        integer :: ode_dim, q, i, j, thread_id
        type(imexpbm_F_coefficients) :: propagator, iterator

        t0 = tspan(1)
        h  = (tspan(2)-tspan(1))/(real(Nt,dp))
        r  = h / settings%propagator%alpha

        q = size(settings%propagator%z)
        ode_dim = size(L)
        allocate(Y(ode_dim, q), YN(ode_dim, q))

        ! == Initialize coefficients ==========================================
        call tic()
        propagator = settings2Coefficient(settings%propagator, r, L)
        iterator   = settings2Coefficient(settings%iterator, r, L)
        call toc(times(3:4))

        ! == compute initial solution vector ==================================
        call initY(N, y_in, t0, r, iterator, settings%kappa, Y, YN)
        
        call tic()

        ! timestepping loop
        !$ call omp_set_num_threads(q)
        !$omp parallel private(thread_id, t0)
        !$ thread_id = omp_get_thread_num() + 1        
        do i = 1, Nt
    
            ! propagator
            call step(Y, YN, N, r, t0, propagator, thread_id)
        
            ! iterator
            do j = 1, settings%kappa
                call step(Y, YN, N, r, t0, iterator, thread_id)
            enddo

        enddo

        !$omp end parallel

        call toc(times(1:2))

        y_out = Y(:, 1)
        
    end subroutine imexficpbm_F_omp

    ! ! =========================================================================
    ! ! INITY   Initializes starting values for IMEX-BDF integrator using a 
    ! !         fully-implicit polynomial iterator.
    ! !
    ! ! Arguments
    ! !
    ! !   L       (input) COMPLEX*16 array, dimension(n)
    ! !           vector cooresponding to PDE linear operator
    ! !
    ! !   N       (input) SUBROUTINE:
    ! !           cooresponds to PDE nonlinear operator. Must be of the form:
    ! !               subroutine N(t,y_in,N_out)
    ! !                   real(kind=8),    intent(in)  :: t
    ! !                   complex(kind=8), intent(in)  :: y_in(:)
    ! !                   complex(kind=8), intent(out) :: N_out(:)
    ! !               end subroutine N
    ! !
    ! !   h   (input) DOUBLE array, dimension(2)
    ! !           contains left and right integration bounds
    ! !
    ! !   y_in    (input) COMPLEX*16 array, dimension (N)
    ! !           Initial condition.
    ! !
    ! !   settings (input) IMEXDIRK_TABLEAU
    ! !           derived data type that contains the IMEX RK tableau
    ! !
    ! !   y_out   (output) COMPLEX*16 array, dimension (N)
    ! !           Initial condition.
    ! !
    ! !   times   (output) DOUBLE array, dimension(4)
    ! !           times(1:2) = [cputime,clocktime] for timestepping loop
    ! ! =========================================================================

    subroutine initY(N, y0, t0, r, iterator, kappa, Y, YN)

        ! Arguments
        complex(dp), intent(inout) :: t0
        real(dp),    intent(in) :: r
        type(imexpbm_F_coefficients), intent(in) :: iterator
        complex(dp), intent(inout) :: y0(:)
        complex(dp), intent(inout) :: Y(:,:), YN(:,:)
        integer, intent(in) :: kappa
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
        integer :: i, q, thread_id
        q = size(iterator%z)

        ! verify valid nodes
        if(iterator%z(1) .ne. -1.0_dp) then
            write (*,*) "Invalid nodes provided to IMEX-PBDF; translate nodes so that z(1) = -1."
            call exit(1)
        endif

        ! zeroth-order guess
        Y(:,1) = y0
        call N(t0, y0, YN(:,1))

        do i = 2, q
            Y(:, i)  = Y(:,1);
            YN(:, i) = YN(:,1)
        enddo

        ! iterator

        !$ call omp_set_num_threads(q)
        !$omp parallel private(thread_id, t0)
        !$ thread_id = omp_get_thread_num() + 1

        call step(Y, YN, N, r, t0, iterator, thread_id, .False.)
        do i = 1, (q + kappa)
            call step(Y, YN, N, r, t0, iterator, thread_id)
        enddo

        !$omp end parallel

    end subroutine initY

    subroutine step(Y, YN, N, r, t0, method_coeffs, thread_id, evalN_in)
    ! =========================================================================
    ! STEP computes a single step of a pbm
    !
    !   Y^[n+1] = I \kron Y_c + h * B_im LY^[n+1] + h * B_ex N(y^[n+1])
    !
    ! Arguments
    !
    !   Y       (inout) COMPLEX*16, dimension(nL,q)
    !           Solution at each node Y(:,j) = y^{[n]}_j
    !
    !   YN      (inout) COMPLEX*16, dimension(nL,q)
    !           Storage vector for nonlinear component. It should be 
    !           pre-evaluated only if evalN_in = false
    !
    !   N       (input) SUBROUTINE:
    !           cooresponds to PDE nonlinear operator. Must be of the form:
    !               subroutine N(t,y_in,N_out)
    !                   real(kind=8),    intent(in)  :: t
    !                   complex(kind=8), intent(in)  :: y_in(:)
    !                   complex(kind=8), intent(out) :: N_out(:)
    !               end subroutine N
    !
    !   r      (input) DOUBLE
    !          node radius
    !
    !   t      (input) complex
    !          timestep center
    !
    !   method_coeffs (input) imexpbm_F_coefficients
    !          type defined in this module that contains  method coefficients
    !
    !   evalN_in (input) LOGICAL
    !          if true, YN will be populated with YN(:,j) = N(y^{[n]}_j)
    !          if false, subroutine will assume that YN already contains the
    !          correct values
    !    
    ! =========================================================================

        ! Arguments
        complex(dp), intent(inout) :: Y(:, :)
        complex(dp), intent(inout) :: YN(:, :)
        complex(dp), intent(inout) :: t0
        real(dp), intent(in) :: r
        integer :: thread_id
        type(imexpbm_F_coefficients) :: method_coeffs
        logical, optional :: evalN_in
        interface ! interface for subroutine N
            subroutine N(t, yh_in, N_out, thread_id)
            import :: dp
            complex(dp), intent(in)  :: t
            complex(dp), intent(in)  :: yh_in(:)
            complex(dp), intent(out) :: N_out(:)
            integer, intent(in), optional :: thread_id
            end subroutine N
        end interface

        ! Local variables
        integer :: q, nL, ind, nAI_ex
        logical :: evalN

        q = size(method_coeffs%z)
        nL = size(Y, 1)
        nAI_ex = size(method_coeffs%AI_ex)

        if(present(evalN_in)) then
            evalN = evalN_in
        else
            evalN = .True.
        end if

        ! ---> eval nonlinear term
        if ( (evalN .eqv. .TRUE.) .and. (thread_id <= nAI_ex) ) then
            ind = method_coeffs%AI_ex(thread_id)    
            call N(t0 + r * method_coeffs%z(ind), Y(:,ind), YN(:,ind), thread_id)
        endif
        !$omp barrier

        ! ---> form explicit right-hand-side
        Y(:, thread_id) = Y(:, method_coeffs%c)
        !$omp barrier       
        if ( thread_id <= nAI_ex ) then
            ind = method_coeffs%AI_ex(thread_id)  
            Y(:, ind) = Y(:, ind) + matmul(YN, r * method_coeffs%B_ex(:, ind))
        endif
        !$omp barrier

        ! ---> invert implicit component
        YN(:, thread_id) = sum(Y * method_coeffs%ImhLB_im(:, :, thread_id), 2) ! store result in YN
        !$omp barrier
        Y(:,thread_id) = YN(:,thread_id) ! copy back to Y
        !$omp barrier

        t0 = t0 + r * method_coeffs%alpha

    end subroutine

end module imexficpbm_F_omp_mod