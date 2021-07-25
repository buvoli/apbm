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

module imexficpbm_F_mod

    ! Module Parameters
    use tools_mod, only: dp, tic, toc, linspace, constvec_r, zeye, quadw, interpw, Zinv
    implicit none
    
    type imexfipbm_F_settings
        real(dp), allocatable :: z(:)       ! nodes
        real(dp) :: alpha                   ! extrapolation factor
        integer,  allocatable :: AO_im(:)   ! Active output index set for propagator's implicit component (output nodes in this set will be used for polynomial)
        integer,  allocatable :: AI_ex(:)   ! Active input index set for propagator's explicit component (output nodes in this set will be used for polynomial)
        integer :: c                        ! expansion point (or left endpoint) for propagator's Adams polynomial
    end type imexfipbm_F_settings
    
    type imexficpbm_F_settings
        integer :: kappa  ! number of iterator applications
        type(imexfipbm_F_settings) :: propagator
        type(imexfipbm_F_settings) :: iterator
    end type imexficpbm_F_settings

    ! =========================================================================
    ! IMEX PBM coefficient type (Used Internally)
    ! -------------------------------------------------------------------------
    !
    !   c       integer
    !           index of input that is used as initial condition for the
    !           Adams ODE polynomial
    !
    !   ImhB_im DOUBLE array, dimension(q,q,size(L))
    !           precomputed inverses ImhB_im(:,:,j) = (I - h L(j) B_im)^{-1}
    !           where B_im is the coefficient matrix for the implicit component
    !
    !   B_ex    DOUBLE array, dimension(q,q)
    !           explicit coefficient matrix
    !
    !   AIS_im  INTEGER array
    !           indices of nodes used for implicit polynomial approximation for
    !           the nonlinear component.
    ! 
    !   AIS_ex  INTEGER array
    !           indices of nodes used for explicit polynomial approximation for
    !           the nonlinear component.
    !   nodes  
    !
    ! =========================================================================
    type imexpbm_F_coefficients
        real(dp), allocatable :: z(:)      ! nodes
        complex(dp), allocatable :: imhLB_im(:,:,:)    
        real(dp), allocatable :: B_ex(:,:)       
        integer, allocatable  :: AO_im(:)  ! active input index set for implicit component
        integer, allocatable  :: AI_ex(:)  ! active input index set for explicit component
        real(dp) :: alpha
        integer  :: c         ! expansion point for Adams ODE p
    end type imexpbm_F_coefficients

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

    subroutine imexficpbm_F(L, N, tspan, y_in, Nt, settings, y_out, times)

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
        integer :: ode_dim, q, i, j
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
        call tic()
        call initY(N, y_in, t0, r, iterator, settings%kappa, Y, YN)
        
        ! timestepping loop
        do i = 1, Nt
    
            ! propagator
            call step(Y, YN, N, r, t0, propagator)
        
            ! iterator
            do j = 1, settings%kappa
                call step(Y, YN, N, r, t0, iterator)
            enddo

        enddo

        call toc(times(1:2))

        y_out = Y(:, 1)
        
    end subroutine imexficpbm_F

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
        integer :: i, q
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
        call step(Y, YN, N, r, t0, iterator, .False.) ! No need to eval N on first step
        do i = 1, (q + kappa)
            call step(Y, YN, N, r, t0, iterator)
        enddo

    end subroutine initY

    subroutine step(Y, YN, N, r, t0, method_coeffs, evalN_in)
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
        integer :: j, q, nL, ind, nAI_ex
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
        if ( evalN .eqv. .TRUE.) then ! allow skipping for first step when computing Y^{[0]}
            do j = 1, nAI_ex
                ind = method_coeffs%AI_ex(j)
                call N(t0 + r * method_coeffs%z(ind), Y(:,ind), YN(:,ind))
            enddo
        endif

        ! store explicit right-hand-side into YN array
        YN(:, method_coeffs%AI_ex) = matmul(YN, r * method_coeffs%B_ex(:, method_coeffs%AI_ex)) ! overwrite YN (since values are no longer needed after)
        do j = 1, q
            Y(:, j) = Y(:,method_coeffs%c)
        enddo
        do j = 1, nAI_ex
            ind = method_coeffs%AI_ex(j)
            Y(:, ind) = Y(:, ind) + YN(:, ind)
        enddo
        
        ! ---> invert implicit component
        ! Method 1: Non-permuted form ( slower but simpler to understand - requires ImhLB_im from line 391 )
        ! do j = 1, Nl
        !     Y(j, :) = matmul(Y(j, :), method_coeffs%ImhLB_im(:, :, j))
        ! enddo
        ! Method 2: Permuted form
        do j = 1, q
            YN(:, j) = sum(Y * method_coeffs%ImhLB_im(:, :, j), 2)
        enddo
        Y = YN

        t0 = t0 + r * method_coeffs%alpha

    end subroutine

    function settings2Coefficient(settings, r, L) result(coefficients)

        ! Arguments
        type(imexfipbm_F_settings), intent(in) :: settings
        complex(dp), intent(in) :: L(:)
        real(dp), intent(in) :: r
        type(imexpbm_F_coefficients) :: coefficients

        ! Local Variables
        real(dp), allocatable :: B_im(:,:), B_ex(:,:)
        complex(dp), allocatable :: ImhLB_im(:,:,:)
        integer  :: i, q, nL
        
        q  = size(settings%z)
        nL = size(L)

        ! == Propagator =======================================================
        
        ! ----> explicit component
        allocate(B_ex(q,q))
        B_ex = 0.0_dp
        B_ex(settings%AI_ex, :) = transpose( &
            quadW( & 
                settings%z(settings%AI_ex), &
                constvec_r(q, settings%z(settings%c)), &
                settings%z + settings%alpha &
            ) &
        )
        
        ! ----> implicit component
        allocate(B_im(q,q))
        B_im = 0.0_dp
        B_im(settings%AO_im, :) = transpose( &
            quadW( & 
                settings%z(settings%AO_im), &
                constvec_r(q, settings%z(settings%c) - settings%alpha), &
                settings%z &
            ) &
        )
        
        ! allocate(ImhLB_im(q, q, nL))
        ! do i = 1, nL
        !     ImhLB_im(:, :, i) = Zinv( zeye(q) - r * L(i) * B_im )
        ! enddo

        allocate(ImhLB_im(NL, q, q))
        do i = 1, nL
            ImhLB_im(i, :, :) = Zinv( zeye(q) - r * L(i) * B_im)
        enddo

        ! == Add values to coefficient type
        coefficients%z = settings%z
        coefficients%alpha = settings%alpha
        coefficients%ImhLB_im = ImhLB_im
        coefficients%B_ex = B_ex
        coefficients%AO_im = settings%AO_im
        coefficients%AI_ex = settings%AI_ex
        coefficients%c = settings%c

    end function

end module imexficpbm_F_mod

module imexficpbm_F_methods_mod

    ! Module Parameters
    use tools_mod, only: dp, tic, toc, weights, radaupts
    use imexficpbm_F_mod, only: imexficpbm_F_settings, imexfipbm_F_settings

    implicit none

    contains
    
    function imexradau(q, kappa) result(s)

        ! Arguments
        integer, intent(in) :: q
        integer, intent(in) :: kappa
        type(imexficpbm_F_settings) :: s

        ! Local Variables
        real(dp), allocatable :: z(:)
        integer, allocatable :: IS(:)
        integer :: i
        type(imexfipbm_F_settings) :: propagator, iterator
        
        ! nodes
        allocate(z(q))
        z(1) = -1.0_dp
        z(q:2:-1) = -1.0_dp * radaupts(q - 1)

        ! Active Index Sets
        allocate(IS(q - 1))
        do i = 1, q - 1
            IS(i) = i + 1
        end do
        
        ! propagator
        propagator%z     = z       ! nodes
        propagator%alpha = 2.0_dp  ! extrapolation factor
        propagator%c     = q
        propagator%AI_ex = IS
        propagator%AO_im = IS
        
        ! iterator
        iterator%z     = z       ! nodes
        iterator%alpha = 0.0_dp  ! extrapolation factor
        iterator%c     = 1
        iterator%AI_ex = IS
        iterator%AO_im = IS

        s%propagator = propagator   !
        s%iterator = iterator       ! 
        s%kappa = kappa             ! number of iterator applications

    end function

    function imexradauS(q, kappa) result(s)

        integer, intent(in) :: q
        integer, intent(in) :: kappa
        integer :: i
        type(imexficpbm_F_settings) :: s

        integer, allocatable :: IS(:)
        
        ! Active Index Sets
        allocate(IS(q))
        do i = 1, q
            IS(i) = i
        end do

        s = imexradau(q, kappa)
        s%propagator%AI_ex = IS

    end function

end module imexficpbm_F_methods_mod