! ============================ Module Description =============================
! imexrk_mod Implementation of IMEX-RK method. Contains the subroutine:
!   1. imexrk  - IMEXRK timestepping code
! =============================================================================

module imexdirk_mod

    ! Module Parameters
    use tools_mod, only: dp, tic, toc, weights
    implicit none
    type imexdirk_tableau
        ! Implicit Method
        real(dp), allocatable :: A_im(:,:) ! stage coefficients
        real(dp), allocatable :: b_im(:)   ! output coefficients
        real(dp), allocatable :: c_im(:)   ! stage nodes
        ! Explicit Method
        real(dp), allocatable :: A_ex(:,:) ! stage coefficients
        real(dp), allocatable :: b_ex(:)   ! output coefficients
        real(dp), allocatable :: c_ex(:)   ! stage nodes
        ! Output Flag
        LOGICAL :: output_is_final_stage
    end type imexdirk_tableau

    contains

    ! =========================================================================
    ! IMEXDIRK   Implicit-Explicit Diagonally Implicit RK integrator
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
    !   tableau (input) IMEXDIRK_TABLEAU
    !           derived data type that contains the IMEX RK tableau
    !
    !   y_out   (output) COMPLEX*16 array, dimension (N)
    !           Initial condition.
    !
    !   times   (output) DOUBLE array, dimension(4)
    !           times(1:2) = [cputime,clocktime] for timestepping loop
    ! =========================================================================

    subroutine imexdirk(L, N, tspan, y_in, Nt, tableau, y_out, times)

        ! Arguments
        complex(dp), intent(in)         :: L(:)
        real(dp),    intent(in)         :: tspan(2)
        complex(dp), intent(inout)      :: y_in(:)
        type(imexdirk_tableau), intent(in) :: tableau
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
        complex(dp), allocatable :: Y(:,:), YL(:,:), YN(:,:), y0(:), t0
        real(dp) :: h
        integer :: ode_dim, nstages, i, j, k
        logical, allocatable :: evalL_flags(:), evalN_flags(:)
        
        nstages = size(tableau%A_im, 1)
        ode_dim = size(L)
        allocate(Y(ode_dim, nstages), YL(ode_dim, nstages) , YN(ode_dim, nstages), y0(ode_dim))
        YL = 0 ! set values to zero to avoid any Nans
        YN = 0 ! set values to zero to avoid any Nans
        
        y0 = y_in
        t0 = tspan(1)
        h  = (tspan(2)-tspan(1))/(real(Nt,8))

        evalL_flags = evalFlags(tableau%A_im, tableau%b_im, tableau%output_is_final_stage)
        evalN_flags = evalFlags(tableau%A_ex, tableau%b_ex, tableau%output_is_final_stage)

        call tic();
        do k = 1, Nt
            ! -- compute stage values
            do i = 1, nstages
                Y(:, i) = y0;
                do j = 1, (i - 1)
                    Y(:, i) = Y(:, i) + (h * tableau%A_im(i,j)) * (YL(:, j)) + (h * tableau%A_ex(i,j)) * YN(:, j)
                enddo
                Y(:, i) = Y(:, i) / ( 1 - (h * tableau%A_im(i,i)) * L)
                ! -- evaluate stage derivatives
                if(evalL_flags(i)) then
                    YL(:,i) = L * Y(:, i)
                endif
                if(evalN_flags(i)) then
                    call N(t0 + h * tableau%c_ex(i), Y(:,i), YN(:,i))
                endif
            enddo
            ! -- compute outputs
            if ( tableau%output_is_final_stage ) then
                y0 = Y(:, nstages)
            else            
                do i = 1, nstages
                    y0 = y0 + (h * tableau%b_im(i)) * YL(:, i) + (h * tableau%b_ex(i)) * YN(:, i)
                enddo
            endif
            t0 = t0 + h   
        enddo
        call toc(times(1:2))

        y_out = y0
        
    end subroutine imexdirk

    function evalFlags(A, b, output_is_final_stage) result(flags)
    ! Analyzes sparsity of RK tablue and returns array flags where:
    !   flags(i) = true \implies it is necessary to evaluate ith stage derivative
    !   flags(i) = false \implies it is not necessary to evaluate ith stage derivative
        
        ! Arguments
        real(dp), intent(in) :: A(:,:), b(:)
        logical, intent(in)  :: output_is_final_stage
        logical, allocatable :: flags(:)

        ! Local variables
        integer :: nstages, i
        
        nstages = size(A, 1)
        allocate(flags(nstages))
        flags = .TRUE.
        
        if(output_is_final_stage) then ! do not consider b if output_is_final_stage = true
            do i = 1, nstages - 1
                if ( maxval(abs(A(i+1:nstages,i))) == 0 ) then
                    flags(i) = .False.
                end if
                flags(nstages) = .False.
            end do
        else ! consider b if output_is_final_stage = false
            do i = 1, nstages - 1
                if ( max( abs(b(i)), maxval(abs(A(i+1:nstages,i)))) == 0) then
                    flags(i) = .False.
                end if
            end do
            if( abs(b(nstages)) == 0 ) then
                flags(nstages) = .False.
            endif
        endif

    end function evalFlags

end module imexdirk_mod

! ============================ Module Description =============================
! imexrk_mod Implementation of IMEX-RK method. Contains the subroutine:
!   1. imexrk  - IMEXRK timestepping code
! =============================================================================

module imexdirk_tableau_mod

    use tools_mod, only: dp, tic, toc, weights    
    use imexdirk_mod, only: imexdirk_tableau

    contains

    function imexdirk1() result(T)
    ! Returns IMEX-RK Tablue of forward-backward Euler

        type(imexdirk_tableau) :: T

        ! implicit component
        T%A_im = reshape( [ 0, 0,   &
                            0, 1 ], &
                        [2, 2], order=[2,1] )
        T%b_im = [0, 1]
        T%c_im = [0, 1]
    
        ! explicit component
        T%A_ex = reshape( [ 0, 0,   &
                            1, 0 ], &
                        [2, 2], order=[2,1] )
        T%b_ex = [1, 0]
        T%c_ex = [0, 1]

        T%output_is_final_stage = .TRUE.
 
    end function imexdirk1

    function imexdirk2() result(T)
    !IMEXDIRK2 Second-Order IMEX RK method (2,3,2) from:
    ! U. M. Ascher, S. J. Ruuth, and R. J. Spiteri, Implicit-Explicit 
    ! Runge-Kutta methods for time-dependent partial differential equations, 
    ! Appl. Numer. Math., 25 (1997), pp. 151-167.   

        real(dp) :: gamma  = (2.0_dp - sqrt(2.0_dp)) / 2;
        real(dp) :: delta  = -2.0_dp * sqrt(2.0_dp) / 3.0_dp;
        real(dp) :: A_im(3,3), b_im(3), c_im(3), A_ex(3,3), b_ex(3), c_ex(3)
        type(imexdirk_tableau) :: T

        ! -- implicit coefficients
        A_im = 0.0_dp
        A_im(2,2) = gamma
        A_im(3,2) = 1 - gamma
        A_im(3,3) = gamma
        b_im = [0.0_dp, 1.0_dp - gamma, gamma]
        c_im = [0.0_dp, gamma, 1.0_dp]

        ! -- explicit coefficients
        A_ex = 0.0_dp
        A_ex(2,1) = gamma
        A_ex(3,1) = delta
        A_ex(3,2) = 1 - delta
        b_ex = [0.0_dp, 1.0_dp - gamma, gamma]
        c_ex = [0.0_dp, gamma, 1.0_dp]

        ! Save values to tablaeu
        T%A_im = A_im
        T%b_im = b_im
        T%c_im = c_im

        T%A_ex = A_ex
        T%b_ex = b_ex
        T%c_ex = c_ex

        T%output_is_final_stage = .FALSE.

    end function imexdirk2

    function imexdirk3() result(T)
    !IMEXDIRK3 Third-order ARK3(2)4L[2]SA
    !   C. A. Kennedy and M. H. Carpenter, Additive Runge-Kutta schemes for 
    !   convection-diffusion-reaction equations, Appl. Numer. Math., 44
    !   (2003), pp. 139-181.

    real(dp) :: A_im(4,4), b_im(4), c_im(4), A_ex(4,4), b_ex(4), c_ex(4)
    type(imexdirk_tableau) :: T

    ! -- implicit coefficients -----------------------------------------------------------
    A_im = 0.0_dp;
    A_im(2,1) =   1767732205903.0_dp  / 4055673282236.0_dp
    A_im(2,2) =   1767732205903.0_dp  / 4055673282236.0_dp
    A_im(3,1) =   2746238789719.0_dp  / 10658868560708.0_dp
    A_im(3,2) = - 640167445237.0_dp   / 6845629431997.0_dp
    A_im(3,3) =   1767732205903.0_dp  / 4055673282236.0_dp
    A_im(4,1) =   1471266399579.0_dp  / 7840856788654.0_dp
    A_im(4,2) = - 4482444167858.0_dp  / 7529755066697.0_dp
    A_im(4,3) =   11266239266428.0_dp / 11593286722821.0_dp
    A_im(4,4) =   1767732205903.0_dp  / 4055673282236.0_dp

    b_im = [  1471266399579.0_dp  / 7840856788654.0_dp    , &
            - 4482444167858.0_dp / 7529755066697.0_dp     , &
              11266239266428.0_dp / 11593286722821.0_dp   , &
              1767732205903.0_dp / 4055673282236.0_dp       &   
    ]

    c_im = [ 0.0_dp, 1767732205903.0_dp / 2027836641118.0_dp, 3.0_dp / 5.0_dp, 1.0_dp ];

    ! -- explicit coefficients -----------------------------------------------------------
    A_ex = 0.0_dp;
    A_ex(2,1) =   1767732205903.0_dp  / 2027836641118.0_dp
    A_ex(3,1) =   5535828885825.0_dp  / 10492691773637.0_dp
    A_ex(3,2) =   788022342437.0_dp   / 10882634858940.0_dp
    A_ex(4,1) =   6485989280629.0_dp  / 16251701735622.0_dp
    A_ex(4,2) = - 4246266847089.0_dp  / 9704473918619.0_dp
    A_ex(4,3) =   10755448449292.0_dp / 10357097424841.0_dp

    b_ex = b_im
    c_ex = c_im

    ! Save values to tablaeu
    T%A_im = A_im
    T%b_im = b_im
    T%c_im = c_im

    T%A_ex = A_ex
    T%b_ex = b_ex
    T%c_ex = c_ex

    T%output_is_final_stage = .FALSE.

    end function imexdirk3

    function imexdirk4() result(T)
    ! IMEXDIRK4 Fourth-order ARK4(3)6L[2]SA from
    !   C. A. Kennedy and M. H. Carpenter, Additive Runge-Kutta schemes for 
    !   convection-diffusion-reaction equations, Appl. Numer. Math., 44
    !   (2003), pp. 139-181.

    real(dp) :: A_im(6,6), b_im(6), c_im(6), A_ex(6,6), b_ex(6), c_ex(6), gamma
    type(imexdirk_tableau) :: T

    ! -- implicit coefficients
    gamma = 0.25_dp
    A_im  = 0.0_dp

    A_im(2,1) =   gamma
    A_im(2,2) =   gamma
    A_im(3,1) =   8611.0_dp        / 62500.0_dp
    A_im(3,2) = - 1743.0_dp        / 31250.0_dp
    A_im(3,3) =   gamma
    A_im(4,1) =   5012029.0_dp     / 34652500.0_dp
    A_im(4,2) = - 654441.0_dp      / 2922500.0_dp
    A_im(4,3) =   174375.0_dp      / 388108.0_dp
    A_im(4,4) =   gamma
    A_im(5,1) =   15267082809.0_dp / 155376265600.0_dp
    A_im(5,2) = - 71443401.0_dp    / 120774400.0_dp
    A_im(5,3) =   730878875.0_dp   / 902184768.0_dp
    A_im(5,4) =   2285395.0_dp     / 8070912.0_dp
    A_im(5,5) =   gamma     
    A_im(6,1) =   82889.0_dp       / 524892.0_dp
    A_im(6,2) =   0.0_dp
    A_im(6,3) =   15625.0_dp       / 83664.0_dp
    A_im(6,4) =   69875.0_dp       / 102672.0_dp
    A_im(6,5) = - 2260.0_dp        / 8211.0_dp
    A_im(6,6) =   gamma
      
    b_im = [  82889.0_dp / 524892.0_dp, &
              0.0_dp,                   &
              15625.0_dp /  83664.0_dp, &
              69875.0_dp / 102672.0_dp, &
            - 2260.0_dp / 8211.0_dp,    &
              0.25_dp                   &
    ];

    c_im = [ 0.0_dp, 0.5_dp, 83.0_dp / 250.0_dp, 31.0_dp / 50.0_dp, 17.0_dp / 20.0_dp, 1.0_dp ];

    ! -- explicit coefficients
    A_ex = 0.0_dp
    A_ex(2,1) =   0.5_dp;
    A_ex(3,1) =   13861.0_dp          / 62500.0_dp
    A_ex(3,2) =   6889.0_dp           / 62500.0_dp
    A_ex(4,1) = - 116923316275.0_dp   / 2393684061468.0_dp
    A_ex(4,2) = - 2731218467317.0_dp  / 15368042101831.0_dp
    A_ex(4,3) =   9408046702089.0_dp  / 11113171139209.0_dp
    A_ex(5,1) = - 451086348788.0_dp   / 2902428689909.0_dp
    A_ex(5,2) = - 2682348792572.0_dp  / 7519795681897.0_dp
    A_ex(5,3) =   12662868775082.0_dp / 11960479115383.0_dp
    A_ex(5,4) =   3355817975965.0_dp  / 11060851509271.0_dp
    A_ex(6,1) =   647845179188.0_dp   / 3216320057751.0_dp
    A_ex(6,2) =   73281519250.0_dp    / 8382639484533.0_dp
    A_ex(6,3) =   552539513391.0_dp   / 3454668386233.0_dp
    A_ex(6,4) =   3354512671639.0_dp  / 8306763924573.0_dp
    A_ex(6,5) =   4040.0_dp           / 17871.0_dp

    b_ex = b_im;
    c_ex = c_im;

    ! Save values to tablaeu
    T%A_im = A_im
    T%b_im = b_im
    T%c_im = c_im

    T%A_ex = A_ex
    T%b_ex = b_ex
    T%c_ex = c_ex

    T%output_is_final_stage = .FALSE.

    end function imexdirk4
    
    function imexdirk4I() result(T)
    ! IMRK4 Improved Fourth-order ARK4(3)7L[2]SA from
    ! Christopher A., and Mark H. Carpenter. "Higher-order additive Runge–Kutta
    ! schemes for ordinary differential equations." Applied Numerical 
    ! Mathematics 136 (2019): 183-205.

    real(dp) :: A_im(7,7), b_im(7), c_im(7), A_ex(7,7), b_ex(7), c_ex(7), gamma
    type(imexdirk_tableau) :: T

    ! -- implicit coefficients
    gamma = 1235.0_dp / 10000.0_dp;
    A_im  = 0.0_dp

    A_im(2,1) =   gamma;
    A_im(2,2) =   gamma;
    A_im(3,1) =   624185399699.0_dp / 4186980696204.0_dp
    A_im(3,2) =   624185399699.0_dp / 4186980696204.0_dp
    A_im(3,3) =   gamma;
    A_im(4,1) =   1258591069120.0_dp / 10082082980243.0_dp
    A_im(4,2) =   1258591069120.0_dp / 10082082980243.0_dp
    A_im(4,3) = - 322722984531.0_dp  / 8455138723562.0_dp
    A_im(4,4) =   gamma;
    A_im(5,1) = - 436103496990.0_dp  / 5971407786587.0_dp
    A_im(5,2) = - 436103496990.0_dp  / 5971407786587.0_dp
    A_im(5,3) = - 2689175662187.0_dp / 11046760208243.0_dp
    A_im(5,4) =   4431412449334.0_dp / 12995360898505.0_dp
    A_im(5,5) =   gamma;     
    A_im(6,1) = - 2207373168298.0_dp / 14430576638973.0_dp
    A_im(6,2) = - 2207373168298.0_dp / 14430576638973.0_dp
    A_im(6,3) =   242511121179.0_dp  / 3358618340039.0_dp
    A_im(6,4) =   3145666661981.0_dp / 7780404714551.0_dp
    A_im(6,5) =   5882073923981.0_dp / 14490790706663.0_dp
    A_im(6,6) =   gamma;
    A_im(7,1) =   0.0_dp;
    A_im(7,2) =   0.0_dp;
    A_im(7,3) =   9164257142617.0_dp  / 17756377923965.0_dp
    A_im(7,4) = - 10812980402763.0_dp / 74029279521829.0_dp
    A_im(7,5) =   1335994250573.0_dp  / 5691609445217.0_dp
    A_im(7,6) =   2273837961795.0_dp  / 8368240463276.0_dp
    A_im(7,7) =   gamma;

    b_im = [  0.0_dp, &
              0.0_dp, &
              9164257142617.0_dp / 17756377923965.0_dp, &
            - 10812980402763.0_dp / 74029279521829.0_dp, &
              1335994250573.0_dp / 5691609445217.0_dp, &
              2273837961795.0_dp / 8368240463276.0_dp, &
              gamma &
    ];

    c_im = [ 0.0_dp, &
             247.0_dp / 1000.0_dp, &
             4276536705230.0_dp / 10142255878289.0_dp, &
             67.0_dp / 200.0_dp, &
             3.0_dp / 40.0_dp, &
             7.0_dp / 10.0_dp, &
             1.0_dp &
    ];

    ! -- explicit coefficients
    A_ex = 0.0_dp
    A_ex(2,1) =   247.0_dp / 1000.0_dp;
    A_ex(3,1) =   247.0_dp / 4000.0_dp;
    A_ex(3,2) =   2694949928731.0_dp / 7487940209513.0_dp
    A_ex(4,1) =   464650059369.0_dp  / 8764239774964.0_dp
    A_ex(4,2) =   878889893998.0_dp  / 2444806327765.0_dp
    A_ex(4,3) = - 952945855348.0_dp  / 12294611323341.0_dp
    A_ex(5,1) =   476636172619.0_dp  / 8159180917465.0_dp
    A_ex(5,2) = - 1271469283451.0_dp / 7793814740893.0_dp
    A_ex(5,3) = - 859560642026.0_dp  / 4356155882851.0_dp
    A_ex(5,4) =   1723805262919.0_dp / 4571918432560.0_dp
    A_ex(6,1) =   6338158500785.0_dp / 11769362343261.0_dp
    A_ex(6,2) = - 4970555480458.0_dp / 10924838743837.0_dp
    A_ex(6,3) =   3326578051521.0_dp / 2647936831840.0_dp
    A_ex(6,4) = - 880713585975.0_dp  / 1841400956686.0_dp
    A_ex(6,5) = - 1428733748635.0_dp / 8843423958496.0_dp
    A_ex(7,1) =   760814592956.0_dp  / 3276306540349.0_dp
    A_ex(7,2) =   760814592956.0_dp / 3276306540349.0_dp
    A_ex(7,3) = - 47223648122716.0_dp / 6934462133451.0_dp
    A_ex(7,4) =   71187472546993.0_dp / 9669769126921.0_dp
    A_ex(7,5) = - 13330509492149.0_dp / 9695768672337.0_dp
    A_ex(7,6) =   11565764226357.0_dp / 8513123442827.0_dp

    b_ex = b_im;
    c_ex = c_im;

    ! Save values to tablaeu
    T%A_im = A_im
    T%b_im = b_im
    T%c_im = c_im

    T%A_ex = A_ex
    T%b_ex = b_ex
    T%c_ex = c_ex

    T%output_is_final_stage = .FALSE.

    end function imexdirk4I

end module imexdirk_tableau_mod