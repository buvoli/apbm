! ============================ Program Description ==========================
! SISC-BDF-Experiments  Generates BDF timing results for:
! 	Buvoli T, "Additive Polynomial Time Integrators"
! ===========================================================================

program SiscExperiment

    ! ===================== Specify Equation Module Here ====================
    use kdv_mod, only: L,N,init,y0,Np,tspan,Fs,reference_methods,eqn_name,error_filter ! quasigeostrophic, nls, nikolaevskiy, kuramoto, kdv
    ! =======================================================================

    use tools_mod,   only: dp, linspace, chebpts, legpts, print_cvector, isfinite, relerror_c, &
        constvec_r, linspace, write_rmatrix, write_rvector, write_cvector, read_cvector
    use etdrk4_mod,  only: etdrk4
    use etdsdc_mod,  only: etdsdc,  etdsdc_settings
    use imexsdc_mod, only: imexsdc, imexsdc_settings
    ! Test Methods
    use imexdirk_mod, only: imexdirk, imexdirk_tableau
    use imexdirk_tableau_mod
    use imexbdf_mod, only: imexbdf, imexbdf_settings
    use imexpbdf_mod, only: imexpbdf, imexpbdf_settings
    use imexpbdf_omp_mod, only: imexpbdf_omp
    use imexficpbm_F_mod, only: imexficpbm_F, imexficpbm_F_settings
    use imexficpbm_F_methods_mod
    use imexficpbm_F_omp_mod, only: imexficpbm_F_omp

    implicit none

    ! Local Variables
    complex(dp),allocatable     :: y_reference(:), y_out(:),Lambda(:)
    integer                     :: F,relcost,Nt,solution_count,i,j, n_Fs, sdc_reference_n, num_methods
    integer                     :: nn, shift, k, sub_index
    TYPE(etdsdc_settings)       :: etdsdc_s
    TYPE(imexsdc_settings)      :: imexsdc_s
    TYPE(imexdirk_tableau)      :: imexdirk_t
    TYPE(imexbdf_settings)      :: imexbdf_s
    TYPE(imexpbdf_settings)     :: imexpbdf_s
    TYPE(imexficpbm_F_settings) :: imexficpbm_F_s
    integer, allocatable        :: IMEX_DIRK_orders(:), IMEX_BDF_orders(:), IMEX_PBDF_orders(:), & 
        IMEX_RADAU_qs(:), IMEX_RADAU_kappa(:)
    real(dp)                    :: time(4)
    real(dp),allocatable        :: times(:,:), errors(:,:), Nts(:,:), hs(:,:)
    character(len=*), parameter :: results_dir = "results/"
    LOGICAL                     :: omp_active, ref_load_flag
    
    omp_active = .FALSE.
    ! check for omp
    !$ omp_active = .TRUE.

    ! IMEX Methods to test
    IMEX_DIRK_orders = [1, 2, 3, 4]
    IMEX_BDF_orders  = [2, 3, 4, 5, 6]
    IMEX_PBDF_orders = [2, 3, 4, 5, 6, 7]
    IMEX_RADAU_qs    = [2, 3, 4, 5, 6]
    IMEX_RADAU_kappa = [0, 1, 2]

    ! Init Equation
    allocate(Lambda(Np))
    call init(8)
    call L(Lambda)

    write (*,"(a,a,a)") "Running experiment for ", eqn_name, " equation"
    
    allocate(y_out(Np), y_reference(Np))
    call read_cvector(y_reference, results_dir//eqn_name//"/reference.txt", ref_load_flag)
    if ( ref_load_flag .eqv. .TRUE. ) then 
        write (*,"(a)") "Loaded reference solution from file "//results_dir//eqn_name//"/reference.txt"
    else

        write (*,"(a,$)") "Computing reference solution... "

        ! === Compute Reference Solution ===
        F = 4*Fs(size(Fs))
        y_reference = (0.d0,0.d0)
        solution_count = 0
        sdc_reference_n = 32

        if(reference_methods(1)) then ! Use ETDSDC_N^N-1
            etdsdc_s%tau = chebpts(sdc_reference_n,(/0.d0, 1.d0 /))
            etdsdc_s%m   = sdc_reference_n-1
            relcost      = size(etdsdc_s%tau) * (etdsdc_s%m + 1)
            Nt = F/relcost + 1
            call etdsdc(Lambda,N,tspan,y0,Nt,etdsdc_s,y_out,time)
            if(isfinite(y_out)) then
                y_reference = y_reference + y_out
                solution_count = solution_count + 1
            endif
        endif

        if(reference_methods(2)) then ! Use IMEXSDC_N^N-1
            imexsdc_s%tau = chebpts(sdc_reference_n,(/0.d0, 1.d0 /))
            imexsdc_s%m   = sdc_reference_n-1
            relcost       = size(imexsdc_s%tau) * (imexsdc_s%m + 1)
            Nt = F/relcost + 1
            call imexsdc(lambda,N,tspan,y0,Nt,imexsdc_s,y_out,time)
            if(isfinite(y_out)) then
                y_reference = y_reference + y_out
                solution_count = solution_count + 1
            endif
        endif

        if(reference_methods(3)) then ! Use ETDRK4
            relcost = 4
            Nt = F/relcost + 1
            call etdrk4(Lambda,N,tspan,y0,Nt,y_out,time)
            if(isfinite(y_out)) then
                y_reference = y_reference + y_out
                solution_count = solution_count + 1
            endif
        endif

        y_reference = y_reference/real(solution_count,8)
        write (*,"(a)") "done."

        call write_cvector(y_reference, results_dir//eqn_name//"/reference.txt")

    end if

    ! === Run Numerical tests ==================================================
    write (*,"(a)") "Running Numerical Tests... "
    num_methods = size(IMEX_DIRK_orders) + size(IMEX_BDF_orders) + 2 * size(IMEX_PBDF_orders) + & 
        4 * size(IMEX_RADAU_qs) * size(IMEX_RADAU_kappa)
    n_Fs = size(Fs)
    allocate(times(n_Fs,num_methods), errors(n_Fs,num_methods), Nts(n_Fs,num_methods), hs(n_Fs,num_methods))

    ! clear arrays
    times  = 0.0_dp
    errors = 0.0_dp
    Nts    = 0.0_dp
    hs     = 0.0_dp

    do i=1,n_Fs
        F = Fs(i)
        shift = 0
        write (*,"(a,I4,a,I4,a,I10)") "     F (",i,"/",n_Fs,") = ",F
        
        ! == Run IMEX-RK Methods ==============================================
        do j=1,size(IMEX_DIRK_orders)
            nn = IMEX_DIRK_orders(j)
            select case (nn)
                case (1)
                    imexdirk_t = imexdirk1()
                    relcost = 1
                case (2)
                    imexdirk_t = imexdirk2()
                    relcost = 2
                case(3)
                    imexdirk_t = imexdirk3()
                    relcost = 3
                case(4)
                    imexdirk_t = imexdirk4()
                    relcost = 5
            end select
            Nt = F / (relcost) + 1
            call imexdirk(Lambda, N, tspan, y0, Nt, imexdirk_t, y_out, time)
            Nts(i,    j + shift) = Nt
            errors(i, j + shift) = error_filter(y_out, y_reference)
            times(i,  j + shift) = time(2)
        enddo
        shift = shift + size(IMEX_DIRK_orders)
        
        ! == Run IMEX-BDF Methods ==================================================
        do j=1,size(IMEX_BDF_orders)
            nn = IMEX_BDF_orders(j)
            imexbdf_s%order = nn
            relcost   = 1
            Nt        = F / (relcost) + 1
            call imexbdf(lambda, N, tspan, y0, Nt, imexbdf_s, y_out, time)
            Nts(i,    j + shift) = Nt
            errors(i, j + shift) = error_filter(y_out, y_reference)
            times(i,  j + shift) = time(2)
        enddo
        shift = shift + size(IMEX_BDF_orders)
        
        ! == Run IMEX-PBDF Methods (Serial) ========================================
        do j=1,size(IMEX_PBDF_orders)
            nn = IMEX_PBDF_orders(j)
            imexpbdf_s%z = linspace(0.0_dp, 2.0_dp, nn)
            imexpbdf_s%alpha = 1.0_dp / (real(nn,dp) - 1.0_dp)
            relcost   = nn
            Nt        = F / (relcost) + 1
            call imexpbdf(lambda, N, tspan, y0, Nt, imexpbdf_s, y_out, time)
            Nts(i,    j + shift) = Nt
            errors(i, j + shift) = error_filter(y_out, y_reference)
            times(i,  j + shift) = time(2)
        enddo
        shift = shift + size(IMEX_PBDF_orders)

        ! == Run IMEX-PBDF Methods (Parallel) ========================================
        if(omp_active .eqv. .TRUE.) then
            do j=1,size(IMEX_PBDF_orders)
                nn = IMEX_PBDF_orders(j)
                imexpbdf_s%z = linspace(0.0_dp, 2.0_dp, nn)
                imexpbdf_s%alpha = 1.0_dp / (real(nn,dp) - 1.0_dp)
                relcost   = 1
                Nt        = F / (relcost) + 1
                call imexpbdf_omp(lambda, N, tspan, y0, Nt, imexpbdf_s, y_out, time)
                Nts(i,    j + shift) = Nt
                errors(i, j + shift) = error_filter(y_out, y_reference)
                times(i,  j + shift) = time(2)
            enddo
        endif
        shift = shift + size(IMEX_PBDF_orders)


        ! == Run IMEX-Radau Methods (Serial) ========================================
        do j=1,size(IMEX_RADAU_qs)
            do k=1,size(IMEX_RADAU_kappa)
                imexficpbm_F_s = imexradau(IMEX_RADAU_qs(j), IMEX_RADAU_kappa(k))
                relcost   = (IMEX_RADAU_qs(j) - 1) * (IMEX_RADAU_kappa(k) + 1)
                Nt        = F / (relcost) + 1
                call imexficpbm_F(lambda, N, tspan, y0, Nt, imexficpbm_F_s, y_out, time)
                sub_index = shift + k + (j-1) * size(IMEX_RADAU_kappa)
                Nts(i, sub_index) = Nt
                errors(i, sub_index) = error_filter(y_out, y_reference)
                times(i, sub_index) = time(2)
            enddo
        enddo
        shift = shift + size(IMEX_RADAU_qs) * size(IMEX_RADAU_kappa)

        ! == Run IMEX-Radau Methods (OMP) ========================================
        if(omp_active .eqv. .TRUE.) then
            do j=1,size(IMEX_RADAU_qs)
                do k=1,size(IMEX_RADAU_kappa)
                    imexficpbm_F_s = imexradau(IMEX_RADAU_qs(j), IMEX_RADAU_kappa(k))
                    relcost   = (IMEX_RADAU_qs(j) - 1) * (IMEX_RADAU_kappa(k) + 1)
                    Nt        = F / (relcost) + 1
                    call imexficpbm_F_omp(lambda, N, tspan, y0, Nt, imexficpbm_F_s, y_out, time)
                    sub_index = shift + k + (j-1) * size(IMEX_RADAU_kappa)
                    Nts(i, sub_index) = Nt
                    errors(i, sub_index) = error_filter(y_out, y_reference)
                    times(i, sub_index) = time(2)
                enddo
            enddo
        endif
        shift = shift + size(IMEX_RADAU_qs) * size(IMEX_RADAU_kappa)

        ! == Run IMEX-Radau* Methods (Serial) ========================================
        do j=1,size(IMEX_RADAU_qs)
            do k=1,size(IMEX_RADAU_kappa)
                imexficpbm_F_s = imexradauS(IMEX_RADAU_qs(j), IMEX_RADAU_kappa(k))
                relcost   = IMEX_RADAU_qs(j) + (IMEX_RADAU_qs(j) - 1) * (IMEX_RADAU_kappa(k))
                Nt        = F / (relcost) + 1
                call imexficpbm_F(lambda, N, tspan, y0, Nt, imexficpbm_F_s, y_out, time)
                sub_index = shift + k + (j-1) * size(IMEX_RADAU_kappa)
                Nts(i, sub_index) = Nt
                errors(i, sub_index) = error_filter(y_out, y_reference)
                times(i, sub_index) = time(2)
            enddo
        enddo
        shift = shift + size(IMEX_RADAU_qs) * size(IMEX_RADAU_kappa)

        ! == Run IMEX-Radau* Methods (OMP) ========================================
        if(omp_active .eqv. .TRUE.) then
            do j=1,size(IMEX_RADAU_qs)
                do k=1,size(IMEX_RADAU_kappa)
                    imexficpbm_F_s = imexradauS(IMEX_RADAU_qs(j), IMEX_RADAU_kappa(k))
                    relcost   = IMEX_RADAU_qs(j) + (IMEX_RADAU_qs(j) - 1) * (IMEX_RADAU_kappa(k))
                    Nt        = F / (relcost) + 1
                    call imexficpbm_F_omp(lambda, N, tspan, y0, Nt, imexficpbm_F_s, y_out, time)
                    sub_index = shift + k + (j-1) * size(IMEX_RADAU_kappa)
                    Nts(i, sub_index) = Nt
                    errors(i, sub_index) = error_filter(y_out, y_reference)
                    times(i, sub_index) = time(2)
                enddo
            enddo
        endif
        shift = shift + size(IMEX_RADAU_qs) * size(IMEX_RADAU_kappa)

        ! Compute stepsizes
        hs = (tspan(2) - tspan(1)) / Nts

        ! Save Partial Results
        call write_rmatrix(Nts,         results_dir//eqn_name//"/Nts.txt")
        call write_rmatrix(hs,          results_dir//eqn_name//"/hs.txt")
        call write_rmatrix(errors,      results_dir//eqn_name//"/errors.txt")
        call write_rmatrix(times,       results_dir//eqn_name//"/times.txt")
        call write_rvector(Fs,          results_dir//eqn_name//"/Fs.txt")

    enddo
    write (*,"(a)") "done."
    
    ! Store stepsizes
    hs = (tspan(2) - tspan(1))/Nts

    ! Save Full Results
    call write_rmatrix(Nts,         results_dir//eqn_name//"/Nts.txt")
    call write_rmatrix(hs,          results_dir//eqn_name//"/hs.txt")
    call write_rmatrix(errors,      results_dir//eqn_name//"/errors.txt")
    call write_rmatrix(times,       results_dir//eqn_name//"/times.txt")
    call write_rvector(Fs,          results_dir//eqn_name//"/Fs.txt")

end program SiscExperiment