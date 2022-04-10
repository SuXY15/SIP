module mod_para
implicit none

save
    !* simulation parameter
    integer,parameter::kdp = 8  ! kind double precision
    integer,parameter::eqn = 4  ! eqn-   num of equations
    integer,parameter::ex  = 6  ! ex-    extra num in x-dir
    
    character(40):: nml_path="src/input.nml"
    character(40):: data_dir="data/"

    !* physical constants
    real(kdp),parameter::h_   = 6.6260693e-34    ! Plank constant
    real(kdp),parameter::hbar = 1.05457266e-34   ! reduced Plank constant
    real(kdp),parameter::c_   = 2.99792458e8     ! velocity of light in vacuum
    real(kdp),parameter::k_   = 1.3806488e-23    ! Boltzmann constant
    real(kdp),parameter::N_   = 500*6.0221413e23 ! Avogadro constant
    
    !* numerical constants
    real(kdp),parameter::ep   = 1e-20            ! ep- epsilon for numerical non-zero division
    real(kdp),parameter::ev   = 1e-20            ! ev- epsilon for non-zero velocity
    real(kdp),parameter::n1_2 = 1./2.
    real(kdp),parameter::n1_3 = 1./3.
    real(kdp),parameter::n2_3 = 2./3.
    real(kdp),parameter::n4_3 = 4./3.
    real(kdp),parameter::n3_2 = 3./2.
    real(kdp),parameter::n3_4 = 3./4.
    real(kdp),parameter::n1_6 = 1./6.
    real(kdp),parameter::pi   = 4.0*ATAN(1.0)
    real(kdp),parameter::pi2  = pi**2.
    
    !* control boolean
    logical::ifsave = .true.           ! ifsave-   save data file or not
    logical::ifdiff = .true.           ! ifdiff-   add diffusion term or not
    logical::ifreac = .true.           ! ifreac-   add reaction term or not
    logical::ifcurv = .false.          ! ifcurv-   add curvature or not
    logical::ifchec = .false.          ! ifchec-   check eigenvalue or not
    logical::ifstop = .false.          ! ifstop-   check if stop

    !* main parameters
    real(kdp)::dx    = 0.02            ! dx-       space step size
    real(kdp)::tstp  = 1.              ! tstp-     timestep, when to save
    real(kdp)::tend  = 1000            ! tend-     time end, when to stop
    real(kdp)::cfl   = 0.25            ! cfl-      Courant-Friedichs-Lewy condition

#ifdef supernova
    real(kdp)::Le       = 1e0          ! Le-       Lewis number
    real(kdp)::gamma    = 4./3.        ! gamma-    heat capacity ratio 
    real(kdp)::lambda   = 0.14541      ! lambda-   heat conductivity           (dimensionless)
    real(kdp)::alpha    = 2./3.        ! alpha-    diffusivity D = alpha / Le  (dimensionless)
    real(kdp)::Q        = 0.904        ! Q-        heat value                  (dimensionless)
    real(kdp)::eos_C1   = 2.60559e-1   ! eos_C1-   C1 of equation of state
    real(kdp)::eos_C2   = 6.82973e-2   ! eos_C2-   C2 of equation of state
    real(kdp)::reac_A   = 4.251e24     ! reac_A-   pre-exponential factor      (dimensionless)
    real(kdp)::reac_Ea  = 57.7224      ! reac_Ea-  activation energy           (dimensionless)
    
    !* case parameters
    real(kdp)::Din = 1.0d0             ! Din-      initial density
    real(kdp)::Vin = 0.0d0             ! Vin-      initial velocity
    real(kdp)::Tin = 0.1d0             ! Tin-      initial temprature
    real(kdp)::Yin = 1.0d0             ! Yin-      initial concentration
    real(kdp)::Tig = 1.0d0             ! Tig-      igniting temprature
    real(kdp)::Yig = 1.0d0             ! Yig-      igniting concentration
    real(kdp)::reac_acc = 1.0d0        ! reac_acc- reaction rate accelerate ratio
#else
    real(kdp)::Le       = 1.2          ! Le-       Lewis number
    real(kdp)::gamma    = 1.3          ! gamma-    heat capacity factor
    real(kdp)::lambda   = 3.9e-4       ! lambda-   heat conductivity           (dimensionless)
    real(kdp)::alpha    = 6.6e-3       ! alpha-    diffusivity D = alpha / Le  (dimensionless)
    real(kdp)::Q        = 36           ! Q-        heat value                  (dimensionless)
    real(kdp)::reac_A   = 21.9         ! reac_A-   pre-exponential factor      (dimensionless)
    real(kdp)::reac_Ea  = 27           ! reac_Ea-  activation energy           (dimensionless)
    
    !* case parameters
    real(kdp)::Din = 1.0d0             ! Din-      initial density
    real(kdp)::Vin = 0.                ! Vin-      initial velocity
    real(kdp)::Tin = 1.0d0             ! Tin-      initial temprature
    real(kdp)::Yin = 1.0d0             ! Yin-      initial concentration
    real(kdp)::Tig = 8.0d0             ! Tig-      igniting temprature
    real(kdp)::Yig = 0.3d0             ! Yig-      igniting concentration
    real(kdp)::reac_acc = 1.0          ! reac_acc- reaction rate accelerate ratio
#endif

    real(kdp)::mdt                     ! mdt-      minimal dt
    
    real(kdp)::bl                      ! bl-       inner bound length
    real(kdp)::nl                      ! nl-       outer bound length
    real(kdp)::tl                      ! tl-       perturb range length
    real(kdp)::pl                      ! pl-       plane range length

    integer::totN = 1000000            ! totN-     total number of timesteps
    integer::bxm                       ! bxm-      x_min_range
    integer::nxm                       ! nxm-      x_max_range

    integer::bx                        ! bx-       start points num in r-dir
    integer::nx                        ! nx-       points num in r-dir
    integer::tx                        ! tx-       perturb num in x-dir
    integer::px                        ! px-       plane num in x-dir
    
    integer::slx
    integer::srx
    integer::tlx
    integer::trx

    integer::direct = 1                ! direct-   propagating direction, =-1 for inward propagation
                                       !                                  = 1 for outward propagation
    integer::bc_l                      ! bc_l-     left boundary type
    integer::bc_r                      ! bc_r-     right boundary type

    logical::final_save = .true.       ! ifstop-   check if stop

    !!** =========================================================
    !* control parameters:
    integer::i,j,k,m                   ! i,j,k,m-  for loop
    integer::nt                        ! nt-       num of time-steps
    integer::count                     ! count-    count for turns
    real(kdp)::dt                      ! dt-       time step (changeable)
    real(kdp)::em                      ! em-       for calculating time-step
    real(kdp)::time                    ! time-     in simulation
    real(kdp)::temp                    ! temp-     temp val

    ! parameter of equations
    real(kdp)::gamma1
    
    !* necessary vars for solutions:
    !ALLOCATABLE
    real(kdp),allocatable::x(:)        ! x-        x val in range
    real(kdp),allocatable::rho(:)      ! rho-      density
    real(kdp),allocatable::vex(:)      ! vex-      velocity
    real(kdp),allocatable::pre(:)      ! pre-      pressure
    
    real(kdp),allocatable::E(:)        ! E-        total energy
    real(kdp),allocatable::H(:)        ! H-        enthalpy
    real(kdp),allocatable::Y(:)        ! Y-        reactant mass fraction
    real(kdp),allocatable::T(:)        ! T-        temprature
    real(kdp),allocatable::omega(:)

    real(kdp),allocatable::U(:,:)      ! U-        U in equations, conserved variables
    real(kdp),allocatable::F(:,:)      ! F-        F in equations, convectiono flux
    real(kdp),allocatable::U0(:,:)     ! U0-       for storing U^(n)
    real(kdp),allocatable::LU(:,:)     ! LU-       L(U):differential in space
    real(kdp),allocatable::fffd(:,:)   ! fffd-     

    real(kdp),allocatable::evl(:,:,:)  ! evl-      left eigen matrix (LF)
    real(kdp),allocatable::evr(:,:,:)  ! evr-      right eigen matrix (RF)

    real(kdp)::flamePx                 ! flamePx-  x pos of flame
    real(kdp)::flameLP                 ! flameLP-  last pos of flame
    real(kdp)::flameLT                 ! flameLT-  flame pos time 
    real(kdp)::flameVe                 ! flameVe-  flame velocity

    NAMELIST/PARA/  ifsave, ifdiff, ifreac, ifcurv, ifchec, ifstop,      &
                    direct, bc_l,   bc_r,   bl,     nl,     tl,     pl,  &
                    dx,     tstp,   tend,   cfl,    Le,     reac_acc,    &
                    Tin,    Tig,    Yin,    Yig,    Din,    Vin,         &
                    data_dir

contains

#ifdef supernova

    ! EOS: calculate pre from rho and T
    real(kdp) function EOS_pre_forward(rho_i, T_i)
        real(kdp), intent(in)::rho_i, T_i
        EOS_pre_forward = eos_C1*rho_i**(4./3.) + eos_C2*rho_i**(2./3.)*T_i**2
    end function EOS_pre_forward

    ! EOS: calculate pre from rho, vex, E and Y
    real(kdp) function EOS_pre_backward(rho_i, vex_i,  E_i, Y_i)
        real(kdp),intent(in)::rho_i, vex_i, E_i, Y_i
        EOS_pre_backward = rho_i*gamma1*(E_i - 0.5*vex_i**2 - Q*Y_i)
    end function EOS_pre_backward

    ! EOS: calculate rho from T, pre
    real(kdp) function EOS_rho(T_i, pre_i)
        real(kdp),intent(in)::T_i, pre_i
        real(kdp)::tmp, rho23 ! rho**(2./3.)
        rho23 = (-eos_C2 * T_i**2 + sqrt(eos_C2**2*T_i**4 + 4.*pre_i*eos_C1))/2/eos_C1
        EOS_rho = rho23**(3./2.)
    end function EOS_rho
    
    ! EOS: calculate T from rho, pre
    real(kdp) function EOS_T(rho_i, pre_i)
        real(kdp),intent(in)::rho_i, pre_i
        real(kdp)::tmp
        tmp = (pre_i-eos_C1*rho_i**(4./3.))
        if (tmp < 0.) then
            write(*,*) "NaN in EOS_T, T^2 =", tmp / eos_C2 / rho_i**(2./3.)
            ifstop = .true.
            tmp = abs(tmp)
        endif
        EOS_T = sqrt( tmp / eos_C2 / rho_i**(2./3.) )
    end function EOS_T
    
    ! CKM: calculate
    real(kdp) function reaction_rate(Y_i, T_i)
        real(kdp)::Y_i, T_i
        reaction_rate = reac_acc * reac_A * Y_i * exp(-reac_Ea/(T_i)**(1./3.))
    end function reaction_rate
    
#else
    ! EOS: E = e + u**2 / 2
    !      e = T / (gamma-1) + Q Y
    !      p = rho (gamma-1) (e - Q Y)
    !      h = E + p/rho
    
    ! EOS: calculate pre from rho, T and Y
    real(kdp) function EOS_pre_forward(rho_i, T_i, Y_i)
        real(kdp), intent(in)::rho_i, T_i, Y_i
        ! EOS_pre_forward = rho_i * T_i + rho_i * gamma1 * Q * Y_i
        EOS_pre_forward = rho_i * T_i
    end function EOS_pre_forward

    ! EOS: calculate pre from rho, vex and E
    real(kdp) function EOS_pre_backward(rho_i, vex_i,  E_i, Y_i)
        real(kdp),intent(in)::rho_i, vex_i, E_i, Y_i
        ! EOS_pre_backward = rho_i*gamma1*(E_i - 0.5*vex_i**2)
        EOS_pre_backward = rho_i*gamma1*(E_i - 0.5*vex_i**2 - Q * Y_i)
    end function EOS_pre_backward
    
    ! EOS: calculate E from rho, pre and vex
    real(kdp) function EOS_E(rho_i, pre_i, vex_i, Y_i)
        real(kdp),intent(in)::rho_i, pre_i, vex_i, Y_i
        ! EOS_E = pre_i / rho_i  / gamma1 + vex_i * vex_i / 2.
        EOS_E = pre_i / rho_i  / gamma1 + vex_i * vex_i / 2. + Q * Y_i
    end function EOS_E
    
    ! EOS: calculate rho from pre, T and Y
    real(kdp) function EOS_rho(pre_i, T_i, Y_i)
        real(kdp),intent(in)::pre_i, T_i, Y_i
        ! EOS_rho = pre_i / (T_i / gamma1 + Q * Y_i) / gamma1
        EOS_rho = pre_i / T_i
    end function EOS_rho
    
    ! EOS: calculate T from rho, pre and Y_i
    real(kdp) function EOS_T(rho_i, pre_i, Y_i)
        real(kdp),intent(in)::rho_i, pre_i, Y_i
        EOS_T = pre_i / rho_i
    end function EOS_T
    
    ! CKM: calculate
    real(kdp) function reaction_rate(Y_i, T_i)
        real(kdp)::Y_i, T_i
        reaction_rate = reac_acc * reac_A * Y_i * exp(-reac_Ea/T_i)
    end function reaction_rate

#endif

    !* =============================================
    !* push data from field variables to U
    subroutine push_data(l,r)
    integer::l,r
        do i = l, r
            U(i,1) = rho(i)            ! rho
            U(i,2) = rho(i)*vex(i)     ! rho*vex
            U(i,3) = rho(i)*E(i)       ! rho*E
            U(i,4) = rho(i)*Y(i)       ! rho*Y
        end do
    end subroutine push_data

    !* =============================================
    !* load forward, from U to field variables
    subroutine load_data(l,r)
    integer::l,r
        do i = l, r
            rho(i) = U(i,1)
            vex(i) = U(i,2) / rho(i)
            E(i)   = U(i,3) / rho(i)
            Y(i)   = U(i,4) / rho(i)
            if (Y(i) < 0.) then
                Y(i) = 0.0
            endif
            if (Y(i) > 1.) then
                Y(i) = 1.0
            endif
            pre(i) = EOS_pre_backward(rho(i), vex(i), E(i), Y(i))
            T(i)   = EOS_T(rho(i), pre(i), Y(i))
        end do
    end subroutine load_data
end module mod_para
