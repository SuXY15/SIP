!------------------------------------------
! para1.f90
! set paras and variables
!------------------------------------------
! parameter setting
!!** ======================================================
integer,parameter::acc = 8

!* control boolean
logical::ifdiff = .true.    ! ifdiff- add diffusion term or not
logical::ifreac = .true.    ! ifreac- add reaction term or not
logical::ifcurv = .true.   ! ifcurv- add curvature or not
logical::ifchec = .false.   ! ifchec- check eigenvalue or not
logical::ifstop = .false.   ! ifstop- check if stop
logical::ifsave = .true.    ! ifsave- save data file or not
!* main parameters
real(kind=acc),parameter::Le=1e5         ! Le-   Lewis number
real(kind=acc),parameter::dx=0.05        ! dx-   space step size
real(kind=acc),parameter::tstp=1         ! tstp- timestep, when to save
real(kind=acc),parameter::tend=2000       ! tend- time end, when to stop
real(kind=acc),parameter::mdt=0.5*dx**2. ! mdt-  minimal dt
real(kind=acc),parameter::cfl=0.1        ! cfl-  Courant-Friedichs-Lewy condition
real(kind=acc),parameter::ep=1.e-16      ! ep-   epsilon
real(kind=acc),parameter::scale=1.0d0
!* scale: the velocity fractor
!  scale = 37555.7641d0  for v0 = 4.66e2 m/s
!  scale = 104.17d0      for v0 = 1.68e5 m/s
!  scale = 1.0d0         for v0 = 1.75e7 m/s

real(kind=acc),parameter::Tin=0.2d0   ! Tin- initial temprature
real(kind=acc),parameter::Tig=1.2d0   ! Tig- ignite temprature
real(kind=acc),parameter::Yin=1.0d0    ! Yin- initial concentration
real(kind=acc),parameter::Yig=1.0d0    ! Yig- ignite concentration
real(kind=acc),parameter::Vin=1.1e-5   ! Vin- initial velocity
real(kind=acc),parameter::Ace=20.0d0     ! Ace- accelerate ratio
real(kind=acc),parameter::Ea=84.0d0    ! Ea-  active energy

integer,parameter::totN=1000000000   ! totN-  total number of timesteps
integer,parameter::bx=100/dx           ! bx-    start points num in r-dir
integer,parameter::nx=120/dx          ! nx-    points num in r-dir
integer,parameter::tx=1/dx           ! tx-    turb num in x-dir
integer,parameter::ex=6              ! ex-    extra num in x-dir
integer,parameter::bxm=bx-ex         ! bxm-   x_min_range
integer,parameter::nxm=nx+ex         ! nxm-   x_max_range
integer,parameter::eqN=4             ! eqN-   num of equations

!!** =========================================================
! parameter of equations
real(kind=acc),parameter::h_ = 6.6260693e-34    ! Plank constant
real(kind=acc),parameter::hbar=1.05457266e-34   ! reduced Plank constant
real(kind=acc),parameter::c_ = 2.99792458e8     ! velocity of light in vacuum
real(kind=acc),parameter::k_ = 1.3806488e-23    ! Boltzmann constant
real(kind=acc),parameter::N_ = 500*6.0221413e23 ! Avogadro constant
real(kind=acc),parameter::rho_0 = 3.5e10
real(kind=acc),parameter::rho_b = 3.5e10
real(kind=acc),parameter::S_c = 466.0d0
real(kind=acc),parameter::T_0 = 0.1e9
real(kind=acc),parameter::T_b = 3.2e9
real(kind=acc),parameter::pi  = 3.141592653589793D0

real(kind=acc)::E_0,E_T0,E_Tb
real(kind=acc)::C_1,C_2
real(kind=acc)::A_1,A_2,qcon
real(kind=acc)::C_t

real(kind=acc)::rho_ratio,rho_outer,rho_inner
real(kind=acc)::E_inner, E_outer

!* control parameters:
integer::i,j,m            ! i,j,m- for loop
integer::nt               ! nt-    num of time-steps
integer::count            ! count- count for turns
real(kind=acc)::dt        ! dt-    time step (changeable)
real(kind=acc)::em        ! em-    for cal time-step
real(kind=acc)::time      ! time-  in simulation
real(kind=acc)::temp      ! temp-  temp val

!* necessary vars for solutions:
real(kind=acc),dimension(bxm:nxm)::x      ! x- x val in range
real(kind=acc),dimension(bxm:nxm)::rho    ! rho- density
real(kind=acc),dimension(bxm:nxm)::vex    ! vex- velocity
real(kind=acc),dimension(bxm:nxm)::pre    ! pre- pressure

real(kind=acc),dimension(bxm:nxm)::E      ! E- energy
real(kind=acc),dimension(bxm:nxm)::H      ! H- enthalpy
real(kind=acc),dimension(bxm:nxm)::Y      ! Y- concentration
real(kind=acc),dimension(bxm:nxm)::T      ! T- temprature
real(kind=acc),dimension(bxm:nxm)::omega

real(kind=acc),dimension(bxm:nxm,eqN)::F         ! F- F in equations
real(kind=acc),dimension(bxm:nxm,eqN)::U         ! U- U in equations
real(kind=acc),dimension(bxm:nxm,eqN)::U0        ! U0- for storing U^(n)
real(kind=acc),dimension(bxm:nxm,eqN)::LU        ! LU- L(U):differential in space
real(kind=acc),dimension(bxm:nxm,eqN)::fffd      ! fffd-

real(kind=acc)::flamePx      ! flamePx- x pos of flame
real(kind=acc)::flameLP      ! flameLP- last pos of flame
real(kind=acc)::flameLT      ! flameLT- flame pos time 
real(kind=acc)::flameVe      ! flameVe- flame velocity

!* get variables global
common /boolPara/     ifdiff, ifreac, ifcurv, ifchec, ifstop, ifsave
common /intPara/      nt, count
common /floatPara/    em, dt, time
common /solutePara/   x, F, U, LU, fffd, omega
common /flamePara/    flamePx, flameLP, flameVe, flameLT
common /boundaryPara/ rho_ratio, rho_inner, rho_outer, E_inner, E_outer

