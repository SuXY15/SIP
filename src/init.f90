module mod_init
use mod_para
implicit none

contains
    !* =========================================================
    !* global initialization
    subroutine init
        call init_namelist
        call init_variable
        call init_field
    end subroutine init


    !* =========================================================
    !* read namelist 'src/input.nml'
    subroutine init_namelist
    logical::if_exist
        
        if (iargc()>0) then
            call getarg(1, nml_path)
        endif

        inquire(file=trim(nml_path), exist=if_exist)
        if (if_exist) then
            open(250, file=trim(nml_path), action="read")
            read(250, nml=PARA)
            close(250)
            write(*,*) 'File '//trim(nml_path)//' loaded.'
        else
            write(*,*) 'File '//trim(nml_path)//' not exist!'
            stop
        endif
    end subroutine init_namelist


    !* =========================================================
    !* initialize variables
    subroutine init_variable
        mdt = .5*dx**2      ! mdt-  minimal dt
        bx  = int(bl/dx)    ! bx-    start points num in r-dir
        nx  = int(nl/dx)    ! nx-    points num in r-dir
        tx  = int(tl/dx)    ! tx-    turb num in x-dir
        px  = int(pl/dx)    ! px-    plane num in x-dir
        bxm = bx-ex         ! bxm-   x_min_range
        nxm = nx+ex         ! nxm-   x_max_range
        
        gamma1  = gamma - 1
    end subroutine init_variable


    !* =========================================================
    !* initialize field values
    subroutine init_field
    real(kdp)::t1, t2, Pin
        allocate(  x(bxm:nxm), rho(bxm:nxm), vex(bxm:nxm), pre(bxm:nxm))
        allocate(  E(bxm:nxm),   H(bxm:nxm),   Y(bxm:nxm),   T(bxm:nxm), omega(bxm:nxm))
        allocate(  F(bxm:nxm, eqn),    U(bxm:nxm, eqn), U0(bxm:nxm, eqn), &
                  LU(bxm:nxm, eqn), fffd(bxm:nxm, eqn) )
        allocate(evl(bxm:nxm, eqn, eqn), evr(bxm:nxm, eqn, eqn))
        !* init x in range
        do i = bxm, nxm
            x(i) = dx*i
        end do

        ! px  - the length of high temperature ignition plane
        ! px  - the length of tanh curved transition temperature profile
        ! slx - the left side of non-initial temperature part
        ! srx - the right side of non-initial temperature part
        ! tlx - the left side of tanh curved transition temperature profile
        ! trx - the right side of tanh curved transition temperature profile

        ! - - -                       !
        !       \                     !
        !        \                    ! 
        !         \                   !
        !          -------------------! 
        !|    |    |                  !
        !  px   tx                    !
        !slx      srx                 !
        !bx  tlx  trx                 !nx
        if (direct .eq. 1) then
            slx = bx
            tlx = slx + px
            trx = tlx + tx
            srx = trx
        endif

        !                       - - - !
        !                     /       !
        !                    /        ! 
        !                   /         !
        !-------------------          !
        !                  |    |    |!
        !                    tx   px  !
        !                 slx      srx!
        !bx               tlx  trx    !nx
        if (direct .eq. -1) then
            srx = nx
            trx = srx - px
            tlx = trx - tx
            slx = tlx
        endif

        !* init with const pressure
        Pin = EOS_pre_forward(Din, Tin, Yin)
        do i = bx, nx
            pre(i) = Pin
            vex(i) = Vin
            Y(i)   = Yin
            T(i)   = Tin
        end do
                
        !* burnt/ignition temperature
        do i = slx, tlx
            T(i) = Tig
            Y(i) = Yig
        end do
        do i = trx, srx
            T(i) = Tig
            Y(i) = Yig
        end do

        !* tanh curv
        do i = tlx, trx
            t1 = direct*(10.*(tlx-i)/(trx-tlx)+5.)
            t1 = (exp(t1)-exp(-t1))/(exp(t1)+exp(-t1))+&
                   & (exp(5.)-exp(-5.))/(exp(-5.)+exp(5.))
            T(i) = t1*(Tig-Tin)*(exp(5.)+exp(-5.))/2./(exp(5.)-exp(-5.))+Tin
            Y(i) = t1*(Yig-Yin)*(exp(5.)+exp(-5.))/2./(exp(5.)-exp(-5.))+Yin
        end do
        
        !* cal rho from pre & T, cal E from rho & T
        do i = bx,nx
            rho(i) = EOS_rho( pre(i), T(i), Y(i) )
            ! pre(i) = EOS_pre_forward( rho(i), T(i) )
            E(i) = EOS_E( rho(i), pre(i), vex(i), Y(i) )
        end do
        
        write(*,*) rho(bx)
        
        !* init for U
        call push_data(bx, nx)
        
        !* init setting values
        time = 0.
        count = 1
        flameLP = 0.0d0
        flameLT = 0.0d0

        open(202, file=trim(adjustl(data_dir))//"flame_position.dat", action='WRITE')
        write(202,*) " time           flamePx        flameVe        fluidVe        flameVe'"
        close(202)
        write(*,*) "Initialize finished."
    end subroutine init_field
     
    !* =========================================================
    !* show initialization results
    subroutine init_show
        !* display some important value
        write(*,'(A)')               '-------------------------------------------'
        write(*,'(A)')               'para| function            | value'
        write(*,'(A)')               '-------------------------------------------'
        write(*,'(A,I12)')           'eqn | num of equations    |', eqn
        write(*,'(A,I12)')           'nx  | points num in x-dir |', nx
        write(*,'(A,I12)')           'bx  | bounds num in x-dir |', bx
        write(*,'(A,E12.5)')         'CFL | Courant Fried Lewy  |', cfl
        write(*,'(A,E12.5)')         'tend| final printing time |', tend
        write(*,'(A,E12.5)')         'Q   | heat value          |', Q
#ifdef supernova
        write(*,'(A,E12.5)')         'C1  | eos para C1         |', eos_C1
        write(*,'(A,E12.5)')         'C2  | eos para C2         |', eos_C2
#endif
        !* display some initialize results
        write(*,'(A)')               '-------------------------------------------'
        write(*,'(A)')               'Init| Non perturb area    | perturbed area'
        write(*,'(A,E19.5,A,E12.5)') 'pre | ', pre(bx)       ,' |', pre(nx)
        write(*,'(A,E19.5,A,E12.5)') 'rho | ', rho(bx)       ,' |', rho(nx)
        write(*,'(A,E19.5,A,E12.5)') 'vex | ', vex(bx)       ,' |', vex(nx)
        write(*,'(A,E19.5,A,E12.5)') 'E   | ', E(bx)         ,' |', E(nx)
        write(*,'(A,E19.5,A,E12.5)') 'Y   | ', Y(bx)         ,' |', Y(nx)
        write(*,'(A,E19.5,A,E12.5)') 'T   | ', T(bx)         ,' |', T(nx)
        write(*,'(A,E19.5,A,E12.5)') 'U(1)| ', U(bx,1)       ,' |', U(nx,1)
        write(*,'(A,E19.5,A,E12.5)') 'U(2)| ', U(bx,2)       ,' |', U(nx,2)
        write(*,'(A,E19.5,A,E12.5)') 'U(3)| ', U(bx,3)       ,' |', U(nx,3)
        write(*,'(A,E19.5,A,E12.5)') 'U(4)| ', U(bx,4)       ,' |', U(nx,4)
        write(*,'(A,E19.5,A,E12.5)') '-------------------------------------------'
    end subroutine init_show
end module mod_init
