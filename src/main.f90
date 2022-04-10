program SIP
use mod_init
implicit none

    call init
    call bc
    call init_show
    call save_data
    
    do nt = 1, totN
        call bc       ! set boundary conditions
        call diff     ! calculate the diffusion part
        call enox     ! calculate the convection part
        call enos     ! calculate the source part

        dt = cfl / em * dx   ! cal the time step
        if ((time+dt) .ge. tend) then
            dt = tend - time
            ifstop = .true.
            write(*,*) 'time stop!'
        endif
        time = time + dt
        
        !* now U = U^(n)

        do i = bx, nx
            do j = 1, eqn
                U0(i,j) = U(i,j)
                U(i,j)  = U(i,j) + dt*LU(i,j)
            end do
        end do

        call bc
        call enox
        call enos

        !* now U = U^(1)

        do i = bx, nx
            do j = 1, eqn
                U(i,j) = 0.75*U0(i,j) + 0.25*(U(i,j) + dt*LU(i,j))
            end do
        end do

        call bc
        call enox
        call enos
        
        !* now U = U^(2)
        
        do i = bx, nx
            do j = 1, eqn
                U(i,j) = 1./3.*U0(i,j) + 2./3.*(U(i,j) + dt*LU(i,j))
            end do
        end do
        
        !* now U = U^(n+1)

        !* display calculating status
        write(*,'("nt=",I8,"  dt=",1e10.3,"  t=",1e10.3,"  ",1e10.3)') nt, dt, time

        !* output when sample
        if ((time/tstp) .gt. count) then
            count = count + 1
            call save_data
        endif
        if (ifstop) goto 999
    end do

    !* if 999 stop calculation and finish
    999 continue

    !* output nt and time when finish
    write(*,*) 'Calculation finish!'
    write(*,'("nt=",I8,"  dt=",1e10.3,"  t=",1e10.3,"  ",1e10.3)') nt,dt,time
    if (final_save) then
        call save_data
    endif
    stop

contains

    !* =============================================
    !* boundary conditions
    subroutine bc
    integer::L,R
        !* TYPE 0: symmetric boundary conditions
        if (bc_l.eq.0) then
            do i = 1, ex
                L = bx - i
                R = bx + i - 1
                U(L, 1) = U(R, 1)
                U(L, 2) = U(R, 2)
                U(L, 3) = U(R, 3)
                U(L, 4) = U(R, 4)
            end do
        endif
        if (bc_r.eq.0) then
            do i = 1, ex
                R = nx + i
                L = nx - i + 1
                U(R, 1) = U(L, 1)
                U(R, 2) = U(L, 2)
                U(R, 3) = U(L, 3)
                U(R, 4) = U(L, 4)
            end do
        endif

        !* TYPE 1: opening boundary conditions
        if (bc_l .eq. 1) then
            do i = 1, ex
                L = bx - i
                R = bx
                U(L, 1) = U(R, 1)
                U(L, 2) = U(R, 2)
                U(L, 3) = U(R, 3)
                U(L, 4) = U(R, 4)
            end do
        endif
        if (bc_r .eq. 1) then
            do i = 1, ex
                R = nx + i
                L = nx
                U(R, 1) = U(L, 1)
                U(R, 2) = U(L, 2)
                U(R, 3) = U(L, 3)
                U(R, 4) = U(L, 4)
            end do
        endif

        !* TYPE 2: reflective & adiabatic wall boundary conditions
        if (bc_l.eq.2) then
            do i = 1, ex
                L = bx - i
                R = bx + i
                U(L, 1) = U(R, 1)
                U(L, 2) =-U(R, 2)
                U(L, 3) = U(R, 3)
                U(L, 4) = U(R, 4)
            end do
        endif
        if (bc_r.eq.2) then
            do i=1, ex
                R = nx + i
                L = nx - i
                U(R, 1) = U(L, 1)
                U(R, 2) =-U(L, 2)
                U(R, 3) = U(L, 3)
                U(R, 4) = U(L, 4)
            end do
        endif
    return
    end subroutine bc

    !* =============================================
    !* diff: diffusion term
    subroutine diff
    real(kdp),parameter::aa1=  1./60., aa2= -9./60., aa3= 45./60.
    real(kdp),parameter::aa4=-45./60., aa5=  9./60., aa6= -1./60.
    real(kdp),dimension(bxm:nxm)::ayx       ! ayx- partial Y partial x
    real(kdp),dimension(bxm:nxm)::atx       ! atx- partial T partial x
    real(kdp),dimension(bxm:nxm,eqn)::ll
    real(kdp),dimension(bxm:nxm,eqn)::fff
        
        call load_data(bxm,nxm)

        !* cal aux and atx differential
        do i = bx-ex/2, nx+ex/2
           ayx(i) = aa1*Y(i+3) + aa2*Y(i+2) + aa3*Y(i+1) &
                & + aa4*Y(i-1) + aa5*Y(i-2) + aa6*Y(i-3)
           atx(i) = aa1*T(i+3) + aa2*T(i+2) + aa3*T(i+1) &
                & + aa4*T(i-1) + aa5*T(i-2) + aa6*T(i-3)
        end do

        !* load data for ll the eq3. and eq4. need for conduction
        do i = bx-ex/2, nx+ex/2
            ll(i,1) = 0.0d0
            ll(i,2) = 0.0d0
            ll(i,3) = lambda * atx(i) / dx ! + rho(i) * alpha / Le * Q * ayx(i) / dx
            ll(i,4) = rho(i) * alpha / Le * ayx(i) / dx
        end do

        !* cal fffd
        do i = bx, nx
            do m = 1, eqn
                fffd(i,m) = aa1*ll(i+3,m) + aa2*ll(i+2,m) + aa3*ll(i+1,m) &
                        & + aa4*ll(i-1,m) + aa5*ll(i-2,m) + aa6*ll(i-3,m)
            end do
        end do

        !* 2/r part of spherical coordinate
        if (ifcurv) then
            do i = bx, nx
                fffd(i,3) = fffd(i,3) + lambda * atx(i) * 2./(x(i)+ep)
                ! fffd(i,3) = fffd(i,3) + rho(i) * Q *  alpha / Le * ayx(i) * 2./(x(i)+ep)
                fffd(i,4) = fffd(i,4) + rho(i) * alpha / Le * ayx(i) * 2./(x(i)+ep)
            end do
        end if
    end subroutine diff


    !* =============================================
    !* enox: convection term
    subroutine enox
    real(kdp),dimension(eqn)::am                   ! am-
    real(kdp),dimension(bxm:nxm)::w                ! w- weight
    real(kdp),dimension(bxm:nxm)::ev1,ev2,ev3,ev4  ! eigenvalues
    real(kdp),dimension(bxm:nxm,eqn,eqn)::evl,evr  ! eigen vectors
    real(kdp),dimension(bxm:nxm,eqn)::dF,dU        ! dF-    dU-
    real(kdp),dimension(bxm:nxm,eqn)::ff,fh        ! ff-    fh-
    real(kdp),dimension(bxm:nxm,eqn,2)::gg,hh      ! gg-    hh-
    real(kdp)::w0, w1                              ! w0 w1 weight for interp
    real(kdp)::rho_m, vex_m, pre_m, E_m, Y_m, H_m, Q_m
    real(kdp)::K_1, K_2, K_3, K_4
    real(kdp)::c_m, err, cvel
    real(kdp)::t0,t1,t2,t3,s1,s2,s3
    integer::i, j, k, ip, k0, k1, m1
        
        
        call load_data(bxm,nxm)
        
        !* cal eigenvalues
        do i = bxm, nxm
            ! sound velocity
            cvel = sqrt(gamma*pre(i)/rho(i))

            ! load F
            F(i,1) = rho(i)*vex(i)                        ! rho*u
            F(i,2) = rho(i)*vex(i)**2 + pre(i)            ! rho*u**2 + p
            F(i,3) = rho(i)*vex(i)*E(i) + vex(i)*pre(i)   ! rho*u*E + u*p
            F(i,4) = rho(i)*vex(i)*Y(i)                   ! rho*u*Y

            w(i) = sqrt(abs(rho(i)))  ! weight 
            
            !* eigenvalues
            ev1(i)=abs(vex(i)-cvel)
            ev2(i)=abs(vex(i))
            ev3(i)=abs(vex(i))
            ev4(i)=abs(vex(i)+cvel)
        end do

        !* cal am
        am(1:4) = ep
        do i = bxm, nxm
            am(1) = max( am(1), ev1(i) )
            am(2) = max( am(2), ev2(i) )
            am(3) = max( am(3), ev3(i) )
            am(4) = max( am(4), ev4(i) )
        end do

        !* get em
        em = ep
        am(1) = am(1) * 1.1
        am(2) = am(2) * 1.1
        am(3) = am(3) * 1.1
        am(4) = am(4) * 1.1
        em = max(em, max(am(3),am(4)))

        
        Q_m = Q     ! if p = rho(RT/(gamma-1)+QY)/(gamma-1), Q_m = 0.
                    ! if p = rho*R*T                         Q_m = Q
        
        !* cal evr & evl
        do i = bx-1, nx
            !* Compute e'vectors using Roe's average:
            ip = i+1
            w0 = w(i)/(w(i)+w(ip))
            w1 = 1.-w0
            
            rho_m = w0*rho(i) + w1*rho(ip)
            vex_m = w0*vex(i) + w1*vex(ip)
            E_m   = w0*E(i)   + w1*E(ip)
            Y_m   = w0*Y(i)   + w1*Y(ip)
            pre_m = EOS_pre_backward(rho_m, vex_m, E_m, Y_m)
            H_m   = E_m + pre_m / rho_m
            
            c_m = sqrt(gamma*pre_m/rho_m)
            
            K_1 = gamma1 * vex_m * vex_m / 2
            K_2 = gamma1 * (-vex_m)
            K_3 = gamma1
            K_4 = gamma1 * (-Q_m)
            
            !* evr: right e'vectors of Jaccobi Matrix
            evr(i,1,1) = 1./c_m
            evr(i,1,2) = 1./c_m
            evr(i,1,3) = 0.
            evr(i,1,4) = 1./c_m
            
            evr(i,2,1) = vex_m / c_m - 1.
            evr(i,2,2) = vex_m / c_m
            evr(i,2,3) = 0.
            evr(i,2,4) = vex_m / c_m + 1.
            
            evr(i,3,1) = H_m / c_m - vex_m
            evr(i,3,2) = vex_m * vex_m / c_m / 2.
            evr(i,3,3) = Q_m / c_m
            evr(i,3,4) = H_m / c_m + vex_m
            
            evr(i,4,1) = Y_m / c_m
            evr(i,4,2) = 0.
            evr(i,4,3) = 1. / c_m
            evr(i,4,4) = Y_m / c_m

            !* evl: left e'vectors of Jaccobi Matrix
            evl(i,1,1) =  vex_m / 2. + gamma1 * vex_m * vex_m / c_m / 4.
            evl(i,1,2) =  -1./2. - gamma1 * vex_m / c_m / 2.
            evl(i,1,3) =  gamma1 / c_m / 2.
            evl(i,1,4) = -gamma1 * Q_m / c_m / 2.
            
            evl(i,2,1) =  c_m - gamma1 * vex_m * vex_m / c_m / 2.
            evl(i,2,2) =  gamma1 * vex_m / c_m
            evl(i,2,3) = -gamma1 / c_m
            evl(i,2,4) =  gamma1 * Q_m / c_m

            evl(i,3,1) = -Y_m * gamma1 * vex_m * vex_m / c_m / 2.
            evl(i,3,2) =  Y_m * gamma1 * vex_m / c_m
            evl(i,3,3) = -Y_m * gamma1 / c_m 
            evl(i,3,4) =  c_m + Y_m * gamma1 * Q_m / c_m

            evl(i,4,1) = -vex_m / 2. + gamma1 * vex_m * vex_m / c_m / 4.
            evl(i,4,2) =  1. / 2. - gamma1 * vex_m / c_m / 2.
            evl(i,4,3) =  gamma1 / c_m / 2.
            evl(i,4,4) = -gamma1 * Q_m / c_m / 2.

            !* check eigenvalues
            if (ifchec) then
                call check_eigen(i,evl,evr)
            endif
        end do
    
        !* cal dF & dU
        do i = bxm, nxm-1
            dF(i,:) = F(i+1,:) - F(i,:)  ! derta F
            dU(i,:) = U(i+1,:) - U(i,:)  ! derta U
        end do

        !* init ff
        ff(bxm:nxm,1:eqn) = 0.0d0

        do m = 1, eqn
        !* Project the relevant first undivided differences into the
        !* 'm'th characteristic field
            !* cal gg
            gg(bxm:nxm-1,:,1) = 0.5*(am(m)*dU(bxm:nxm-1,:) + dF(bxm:nxm-1,:))
            gg(bxm:nxm-1,:,2) = 0.5*(am(m)*dU(bxm:nxm-1,:) - dF(bxm:nxm-1,:))

            !* cal hh
            do m1 = 1, 4
                k0 = m1 - 3
                k1 = 3 - m1
                do i = bx-1, nx
                    hh(i,m1,1) = evl(i,m,1)*gg(i+k0,1,1) + evl(i,m,2)*gg(i+k0,2,1) &
                             & + evl(i,m,3)*gg(i+k0,3,1) + evl(i,m,4)*gg(i+k0,4,1)
                    hh(i,m1,2) = evl(i,m,1)*gg(i+k1,1,2) + evl(i,m,2)*gg(i+k1,2,2) &
                             & + evl(i,m,3)*gg(i+k1,3,2) + evl(i,m,4)*gg(i+k1,4,2)
                end do
            end do

            !* Compute numerical flux in each characteristic field:
            do m1 = 1, 2
                do i = bx-1, nx
                    t1 = hh(i,1,m1)-hh(i,2,m1)
                    t2 = hh(i,2,m1)-hh(i,3,m1)
                    t3 = hh(i,3,m1)-hh(i,4,m1)

                    t0 = hh(i,1,m1)-3.*hh(i,2,m1)
                    s1 = 1./(ep+13.*t1*t1+3.*t0*t0)**2.

                    t0 = hh(i,2,m1)+hh(i,3,m1)
                    s2 = 6./(ep+13.*t2*t2+3.*t0*t0)**2.

                    t0 = 3.*hh(i,3,m1)-hh(i,4,m1)
                    s3 = 3./(ep+13.*t3*t3+3.*t0*t0)**2.

                    t0 = 1./(s1+s2+s3)
                    s1 = s1*t0
                    s3 = s3*t0

                    ff(i,m) = ff(i,m) + (s1*(t2-t1)+(0.5*s3-0.25)*(t3-t2))/3.
                end do
            end do
        end do
        !*  end of the big loop

        !* Project the numerical flux back to the physical space:
        do m = 1, eqn 
            do i = bx-1, nx
                fh(i,m) = (-F(i-1,m) + 7.*(F(i,m)+F(i+1,m))-F(i+2,m))/12. &
                      & + evr(i,m,1)*ff(i,1) + evr(i,m,2)*ff(i,2) &
                      & + evr(i,m,3)*ff(i,3) + evr(i,m,4)*ff(i,4)
            end do
        end do

        cal_f: do i = bx, nx
            LU(i,1) = (fh(i-1,1)-fh(i,1))/dx
            LU(i,2) = (fh(i-1,2)-fh(i,2))/dx
            LU(i,3) = (fh(i-1,3)-fh(i,3))/dx
            LU(i,4) = (fh(i-1,4)-fh(i,4))/dx
            
            !* diffusion part
            if (ifdiff) then
                LU(i,3) = LU(i,3) + 0.25d0*fffd(i,3)/dx
                LU(i,4) = LU(i,4) + 0.25d0*fffd(i,4)/dx
            end if

            !* 2/r part of spherical coordinate
            if (ifcurv) then
                LU(i,1) = LU(i,1) - 2.*F(i,1)/(x(i)+ep)
                LU(i,2) = LU(i,2) - 2.*(F(i,2)-pre(i))/(x(i)+ep)
                LU(i,3) = LU(i,3) - 2.*F(i,3)/(x(i)+ep)
                LU(i,4) = LU(i,4) - 2.*F(i,4)/(x(i)+ep)
            end if
        end do cal_f
    return
    end subroutine enox


    !* =============================================
    !* enos: source term
    subroutine enos
    real(kdp)::theta,s,t0
    
        call load_data(bx,nx)

        ! cal reaction
        do i=bx,nx
            if (ifreac) then
                omega(i) = reaction_rate(Y(i), T(i))
                
                ! check if NaN arouse, and output wrong message
                call check_nan(i)
                
                LU(i,4) = LU(i,4) - rho(i) * omega(i)
            end if
        end do
    return
    end subroutine enos


    !* =============================================
    !* save data
    subroutine save_data
    integer::number=0
    character(len=15)::filename
    real(kdp)::t1
    real(kdp)::fluidVe

        call load_data(bx,nx)
        
        flamePx = x(bx)

        !* cal flame position with half height
        t1 = (T(bx)+T(nx))/2.
        do i = bx, nx
            if ((T(i)-t1)*(T(i-1)-t1).lt.0) then
                flamePx = (t1-T(i-1))/(T(i)-T(i-1))*(x(i)-x(i-1))+x(i-1)
                fluidVe = (flamePx-x(i-1))/(x(i)-x(i-1))*(vex(i)-vex(i-1))+vex(i-1)
                if (count .gt. 1) then
                    flameVe = (flamePx-flameLP)/tstp
                end if
                flameLP = flamePx
            endif
        end do

        !* every sample time, record
        open(202,file=trim(adjustl(data_dir))//"flame_position.dat",access='APPEND')
        write(202,'(5e15.7)') time, flamePx, flameVe, fluidVe, flameVe-fluidVe
        close(202)

        if (ifsave) then
            write(filename,'(I4.4)') number
            open(201,file=trim(adjustl(data_dir))//'data'//trim(adjustl(filename))//'.dat')
                write(201,*) 'TITLE="t =', number*tstp, '"'
                write(201,301)
                write(201,*) 'ZONE SOLUTIONTIME=',number*tstp
                do i = bxm, nxm
                    write(201,201) x(i), output(pre(i)), output(rho(i)), output(E(i)), &
                                output(T(i)), output(Y(i)), output(vex(i)), output(omega(i))
                end do
            close(201)
        endif
    301 format(' VARIABLES=',1x,'"x"',10x,'"pre"',10x,'"rho"',12x,'"E"',12x,'"T"',12x,'"Y"',10x,'"vex"',8x,'"omega"')
    201 format(8e15.7)
        number=number+1
    return
    end subroutine save_data

    !* =============================================
    !* handle non-numerical conditions
    real(kdp) function output(input)
    real(kdp),intent(in)::input
        output = input
        if (abs(input).le.(.1e-36)) then
            output = 0.0d0
        else if (input > .1e30) then
            write(*,*) "Large Number arouse!"
            output = .1e30
        else if (isnan(input)) then
            ! output = 0.0d0
            if (ifstop .eqv. .false.) then
                write(*,*) "NaN in save_data:output!"
                ifstop = .true.
                final_save = .false.
            endif
        endif
    end function output

    !* =============================================
    !* check Not A Number case
    subroutine check_nan(i)
    integer::i
    real(kdp)::temp
        if (isnan(omega(i))) then
            write(*,'(A)')               'ERROR: NaN in omega(i)   |----------------------'
            write(*,'(A,I16,  A,E12.5)') 'Old!! i:', i,          ' |omega:',omega(i)
            write(*,'(A,E16.5,A,E12.5)') '      T:', T(i),       ' |  rho:',rho(i)
            write(*,'(A,E16.5,A,E12.5)') '    pre:', pre(i),     ' |  vex:',vex(i)
            write(*,'(A,E16.5,A,E12.5)') '      Y:', Y(i),       ' |    E:',E(i)
            write(*,'(A,E16.5,A,E12.5)') '  rho*E:', rho(i)*E(i),' | U(3):',U(i,3)
            write(*,'(A)')               '-------------------------|----------------------'
            ifsave = .true.
            call save_data
            stop
        endif
    end subroutine check_nan

    !* =============================================
    !* check Eigen
    subroutine check_eigen(i,evl,evr)
    integer::i
    real(kdp)::temp, err
    real(kdp),dimension(bxm:nxm,eqn,eqn)::evl,evr  ! eigen vectors
        !* check eigenvalues
        !write(*,*) "i:",i,"v(i)",vex(i),"v(ip):",vex(i+1)
        !do m = 1, 4
        !    do k = 1, 4
        !        write(*,'("row:",I4,"  col:",I4,"  evr:",1e15.7, "  evl:",e15.7)') &
        !                    m, k, evr(i,m,k), evl(i,m,k)
        !    enddo
        !enddo

        do m = 1, 4
            do k = 1, 4
                temp = 0.0d0
                do j = 1, 4
                    temp = temp + evr(i,m,j)*evl(i,j,k)
                    !write(*,'("m:",I4,"  k:",I4,"  j:",I4,"  val",1e15.7)') &
                    !        m, k, j, evl(i,m,j)*evr(i,j,k)
                end do
                if (m .eq. k) then
                    err = temp - 1.
                else
                    err = temp - 0.
                endif
                if (abs(err) > 1e-6) then
                    write(*,'("row:",I4,"  col:",I4,"  I:",1e15.7,"  err:",1e15.7)') &
                            m, k, temp, err
                    ifstop = .true.
                endif
            enddo
        enddo
    end subroutine check_eigen

end program SIP
