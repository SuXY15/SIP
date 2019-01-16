program main
implicit none
!***************************************************************
!* main program to solve compressible SN Ia supernova problem
!***************************************************************
    include 'para1.f90'
    include 'para2.f90'

    !* initialization and storage
    call init     ! initialization
    call bc       ! set boundary conditions
    call setup    ! display setup message
    call res_data ! restore initial field

    !* Runge-Kutta time stepping:
    main_do: do nt=1,totN
        call check_continuity
        call bc       ! set boundary conditions
        call vis_fx_L ! calculate the diffusion part
        call enox     ! calculate the convection part
        call enos     ! calculate the source part

        dt = cfl/em*dx   ! cal the time step
        if((time+dt).ge.tend) then
            dt = tend-time
            ifstop = .true.
            write(*,*) 'time stop!'
        endif
        time=time+dt
        
        !* now U = U^(n)
        
        cal_u0u: do i=bx,nx
            do j=1,eqN
                U0(i,j)= U(i,j)
                U(i,j) = U(i,j)+dt*LU(i,j)
            end do
        end do cal_u0u

        call check_continuity
        call bc
        call enox
        call enos
        
        !* now U = U^(1)
        
        cal_u1: do i=bx,nx
            do j=1,eqN
                U(i,j)=0.75*U0(i,j)+0.25*(U(i,j)+dt*LU(i,j))
            end do
        end do cal_u1

        call check_continuity
        call bc
        call enox
        call enos
        
        !* now U = U^(2)
        
        cal_u2: do i=bx,nx
            do j=1,eqN
                U(i,j)=(2./3.)*(0.5*U0(i,j)+U(i,j)+dt*LU(i,j))
            end do
        end do cal_u2
        
        !* now U = U^(n+1)

        !* display calculating status
        write(*,'("nt=",I8,"  dt=",1e10.3,"  t=",1e10.3,"  ",1e10.3)') nt,dt,time

        !* output when sample
        if((time*1./tstp).gt.count) then
            count = count + 1
            call res_data
        endif
        if(ifstop) goto 999
    end do main_do

    !* if 999 stop calculation and finish
    999 continue

    !* output nt and time when finish
    write(*,*) 'Calculation finish!'
    write(*,'("nt=",I8,"  dt=",1e10.3,"  t=",1e10.3,"  ",1e10.3)') nt,dt,time
    call res_data
    stop
end


!!** setup ----------------------------------
subroutine setup
include 'para1.f90'
include 'para2.f90'
!---------------------------
! functions:  read and show parameters
!---------------------------
    !* display some important value
    write(*,*) ' Message you want to show         '
    write(*,*) '-------------------------------------------'
    write(*,*) 'para|-function------------|------value-----'
    write(*,*) '-------------------------------------------'
    write(*,*) 'eqN | num of equations:   ', eqN
    write(*,*) 'nx  | points num in x-dir:', nx
    write(*,*) 'CFL | CFL =               ', cfl
    write(*,*) 'tend| final printing time:', tend
    write(*,*) 'totN| max num of timestep:', totN
    write(*,*) 'E_0 |                     ', E_0
    write(*,*) 'E_T0|                     ', E_T0
    write(*,*) 'E_Tb|                     ', E_Tb
    write(*,*) 'A_1 |                     ', A_1
    write(*,*) 'A_2 |                     ', A_2
    write(*,*) 'qcon|                     ', qcon
    write(*,*) '-------------------------------------------'
    !* display some initilize results
    write(*,*) ' Initialize results         '
    write(*,*) '----------------------------------------------------'
    write(*,*) 'Init| Non disturbulence area     | disturbulence area'
    write(*,*) 'rho | ', rho(bx)             ,'  |', rho(nx-tx)
    write(*,*) 'vex | ', vex(bx)             ,'  |', vex(nx-tx)
    write(*,*) 'Y   | ', Y(bx)               ,'  |', Y(nx-tx)
    write(*,*) 'E   | ', E(bx)               ,'  |', E(nx-tx)
    write(*,*) 'T   | ', T(bx)               ,'  |', T(nx-tx)
    write(*,*) 'pre | ', pre(bx)             ,'  |', pre(nx-tx)
    write(*,*) 'U(2)| ', U(bx,2)             ,'  |', U(nx-tx,2)
    write(*,*) 'U(3)| ', U(bx,3)             ,'  |', U(nx-tx,3)
    write(*,*) 'U(4)| ', U(bx,4)             ,'  |', U(nx-tx,4)
    write(*,*) '----------------------------------------------------'
return
end subroutine setup


!!** check_continuity -------------------------
subroutine check_continuity
    include 'para1.f90'
    include 'para2.f90'
    ! check continutiy:
    do i=bx,nx
        if((U(i,1)-U(i-1,1))*(U(i,1)-U(i+1,1))/U(i,1)/U(i,1)>1e-6 .and.&
          (U(i+1,1)-U(i,1))*(U(i+1,1)-U(i+2,1))/U(i+1,1)/U(i+1,1)>1e-6) then
            write(*,*) "Revising rho:", i
            U(i,1) = U(i-1,1)*2./3.+U(i+2,1)/3.
            U(i+1,1) = U(i-1,1)/3.+U(i+2,1)*2./3.
        end if
    end do
return
end subroutine check_continuity


!!** init -------------------------------------
subroutine init
    include 'para1.f90'
    include 'para2.f90'

    !* init x in range
    init_x: do i=bxm,nxm
        x(i)=dx*i
    end do init_x

    !* init with Tin and Tig - pre = const
    do i=bx,nx
        pre(i) = E_T0/E_Tb
        vex(i) = Vin
        Y(i)   = Yin
        T(i)   = Tin*(T_b-T_0)+T_0
    end do
    
    !* init with turbulence area
    turb: do i=nx-tx,nx        ! tx- turbulence x range
        Y(i) = Yig
        T(i) = Tig*(T_b-T_0)+T_0
    end do turb

    !* tanh curv smooth for gap continutiy
    curv: do i=(nx-2*tx),nx-tx
        temp = 10./((2.-1.)*tx)*(i-nx+2.*tx)-5.
        temp = (exp(temp)-exp(-temp))/(exp(temp)+exp(-temp))+&
               & (exp(5.)-exp(-5.))/(exp(-5.)+exp(5.))
        tmp2 = temp*(Tig-Tin)*(exp(5.)+exp(-5.))/2./(exp(5.)-exp(-5.))+Tin
        T(i) = tmp2*(T_b-T_0)+T_0
        Y(i) = temp*(Yig-Yin)*(exp(5.)+exp(-5.))/2./(exp(5.)-exp(-5.))+Yin
    end do curv
    
    !* cal rho from pre & T, cal E from rho & T
    do i=bx,nx
        temp = C_2*T(i)**2./2./C_1
        rho(i) = ((pre(i)*E_Tb/C_1+temp*temp)**(1./2.)-temp)**(3./2.)
        E(i) = pre(i)/rho(i)+qcon*Y(i)+0.5/A_1*vex(i)**2.
    end do

    rho_inner = rho(bx)
    rho_outer = rho(nx)
    E_inner = E(bx)
    E_outer = E(nx)

    !* init for U
    init_u: do i=bx,nx
        U(i,1) = rho(i)            ! rho
        U(i,2) = rho(i)*vex(i)/A_1 ! rho*vex
        U(i,3) = rho(i)*E(i)       ! rho*E
        U(i,4) = rho(i)*Y(i)       ! rho*Y
    end do init_u

    !* init setting values
    time = 0.
    count = 1
    flameLP = 0.0d0
    flameLT = 0.0d0

    open(202, file="flame_position.txt",action='WRITE')
    write(202,*) " time           flamePx        flameVe"
    close(202)
return
end subroutine init


!!** bc -------------------------------
subroutine bc
include 'para1.f90'
include 'para2.f90'
!--------------------------------------
!     ex |           nx            | ex
! rho +++|+++*******************+++|+++
! vex +++|+++*******************+++|+++
! E   +++|+++*******************+++|+++
! Y   +++|+++*******************+++|+++
!--------------------------------------
    if(ifcurv) then
        !* spherical boundary conditions
        ! rho_ratio = 1-U(bx,2)/U(bx,1)/(bx*dx)*3*dt
        ! bc_inner_s: do i=1,ex
        !     U(bx-i,1) = U(bx+i,1)*rho_ratio
        !     U(bx-i,2) = U(bx+i,2)*rho_ratio
        !     U(bx-i,3) = U(bx+i,3)*rho_ratio
        !     U(bx-i,4) = U(bx+i,4)*rho_ratio
        ! end do bc_inner_s
        bc_inner_s: do i=1,ex
            U(bx-i,1) = rho_inner
            U(bx-i,2) = 2.*U(bx,2)-U(bx+i,2)
            U(bx-i,3) = rho_inner*E_inner
            U(bx-i,4) = rho_inner*1.0
        end do bc_inner_s
        bc_outer_s: do i=1,ex
            U(nx+i,1) = rho_outer
            U(nx+i,2) = 2.*U(nx,2)-U(nx-i,2)
            U(nx+i,3) = rho_outer*E_outer
            U(nx+i,4) = rho_outer*1.0
        end do bc_outer_s
        ! bc_outer_s: do i=1,ex
        !     U(nx+i,1) = 2*U(nx+i-1,1)-U(nx+i-2,1)
        !     U(nx+i,2) = 2*U(nx+i-1,2)-U(nx+i-2,2)
        !     U(nx+i,3) = 2*U(nx+i-1,3)-U(nx+i-2,3)
        !     U(nx+i,4) = 2*U(nx+i-1,4)-U(nx+i-2,4)
        ! end do bc_outer_s
    else
        !* liner boundary conditions
        bc_inner_l: do i=1,ex
            U(bx-i,1) = U(bx+i-1,1)
            U(bx-i,2) = U(bx+i-1,2)
            U(bx-i,3) = U(bx+i-1,3)
            U(bx-i,4) = U(bx+i-1,4)
        end do bc_inner_l
        
        bc_outer_l: do i=1,ex
            U(nx+i,1) = U(nx-i+1,1)
            U(nx+i,2) = U(nx-i+1,2)
            U(nx+i,3) = U(nx-i+1,3)
            U(nx+i,4) = U(nx-i+1,4)
        end do bc_outer_l
    end if 
return
end subroutine bc


!!** res_data -----------------------------------
subroutine res_data
    include 'para1.f90'
    ! interface for output format
    interface
        real(kind=8) function output(input)
            real(kind=8)::input
        end function output
        real(kind=16) function output16(input)
            real(kind=16)::input
        end function output16
    end interface

    integer::number=0
    character(len=15)::filename
    real(kind=acc)::tmp1,tmp2   ! temp values
    include 'para2.f90'

    !* load data
    do i=bxm,nxm
        rho(i) = U(i,1)
        vex(i) = U(i,2)/rho(i)*A_1
        E(i)   = U(i,3)/rho(i)
        Y(i)   = U(i,4)/rho(i)

        pre(i)=rho(i)*(E(i)-0.5/A_1*vex(i)**2.-qcon*Y(i))
        T(i) = ((pre(i)*E_Tb-C_1*rho(i)**(4./3.))/C_2/rho(i)**(2./3.))**(1./2.)
    end do
    
    flamePx = x(bx)
    tmp1 = 0.0d0
    tmp2 = 0.0d0

    !* cal flame position with half height
    tmp1 = (T(bx)+T(nx))/2.
    cal_flame: do i=bx,nx
        if ((T(i)-tmp1)*(T(i-1)-tmp1).lt.0) then
            flamePx = (tmp1-T(i-1))/(T(i)-T(i-1))*(x(i)-x(i-1))+x(i-1)
            if (flameLP.ne.0) then
                flameVe = (flamePx-flameLP)/tstp
            end if
            flameLP = flamePx
        endif
    end do cal_flame

    !* every sample time, record
    open(202,file='flame_position.txt',access='APPEND')
    write(202,'(3e15.7)') time,flamePx,flameVe
    close(202)

    if (ifsave) then
        write(filename,'(I4.4)') number
        open(201,file='data'//trim(adjustl(filename))//'.txt')
        write(201,301)
        open(201,file='data'//trim(adjustl(filename))//'.txt')
        !* write data bx:nx
        write_data: do i=bx,nx
            write(201,201) x(i),pre(i),rho(i),E(i),T(i),output(Y(i)),output(vex(i)),output(omega(i))
            ! write(201,201) x(i),pre(i),rho(i),E(i),T(i),output16(Y(i)),output16(vex(i)),output16(omega(i))
        end do write_data

        close(201)
    endif
301 format(14x,'x',12x,'pre',12x,'rho',14x,'E',14x,'T',14x,'Y',12x,'vex',10x,'omega')
201 format(8e15.7)
    number=number+1
return
end subroutine res_data


!!** vis_fx_L -------------------------
subroutine vis_fx_L
    include 'para1.f90'
    real(kind=acc),parameter::aa1=  1./60., aa2= -9./60., aa3= 45./60.
    real(kind=acc),parameter::aa4=-45./60., aa5=  9./60., aa6= -1./60.
    real(kind=acc),dimension(bxm:nxm)::ayx       ! ayx- partial Y partial x
    real(kind=acc),dimension(bxm:nxm)::atx       ! atx- partial T partial x
    real(kind=acc),dimension(bxm:nxm,eqN)::ll
    real(kind=acc),dimension(bxm:nxm,eqN)::fff
    include 'para2.f90'
    
    !* load data
    do i=bxm,nxm
        rho(i) = U(i,1)
        vex(i) = U(i,2)/rho(i)*A_1
        E(i)   = U(i,3)/rho(i)
        Y(i)   = U(i,4)/rho(i)
    end do
    
    !* cal pressuer and temperature from rho/vex/E/Y
    load_pre: do i=bxm, nxm
        pre(i)=rho(i)*(E(i)-0.5/A_1*vex(i)**2.-qcon*Y(i))
        T(i) = ((pre(i)*E_Tb-C_1*rho(i)**(4./3.))/C_2/rho(i)**(2./3.))**(1./2.)
        T(i) = (T(i)-T_0)/(T_b-T_0)   ! temperature -dimensionless
    end do load_pre

    !* cal aux and atx differential
    cal_ytx: do i=bx-ex/2,nx+ex/2
       ayx(i)=aa1*Y(i+3)+aa2*Y(i+2)+aa3*Y(i+1)&
           & +aa4*Y(i-1)+aa5*Y(i-2)+aa6*Y(i-3)
       atx(i)=aa1*T(i+3)+aa2*T(i+2)+aa3*T(i+1)&
           & +aa4*T(i-1)+aa5*T(i-2)+aa6*T(i-3)
    end do cal_ytx

    !* load data for ll the eq3. and eq4. need for conduction
    init_ll: do i=bx-ex/2, nx+ex/2
        ll(i,1:2)=0.0d0
        ll(i,3)=A_2*atx(i)/dx
        ll(i,4)=rho(i)/Le/C_t*ayx(i)/dx
        !* 2/r part of spherical coordinate
        if (ifcurv) then 
            ll(i,3)=ll(i,3)+A_2*T(i)*2./x(i)
            ll(i,4)=ll(i,4)+rho(i)/Le/C_t*Y(i)*2./x(i)
        end if
    end do init_ll

    !* cal fffd
    fffd_i: do i=bx, nx
        fffd_m: do m=1,eqN
            fffd(i,m)=aa1*ll(i+3,m)+aa2*ll(i+2,m)+aa3*ll(i+1,m)&
                   & +aa4*ll(i-1,m)+aa5*ll(i-2,m)+aa6*ll(i-3,m)
        end do fffd_m
    end do fffd_i
end subroutine vis_fx_L


!!** enox ----------------------------
subroutine enox
    include 'para1.f90'
    real(kind=acc),dimension(eqN)::am                   ! am-
    real(kind=acc),dimension(bxm:nxm)::w                ! w- weight
    real(kind=acc),dimension(bxm:nxm)::ev1,ev2,ev3,ev4  ! eigenvalues
    real(kind=acc),dimension(bxm:nxm,eqN,eqN)::evl,evr  ! eigen vectors
    real(kind=acc),dimension(bxm:nxm,eqN)::dF,dU        ! dF-    dU-
    real(kind=acc),dimension(bxm:nxm,eqN)::ff,fh        ! ff-    fh-
    real(kind=acc),dimension(bxm:nxm,eqN,2)::gg,hh      ! gg-    hh-
    real(kind=acc)::w0,w1        ! w0 w1 weight for interp
    real(kind=acc)::vex_m,c_m    ! vecm-vel in x-dir;  cm- sound vel
    real(kind=acc)::H_m,Y_m,E_m  ! middle H Y
    real(kind=acc)::t0,t1,t2,t3,s1,s2,s3
    include 'para2.f90'
    
    em=1.e-14
    
    !* load data
    do i = bxm, nxm
        rho(i) = U(i,1)
        vex(i) = U(i,2)/rho(i)*A_1
        E(i)   = U(i,3)/rho(i)
        Y(i)   = U(i,4)/rho(i)
    end do
    
    cal_ev: do i=bxm, nxm
        am(1:4)=ep
        
        ! pressure
        pre(i)=rho(i)*(E(i)-0.5/A_1*vex(i)**2.-qcon*Y(i))

        ! sound velocity
        cvel = 2./3.*sqrt(abs(A_1*pre(i)/rho(i)))
        
        ! load F
        F(i,1)=rho(i)*vex(i)                        ! Rho*v
        F(i,2)=rho(i)*vex(i)**2./A_1+pre(i)/3.      ! Rho*v**2+A_1*P/3
        F(i,3)=rho(i)*vex(i)*E(i)+vex(i)*pre(i)/3.  ! Rho*v*E+v*P
        F(i,4)=rho(i)*vex(i)*Y(i)                   ! Rho*v*Y

        w(i)=sqrt(abs(rho(i)))     ! weight 
        
        ! eigenvalues
        ev1(i)=abs(vex(i))
        ev2(i)=abs(vex(i))
        ev3(i)=abs(vex(i)-cvel)
        ev4(i)=abs(vex(i)+cvel)
    end do cal_ev

    cal_am: do i = bxm, nxm
        am(1) = max( am(1), ev1(i) )
        am(2) = max( am(2), ev2(i) )
        am(3) = max( am(3), ev3(i) )
        am(4) = max( am(4), ev4(i) )
    end do cal_am

    !* get em
    ! am(1) = am(1) * 1.1
    ! am(2) = am(2) * 1.1
    ! am(3) = am(3) * 1.1
    ! am(4) = am(4) * 1.1
    em = max(em,max(am(3),am(4)))

    !* cal evr & evl
    cal_evrl: do i=bx-1,nx
        !* Compute e'vectors using Roe's average:
        ip=i+1
        w0=w(i)/(w(i)+w(ip)+ep)
        w1=1.-w0
        
        rho_m= w0*rho(i) + w1*rho(ip)
        vex_m= w0*vex(i) + w1*vex(ip)
        E_m  = w0*E(i)   + w1*E(ip)
        Y_m  = w0*Y(i)   + w1*Y(ip)

        !* check the min
        if(abs(vex_m).le.ep) then
            vex_m = abs(vex_m)/vex_m*ep
        endif

        c_m  = sqrt(abs(2.*A_1*(E_m-0.5/A_1*vex_m**2.-qcon*Y_m)))

        !* evr: right e'vectors of Jaccobi Matrixs
        evr(i,1,1) = 2.*A_1*rho_m/vex_m
        evr(i,1,2) =-2.*A_1*qcon*rho_m/vex_m
        evr(i,1,3) = 1./Y_m
        evr(i,1,4) = 1./Y_m
        
        evr(i,2,1) = 2.*A_1*rho_m
        evr(i,2,2) =-2.*A_1*qcon*rho_m
        evr(i,2,3) = (3.*vex_m+2.**(1/2.)*c_m)/(3.*Y_m)
        evr(i,2,4) = (3.*vex_m-2.**(1/2.)*c_m)/(3.*Y_m)
        
        evr(i,3,1) = rho_m*vex_m
        evr(i,3,2) = 0.0d0
        evr(i,3,3) =-(vex_m**2.-8.*A_1*E_m+2.*A_1*Y_m*qcon-2.*2.**(1/2.)*vex_m*c_m)/(6.*A_1*Y_m)
        evr(i,3,4) =-(vex_m**2.-8.*A_1*E_m+2.*A_1*Y_m*qcon+2.*2.**(1/2.)*vex_m*c_m)/(6.*A_1*Y_m)
        
        evr(i,4,1) = 0.0d0
        evr(i,4,2) = rho_m*vex_m
        evr(i,4,3) = 1.0d0
        evr(i,4,4) = 1.0d0

        !* evl: inv(evr)
        evl(i,1,1) = vex_m*(7.*vex_m**2.-8.*A_1*E_m+14.*A_1*Y_m*qcon)/(8.*A_1*rho_m*(-c_m**2.))
        evl(i,1,2) =-(3.*(vex_m**2.+2.*A_1*Y_m*qcon))/(4.*A_1*rho_m*(-c_m**2.))
        evl(i,1,3) = (3.*(vex_m**2.+2.*A_1*Y_m*qcon))/(4.*rho_m*vex_m*(-c_m**2.))
        evl(i,1,4) = (qcon*(vex_m**2.-8.*A_1*E_m+2.*A_1*Y_m*qcon))/(4.*rho_m*vex_m*(-c_m**2.))
        
        evl(i,2,1) = (3.*Y_m*vex_m)/(4.*rho_m*(-c_m**2.))
        evl(i,2,2) =-(3.*Y_m)/(2.*rho_m*(-c_m**2.))
        evl(i,2,3) = (3.*A_1*Y_m)/(2.*rho_m*vex_m*(-c_m**2.))
        evl(i,2,4) = (2.*vex_m**2.-4.*A_1*E_m+A_1*Y_m*qcon)/(2.*rho_m*vex_m*(-c_m**2.))

        evl(i,3,1) = (3.*2.**(1/2.)*Y_m*vex_m*(4.*rho_m*(-c_m**2.)+2.**(1/2.)*rho_m*vex_m*c_m)) &
                    /(16.*rho_m*c_m**3.)
        evl(i,3,2) =-(3.*2.**(1/2.)*Y_m*(-2.*c_m**2+2.**(1/2.)*vex_m*c_m)) &
                    /(8.*c_m**3.)
        evl(i,3,3) =-(3.*A_1*Y_m)/(-4.*c_m**2.)
        evl(i,3,4) = (3.*A_1*Y_m*qcon)/(-4.*c_m**2.)

        evl(i,4,1) =-(3.*2.**(1/2.)*Y_m*vex_m*(4.*rho_m*(-c_m**2.)-2.**(1/2.)*rho_m*vex_m*c_m)) &
                    /(16.*rho_m*c_m**3.)
        evl(i,4,2) = (3.*2.**(1/2.)*Y_m*(-2.*c_m**2-2.**(1/2.)*vex_m*c_m)) &
                    /(8.*c_m**3.)
        evl(i,4,3) =-(3.*A_1*Y_m)/(-4.*c_m**2.)
        evl(i,4,4) = (3.*A_1*Y_m*qcon)/(-4.*c_m**2.)
        
        !* check eigenvalues
        if (ifchec) then
            write(*,*) "rho_m:",rho_m,"vex_m:",vex_m,"c_m:",c_m
            write(*,*) "Y_m:",Y_m,"E_m:",E_m,"A_1:",A_1
            write(*,*) (3.*(vex_m**2.+2.*A_1*Y_m*qcon))
            write(*,*) (4.*rho_m**2.*vex_m**2.*(-c_m**2.))
            do m = 1,4
            do k = 1,4
                write(*,'("row:",I4,"  col:",I4,"  evr:",1e15.7, "  evl:",e15.7)') m,k,evr(i,m,k),evl(i,m,k)
            enddo
            enddo

            do m = 1,4
            do k = 1,4
                temp = 0.0d0
                do j = 1,4
                    temp = temp + evl(i,m,j)*evr(i,j,k)
                enddo
                if(m.eq.k) then
                    err = temp-1.
                else
                    err = temp
                endif
                write(*,'("row:",I4,"  col:",I4,"  I:",1e15.7,"  err:",1e15.7)') m,k,temp,err
            enddo
            enddo
            stop
        endif
    end do cal_evrl
    
    !* cal dF & dU
    cal_dFU: do i=bxm,nxm-1
        dF(i,:)=F(i+1,:)-F(i,:)  ! derta F
        dU(i,:)=U(i+1,:)-U(i,:)  ! derta U
    end do cal_dFU

    !* init ff
    ff(bxm:nxm,1:eqN)=0.0d0

    cal_ff: do m=1,eqN
    !* Project the relevant first undivided differences into the
    !* 'm'th characteristic field
        !* cal gg
        gg(bxm:nxm-1,:,1)=0.5*(am(m)*dU(bxm:nxm-1,:)+dF(bxm:nxm-1,:))
        gg(bxm:nxm-1,:,2)=0.5*(am(m)*dU(bxm:nxm-1,:)-dF(bxm:nxm-1,:))

        !* cal hh
        hh_m1: do m1=1,4
            k0=m1-3
            k1=3-m1
            hh_i: do i=bx-1,nx
                hh(i,m1,1)=evl(i,m,1)*gg(i+k0,1,1)+evl(i,m,2)*gg(i+k0,2,1)&
                          +evl(i,m,3)*gg(i+k0,3,1)+evl(i,m,4)*gg(i+k0,4,1)
                hh(i,m1,2)=evl(i,m,1)*gg(i+k1,1,2)+evl(i,m,2)*gg(i+k1,2,2)&
                          +evl(i,m,3)*gg(i+k1,3,2)+evl(i,m,4)*gg(i+k1,4,2)
            end do hh_i
        end do hh_m1

        !* Compute numerical flux in each characteristic field:
        cal_s: do m1=1,2
            cal_t: do i=bx-1,nx
                t1=hh(i,1,m1)-hh(i,2,m1)
                t2=hh(i,2,m1)-hh(i,3,m1)
                t3=hh(i,3,m1)-hh(i,4,m1)

                t0=hh(i,1,m1)-3.*hh(i,2,m1)
                s1=1./(ep+13.*t1*t1+3.*t0*t0)**2.

                t0=hh(i,2,m1)+hh(i,3,m1)
                s2=6./(ep+13.*t2*t2+3.*t0*t0)**2.

                t0=3*hh(i,3,m1)-hh(i,4,m1)
                s3=3./(ep+13.*t3*t3+3.*t0*t0)**2.

                t0=1./(s1+s2+s3)
                s1=s1*t0
                s3=s3*t0

                ff(i,m)=ff(i,m)+(s1*(t2-t1)+(0.5*s3-0.25)*(t3-t2))/3.
            end do cal_t
        end do cal_s
    end do cal_ff
    !*  end of the big loop

    !* Project the numerical flux to the physical space:
    fh_m: do m=1,eqN 
        fh_i: do i=bx-1,nx
            fh(i,m)=(-F(i-1,m)+7.*(F(i,m)+F(i+1,m))-F(i+2,m))/12.&
                   +evr(i,m,1)*ff(i,1)+evr(i,m,2)*ff(i,2)&
                   +evr(i,m,3)*ff(i,3)+evr(i,m,4)*ff(i,4)
        end do fh_i
    end do fh_m

    cal_f: do i=bx,nx
        LU(i,1)=(fh(i-1,1)-fh(i,1))/dx
        LU(i,2)=(fh(i-1,2)-fh(i,2))/dx
        LU(i,3)=(fh(i-1,3)-fh(i,3))/dx
        LU(i,4)=(fh(i-1,4)-fh(i,4))/dx
        
        !* diffusion part
        if (ifdiff) then
            LU(i,3)=LU(i,3)+0.25d0*fffd(i,3)/dx
            LU(i,4)=LU(i,4)+0.25d0*fffd(i,4)/dx
        end if

        !* 2/r part of spherical coordinate
        if (ifcurv) then
            LU(i,1)=LU(i,1) - 2.*F(i,1)/x(i)
            LU(i,2)=LU(i,2) - 2.*(F(i,2)-pre(i)/3.)/x(i)
            LU(i,3)=LU(i,3) - 2.*F(i,3)/x(i)
            LU(i,4)=LU(i,4) - 2.*F(i,4)/x(i)
        end if
    end do cal_f
return
end subroutine enox


!!** enos ---------------------------------
subroutine enos
    include 'para1.f90'
    real(kind=acc)::theta,s,t0
    include 'para2.f90'
    
    !* load data
    do i=bxm,nxm
        rho(i) = U(i,1)
        vex(i) = U(i,2)/rho(i)*A_1
        E(i)   = U(i,3)/rho(i)
        Y(i)   = U(i,4)/rho(i)
    end do
    
    ! cal reaction
    do i=bx,nx
        ! pressure    - from rho & E & vex & Y
        pre(i)=rho(i)*(E(i)-0.5/A_1*vex(i)**2.-qcon*Y(i))
        ! temperature - from rho & pre
        T(i) = ((pre(i)*E_Tb-C_1*rho(i)**(4./3.))/C_2/rho(i)**(2./3.))**(1./2.)
        omega(i)= Ace*Y(i)*rho(i)*exp(-Ea/((T(i)/1.e9)**(1./3.))+Ea/((T_b/1.e9)**(1./3.)))
        ! check if NaN arouse, and output wrong message
        if (isnan(omega(i))) then
            write(*,*) "NaN!!===================================|============================"
            write(*,*) "Old!! i:",i,        "              |omega:",omega(i)
            write(*,*) "      T:",T(i),                   "|  rho:",rho(i)
            write(*,*) "    pre:",pre(i),                 "|  vex:",vex(i)
            write(*,*) "      Y:",Y(i),                   "|    E:",E(i)
            write(*,*) "  rho*E:",rho(i)*E(i),            "| U(3):",U(i,3)
            write(*,*) "-----------------------------------|--------------------------------"
            temp = (pre(i)*E_Tb/E_0-rho(i)**(4./3.))*E_0/(E_T0-E_0)/rho(i)**(2./3.)*T_0*T_0
            write(*,*) "T(i)^2:",temp
            T(i) = 2.*T(i-1)-T(i-2)
            temp = E(i)-0.5*(vex(i)**2.)/A_1-qcon*Y(i)
            rho(i) = (temp*E_Tb/2./E_0+sqrt((temp*E_Tb/2./E_0)**2.-(E_T0/E_0-1.)*(T(i)/T_0)**2.))**3
            omega(i)= Y(i)*rho(i)*exp(-Ea/((T(i)/1.e9)**(1./3.))+Ea/((T_b/1.e9)**(1./3.)))
            write(*,*) "Fix!! i:",i,        "              |omega:",omega(i)
            write(*,*) "      T:",T(i),                   "|  rho:",rho(i)
            write(*,*) "===================================|============================NaN!!"
            s = s+1
            ifsave = .true.
            call res_data
            stop
        endif
        if (ifreac) then
            LU(i,4)=LU(i,4)-omega(i)/C_t
        end if
    end do
return
end subroutine enos

!!**  Output: check Nan & min val =======
real(kind=8) function output(input)
    implicit none
    real(kind=8)::input
    output = input
    if (abs(input).le.(.1e-36)) then
        output = 0.0d0
    else if (input > .1e30) then
        write(*,*) "Large Number arouse!"
        output = .1e30
    else if (isnan(input)) then
        write(*,*) "NaN arouse!"
        output = 0.0d0
    endif
end function output

real(kind=16) function output16(input)
    implicit none
    real(kind=16)::input
    if (input.le.(.1e-36)) then
        output16 = 0.0d0
    else if (input > .1e30) then
        write(*,*) "Large Number arouse!"
        output16 = .1e30
    else if (isnan(input)) then
        write(*,*) "NaN arouse!"
        output16 = 0.0d0
    endif
end function output16