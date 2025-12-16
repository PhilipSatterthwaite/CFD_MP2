module compressible
    use precision
    implicit none
    real(DP), dimension(:,:), pointer :: dens
    real(DP), dimension(:,:), pointer :: xvel
    real(DP), dimension(:,:), pointer :: yvel
    real(DP), dimension(:,:), pointer :: ener
    real(DP), dimension(:,:), pointer :: pressure
    real(DP), dimension(:,:), pointer :: temperature

    real(DP), parameter :: wall_temp = 300.0_DP
    real(DP), parameter :: gamma = 1.4_DP
    real(DP), parameter :: R = 287.0_DP !J/kgK
    real(DP), parameter :: Re = 100.0_DP
    real(DP), parameter :: Pr = 0.7_DP

    real(DP), parameter :: pressure_init = 100000.0_DP !pascals
    real(DP), parameter :: temp_init = 300.0_DP !K

    real(DP) :: cp, omega, a_0, vel_wall
    real(DP) :: nu, lam_rhocp

    real(DP) :: delx, dely, delt
    real(DP) :: time

    ! Inputs
    integer :: NP
    real(DP) :: Ma
    real(DP) :: len

    contains

    subroutine compressible_initialize(NP_in, Ma_in, len_in)
        use precision
        implicit none
        integer, intent(in) :: NP_in
        real(DP), intent(in) :: Ma_in
        real(DP), intent(in) :: len_in

        real(DP), dimension(0:NP+1,0:NP+1):: u_vec, u_new

        NP = NP_in
        Ma = Ma_in
        len = len_in

        ! Allocate and initialize variables
        allocate(pressure(0:NP+1,0:NP+1))
        pressure(:,:) = pressure_init
        allocate(temperature(0:NP+1,0:NP+1))
        temperature(:,:) = temp_init
        allocate(dens(0:NP+1,0:NP+1))
        dens(:,:) = pressure_init/(R*temp_init)
        allocate(xvel(0:NP+1,0:NP+1))
        xvel(:,:) = 0.0_DP
        allocate(yvel(0:NP+1,0:NP+1))
        yvel(:,:) = 0.0_DP
        allocate(ener(0:NP+1,0:NP+1))
        ener(:,:) = R*temp_init/(gamma - 1.0_DP)

        ! Compute global constants
        a_0 = sqrt(gamma*R*temp_init)
        vel_wall = Ma*a_0
        nu = 0.025_DP*a_0/100.0_DP
        !print *, nu
        !stop
        lam_rhocp = nu/0.7_DP
        cp = R*gamma/(gamma - 1.0_DP)
        omega = 2.0_DP/len**2*nu


        ! Discretizations
        delx = len/real(NP, DP)
        dely = len/real(NP, DP)

        time = 0.0_DP

    end subroutine compressible_initialize


    subroutine compressible_run(delt, t_end)
        use precision
        implicit none
        real(DP), intent(in) :: delt
        real(DP), intent(in) :: t_end

        integer :: NT, time_int
        integer :: p,q

        real(DP) :: CFL

        CFL = (abs(vel_wall) + abs(a_0))/delx*delt
        print *, "CFL", CFL
        if (CFL.ge.1.0_DP) then
            print *, "CFL too large"
            stop
        end if

        !t_end = 10.0_DP*omega

        if (t_end.gt.0.0_DP) then
            NT = int(t_end/delt,8)
        else
            NT = int(10.0_DP*omega/delt,8)
        end if

        do time_int = 1,(NT-1)
            time = time + delt
            print *, "time: ", time, 10.0_DP*omega, CFL
            ! Update parameters
            do q = 1,NP
                do p = 1,NP
                    temperature(p,q) = ener(p,q)*(gamma - 1.0_DP)/R
                    pressure(p,q) = dens(p,q)*R*temperature(p,q)
                end do
            end do


            ! Update Boundary Conditions
            temperature(0,1:NP) = 2.0_DP*wall_temp - temperature(1,1:NP)
            temperature(NP+1,1:NP) = 2.0_DP*wall_temp - temperature(NP,1:NP)
            temperature(1:NP,0) = 2.0_DP*wall_temp - temperature(1:NP,1)
            temperature(1:NP,NP+1) = 2.0_DP*wall_temp - temperature(1:NP,NP)

            pressure(0,1:NP) = pressure(1,1:NP)
            pressure(NP+1,1:NP) = pressure(NP,1:NP)
            pressure(1:NP,0) = pressure(1:NP,1)
            pressure(1:NP,NP+1) = pressure(1:NP,NP)

            dens(0,1:NP) = pressure(0,1:NP)/(R*temperature(0,1:NP))
            dens(NP+1,1:NP) = pressure(NP+1,1:NP)/(R*temperature(NP+1,1:NP))
            dens(1:NP,0) = pressure(1:NP,0)/(R*temperature(1:NP,0))
            dens(1:NP,NP+1) = pressure(1:NP,NP+1)/(R*temperature(1:NP,NP+1))

            xvel(0,1:NP)  = - xvel(1,1:NP)
            xvel(NP+1,1:NP) = - xvel(NP,1:NP)
            xvel(1:NP,0) = 2.0_DP*vel_wall*sin(omega*time) - xvel(1:NP,1)
            xvel(1:NP,NP+1) = - xvel(1:NP,NP)

            yvel(0,1:NP) = - yvel(1,1:NP)
            yvel(NP+1,1:NP) = - yvel(NP,1:NP)
            yvel(1:NP,0) = - yvel(1:NP,1)
            yvel(1:NP,NP+1) = - yvel(1:NP,NP)

            ener(0,1:NP) = temperature(0,1:NP)*R/(gamma - 1.0_DP)
            ener(NP+1,1:NP) = temperature(NP+1,1:NP)*R/(gamma - 1.0_DP)
            ener(1:NP,0) = temperature(1:NP,0)*R/(gamma - 1.0_DP)
            ener(1:NP,NP+1) = temperature(1:NP,NP+1)*R/(gamma - 1.0_DP)

            ! Corner velocities

            ! Bottom right
            ! parallel to wall, du/dx = 0, dv/dy = 0
            !xvel(NP+1,NP+1) - xvel(NP-1,NP+1) + xvel(NP+1,NP) - xvel(NP-1,NP) = 0
            xvel(NP+1,NP+1) =  xvel(NP-1,NP+1) - xvel(NP+1,NP) + xvel(NP-1,NP)
            !yvel(NP+1,NP+1) - yvel(NP+1,NP-1) + yvel(NP,NP+1) - yvel(NP,NP-1) = 0
            yvel(NP+1,NP+1) = yvel(NP+1,NP-1) - yvel(NP,NP+1) + yvel(NP,NP-1)

            ! Bottom left
            !xvel(2,NP+1) - xvel(0,NP+1) + xvel(2,NP) - xvel(0,NP) = 0
            xvel(0,NP+1) = xvel(2,NP+1) + xvel(2,NP) - xvel(0,NP)
            !yvel(1,NP+1) - yvel(1,NP-1) + yvel(0,NP+1) - yvel(0,NP-1) = 0
            yvel(0,NP+1) = - yvel(1,NP+1) + yvel(1,NP-1) + yvel(0,NP-1)

            ! Top left
            !xvel(2,0) - xvel(0,0) + xvel(2,1) - xvel(0,1) = 0
            xvel(0,0) = xvel(2,1) - xvel(0,1) + xvel(2,0)
            !yvel(0,2) - yvel(0,0) + yvel(1,2) - yvel(1,0) = 0
            yvel(0,0) = yvel(1,2) - yvel(1,0) + yvel(0,2)

            ! Top right
            !xvel(NP+1,0) - xvel(NP-1,0) + xvel(NP+1,1) - xvel(NP-1,1) = 0
            xvel(NP+1,0) = - xvel(NP+1,1) + xvel(NP-1,1) + xvel(NP-1,0)
            !yvel(NP+1,2) - yvel(NP+1,0) + yvel(NP,2) - yvel(NP,0) = 0
            yvel(NP+1,0) = yvel(NP+1,2) + yvel(NP,2) - yvel(NP,0)



            
            !print *,xvel(1,1:NP)
            !print *, xvel(0,1:NP)

            call compressible_updatevecs()
            
            if (.false.) then
                print *, "dens"
                do p = 0,NP+1
                    print *, dens(:,p)
                end do
                print *, "X"
                do p = 0,NP+1
                    print *, xvel(:,p)
                end do
                print *, "Y"
                do p = 0,NP+1
                    print *, yvel(:,p)
                end do
                print *, "ener"
                do p = 0,NP+1
                    print *, ener(:,p)
                end do
                print *, "temp"
                do p = 0,NP+1
                    print *, temperature(:,p)
                end do
                print *, "press"
                do p = 0,NP+1
                    print *, pressure(:,p)
                end do
            end if
            !stop

            !print *, yvel(1:NP,1)
            !print *, ener(1:NP,1)
            !if (time_int.eq.20) stop
            if (any(dens.ne.dens)) stop
            if (any(xvel.ne.xvel)) stop
            if (any(yvel.ne.yvel)) stop
            if (any(ener.ne.ener)) stop

        end do


        return
        contains 
        subroutine compressible_updatevecs()
            implicit none

            real(DP), dimension(4,NP*NP) :: rhs

            integer :: p,q
            real(DP) :: usign, vsign
            real(DP) :: diff_minus, diff_plus, dissip
            real(DP) :: Flux_minus, Flux_plus
             
            real(DP), dimension(NP,NP) :: dens_new, xvel_new, yvel_new, ener_new

            rhs = 0.0_DP

            !$OMP PARALLEL DO PRIVATE(p,q,usign,vsign,Flux_plus,Flux_minus,diff_plus,diff_minus,dissip)
            do q = 1,NP
                do p = 1,NP
                    ! Continuity
                    !if (abs(xvel(p,q)) > 0.0_DP) then
                    !    usign = xvel(p,q)/abs(xvel(p,q))
                    !else
                    !    usign = 0.0_DP
                    !end if
                    !if (abs(yvel(p,q)) > 0.0_DP) then
                    !    vsign = yvel(p,q)/abs(yvel(p,q))
                    !else
                    !    vsign = 0.0_DP
                    !end if
                    !rhs(1,(p-1)*NP + q) = (dens_new(p,q) - dens(p,q))/delt
                    !rhs(1,(p-1)*NP + q) = rhs(1,(p-1)*NP + q) &
                    !    + max(usign,0.0_DP)*(dens(p,q)*xvel(p,q) - dens(p-1,q)*xvel(p-1,q))/(delx) &
                    !    - min(usign,0.0_DP)*(dens(p+1,q)*xvel(p+1,q) - dens(p,q)*xvel(p,q))/(delx) &
                    !    + max(vsign,0.0_DP)*(dens(p,q)*yvel(p,q) - dens(p,q-1)*yvel(p,q-1))/(delx) &
                    !    - min(vsign,0.0_DP)*(dens(p,q+1)*yvel(p,q+1) - dens(p,q)*yvel(p,q))/(delx)
                    !Flux_plus = (dens(p,q)*xvel(p,q) + dens(p+1,q)*xvel(p+1,q))/2.0_DP
                    !Flux_minus = (dens(p,q)*xvel(p,q) + dens(p-1,q)*xvel(p-1,q))/2.0_DP
                    !rhs(1,(p-1)*NP + q) = rhs(1,(p-1)*NP + q) + (Flux_plus - Flux_minus)/delx
                    !Flux_plus = (dens(p,q)*yvel(p,q) + dens(p,q+1)*yvel(p,q+1))/2.0_DP
                    !Flux_minus = (dens(p,q)*yvel(p,q) + dens(p,q-1)*yvel(p,q-1))/2.0_DP
                    !rhs(1,(p-1)*NP + q) = rhs(1,(p-1)*NP + q) + (Flux_plus - Flux_minus)/delx
                    usign = (xvel(p,q)+xvel(p+1,q)+1.0e-10_DP)/abs((xvel(p,q)+xvel(p+1,q))+1.0e-10_DP)
                    Flux_plus = max(usign,0.0_DP)*(dens(p,q)*0.5_DP*(xvel(p,q)+xvel(p+1,q))) &
                        - min(usign,0.0_DP)*(dens(p+1,q)*0.5_DP*(xvel(p,q)+xvel(p+1,q)))
                    usign = (xvel(p,q)+xvel(p-1,q)+1.0e-10_DP)/abs((xvel(p,q)+xvel(p-1,q))+1.0e-10_DP)
                    Flux_minus = max(usign,0.0_DP)*(dens(p-1,q)*0.5_DP*(xvel(p,q)+xvel(p-1,q))) &
                        - min(usign,0.0_DP)*(dens(p,q)*0.5_DP*(xvel(p,q)+xvel(p-1,q)))
                    rhs(1,(p-1)*NP + q) = rhs(1,(p-1)*NP + q) - (Flux_plus - Flux_minus) / delx
                    
                    vsign = (yvel(p,q)+yvel(p,q+1)+1.0e-10_DP)/abs((yvel(p,q)+yvel(p,q+1))+1.0e-10_DP)
                    Flux_plus = max(vsign,0.0_DP)*(dens(p,q)*0.5_DP*(yvel(p,q)+yvel(p,q+1))) &
                        - min(vsign,0.0_DP)*(dens(p,q+1)*0.5_DP*(yvel(p,q)+yvel(p,q+1)))
                    vsign = (yvel(p,q)+yvel(p,q-1)+1.0e-10_DP)/abs((yvel(p,q)+yvel(p,q-1))+1.0e-10_DP)
                    Flux_minus = max(vsign,0.0_DP)*(dens(p,q-1)*0.5_DP*(yvel(p,q)+yvel(p,q-1))) &
                        - min(vsign,0.0_DP)*(dens(p,q)*0.5_DP*(yvel(p,q)+yvel(p,q-1)))
                    rhs(1,(p-1)*NP + q) = rhs(1,(p-1)*NP + q) - (Flux_plus - Flux_minus) / dely

                    ! X Momentum
                    !rhs(2,(p-1)*NP + q) = (dens_new(p,q)*xvel_new(p,q) - dens(p,q)*xvel(p,q))/delt
                    !rhs(2,(p-1)*NP + q) = rhs(2,(p-1)*NP + q) &
                    !    + max(usign,0.0_DP)*(dens(p,q)*xvel(p,q)*xvel(p,q) - dens(p-1,q)*xvel(p-1,q)*xvel(p-1,q))/(delx) &
                    !    - min(usign,0.0_DP)*(dens(p+1,q)*xvel(p+1,q)*xvel(p+1,q) - dens(p,q)*xvel(p,q)*xvel(p,q))/(delx) &
                    !    + max(vsign,0.0_DP)*(dens(p,q)*yvel(p,q)*xvel(p,q) - dens(p,q-1)*yvel(p,q-1)*xvel(p,q-1))/(delx) &
                    !    - min(vsign,0.0_DP)*(dens(p,q+1)*yvel(p,q+1)*xvel(p,q+1) - dens(p,q)*yvel(p,q)*xvel(p,q))/(delx)

                    !Flux_plus = (dens(p,q)*xvel(p,q)*xvel(p,q) + dens(p+1,q)*xvel(p+1,q)*xvel(p+1,q))/2.0_DP
                    !Flux_minus = (dens(p,q)*xvel(p,q)*xvel(p,q) + dens(p-1,q)*xvel(p-1,q)*xvel(p-1,q))/2.0_DP
                    !rhs(2,(p-1)*NP + q) = rhs(2,(p-1)*NP + q) + (Flux_plus - Flux_minus)/delx
                    !Flux_plus = (dens(p,q)*xvel(p,q)*yvel(p,q) + dens(p,q+1)*xvel(p,q+1)*yvel(p,q+1))/2.0_DP
                    !Flux_minus = (dens(p,q)*xvel(p,q)*yvel(p,q) + dens(p,q-1)*xvel(p,q-1)*yvel(p,q-1))/2.0_DP
                    !rhs(2,(p-1)*NP + q) = rhs(2,(p-1)*NP + q) + (Flux_plus - Flux_minus)/delx

                    usign = (xvel(p,q)+xvel(p+1,q)-1.0e-10_DP)/abs((xvel(p,q)+xvel(p+1,q))-1.0e-10_DP)
                    Flux_plus = max(usign,0.0_DP)*(dens(p,q)*xvel(p,q)*0.5_DP*(xvel(p,q)+xvel(p+1,q))) &
                        - min(usign,0.0_DP)*(dens(p+1,q)*xvel(p+1,q)*0.5_DP*(xvel(p,q)+xvel(p+1,q)))
                    usign = (xvel(p,q)+xvel(p-1,q)-1.0e-10_DP)/abs((xvel(p,q)+xvel(p-1,q))-1.0e-10_DP)
                    Flux_minus = max(usign,0.0_DP)*(dens(p-1,q)*xvel(p-1,q)*0.5_DP*(xvel(p,q)+xvel(p-1,q))) &
                        - min(usign,0.0_DP)*(dens(p,q)*xvel(p,q)*0.5_DP*(xvel(p,q)+xvel(p-1,q)))
                    rhs(2,(p-1)*NP + q) = rhs(2,(p-1)*NP + q) - (Flux_plus - Flux_minus) / delx
                    
                    vsign = (yvel(p,q)+yvel(p,q+1)-1.0e-10_DP)/abs((yvel(p,q)+yvel(p,q+1))-1.0e-10_DP)
                    Flux_plus = max(vsign,0.0_DP)*(dens(p,q)*xvel(p,q)*0.5_DP*(yvel(p,q)+yvel(p,q+1))) &
                        - min(vsign,0.0_DP)*(dens(p,q+1)*xvel(p,q+1)*0.5_DP*(yvel(p,q)+yvel(p,q+1)))
                    vsign = (yvel(p,q)+yvel(p,q-1)-1.0e-10_DP)/abs((yvel(p,q)+yvel(p,q-1))-1.0e-10_DP)
                    Flux_minus = max(vsign,0.0_DP)*(dens(p,q-1)*xvel(p,q-1)*0.5_DP*(yvel(p,q)+yvel(p,q-1))) &
                        - min(vsign,0.0_DP)*(dens(p,q)*xvel(p,q)*0.5_DP*(yvel(p,q)+yvel(p,q-1)))
                    rhs(2,(p-1)*NP + q) = rhs(2,(p-1)*NP + q) - (Flux_plus - Flux_minus) / dely
                    
                    rhs(2,(p-1)*NP + q) = rhs(2,(p-1)*NP + q) - (pressure(p+1,q) - pressure(p-1,q))/(2.0_DP*delx)

                    diff_plus = 2.0_DP*(xvel(p+1,q) - xvel(p,q))/delx &
                        - 2.0_DP/3.0_DP*((xvel(p+1,q) - xvel(p,q))/delx &
                        + ((yvel(p+1,q+1) - yvel(p+1,q-1))/(2.0_DP*dely) + (yvel(p,q+1) - yvel(p,q-1))/(2.0_DP*dely))/2.0_DP)
                    diff_plus = nu*(dens(p,q) + dens(p+1,q))/2.0_DP*diff_plus
                    diff_minus = 2.0_DP*(xvel(p,q) - xvel(p-1,q))/delx &
                        - 2.0_DP/3.0_DP*((xvel(p,q) - xvel(p-1,q))/delx &
                        + ((yvel(p-1,q+1) - yvel(p-1,q-1))/(2.0_DP*dely) + (yvel(p,q+1) - yvel(p,q-1))/(2.0_DP*dely))/2.0_DP)
                    diff_minus = nu*(dens(p-1,q) + dens(p,q))/2.0_DP*diff_minus
                    rhs(2,(p-1)*NP + q) = rhs(2,(p-1)*NP + q) + (diff_plus - diff_minus)/delx

                    diff_plus = (xvel(p,q+1) - xvel(p,q))/dely &
                        + ((yvel(p+1,q) - yvel(p-1,q))/(2.0_DP*delx) + (yvel(p+1,q+1) - yvel(p-1,q+1))/(2.0_DP*delx))/2.0_DP
                    diff_plus = nu*(dens(p,q) + dens(p,q+1))/2.0_DP*diff_plus
                    diff_minus = (xvel(p,q) - xvel(p,q-1))/dely &
                        + ((yvel(p+1,q) - yvel(p-1,q))/(2.0_DP*delx) + (yvel(p+1,q-1) - yvel(p-1,q-1))/(2.0_DP*delx))/2.0_DP
                    diff_minus = nu*(dens(p,q-1) + dens(p,q))/2.0_DP*diff_minus
                    rhs(2,(p-1)*NP + q) = rhs(2,(p-1)*NP + q) + (diff_plus - diff_minus)/dely


                    ! Y Momentum
                    !rhs(3,(p-1)*NP + q) = (dens_new(p,q)*yvel_new(p,q) - dens(p,q)*yvel(p,q))/delt
                    !rhs(3,(p-1)*NP + q) = rhs(3,(p-1)*NP + q) &
                    !    + max(usign,0.0_DP)*(dens(p,q)*xvel(p,q)*yvel(p,q) - dens(p-1,q)*xvel(p-1,q)*yvel(p-1,q))/(delx) &
                    !    - min(usign,0.0_DP)*(dens(p+1,q)*xvel(p+1,q)*yvel(p+1,q) - dens(p,q)*xvel(p,q)*yvel(p,q))/(delx) &
                    !    + max(vsign,0.0_DP)*(dens(p,q)*yvel(p,q)*yvel(p,q) - dens(p,q-1)*yvel(p,q-1)*yvel(p,q-1))/(delx) &
                    !    - min(vsign,0.0_DP)*(dens(p,q+1)*yvel(p,q+1)*yvel(p,q+1) - dens(p,q)*yvel(p,q)*yvel(p,q))/(delx)
                    !Flux_plus = (dens(p,q)*yvel(p,q)*xvel(p,q) + dens(p+1,q)*yvel(p+1,q)*xvel(p+1,q))/2.0_DP
                    !Flux_minus = (dens(p,q)*yvel(p,q)*xvel(p,q) + dens(p-1,q)*yvel(p-1,q)*xvel(p-1,q))/2.0_DP
                    !rhs(3,(p-1)*NP + q) = rhs(3,(p-1)*NP + q) + (Flux_plus - Flux_minus)/delx
                    !Flux_plus = (dens(p,q)*yvel(p,q)*yvel(p,q) + dens(p,q+1)*yvel(p,q+1)*yvel(p,q+1))/2.0_DP
                    !Flux_minus = (dens(p,q)*yvel(p,q)*yvel(p,q) + dens(p,q-1)*yvel(p,q-1)*yvel(p,q-1))/2.0_DP
                    !rhs(3,(p-1)*NP + q) = rhs(3,(p-1)*NP + q) + (Flux_plus - Flux_minus)/delx

                    usign = (xvel(p,q)+xvel(p+1,q)+1.0e-10_DP)/abs((xvel(p,q)+xvel(p+1,q))+1.0e-10_DP)
                    Flux_plus = max(usign,0.0_DP)*(dens(p,q)*yvel(p,q)*0.5_DP*(xvel(p,q)+xvel(p+1,q))) &
                        - min(usign,0.0_DP)*(dens(p+1,q)*yvel(p+1,q)*0.5_DP*(xvel(p,q)+xvel(p+1,q)))
                    usign = (xvel(p,q)+xvel(p-1,q)+1.0e-10_DP)/abs((xvel(p,q)+xvel(p-1,q))+1.0e-10_DP)
                    Flux_minus = max(usign,0.0_DP)*(dens(p-1,q)*yvel(p-1,q)*0.5_DP*(xvel(p,q)+xvel(p-1,q))) &
                        - min(usign,0.0_DP)*(dens(p,q)*yvel(p,q)*0.5_DP*(xvel(p,q)+xvel(p-1,q)))
                    rhs(3,(p-1)*NP + q) = rhs(3,(p-1)*NP + q) - (Flux_plus - Flux_minus) / delx
                    
                    vsign = (yvel(p,q)+yvel(p,q+1)+1.0e-10_DP)/abs((yvel(p,q)+yvel(p,q+1))+1.0e-10_DP)
                    Flux_plus = max(vsign,0.0_DP)*(dens(p,q)*yvel(p,q)*0.5_DP*(yvel(p,q)+yvel(p,q+1))) &
                        - min(vsign,0.0_DP)*(dens(p,q+1)*yvel(p,q+1)*0.5_DP*(yvel(p,q)+yvel(p,q+1)))
                    vsign = (yvel(p,q)+yvel(p,q-1)+1.0e-10_DP)/abs((yvel(p,q)+yvel(p,q-1))+1.0e-10_DP)
                    Flux_minus = max(vsign,0.0_DP)*(dens(p,q-1)*yvel(p,q-1)*0.5_DP*(yvel(p,q)+yvel(p,q-1))) &
                        - min(vsign,0.0_DP)*(dens(p,q)*yvel(p,q)*0.5_DP*(yvel(p,q)+yvel(p,q-1)))
                    rhs(3,(p-1)*NP + q) = rhs(3,(p-1)*NP + q) - (Flux_plus - Flux_minus) / dely

                    rhs(3,(p-1)*NP + q) = rhs(3,(p-1)*NP + q) - (pressure(p,q+1) - pressure(p,q-1))/(2.0_DP*dely)

                    diff_plus = (yvel(p,q+1) - yvel(p,q))/delx &
                     + ((xvel(p,q+1) - xvel(p,q-1))/(2.0_DP*dely) + (xvel(p+1,q+1) - xvel(p+1,q-1))/(2.0_DP*dely))/2.0_DP
                    diff_plus = nu*(dens(p,q) + dens(p+1,q))/2.0_DP*diff_plus
                    diff_minus = (yvel(p,q) - yvel(p,q-1))/delx &
                     + ((xvel(p,q+1) - xvel(p,q-1))/(2.0_DP*dely) + (xvel(p-1,q+1) - xvel(p-1,q-1))/(2.0_DP*dely))/2.0_DP
                    diff_minus = nu*(dens(p,q) + dens(p-1,q))/2.0_DP*diff_minus
                    rhs(3,(p-1)*NP + q) = rhs(3,(p-1)*NP + q) + (diff_plus - diff_minus)/delx

                    diff_plus = 2.0_DP*(yvel(p,q+1) - yvel(p,q))/dely &
                        - 2.0_DP/3.0_DP*((yvel(p,q+1) - yvel(p,q))/dely &
                        + ((xvel(p+1,q+1) - xvel(p-1,q+1))/(2.0_DP*delx) + (xvel(p+1,q) - xvel(p-1,q))/(2.0_DP*delx))/2.0_DP)
                    diff_plus = nu*(dens(p,q) + dens(p,q+1))/2.0_DP*diff_plus
                    diff_minus = 2.0_DP*(yvel(p,q) - yvel(p,q-1))/dely &
                        - 2.0_DP/3.0_DP*((yvel(p,q) - yvel(p,q-1))/dely &
                        + ((xvel(p+1,q-1) - xvel(p-1,q-1))/(2.0_DP*delx) + (xvel(p+1,q) - xvel(p-1,q))/(2.0_DP*delx))/2.0_DP)
                    diff_minus = nu*(dens(p,q) + dens(p,q-1))/2.0_DP*diff_minus
                    rhs(3,(p-1)*NP + q) = rhs(3,(p-1)*NP + q) + (diff_plus - diff_minus)/dely

                    ! Energy Equation
                    !rhs(4,(p-1)*NP + q) = (dens_new(p,q)*ener_new(p,q) - dens(p,q)*ener(p,q))/delt
                    !rhs(4,(p-1)*NP + q) = rhs(4,(p-1)*NP + q) &
                    !    + max(usign,0.0_DP)*(dens(p,q)*xvel(p,q)*ener(p,q) - dens(p-1,q)*xvel(p-1,q)*ener(p-1,q))/(delx) &
                    !    - min(usign,0.0_DP)*(dens(p+1,q)*xvel(p+1,q)*ener(p+1,q) - dens(p,q)*xvel(p,q)*ener(p,q))/(delx) &
                    !    + max(vsign,0.0_DP)*(dens(p,q)*yvel(p,q)*ener(p,q) - dens(p,q-1)*yvel(p,q-1)*ener(p,q-1))/(delx) &
                    !    - min(vsign,0.0_DP)*(dens(p,q+1)*yvel(p,q+1)*ener(p,q+1) - dens(p,q)*yvel(p,q)*ener(p,q))/(delx)
                    !Flux_plus = (dens(p,q)*ener(p,q)*xvel(p,q) + dens(p+1,q)*ener(p+1,q)*xvel(p+1,q))/2.0_DP
                    !Flux_minus = (dens(p,q)*ener(p,q)*xvel(p,q) + dens(p-1,q)*ener(p-1,q)*xvel(p-1,q))/2.0_DP
                    !rhs(4,(p-1)*NP + q) = rhs(4,(p-1)*NP + q) + (Flux_plus - Flux_minus)/delx
                    !Flux_plus = (dens(p,q)*ener(p,q)*yvel(p,q) + dens(p,q+1)*ener(p,q+1)*yvel(p,q+1))/2.0_DP
                    !Flux_minus = (dens(p,q)*ener(p,q)*yvel(p,q) + dens(p,q-1)*ener(p,q-1)*yvel(p,q-1))/2.0_DP
                    !rhs(4,(p-1)*NP + q) = rhs(4,(p-1)*NP + q) + (Flux_plus - Flux_minus)/delx

                    usign = (xvel(p,q)+xvel(p+1,q)+1.0e-10_DP)/abs((xvel(p,q)+xvel(p+1,q))+1.0e-10_DP)
                    Flux_plus = max(usign,0.0_DP)*(dens(p,q)*ener(p,q)*0.5_DP*(xvel(p,q)+xvel(p+1,q))) &
                        - min(usign,0.0_DP)*(dens(p+1,q)*ener(p+1,q)*0.5_DP*(xvel(p,q)+xvel(p+1,q)))
                    usign = (xvel(p,q)+xvel(p-1,q)+1.0e-10_DP)/abs((xvel(p,q)+xvel(p-1,q))+1.0e-10_DP)
                    Flux_minus = max(usign,0.0_DP)*(dens(p-1,q)*ener(p-1,q)*0.5_DP*(xvel(p,q)+xvel(p-1,q))) &
                        - min(usign,0.0_DP)*(dens(p,q)*ener(p,q)*0.5_DP*(xvel(p,q)+xvel(p-1,q)))
                    rhs(4,(p-1)*NP + q) = rhs(4,(p-1)*NP + q) - (Flux_plus - Flux_minus) / delx
                    
                    vsign = (yvel(p,q)+yvel(p,q+1)+1.0e-10_DP)/abs((yvel(p,q)+yvel(p,q+1))+1.0e-10_DP)
                    Flux_plus = max(vsign,0.0_DP)*(dens(p,q)*ener(p,q)*0.5_DP*(yvel(p,q)+yvel(p,q+1))) &
                        - min(vsign,0.0_DP)*(dens(p,q+1)*ener(p,q+1)*0.5_DP*(yvel(p,q)+yvel(p,q+1)))
                    vsign = (yvel(p,q)+yvel(p,q-1)+1.0e-10_DP)/abs((yvel(p,q)+yvel(p,q-1))+1.0e-10_DP)
                    Flux_minus = max(vsign,0.0_DP)*(dens(p,q-1)*ener(p,q-1)*0.5_DP*(yvel(p,q)+yvel(p,q-1))) &
                        - min(vsign,0.0_DP)*(dens(p,q)*ener(p,q)*0.5_DP*(yvel(p,q)+yvel(p,q-1)))
                    rhs(4,(p-1)*NP + q) = rhs(4,(p-1)*NP + q) - (Flux_plus - Flux_minus) / dely
                    
                    rhs(4,(p-1)*NP + q) = rhs(4,(p-1)*NP + q) - pressure(p,q)*(xvel(p+1,q)-xvel(p-1,q))/(2.0_DP*delx) &
                        - pressure(p,q)*(yvel(p,q+1)-yvel(p,q-1))/(2.0_DP*dely)


                    diff_plus = lam_rhocp*cp*(dens(p,q) + dens(p+1,q))/2.0_DP*(temperature(p+1,q)-temperature(p,q))/delx
                    diff_minus = lam_rhocp*cp*(dens(p,q) + dens(p-1,q))/2.0_DP*(temperature(p,q)-temperature(p-1,q))/delx
                    rhs(4,(p-1)*NP + q) = rhs(4,(p-1)*NP + q) + (diff_plus - diff_minus)/delx

                    diff_plus = lam_rhocp*cp*(dens(p,q) + dens(p,q+1))/2.0_DP*(temperature(p,q+1)-temperature(p,q))/dely
                    diff_minus = lam_rhocp*cp*(dens(p,q) + dens(p,q-1))/2.0_DP*(temperature(p,q)-temperature(p,q-1))/dely
                    rhs(4,(p-1)*NP + q) = rhs(4,(p-1)*NP + q) + (diff_plus - diff_minus)/dely

                    !dissip = (5.0_DP/9.0_DP*(((xvel(p+1,q) - xvel(p-1,q))/(2.0_DP*delx))**2+((yvel(p,q+1) - yvel(p,q-1))/(2.0_DP*delx))**2) &
                    !   - 8.0_DP/9.0_DP*((xvel(p+1,q) - xvel(p-1,q))/(2.0_DP*delx))*((yvel(p,q+1) - yvel(p,q-1))/(2.0_DP*delx)) &
                    !   + 1.0_DP/2.0_DP*(((xvel(p,q+1) - xvel(p,q-1))/(2.0_DP*delx))**2+((yvel(p+1,q) - yvel(p-1,q))/(2.0_DP*delx))**2) &
                    !   + ((xvel(p,q+1) - xvel(p,q-1))/(2.0_DP*delx))*((yvel(p+1,q) - yvel(p-1,q))/(2.0_DP*delx)) )

                    !dissip = ((xvel(p+1,q) - xvel(p-1,q))/(2.0_DP*delx))**2 &
                    !    + ((yvel(p,q+1) - yvel(p,q-1))/(2.0_DP*dely))**2 &
                    !    + 0.5_DP*((xvel(p,q+1) - xvel(p,q-1))/(2.0_DP*dely) + (yvel(p+1,q) - yvel(p-1,q))/(2.0_DP*delx))**2 &
                    !    - 1.0_DP/3.0_DP*((xvel(p+1,q) - xvel(p-1,q))/(2.0_DP*delx) + (yvel(p,q+1) - yvel(p,q-1))/(2.0_DP*dely))**2
                    dissip = 2.0_DP/3.0_DP*(((xvel(p+1,q) - xvel(p-1,q))/(2.0_DP*delx))**2 &
                        + ((yvel(p,q+1) - yvel(p,q-1))/(2.0_DP*dely))**2 &
                        - ((xvel(p+1,q) - xvel(p-1,q))/(2.0_DP*delx))*((yvel(p,q+1) - yvel(p,q-1))/(2.0_DP*dely))) &
                        + 0.5_DP*(((xvel(p+1,q) - xvel(p-1,q))/(2.0_DP*dely))**2 &
                        + ((yvel(p,q+1) - yvel(p,q-1))/(2.0_DP*delx))**2)
                    rhs(4,(p-1)*NP + q) = rhs(4,(p-1)*NP + q) + 2.0_DP*nu*dens(p,q)*dissip

                    ! Update variables
                    dens_new(p,q) = dens(p,q) + rhs(1,(p-1)*NP + q)*delt
                    xvel_new(p,q) = (dens(p,q)*xvel(p,q) + rhs(2,(p-1)*NP + q)*delt)/dens_new(p,q)
                    yvel_new(p,q) = (dens(p,q)*yvel(p,q) + rhs(3,(p-1)*NP + q)*delt)/dens_new(p,q)
                    ener_new(p,q) = (dens(p,q)*ener(p,q) + rhs(4,(p-1)*NP + q)*delt)/dens_new(p,q)
                    !print *, "RHS", p,q, rhs(2,(p-1)*NP + q), xvel_new(p,q), rhs(3,(p-1)*NP + q), yvel_new(p,q)
                end do
            end do
            !$OMP END PARALLEL DO

            dens(1:NP,1:NP) = dens_new
            xvel(1:NP,1:NP) = xvel_new
            yvel(1:NP,1:NP) = yvel_new
            ener(1:NP,1:NP) = ener_new

        end subroutine compressible_updatevecs

    end subroutine compressible_run

subroutine compressible_output(output_file)
    use precision
    implicit none
    character(len=*), intent(in) :: output_file
    integer :: unit, ios
    character(len=200) :: cmd, dirpath
    integer :: i, j, p

    ! === Extract directory path from output_file ===
    p = len_trim(output_file)
    do while (p > 0 .and. output_file(p:p) /= '/')
        p = p - 1
    end do

    if (p > 0) then
        dirpath = output_file(1:p-1)

        ! Create directory if it doesn't exist
        write(cmd,'("mkdir -p ",A)') trim(dirpath)
        call execute_command_line(cmd)
    end if

    ! === Open file ===
    unit = 99
    open(unit=unit, file=trim(output_file), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        print *, "Error: cannot open output file ", trim(output_file)
        return
    end if

    ! === Write arrays ===
    write(unit,*) "dens:"
    do i = 1, NP
        write(unit,'(1000ES15.7)') dens(1:NP,i)
    end do

    write(unit,*) "xvel:"
    do i = 1, NP
        write(unit,'(1000ES15.7)') xvel(1:NP,i)
    end do

    write(unit,*) "yvel:"
    do i = 1, NP
        write(unit,'(1000ES15.7)') yvel(1:NP,i)
    end do

    write(unit,*) "ener:"
    do i = 1, NP
        write(unit,'(1000ES15.7)') ener(1:NP,i)
    end do

    write(unit,*) "temperature:"
    do i = 1, NP
        write(unit,'(1000ES15.7)') temperature(1:NP,i)
    end do

    write(unit,*) "pressure:"
    do i = 1, NP
        write(unit,'(1000ES15.7)') pressure(1:NP,i)
    end do

    close(unit)

end subroutine compressible_output


end module