module scheme2
    use precision
    implicit none
    real(DP), dimension(:,:), pointer :: sol_array

    contains
    subroutine scheme2_run(NX, delt, nu, t_end)
        use precision
        implicit none
        integer, intent(in) :: NX
        real(DP), intent(in) :: delt
        real(DP), intent(in) :: nu
        real(DP), intent(in) :: t_end

        real(DP), dimension(0:NX+1):: u_vec, u_new

        real(DP), parameter :: pi = 3.141592653589793_DP

        real(DP) :: delx
        real(DP) :: time
        real(DP), dimension(NX) :: x_vec
        integer :: p
        integer :: time_int,NT

        delx = (2.0_DP*pi)/real(NX,DP)
        do p = 1,NX
            x_vec(p) = 0.0_DP + delx*real(p-1,DP)
            u_vec(p) = sin(x_vec(p))*exp(-(x_vec(p)-pi)**2)
        end do
        u_vec(0) = u_vec(NX)
        u_vec(NX+1) = u_vec(1)
        time = 0.0_DP
        NT = 2 + int(t_end/delt,8)
        allocate(sol_array(0:NT,0:NX))
        sol_array = -1.0_DP
        sol_array(0,1:) = x_vec
        sol_array(1,0) = time
        sol_array(1,1:) = u_vec(1:NX)

        do time_int = 1,(NT-1)
            time = time + delt
            call scheme2_advancetime(u_vec, u_new)
            sol_array(1+time_int,0) = time
            sol_array(1+time_int,1:) = u_new(1:NX)
            u_vec(:) = u_new(:)
        end do

        return
        contains 
        subroutine scheme2_advancetime(u_old, u_vec_next)
            implicit none
            real(DP), dimension(0:NX+1), intent(in) :: u_old
            real(DP), dimension(0:NX+1), intent(out) :: u_vec_next

            real(DP), dimension(NX,NX) :: jac, jac_approx
            real(DP), dimension(NX) :: F_vec
            real(DP), dimension(NX) :: diff
            real(DP), dimension(NX) :: adv
            integer :: p
            integer :: iter, max_iter, info
            integer, dimension(NX) :: ipiv
            real(DP) :: usign
            real(DP) :: tol = 1.0e-8_DP
            max_iter = 1000

            do iter = 1, max_iter
                jac = 0.0_DP
                do p = 1,NX
                    if (abs(u_vec_next(p)) > 0.0_DP) then
                        usign = u_vec_next(p)/abs(u_vec_next(p))
                    else
                        usign = 0.0_DP
                    end if
                    adv(p) = max(usign,0.0_DP) * u_vec_next(p)*(u_vec_next(p) - u_vec_next(p-1))/(delx) &
                          -  min(usign,0.0_DP) * u_vec_next(p)*(u_vec_next(p+1) - u_vec_next(p))/(delx)
                    diff(p) = nu*(u_vec_next(p+1) - 2.0_DP*u_vec_next(p) + u_vec_next(p-1))/(delx**2)
                    F_vec(p) = u_vec_next(p) - u_old(p) + delt*(adv(p) - diff(p))

                    jac(p,p) = 1.0_DP &
                        + delt/delx*max(usign,0.0_DP)*(2.0_DP*u_vec_next(p) - u_vec_next(p-1)) &
                        - delt/delx*min(usign,0.0_DP)*(u_vec_next(p+1) - 2.0_DP*u_vec_next(p)) &
                        + 2.0_DP*nu*delt/delx**2
                    if (p.eq.1) then
                        jac(p,p+1) = -delt/delx*min(usign,0.0_DP)*u_vec_next(p) - nu*delt/delx**2
                        jac(p,NX) = -delt/delx*max(usign,0.0_DP)*u_vec_next(p) - nu*delt/delx**2
                    else if (p.eq.NX) then
                        jac(p,1) = -delt/delx*min(usign,0.0_DP)*u_vec_next(p) - nu*delt/delx**2
                        jac(p,p-1) = -delt/delx*max(usign,0.0_DP)*u_vec_next(p) - nu*delt/delx**2
                    else
                        jac(p,p+1) = -delt/delx*min(usign,0.0_DP)*u_vec_next(p) - nu*delt/delx**2
                        jac(p,p-1) = -delt/delx*max(usign,0.0_DP)*u_vec_next(p) - nu*delt/delx**2
                    end if
                end do
                F_vec = -F_vec
                call dgesv(NX, 1, jac, NX, ipiv, F_vec, NX, info)
                if (info .ne. 0) then
                    print *, "DGESV failed:", info
                    stop
                end if

                u_vec_next(1:NX) = u_vec_next(1:NX) + F_vec
                u_vec_next(0) = u_vec_next(NX)
                u_vec_next(NX+1) = u_vec_next(1)
                if (sqrt(sum(F_vec**2)) .le. tol) then
                    exit
                end if
            end do

            return
        end subroutine scheme2_advancetime
    end subroutine scheme2_run

    subroutine scheme2_output(output_file)
        use precision
        implicit none
        character(len=100), intent(in) :: output_file
        
        integer :: nrow, ncol, unit, ierr
        integer :: p,t

        nrow = size(sol_array,1)
        ncol = size(sol_array,2)

        unit = 50
        open(unit, file=output_file, form="formatted", iostat = ierr)
        do t = 0,nrow-1
            do p = 0,ncol-1
                if (abs(sol_array(t,p)).ge.1.0E100_DP) sol_array(t,p) = sign(1.0E99_DP, sol_array(t,p))
                if (isnan(sol_array(t,p))) sol_array(t,p) = 1.0E99_DP
                write(unit,'(es20.12)',advance='no') sol_array(t,p)
            end do
            write(unit,*)
        end do
        close(unit)
        deallocate(sol_array)

    end subroutine scheme2_output

end module