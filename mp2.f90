program mp2
    use compressible
    use precision
    implicit none

    integer :: num_x
    real(DP) :: del_time
    real(DP) :: end_time
    real(DP) :: del_x
    character(len=100) :: output_dir
    character(len=100) :: arg

    integer :: s


    call get_command_argument(1,arg)
    read(arg,*) num_x

    call get_command_argument(2,arg)
    read(arg,*) del_time

    call get_command_argument(3,arg)
    read(arg,*) end_time

    call get_command_argument(4,output_dir)
    output_dir = trim(output_dir)

    call compressible_initialize(num_x, 0.025_DP, 1.0_DP)
    call compressible_run(del_time, end_time)
    call compressible_output(output_dir)

end program mp2