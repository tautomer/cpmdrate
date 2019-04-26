module stp
    implicit none

    contains 
    subroutine stopgm(msgp1, msgp2)
        implicit none
    
        character(len=*), intent(in) :: msgp1
        character(len=*), optional, intent(in) :: msgp2
        character(len=50) :: msg
    
        if(present(msgp2)) then
            write(msg, '(3a)') msgp1, msgp2, '.'
        else
            write(msg, '(2a)') msgp1, '.'
        end if
        ! for some unknown reasons, ifort does not stop + arg
#if defined(__INTEL_COMPILER)
        write(*, '(a)') msg
        stop
#else
        stop msg
#endif
    end subroutine
end module

module init
    implicit none

    contains
    subroutine readconf()
        implicit none

        integer :: line = 0, ioerr = 0, pos
        character(len=200) :: buffer, label*20, config*10
    
        call get_command_argument(1, config)
        open(unit=20, file=config, action="read")
        do while(ioerr == 0)
            read(20, '(A)', iostat=ioerr) buffer
            if(ioerr == 0) then
                line = line + 1
                pos = scan(buffer, ' ')
                label = buffer(1: pos)
                buffer = buffer(pos+1: )
                call read_buffer(label, buffer, line)
            end if
        end do
        close(20)
    end subroutine

    subroutine read_buffer(label, buffer, line)
        use pmf
        use kappa
        use stp
        implicit none
    
        integer :: ioerr
        integer, intent(in) :: line
        real(kind=16) :: kb = 3.16679d-6
        character(len=*), intent(in) :: buffer, label

        select case (label)
        ! general settings
        ! need fix twham & tbl & tdiff & tproj
        case ("f")
            read(buffer, *, iostat=ioerr) flnmpmf
            if (ioerr > 0) then
                call stopgm("PMF filename error.")
            else if (ioerr < 0) then
                call stopgm("PMF filename not sepcified")
            end if
        case ("gr")
            read(buffer, *, iostat=ioerr) gr
            if (ioerr > 0) then
                call stopgm("Mass should be a real number in atomic units.")
            else if (ioerr < 0) then
                call stopgm("Mass not sepcified.")
            end if
        case ("temperature")
            read(buffer, *, iostat=ioerr) temp
            if (ioerr > 0) then
                call stopgm("Temperature should be a real number")
            else if (ioerr < 0) then
                call stopgm("Temperature not sepcified")
            end if
            beta = 1 / (temp * kb)
        case ("side")
            read(buffer, *, iostat=ioerr) int0
            if (ioerr > 0) then
                call stopgm("Provide either 0 (left) or 1 (right).")
            else if (ioerr < 0) then
                call stopgm("Reactant side not sepcified")
            end if
        case ("k")
            read(buffer, *, iostat=ioerr) k_plateau
            if (ioerr > 0) then
                call stopgm("kappa should be a real number")
            else if (ioerr < 0) then
                call stopgm("kappa value not sepcified")
            end if
        case default
            write(*, "(a, i0)") "Skipping invalid label at line ", line
        end select
    
    end subroutine
end module
