module pmf
    implicit none

    integer, private :: i, ioerr
    real*8, private :: tmp, tmp1

    integer, public :: int0
    real*8, public :: temp, beta, rate
    real*8, public, allocatable :: xi(:), fe(:)
    character(len=30), public :: flnmpmf
    contains
    subroutine read_free_ener()
        implicit none

        real*8 :: x

        allocate(xi(0), fe(0))
        open(unit=20, file=flnmpmf, action='read', iostat=ioerr)
        if (ioerr /= 0) then
            write(*, "(2a)") flnmpmf, ' not found.'
            stop
        end if

        do while (ioerr == 0)
            read(20, *, iostat=ioerr) x, tmp
            if (ioerr > 0 ) then
                write(*, "(2a)") 'Error in file ', flnmpmf
                stop
            else if (ioerr == 0 ) then
                xi = [xi, x]
                fe = [fe, tmp]
            end if
        end do

    end subroutine

    subroutine tst_rate()
        implicit none

        integer :: idx(3), size_reac
        real*8 :: numer, denom
        real*8, allocatable :: xi_reac(:), fe_reac(:)

        call find_idx(idx)
        size_reac = idx(3) - idx(2) + 1
        allocate(xi_reac(size_reac), fe_reac(size_reac))
        xi_reac = xi(idx(2):idx(3))
        !TODO: parllalize this with omp
        fe_reac = dexp(fe(idx(2):idx(3))) * 0.5
        numer = dexp(-beta * fe(idx(1)))
        denom = 0
        do i = 2, size_reac
            tmp = xi_reac(i) - xi_reac(i-1)
            tmp1 = fe_reac(i) + fe_reac(i-1)
            denom = denom + tmp * tmp1
        end do
        rate = numer / denom

    end subroutine

    subroutine find_idx(idx)
        implicit none

        integer, intent(out) :: idx(3)

        idx(1) = maxloc(fe, dim=1)
        if (int0 == 0) then
            idx(2) = 1
            idx(3) = idx(1)
        else
            idx(2) = idx(1)
            idx(3) = size(fe)
        end if

    end subroutine
end module
