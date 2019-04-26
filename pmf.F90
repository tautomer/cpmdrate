module pmf
    implicit none

    integer, private :: i, ioerr
    real(kind=16), private :: tmp, tmp1

    integer  :: int0
    real(kind=16) :: gr, temp, beta, rate
    real(kind=16), private, allocatable :: xi(:), fe(:)
    character(len=30) :: flnmpmf
    contains
    subroutine read_free_ener()
        use stp
        implicit none

        real(kind=16) :: x, a2au = 1.889725989, kcal2au = 1 / 6.27509d2
        character(len=50) :: buffer

        allocate(xi(0), fe(0))
        open(unit=20, file=flnmpmf, action='read', iostat=ioerr)
        if (ioerr /= 0) then
            write(*, "(2a)") flnmpmf, ' not found.'
            stop
        end if

        do while (ioerr == 0)
            read(20, '(a)', iostat=ioerr) buffer
            if (ioerr > 0 ) then
                call stopgm('Error in file ', flnmpmf)
            else if (ioerr == 0 ) then
                if (index(buffer(1:3), '#') /= 0) cycle
                read(buffer, *, iostat=ioerr) x, tmp
                if (ioerr /= 0 ) call stopgm('Error in file ', flnmpmf)
                xi = [xi, x * a2au]
                fe = [fe, tmp * kcal2au]
            end if
        end do

    end subroutine

    subroutine tst_rate()
        implicit none

        integer :: idx(2), size_reac
        real(kind=16) :: numer, denom, minv, prefac
        real(kind=16), parameter :: pi = 3.14159265
        real(kind=16), allocatable :: xi_reac(:), fe_reac(:)

        call find_idx(idx)
        !idx(1) = size(fe)/2 + 1
        !idx(2) = size(fe)
        size_reac = idx(2) - idx(1) + 1
        allocate(xi_reac(size_reac), fe_reac(size_reac))
        fe_reac = fe(idx(1):idx(2))
        !write(*, '(f6.3,2i5)') fe(idx(1))*627.509, idx(1), idx(2)
        write(*, '(f6.3,2i5)') fe(idx(1)), idx(1), idx(2)
        xi_reac = xi(idx(1):idx(2)) * 0.5
        minv = minval(fe_reac, dim=1)
        fe_reac = fe_reac - minv
        write(*, '(2f9.4)') beta, fe_reac(1)
        fe_reac = exp(-beta * fe_reac)
        if (int0 == 0) then
            numer = fe_reac(size_reac)
        else
            numer = fe_reac(1)
        end if
        denom = 0
        do i = 2, size_reac
            tmp = xi_reac(i) - xi_reac(i-1)
            tmp1 = fe_reac(i) + fe_reac(i-1)
            denom = denom + tmp * tmp1
        end do
        write(*, '(2e14.6)') numer, denom
        prefac = 1 / sqrt(2*pi*beta) * gr
        rate = prefac * numer / denom

    end subroutine

    subroutine find_idx(idx)
        implicit none

        integer, intent(out) :: idx(2)
        integer :: i, j, k, incre, tmp
        real(kind=16) :: fo, fn

        i = minloc(fe, dim=1)
        if (int0 == 0) then
            incre = -1
        else
            incre = 1
        end if
        fo = fe(i)
        do while(.true.)
            fn = fe(i)
            if (fn >= fo) then
                fo = fn
            else
                tmp = i + incre
                if (fe(tmp) < fn) then
                    i = i - incre
                    exit
                else
                    i = tmp
                    fo = fe(i)
                end if
            end if
            i = i + incre
        end do
        if (int0 == 0) then
            idx(1) = 1
            idx(2) = i
        else
            idx(1) = i
            idx(2) = size(fe)
        end if

    end subroutine
end module
