program test
    use init
    use pmf
    use kappa
    implicit none

    call readconf()
    call read_free_ener()
    call tst_rate()
    write(*, '(a,f5.1,a,f11.7)') 'Full rate constant at', temp, 'K is ', &
    k_plateau * rate
end program
