program test
    use init
    use pmf
    use kappa
    implicit none

    call readconf()
    call read_free_ener()
    call tst_rate()
    write(*, '(a,f6.1,a,e13.6)') 'Full rate constant at', temp, 'K is ', &
    k_plateau * rate
end program
