program test_code
    use iso_fortran_env, only: int32, real64, real64
    use mod_particle
    use mod_ib
    implicit none

    type(ib) :: circle
    circle = ib('circle',2)

    circle%boundary(2)%Fx = 2.0d0
    print *, circle%boundary(2)%Fx
end program test_code