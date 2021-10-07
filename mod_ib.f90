module mod_ib
    use iso_fortran_env, only: int32, real32, real64
    implicit none

    private
    public :: particle

    ! Immersed boundary
    type :: particle
        integer(int32)                  :: n
        real(real64), allocatable       :: x, y
        real(real64), allocatable       :: Fx, Fy
    end type particle

    interface particle
        module procedure :: ib_constructor
    end interface particle

contains

    pure type(particle) function ib_constructor(n) result(self)
            integer(int32), intent(in)      :: n

            self%n  = n ! Particle tag
            self%x  = 0.0d0
            self%y  = 0.0d0
            self%Fx = 0.0d0
            self%Fy = 0.0d0
    end function ib_constructor
    
end module mod_ib