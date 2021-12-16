module mod_particle
    use iso_fortran_env, only: int32, real64, real64
    implicit none

    private
    public :: particle

    ! Immersed boundary
    type :: particle
        real(real64), allocatable       :: x, y
        real(real64), allocatable       :: xo, yo
        real(real64), allocatable       :: Fx, Fy
        real(real64), allocatable       :: Ux, Uy 
    end type particle

    interface particle
        module procedure :: ib_constructor
    end interface particle

contains

    pure type(particle) function ib_constructor() result(self)
            self%x  = 0.0d0
            self%y  = 0.0d0
            self%xo  = 0.0d0
            self%yo  = 0.0d0
            self%Fx = 0.0d0
            self%Fy = 0.0d0
            self%Ux = 0.0d0
            self%Uy = 0.0d0
    end function ib_constructor
    
end module mod_particle