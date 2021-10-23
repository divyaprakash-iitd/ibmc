module mod_ib
    use iso_fortran_env, only: int32, real64
    use mod_particle
    implicit none

    private
    public :: ib

    type :: ib
        character(len=:), allocatable :: name
        integer(int32) :: np
        type(particle), allocatable :: boundary(:)
        ! To-do: Add neighbours information
    end type ib

    interface ib
        module procedure :: ib_constructor
    end interface ib

contains

    pure type(ib) function ib_constructor(name,np) result(self)
        ! Creates an immersed boundary of particles
        character(len=*), intent(in) :: name
        integer(int32), intent(in) :: np

        ! Indices
        integer(int32) :: i

        ! Allocate the boundary
        allocate(self%boundary(np))

        self%np = np
        do concurrent (i = 1:np)
            self%boundary(i) = particle(i)
        end do

    end function ib_constructor
    
end module mod_ib
