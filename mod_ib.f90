module mod_ib
    use iso_fortran_env, only: int32, real64
    use mod_particle
    implicit none

    private
    public :: ib

    type :: ib
        integer(int32) :: np
        type(particle), allocatable :: boundary(:)
        character(len=1) :: t ! Open or closed
        ! To-do: Add neighbours information
    end type ib

    interface ib
        module procedure :: ib_constructor
    end interface ib

contains

    pure type(ib) function ib_constructor(np) result(self)
        ! Creates an immersed boundary of particles
        integer(int32), intent(in) :: np

        ! Indices
        integer(int32) :: i

        ! Allocate the boundary
        allocate(self%boundary(np))

        ! Initializes particles in each Immersed Boundary(IB) layer
        self%np = np
        self%t = 'o' ! Open by default
        do i = 1,np
            self%boundary(i) = particle()
        end do

    end function ib_constructor
    
end module mod_ib
