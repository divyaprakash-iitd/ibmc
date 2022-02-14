module mod_ib
!
!   Purpose:
!       To define the derived type to represent an immersed boundary made up of
!       an array of the particles derived data type.
!
!       o__o__o__o__o__o       
!       
!       The structure in the above illustration represents the immersed boundary,
!       which is basically particles connected by links.
!    

    use iso_fortran_env, only: int32, real64
    use mod_particle
    implicit none

    private
    public :: ib

    type :: ib
        integer(int32)              :: np           ! The number of particles in an immersed boundary
        type(particle), allocatable :: boundary(:)  ! An array of particles
        character(len=1)            :: t            ! Describes the type of immersed boundary (Open or closed)
        ! To-do: Add neighbours information
    end type ib

    interface ib
        module procedure :: ib_constructor
    end interface ib

contains

    pure type(ib) function ib_constructor(np) result(self)
        ! 
        !   Purpose:
        !       To create an immersed boundary of particles and initialize it's attributes with zeros.
        !

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
