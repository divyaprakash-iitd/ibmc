module mod_particle
! 
! Purpose:
!   To define the derived data type, particle, used to represent the node
!   of an immersed boundary.
!   o__o__o__o__o__o
!   The circles represent a node/particle in the above illustration.
!   The entire structure is an immersed boundary.
!    
    use iso_fortran_env, only: int32, real64, real64
    implicit none

    private
    public :: particle

    type :: particle
        real(real64), allocatable :: x, y     ! Location of the particle
        real(real64), allocatable :: xo, yo   ! Original location of the particle
        real(real64), allocatable :: Fx, Fy   ! The componenets of force acting on a particle
        real(real64), allocatable :: Ux, Uy   ! The velocity components of a particle
    end type particle

    interface particle
        module procedure :: ib_constructor
    end interface particle

contains

    pure type(particle) function ib_constructor() result(self)
    !
    ! Purpose:
    !   To initialize all the values of the attributes of to zero.
    !
            self%x  = 0.0d0
            self%y  = 0.0d0
            self%xo = 0.0d0
            self%yo = 0.0d0
            self%Fx = 0.0d0
            self%Fy = 0.0d0
            self%Ux = 0.0d0
            self%Uy = 0.0d0
    end function ib_constructor
    
end module mod_particle