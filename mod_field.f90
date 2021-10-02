module mod_field
    use iso_fortran_env, only: int32, real32, int64, real64
    implicit none

    type :: dimension
        character(len=5) :: name
        integer(int32) :: lb, ub
    end type dimension

    ! type :: Field

    !     character(len=:), allocatable :: name
    !     integer(int32):: lb(2), ub(2)
    !     real(real64), allocatable :: data(:,:)

    ! contains
    !     procedure, private, pass(self) :: diffx
    ! end type Field

contains
    
end module mod_field