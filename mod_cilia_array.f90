module mod_cilia_array
    use iso_fortran_env, only: int32, int32, real64, real64
    use mod_cilia
    implicit none
   
    private
    public :: cilia_array

    type :: cilia_array
        type(cilia), allocatable :: array(:)
        integer(int32) :: nc 
    end type cilia_array

    interface cilia_array
        module procedure :: cilia_array_constructor
    end interface cilia_array

contains
    
    pure type(cilia_array) function cilia_array_constructor(nc,nl,np) result(self)
        integer(int32), intent(in) :: nl
        integer(int32), intent(in) :: np
        integer(int32), intent(in) :: nc

        integer(int32) :: ic

        self%nc = nc
        allocate(self%array(nc))

        ! Initializes each cilia in the array
        do ic = 1,nc
            self%array(ic) = cilia(nl,np)
        end do
    end function cilia_array_constructor
    
end module mod_cilia_array
