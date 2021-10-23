module mod_cilia
    use iso_fortran_env, only: int32, real64, real64
    use mod_particle
    use mod_ib
    implicit none
    
    private
    public :: cilia

    type :: cilia
        integer(int32) :: nl ! Number of layers
        integer(int32) :: np ! Number of particles per layer 
        type(ib), allocatable :: layers(:) 
    end type cilia

    interface cilia
        module procedure ::cilia_constructor
    end interface cilia

contains

    pure type(cilia) function cilia_constructor(nl,np) result(self)
        integer(int32), intent(in) :: nl
        integer(int32), intent(in) :: np
        
        integer(int32) :: il, ip 

        self%nl = nl
        self%np = np

        allocate(self%layers(nl))
        ! Create immersed boundaries (ib) for each layer
        do il = 1,nl
            allocate(self%layers(il)%boundary(np))
        end do 
    end function cilia_constructor
    
end module mod_cilia