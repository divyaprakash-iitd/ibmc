module mod_closed_cilia
    use iso_fortran_env, only: int32, real64, real64
    use mod_mesh
    use mod_ib
    use mod_vec
    use mod_cilia
    implicit none

    real(real64), parameter :: PI = 3.141592653589793

    private
    public:: create_closed_loop
    
contains

    subroutine create_closed_loop(C,W,r,origin)
        class(cilia),   intent(in out)    :: C
        class(vec),     intent(in)        :: origin
        real(real64),   intent(in)        :: W          ! Width of cilia
        real(real64),   intent(in)        :: r          ! Radius of the loop

        real(real64)    :: dtheta   ! Spacing of particles in the closed loop
        integer(int32)  :: nl       ! No. of layers
        integer(int32)  :: np       ! No. of particle

        integer(int32) :: il,ip

        ! Assign values
        nl = C%nl
        np = C%np

        ! Calculate angular spacing
        dtheta = 2*PI/np

        ! Assign locations to the particles of the cilia
        do il = 1,nl
            do ip = 1,np 
                c%layers(il)%boundary(ip)%x = origin%x + (r - W*(il-1)) * cos((ip-1)*dtheta)
                c%layers(il)%boundary(ip)%y = origin%y + (r - W*(il-1)) * sin((ip-1)*dtheta)
            end do
        end do




    end subroutine create_closed_loop
    
end module mod_closed_cilia