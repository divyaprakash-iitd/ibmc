module mod_closed_cilia
    use iso_fortran_env, only: int32, real64, real64
    use mod_mesh
    use mod_ib
    use mod_vec
    use mod_cilia
    use mod_cilia_array
    use mod_ibm
    implicit none

    real(real64), parameter :: PI = 3.141592653589793

    private
    public:: create_closed_loop, create_closed_loop_array, calculate_closed_loop_array_force
    
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
            c%layers(il)%t = 'c' ! Closed boundary type
            do ip = 1,np 
                c%layers(il)%boundary(ip)%x = origin%x + (r - W*(il-1)) * cos((ip-1)*dtheta)
                c%layers(il)%boundary(ip)%y = origin%y + (r - W*(il-1)) * sin((ip-1)*dtheta)
            end do
        end do

    end subroutine create_closed_loop
   
    subroutine create_closed_loop_array(CA,W,r,origin)
        class(cilia_array),   intent(in out)    :: CA
        class(vec),     intent(in)        :: origin
        real(real64),   intent(in)        :: W          ! Width of cilia
        real(real64),   intent(in)        :: r          ! Radius of the loop

        integer(int32) :: ic

        ! Assign locations to the particles of the cilia
        do ic = 1,CA%nc
            call create_closed_loop(CA%array(ic),W,r,origin)
        end do
    end subroutine create_closed_loop_array

    subroutine calculate_closed_loop_array_force(CA,ko,kd,Rlv,Rlh,Rld)
        class(cilia_array), intent(in out)  :: CA      ! Cilia array
        real(real64),       intent(in)      :: ko      ! Horizontal spring stiffness
        real(real64),       intent(in)      :: kd      ! Diagonal spring stiffness
        real(real64),       intent(in)      :: Rlv,Rlh,Rld      ! Resting length

        integer(int32) :: ic

        ! Calculate forces for all the cilia within the array
        do ic = 1,CA%nc
            call calculate_closed_loop_force(CA%array(ic),ko,kd,Rlv,Rlh,Rld)
        end do
    end subroutine calculate_closed_loop_array_force

    subroutine calculate_closed_loop_force(C,ko,kd,Rlv,Rlh,Rld)
        class(cilia), intent(in out)    :: C         ! Cilia structure
        real(real64), intent(in)        :: ko        ! Horizontal spring stiffness
        real(real64), intent(in)        :: kd        ! Diagonal spring stiffness
        real(real64), intent(in)        :: Rlv,Rlh,Rld        ! Resting length

        integer(int32)  :: il, ip

        ! Initialize the forces to zero on all the nodes at every time step
        do il = 1,C%nl 
            do ip = 1,C%np
                C%layers(il)%boundary(ip)%Fx = 0.0d0
                C%layers(il)%boundary(ip)%Fy = 0.0d0
            end do
        end do

        ! Calculate forces on the layers
        call calculate_spring_force(C%layers(1),ko,kd,Rlv,C%layers(1)%t)
        call calculate_spring_force(C%layers(2),ko,kd,0.5d0*Rlv,C%layers(2)%t)

        ! do il = 1,C%nl
        !     call calculate_spring_force(C%layers(il),ko,kd,Rlv,C%layers(il)%t)
        ! end do

        ! Calculate forces on the horizontal links
        call calculate_horizontal_link_force(C%layers(1),C%layers(2),ko,Rlh)

        ! Calculate forces on the diagonal links (negative slope)
        call calculate_diagonal_link_force(C%layers(1),C%layers(2),kd,Rld)
        call calculate_diagonal_link_force_pos(C%layers(2),C%layers(1),kd,Rld)
    end subroutine calculate_closed_loop_force



end module mod_closed_cilia