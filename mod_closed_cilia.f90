module mod_closed_cilia
!
!   Purpose:
!       To provide subroutines to create closed cilia loops structures
!       
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
    !
    !   Purpose:
    !       To create closed loop structures consisting of 2 immersed boundaries.
    !       A circular structure is created with support of 2 layers and np particles.
    !       The particles are equally spaced along the circumference.
    !       The particles are assumed to be connected by individual springs.
    !       The boundary type attribute 't' of the structure is specified as closed.
    !
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
   
    subroutine create_closed_loop_ellipse(C,W,r,ar,origin)
    !
    !   Purpose:
    !       To create closed loop structures consisting of 2 immersed boundaries.
    !       An elliptical structure is created with support of 2 layers and np particles.
    !       The particles are equally spaced along the circumference.
    !       The particles are assumed to be connected by individual springs.
    !       The boundary type attribute 't' of the structure is specified as closed.
    !
        class(cilia),   intent(in out)    :: C
        class(vec),     intent(in)        :: origin
        real(real64),   intent(in)        :: W          ! Width of cilia
        real(real64),   intent(in)        :: r          ! Radius of the loop
        real(real64),   intent(in)        :: ar         ! Aspect Ratio

        real(real64)    :: dtheta   ! Spacing of particles in the closed loop
        integer(int32)  :: nl       ! No. of layers
        integer(int32)  :: np       ! No. of particle
        real(real64) :: particle_loc(C%np,2)
        integer(int32) :: il,ip

        real(real64) :: a, b

        ! Semi-minor and semi-major axis
        a = r    
        b = a*ar

        ! Assign values
        nl = C%nl
        np = C%np

        ! Calculate angular spacing
        dtheta = 2*PI/np

    
        ! Assign locations to the particles of the cilia
        do il = 1,nl
            c%layers(il)%t = 'c' ! Closed boundary type
            ! Generate the cooridnates of the ilth layer
            a = a - W*(il-1)
            b = b - W*(il-1)
            particle_loc = ellipse_points(a,b,np)
            ! print *, particle_loc 
            do ip = 1,np 
                c%layers(il)%boundary(ip)%x = origin%x + particle_loc(ip,1)
                c%layers(il)%boundary(ip)%y = origin%y + particle_loc(ip,2)
            end do
        end do

        ! ! Assign locations to the particles of the cilia
        ! do il = 1,nl
        !     c%layers(il)%t = 'c' ! Closed boundary type
        !     do ip = 1,np 
        !         c%layers(il)%boundary(ip)%x = origin%x + (a - W*(il-1)) * cos((ip-1)*dtheta)
        !         c%layers(il)%boundary(ip)%y = origin%y + (b - W*(il-1)) * sin((ip-1)*dtheta)
        !     end do
        ! end do

        contains

        function ellipse_points(r1,r2,n) result(coordinates)
        !
        !   Purpose:
        !       To provide equally distributed points along the circumference of an ellipse
        !
            integer(int32), intent(in) :: n
            real(real64), intent(in) :: r1, r2


            real(real64) :: coordinates(n,2)
            real(real64) :: theta, circ, dpt, run
            real(real64) :: deltaTheta, subIntegral
            integer(int32) :: numIntegrals, nextPoint

            integer(int32) :: i

            nextPoint = 0
            run = 0.0d0
            theta = 0.0d0
            deltaTheta = 0.0001
            numIntegrals = nint(2*PI/deltaTheta)
            circ = 0.0d0
            coordinates = 0.0d0

            ! Calculate the circumference
            do i = 1,numIntegrals
                theta = theta + deltaTheta
                dpt = computeDpt(r1,r2,theta)
                circ = circ + dpt
            end do

            ! print *, 'circumference = ', circ
            ! print *, 'numIntegrals = ', numIntegrals

            theta = 0.0d0
            ! Generate coordinates
            do i = 1,numIntegrals
                theta = theta + deltaTheta
                subIntegral = n*run/circ
                if (subIntegral.ge.nextPoint) then
                    ! print *, 'hi'
                    coordinates(nextPoint+1,1) = r1 * cos(theta)
                    coordinates(nextPoint+1,2) = r2 * sin(theta)
                    nextPoint = nextPoint + 1
                end if
                run = run + computeDpt(r1,r2,theta)
            end do

        end function ellipse_points

        function computeDpt(r1,r2,theta)
            real(real64), intent(in) :: r1, r2, theta
            real(real64) :: computeDpt

            computeDpt = sqrt((r1*sin(theta))**2+(r2*cos(theta))**2)

        end function computeDpt

    end subroutine create_closed_loop_ellipse

    subroutine create_closed_loop_array(CA,W,r,ar,origin)
    !
    !   Purpose:
    !       To create an array of closed loop structures.
    !       Currently an array consisting of only 1 closed loop structure is supported.
    !
        class(cilia_array),   intent(in out)    :: CA
        class(vec),     intent(in)              :: origin
        real(real64),   intent(in)              :: W          ! Width of cilia
        real(real64),   intent(in)              :: r          ! Radius of the loop
        real(real64),   intent(in)              :: ar         ! Aspect Ratio

        integer(int32) :: ic

        ! Assign locations to the particles of the cilia
        do ic = 1,CA%nc
            ! call create_closed_loop(CA%array(ic),W,r,origin)
            call create_closed_loop_ellipse(CA%array(ic),W,r,ar,origin)
        end do
    end subroutine create_closed_loop_array

    subroutine calculate_closed_loop_array_force(CA,ko,kd,Rlv,Rlh,Rld,t)
    !
    !   Purpose:
    !       To calculate the forces in the links in a each closed loop of the array of closed loops
    !
        class(cilia_array), intent(in out)  :: CA      ! Cilia array
        real(real64),       intent(in)      :: ko      ! Horizontal spring stiffness
        real(real64),       intent(in)      :: kd      ! Diagonal spring stiffness
        real(real64),       intent(in)      :: Rlv,Rlh,Rld      ! Resting length
        real(real64),       intent(in)      :: t       ! Time instant

        integer(int32) :: ic

        ! Calculate forces for all the cilia within the array
        do ic = 1,CA%nc
            call calculate_closed_loop_force(CA%array(ic),ko,kd,Rlv,Rlh,Rld,t)
        end do
    end subroutine calculate_closed_loop_array_force

    subroutine calculate_closed_loop_force(C,ko,kd,Rlv,Rlh,Rld,t)
    !
    !   Purpose:
    !       To calculate the forces in the links in a each immersed boundary of a closed loop
    !
        class(cilia), intent(in out)    :: C         ! Cilia structure
        real(real64), intent(in)        :: ko        ! Horizontal spring stiffness
        real(real64), intent(in)        :: kd        ! Diagonal spring stiffness
        real(real64), intent(in)        :: Rlv,Rlh,Rld        ! Resting length
        real(real64), intent(in)        :: t       ! Time instant

        integer(int32)  :: il, ip

        ! Initialize the forces to zero on all the nodes at every time step
        do il = 1,C%nl 
            do ip = 1,C%np
                C%layers(il)%boundary(ip)%Fx = 0.0d0
                C%layers(il)%boundary(ip)%Fy = 0.0d0
            end do
        end do

        ! Calculate forces on the layers
        call calculate_spring_force(C%layers(1),ko,kd,Rlv,C%layers(1)%t,t)
        call calculate_spring_force(C%layers(2),ko,kd,0.5d0*Rlv,C%layers(2)%t,t)

        ! do il = 1,C%nl
            ! call calculate_spring_force(C%layers(il),ko,kd,Rlv,C%layers(il)%t)
        ! end do

        ! Calculate forces on the horizontal links
        call calculate_horizontal_link_force(C%layers(1),C%layers(2),ko,Rlh,t)

        ! Calculate forces on the diagonal links (negative slope)
        call calculate_diagonal_link_force(C%layers(1),C%layers(2),kd,Rld,t)
        call calculate_diagonal_link_force_pos(C%layers(2),C%layers(1),kd,Rld,t)
    end subroutine calculate_closed_loop_force



end module mod_closed_cilia