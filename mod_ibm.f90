module mod_ibm
    use iso_fortran_env, only: int32, real32, real64
    use mod_mesh
    use mod_ib
    use mod_vec
    implicit none
   
    real(real64), parameter :: PI = 3.141592653589793

    private
    public :: initialize_ib, update_ib, spread_force, interpolate_velocity, & 
              write_location, calculate_spring_force, calculate_torsional_spring_force, &
              create_structure
contains
   
    subroutine initialize_ib(B)
        ! Creates and Initializes an immersed boundary structure
        class(ib), intent(in out) :: B

        integer(int32) :: np, inp

        np = size(B%boundary)
        do concurrent (inp = 1:np)
            ! To-do: Implement some function defining a structure
            ! Location
            B%boundary(inp)%x = 5.0d0
            B%boundary(inp)%y = 5.0d0
            ! Forces
            B%boundary(inp)%Fx = 0.1d0
            B%boundary(inp)%Fy = 0.0d0
            ! Velocity
            B%boundary(inp)%Ux = 0.0d0
            B%boundary(inp)%Uy = 0.0d0
        end do 
    end subroutine initialize_ib

    subroutine update_ib(B,dt)
        class(ib), intent(in out) :: B
        real(real32), intent(in) :: dt

        integer(int32) :: np, inp

        np = size(B%boundary)
        
        ! Calculate the new position
        do concurrent (inp = 1:np)
            B%boundary(inp)%x = B%boundary(inp)%x + B%boundary(inp)%Ux * dt
            B%boundary(inp)%y = B%boundary(inp)%y + B%boundary(inp)%Uy * dt
        end do

    end subroutine update_ib

    subroutine spread_force(M,B,Fx,Fy)
        class(mesh), intent(in) :: M
        class(ib), intent(in) :: B
        real(real64), intent(in out) :: Fx(:,:), Fy(:,:)

        ! Define variables to store locations
        real(real64) :: Lx, Ly      ! Lagrangian locations
        real(real64) :: Ex, Ey      ! Eulerian locations
        real(real64) :: Flx, Fly    ! Forces at Lagrangian locations
        integer(int32) :: np ! Number of Lagrangian particles

        ! real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)

        ! Indices
        integer(int32) :: i, j, inp
        np = size(B%boundary)


        ! Initialize the forces at every time-step
        Fx = 0.0d0
        Fy = 0.0d0

        ! Calculate the x-direction force on the u-velocity cells
        ! Iterating over all the grid points including the boundary values
        do j = M%yu%lb,M%yu%ub
            do i = M%xu%lb,M%xu%ub
                Ex = M%u_mesh(i,j)%x
                Ey = M%u_mesh(i,j)%y
                do inp = 1,np
                    Lx = B%boundary(inp)%x
                    Ly = B%boundary(inp)%y
                    Flx = B%boundary(inp)%Fx
                
                    Fx(i,j) = Fx(i,j) + Flx * dirac( [(Ex-Lx), (Ey-Ly)], M%dx)
                end do
            end do
        end do

        ! Calculate the y-direction force on the v-velocity cells
        do j = M%yv%lb,M%yv%ub
            do i = M%xv%lb,M%xv%ub
                Ex = M%v_mesh(i,j)%x
                Ey = M%v_mesh(i,j)%y
                do inp = 1,np
                    Lx = B%boundary(inp)%x
                    Ly = B%boundary(inp)%y
                    Fly = B%boundary(inp)%Fy

                    Fy(i,j) = Fy(i,j) + Fly * dirac( [(Ex-Lx), (Ey-Ly)], M%dy)
                end do
            end do
        end do

        contains

        function dirac(x,h)
            ! Defined for a uniform grid
            real(real64), intent(in) :: x(2)
            real(real32), intent(in) :: h
            real(real64) :: dirac

            integer(int32) :: ii
            real(real64) :: phi, r

            dirac = 1.0d0
            do ii = 1,2
                r = x(ii)/h 
                if (abs(r).le.2) then
                    phi = 0.25d0 * (1 + cos(PI*r/2))
                else 
                    phi = 0.0d0
                end if
            
                dirac = dirac*phi
            end do        

            dirac = (1/h**2) * dirac

        end function dirac

    end subroutine spread_force

    subroutine interpolate_velocity(M,B,u,v)
        class(mesh), intent(in) :: M
        class(ib), intent(in out) :: B
        real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)

        ! Define variables to store locations
        real(real64) :: Lx, Ly      ! Lagrangian locations
        real(real64) :: Ex, Ey      ! Eulerian locations
        real(real64) :: UL, VL      ! Velocity at Lagrangian locations
        real(real64) :: UE, VE      ! Velocity at Eulerian locations
        integer(int32) :: np        ! Number of Lagrangian particles
        
        ! Indices
        integer(int32) :: i, j, inp
        np = size(B%boundary)

        ! Initialize the velocity at every time-step
        do concurrent (inp = 1:np)
            B%boundary(inp)%Ux = 0.0d0
            B%boundary(inp)%Uy = 0.0d0
        end do

        ! Calculate the u velocity of Lagrangian points
        ! Iterating over all the grid points including the boundary values for u-velocity cells
        do j = M%yu%lb,M%yu%ub
            do i = M%xu%lb,M%xu%ub
                Ex = M%u_mesh(i,j)%x ! u-cell x location
                Ey = M%u_mesh(i,j)%y ! u-cell y location
                UE = u(i,j)
                do inp = 1,np
                    Lx = B%boundary(inp)%x
                    Ly = B%boundary(inp)%y
                
                    B%boundary(inp)%Ux = B%boundary(inp)%Ux + UE * dirac( [(Ex-Lx), (Ey-Ly)], M%dx) * M%dx**2
                end do
            end do
        end do

        ! Calculate the v velocity of Lagrangian points
        ! Iterating over all the grid points including the boundary values for v-velocity cells
        do j = M%yv%lb,M%yv%ub
            do i = M%xv%lb,M%xv%ub
                Ex = M%v_mesh(i,j)%x ! v-cell x location
                Ey = M%v_mesh(i,j)%y ! v-cell y location
                VE = v(i,j)
                do inp = 1,np
                    Lx = B%boundary(inp)%x
                    Ly = B%boundary(inp)%y
                
                    B%boundary(inp)%Uy = B%boundary(inp)%Uy + VE * dirac( [(Ex-Lx), (Ey-Ly)], M%dy) * M%dy**2
                end do
            end do
        end do

        contains 

        function dirac(x,h)
            ! Defined for a uniform grid
            real(real64), intent(in) :: x(2)
            real(real32), intent(in) :: h
            real(real64) :: dirac

            integer(int32) :: ii
            real(real64) :: phi, r

            dirac = 1.0d0
            do ii = 1,2
                r = x(ii)/h 
                if (abs(r).le.2) then
                    phi = 0.25d0 * (1 + cos(PI*r/2))
                else 
                    phi = 0.0d0
                end if
            
                dirac = dirac*phi
            end do        

            dirac = (1/h**2) * dirac

        end function dirac

    end subroutine interpolate_velocity

    subroutine write_location(B,timestep)
        class(ib), intent(in) :: B
        integer(int32), intent(in) :: timestep        

        integer(int32) :: fileunit = 8
        character(len=:), allocatable :: filename
        integer(int32) :: inp, np
        character(len=8) :: itnumber

        write(itnumber,"(I8.8)") timestep

        np = size(B%boundary)

        filename = 'ib_loc' // itnumber // '.txt'
        open(unit=fileunit, file=filename, ACTION="write", STATUS="replace")
        do inp = 1,np
            write(fileunit, '((F14.7), (F14.7))') B%boundary(inp)%x, B%boundary(inp)%y
        end do
        close(fileunit)

        ! filename = 'ib_y.txt'
        ! open(unit=fileunit, file=filename, ACTION="write", STATUS="replace")
        ! do j = M%yu%lb,M%yu%ub
        !     write(fileunit, '(*(F14.7))')(M%u_mesh(i,j)%y , i = M%xu%lb,M%xu%ub)
        ! end do
        ! close(fileunit)
    end subroutine write_location

    subroutine calculate_spring_force(B,ks,Rl,t)
        class(ib), intent(in out) :: B    ! Immersed boundary
        real(real64), intent(in)  :: ks   ! Spring stiffness
        real(real64),intent(in)   :: Rl   ! Resting length
        character(1), intent(in)  :: t    ! Boundary type (Open or Closed) 

        integer(int32)  :: np   ! Number of nodes/particle
        integer(int32)  :: lastp ! Last particle while calculating forces
        integer(int32)  :: inp, master, slave  ! Indices
        real(real64)    :: d    ! Distance between two nodes

        real(real64) :: Fmx, Fslx, Fmy, Fsly ! Forces (Master(m) and Slave(sl) node)
        real(real64) :: xm, xsl, ym, ysl ! Location of master and slave nodes

        ! Count the number of nodes/particles in the immersed boundary
        np = size(B%boundary)

        ! Initialize the forces to zero on all the nodes at every time step
        do concurrent (inp = 1:np)
            B%boundary(inp)%Fx = 0.0d0
            B%boundary(inp)%Fy = 0.0d0
        end do

        ! Calculate the forces on each node
        if (t.eq.'o') then
            lastp = np-1
        elseif (t.eq.'c') then 
            lastp = np
        else 
            print *, "Wrong type of Immersed Boundary: Enter 'o' or 'c' "
        end if

        do master = 1,lastp
            if (master == np) then ! This condition will never be satisfied for open case 
                slave = 1
            else 
                slave = master+1
            end if
            ! Master node location
            xm = B%boundary(master)%x
            ym = B%boundary(master)%y

            ! Slave node location
            xsl = B%boundary(slave)%x
            ysl = B%boundary(slave)%y

            ! Calculate distance between master and slave nodes
            d = norm2([(xsl-xm),(ysl-ym)])
            
            ! Calculate forces (Master node)
            Fmx = ks*(1.0d0-Rl/d)*(xsl-xm)
            Fmy = ks*(1.0d0-Rl/d)*(ysl-ym)

            ! Calculate forces (Slave node)
            Fslx = -Fmx
            Fsly = -Fmy

            ! Assign fores to master node
            B%boundary(master)%Fx = B%boundary(master)%Fx + Fmx
            B%boundary(master)%Fy = B%boundary(master)%Fy + Fmy

            ! Assign forces to slave node
            B%boundary(slave)%Fx = B%boundary(slave)%Fx + Fslx
            B%boundary(slave)%Fy = B%boundary(slave)%Fy + Fsly
        end do

    end subroutine calculate_spring_force

    subroutine calculate_torsional_spring_force(B,kb,theta,t)
        class(ib), intent(in out) :: B ! Immersed boundary
        real(real64), intent(in) :: kb ! Torsional spring stiffness 
        real(real64), intent(in) :: theta ! Desired angle
        character(1), intent(in)  :: t    ! Boundary type (Open or Closed) 

        integer(int32)  :: np               ! Number of nodes/particle
        integer(int32)  :: firstp, lastp    ! First and last particle while calculating forces
        real(real64)    :: dLM, dMR         ! Distance between master and the left and right nodes
        real(real64)    :: C                ! Curvature
        integer(int32)  :: inp, master, left, right  ! Indices

        real(real64) :: Fmx, Fmy ! Forces on master nodes only due to bending
        real(real64) :: xm, ym, xl, yl, xr, yr ! Location of master and left and right nodes

        ! Count the number of nodes/particles in the immersed boundary
        np = size(B%boundary)

        !!!!!! IMPORTANT !!!!!!
        ! ! Initialize the forces to zero on all the nodes at every time step
        ! do concurrent (inp = 1:np)
        !     B%boundary(inp)%Fx = 0.0d0
        !     B%boundary(inp)%Fy = 0.0d0
        ! end do
        ! Since this subroutine is called after calculation of spring forces,
        ! there is no need to initialize the forces to zero again
        if (t.eq.'o') then
            firstp = 2
            lastp = np-1
        elseif (t.eq.'c') then 
            firstp = 1
            lastp = np
        else 
            print *, "Wrong type of Immersed Boundary: Enter 'o' or 'c' "
        end if
        do master = firstp,lastp
            if (master == 1) then 
                left    = np
                right   = master + 1
            elseif (master == np) then
                left    = master - 1
                right   = 1
            else
                left    = master - 1
                right   = master + 1
            end if 

            ! Master node location
            xm = B%boundary(master)%x
            ym = B%boundary(master)%y

            ! Left node location
            xl = B%boundary(left)%x
            yl = B%boundary(left)%y

            ! Right node location
            xr = B%boundary(right)%x
            yr = B%boundary(right)%y

            ! Calculate the distance 
            ! Master node and the left node
            dLM = norm2([(xl-xm),(yl-ym)])
            ! Master node and the left node
            dMR = norm2([(xr-xm),(yr-ym)])

            ! Calculate curvature 
            C = dLM*dMR*sin(theta)

            ! Calculate force on the master node
            Fmx = kb*((xr-xm)*(ym-yl) - (yr-ym)*(xm-xl) - C) * &
                  ((ym-yl) + (yr-ym))
            Fmy = kb*((xr-xm)*(ym-yl) - (yr-ym)*(xm-xl) - C) * &
                  (-(xr-xm) - (xm-xl))

            ! Assign force to master node
            B%boundary(master)%Fx = B%boundary(master)%Fx + Fmx
            B%boundary(master)%Fy = B%boundary(master)%Fy + Fmy
        end do
    end subroutine calculate_torsional_spring_force

    subroutine create_structure(B,origin,L,t)
        class(ib), intent(in out)   :: B
        class(vec), intent(in)      :: origin
        real(real64), intent(in)    :: L
        character(1), intent(in)    :: t

        real(real64) :: ds
        integer(int32) :: np, inp
        
        np = size(B%boundary)
        ds = L/(np-1)

        ! Create a line (at 45 degree)
        B%boundary(1)%x = origin%x
        B%boundary(1)%y = origin%y
        do inp = 2,np
            B%boundary(inp)%x = B%boundary(inp-1)%x + ds
            B%boundary(inp)%y = B%boundary(inp-1)%y + ds
        end do

        do inp = 1,np
            print *, B%boundary(inp)%x
            print *, B%boundary(inp)%y
        end do

    end subroutine create_structure

end module mod_ibm
