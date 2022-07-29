module mod_ibm
    use iso_fortran_env, only: int32, real64, real64
    use mod_mesh
    use mod_ib
    use mod_vec
    use mod_cilia
    use mod_cilia_array
    use mod_particle
    implicit none
   
    real(real64), parameter :: PI = 3.141592653589793
    integer(int32) :: switchcounter = 0

    private
    public :: initialize_ib, update_ib, spread_force, interpolate_velocity, & 
              write_location, calculate_spring_force, calculate_torsional_spring_force, &
              create_structure, create_cilia, calculate_horizontal_link_force, &
              calculate_diagonal_link_force, calculate_cilia_force, create_cilia_array, &
              calculate_cilia_array_force, update_cilia, update_cilia_array, spread_force_cilia, &
              spread_force_cilia_array, interpolate_velocity_cilia, &
              interpolate_velocity_cilia_array, initialize_velocity_cilia_array, &
              apply_tip_force_cilia, apply_tip_force_cilia_array, write_location_cilia, &
              write_location_cilia_force, write_location_cilia_velocity, calculate_diagonal_link_force_pos, &
              store_original_locations

contains
   
    subroutine initialize_ib(B)
        ! Creates and Initializes an immersed boundary structure
        class(ib), intent(in out) :: B

        integer(int32) :: np, inp

        np = size(B%boundary)
        do inp = 1,np
            ! To-do: Implement some function defining a structure
            ! Location
            B%boundary(inp)%x = 0.0d0
            B%boundary(inp)%y = 0.0d0
            ! Forces
            B%boundary(inp)%Fx = 0.0d0
            B%boundary(inp)%Fy = 0.0d0
            ! Velocity
            B%boundary(inp)%Ux = 0.0d0
            B%boundary(inp)%Uy = 0.0d0
        end do 
    end subroutine initialize_ib

    subroutine update_ib(B,dt)
        class(ib), intent(in out) :: B
        real(real64), intent(in) :: dt

        integer(int32) :: firstp, np, inp

        np = size(B%boundary)
       
        if (B%t.eq.'o') then
            ! Start from 2 if fixing the first particle in each layer of cilia
            firstp = 2 ! Represents cilia fixed at their base
        else if (B%t.eq.'c') then
            firstp = 1 ! Represents particles which aren't fixed
        end if
        
        ! Calculate the new position
        do inp = firstp,np
            B%boundary(inp)%x = B%boundary(inp)%x + B%boundary(inp)%Ux * dt
            B%boundary(inp)%y = B%boundary(inp)%y + B%boundary(inp)%Uy * dt
        end do

    end subroutine update_ib

    subroutine spread_force_compact(M,B,Fx,Fy)
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
        integer(int32) :: i, j, inp, LeftIndex, RightIndex, BottomIndex, TopIndex, StencilSize

        np = size(B%boundary)
        StencilSize = 4; ! Only on one side

        ! Calculate the x-direction force on the u-velocity cells
        ! Iterating over all the grid points including the boundary values
        do inp = 1,np
            Lx = B%boundary(inp)%x
            Ly = B%boundary(inp)%y
            Flx = B%boundary(inp)%Fx
            ! Make sure that the following indices lie inside the domain            
            LeftIndex   = 0 + M%xu%lb + floor((Lx-StencilSize*M%dx)/M%dx) 
            BottomIndex = 1 + M%yu%lb + floor((Ly-StencilSize*M%dy)/M%dy)   ! Increment by 1 since the first value is negative and not zero

                    ! print *, floor((Ly-StencilSize*M%dy)/M%dy)
            RightIndex  = LeftIndex     +   2*StencilSize
            TopIndex    = BottomIndex   +   2*StencilSize

            ! print *, 'Left index = ', LeftIndex
            ! print *, 'Right index = ', RightIndex
            ! print *, 'Top index = ', TopIndex
            ! print *, 'Bottom index = ', BottomIndex

            do j = BottomIndex,TopIndex
                do i = LeftIndex,RightIndex
                    Ex = M%u_mesh(i,j)%x
                    Ey = M%u_mesh(i,j)%y
                    
                    Fx(i,j) = Fx(i,j) + Flx * dirac( [(Ex-Lx), (Ey-Ly)], M%dx)
                end do
            end do
        end do

        ! Calculate the y-direction force on the v-velocity cells
        do inp = 1,np
            Lx = B%boundary(inp)%x
            Ly = B%boundary(inp)%y
            Fly = B%boundary(inp)%Fy

            LeftIndex   = 1 + M%xv%lb + floor((Lx-StencilSize*M%dx)/M%dx) ! Increment by 1 since the first value is negative and not zero 
            BottomIndex = 0 + M%yv%lb + floor((Ly-StencilSize*M%dy)/M%dy)

            RightIndex  = LeftIndex     +   2*StencilSize
            TopIndex    = BottomIndex   +   2*StencilSize

            do j = BottomIndex,TopIndex
                do i = LeftIndex,RightIndex
                    Ex = M%v_mesh(i,j)%x
                    Ey = M%v_mesh(i,j)%y
                    
                    Fy(i,j) = Fy(i,j) + Fly * dirac( [(Ex-Lx), (Ey-Ly)], M%dy)
                end do
            end do
        end do

        contains

        pure function dirac(x,h)
            ! Defined for a uniform grid
            real(real64), intent(in) :: x(2)
            real(real64), intent(in) :: h
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

    end subroutine spread_force_compact

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
            real(real64), intent(in) :: h
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

    subroutine interpolate_velocity_compact(M,B,u,v)
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
        integer(int32) :: i, j, inp, LeftIndex, RightIndex, BottomIndex, TopIndex, StencilSize

        np = size(B%boundary)
        StencilSize = 4

        ! Calculate the u velocity of Lagrangian points
        ! Iterating over all the grid points including the boundary values for u-velocity cells
        do inp = 1,np
            Lx = B%boundary(inp)%x
            Ly = B%boundary(inp)%y
        
            ! Make sure that the following indices lie inside the domain            
            LeftIndex   = 0 + M%xu%lb + floor((Lx-StencilSize*M%dx)/M%dx) 
            BottomIndex = 1 + M%yu%lb + floor((Ly-StencilSize*M%dy)/M%dy)   ! Increment by 1 since the first value is negative and not zero

            RightIndex  = LeftIndex     +   2*StencilSize
            TopIndex    = BottomIndex   +   2*StencilSize

            do j = BottomIndex,TopIndex
                do i = LeftIndex,RightIndex
                    Ex = M%u_mesh(i,j)%x ! u-cell x location
                    Ey = M%u_mesh(i,j)%y ! u-cell y location
                    UE = u(i,j)
                    B%boundary(inp)%Ux = B%boundary(inp)%Ux + UE * dirac( [(Ex-Lx), (Ey-Ly)], M%dx) * M%dx**2
                end do
            end do
        end do

        ! Calculate the v velocity of Lagrangian points
        ! Iterating over all the grid points including the boundary values for v-velocity cells
        
        do inp = 1,np
            Lx = B%boundary(inp)%x
            Ly = B%boundary(inp)%y
       
            LeftIndex   = 1 + M%xv%lb + floor((Lx-StencilSize*M%dx)/M%dx) ! Increment by 1 since the first value is negative and not zero 
            BottomIndex = 0 + M%yv%lb + floor((Ly-StencilSize*M%dy)/M%dy)

            RightIndex  = LeftIndex     +   2*StencilSize
            TopIndex    = BottomIndex   +   2*StencilSize

            do j = BottomIndex,TopIndex
                do i = LeftIndex,RightIndex
                    Ex = M%v_mesh(i,j)%x ! v-cell x location
                    Ey = M%v_mesh(i,j)%y ! v-cell y location
                    VE = v(i,j)
                    B%boundary(inp)%Uy = B%boundary(inp)%Uy + VE * dirac( [(Ex-Lx), (Ey-Ly)], M%dy) * M%dy**2
                end do
            end do
        end do

        contains 

        function dirac(x,h)
            ! Defined for a uniform grid
            real(real64), intent(in) :: x(2)
            real(real64), intent(in) :: h
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

    end subroutine interpolate_velocity_compact

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
            real(real64), intent(in) :: h
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

    subroutine write_location_cilia_force(CA,timestep,filetype)
        class(cilia_array), intent(in) :: CA
        integer(int32), intent(in) :: timestep        
        character(len=1), intent(in) :: filetype

        integer(int32) :: fileunit = 8
        character(len=:), allocatable :: filename
        integer(int32) :: inp, np, nl, il, ic
        character(len=8) :: itnumber
        character(len=1) :: idl ! Layer id
        character(len=1) :: idc ! Cilia id
        character(len=8) :: file_status
        character(len=3) :: file_advance
        write(itnumber,"(I8.8)") timestep

        file_status = "replace"
        file_advance = "no"

        filename = 'force_ib_loc_' // filetype // itnumber // '.txt'
        
        open(unit=fileunit, file=filename, ACTION="write", Position="Append", & 
                        STATUS=trim(file_status))
        ! Loop through the cilia
        do ic = 1,CA%nc
            nl = CA%array(ic)%nl
            np = CA%array(ic)%np

            ! Loop through the layers of a cilia
            do il = 1,nl
                
                ! Loop through the particles in a cilia layer to write the x position
                do inp = 1,np
                    if (inp.eq.np) then 
                        file_advance='yes'
                    else 
                        file_advance='no'
                    end if
                    write(fileunit, '(F14.7)',ADVANCE=file_advance) CA%array(ic)%layers(il)%boundary(inp)%Fx
                end do
                
                ! Loop through the particles in a cilia layer to write the y position
                do inp = 1,np
                    if (inp.eq.np) then 
                        file_advance='yes'
                    else 
                        file_advance='no'
                    end if
                    write(fileunit, '(F14.7)',ADVANCE=file_advance) CA%array(ic)%layers(il)%boundary(inp)%Fy
                end do
            end do
        end do
        close(fileunit)

    end subroutine write_location_cilia_force

    subroutine write_location_cilia_velocity(CA,timestep,filetype)
        class(cilia_array), intent(in) :: CA
        integer(int32), intent(in) :: timestep        
        character(len=1), intent(in) :: filetype

        integer(int32) :: fileunit = 8
        character(len=:), allocatable :: filename
        integer(int32) :: inp, np, nl, il, ic
        character(len=8) :: itnumber
        character(len=1) :: idl ! Layer id
        character(len=1) :: idc ! Cilia id
        character(len=8) :: file_status
        character(len=3) :: file_advance
        write(itnumber,"(I8.8)") timestep

        file_status = "replace"
        file_advance = "no"

        filename = 'vel_ib_loc_' // filetype // itnumber // '.txt'
        
        open(unit=fileunit, file=filename, ACTION="write", Position="Append", & 
                        STATUS=trim(file_status))
        ! Loop through the cilia
        do ic = 1,CA%nc
            nl = CA%array(ic)%nl
            np = CA%array(ic)%np

            ! Loop through the layers of a cilia
            do il = 1,nl
                
                ! Loop through the particles in a cilia layer to write the x position
                do inp = 1,np
                    if (inp.eq.np) then 
                        file_advance='yes'
                    else 
                        file_advance='no'
                    end if
                    write(fileunit, '(F14.7)',ADVANCE=file_advance) CA%array(ic)%layers(il)%boundary(inp)%Ux
                end do
                
                ! Loop through the particles in a cilia layer to write the y position
                do inp = 1,np
                    if (inp.eq.np) then 
                        file_advance='yes'
                    else 
                        file_advance='no'
                    end if
                    write(fileunit, '(F14.7)',ADVANCE=file_advance) CA%array(ic)%layers(il)%boundary(inp)%Uy
                end do
            end do
        end do
        close(fileunit)

    end subroutine write_location_cilia_velocity

    subroutine write_location_cilia(CA,timestep,filetype)
        class(cilia_array), intent(in) :: CA
        integer(int32), intent(in) :: timestep        
        character(len=1), intent(in) :: filetype

        integer(int32) :: fileunit = 8
        character(len=:), allocatable :: filename
        integer(int32) :: inp, np, nl, il, ic
        character(len=8) :: itnumber
        character(len=1) :: idl ! Layer id
        character(len=1) :: idc ! Cilia id
        character(len=8) :: file_status
        character(len=3) :: file_advance
        write(itnumber,"(I8.8)") timestep

        file_status = 'replace'
        file_advance = "no"

        filename = 'ib_loc_' // filetype // itnumber // '.txt'
        
        open(unit=fileunit, file=filename, ACTION="write", Position="Append", & 
                        STATUS=trim(file_status))
        ! Loop through the cilia
        do ic = 1,CA%nc
            nl = CA%array(ic)%nl
            np = CA%array(ic)%np

            ! Loop through the layers of a cilia
            do il = 1,nl
                
                ! Loop through the particles in a cilia layer to write the x position
                do inp = 1,np
                    if (inp.eq.np) then 
                        file_advance='yes'
                    else 
                        file_advance='no'
                    end if
                    write(fileunit, '(F14.7)',ADVANCE=file_advance) CA%array(ic)%layers(il)%boundary(inp)%x
                end do
                
                ! Loop through the particles in a cilia layer to write the y position
                do inp = 1,np
                    if (inp.eq.np) then 
                        file_advance='yes'
                    else 
                        file_advance='no'
                    end if
                    write(fileunit, '(F14.7)',ADVANCE=file_advance) CA%array(ic)%layers(il)%boundary(inp)%y
                end do
            end do
        end do
        close(fileunit)

    end subroutine write_location_cilia

    subroutine write_location(CA,timestep)
        class(cilia_array), intent(in) :: CA
        integer(int32), intent(in) :: timestep        

        integer(int32) :: fileunit = 8
        character(len=:), allocatable :: filename
        integer(int32) :: inp, np, nl, il, ic
        character(len=8) :: itnumber
        character(len=1) :: idl ! Layer id
        character(len=1) :: idc ! Cilia id
        write(itnumber,"(I8.8)") timestep

        do ic = 1,CA%nc
            nl = CA%array(ic)%nl
            np = CA%array(ic)%np
            write(idc,"(I1.1)") ic
            do il = 1,nl
                write(idl,"(I1.1)") il
                filename = idc // '_' // idl // '_ib_loc' // itnumber // '.txt'
                open(unit=fileunit, file=filename, ACTION="write", STATUS="replace")
                do inp = 1,np
                    write(fileunit, '((F14.7), (F14.7))') CA%array(ic)%layers(il)%boundary(inp)%x, & 
                                    CA%array(ic)%layers(il)%boundary(inp)%y
                end do
                close(fileunit)

                ! filename = 'ib_y.txt'
                ! open(unit=fileunit, file=filename, ACTION="write", STATUS="replace")
                ! do j = M%yu%lb,M%yu%ub
                !     write(fileunit, '(*(F14.7))')(M%u_mesh(i,j)%y , i = M%xu%lb,M%xu%ub)
                ! end do
                ! close(fileunit)
            end do
        end do
    end subroutine write_location

    subroutine calculate_spring_force(B,ko,kd,Rlo,t,ti,ilayer)
        class(ib), intent(in out) :: B      ! Immersed boundary
        real(real64), intent(in)  :: ko     ! Orthogonal Spring stiffness
        real(real64), intent(in)  :: kd     ! Diagonal Spring stiffness
        real(real64),intent(in)   :: Rlo    ! Resting length
        character(1), intent(in)  :: t      ! Boundary type (Open or Closed) 
        real(real64), intent(in)  :: ti     ! Time instant
        integer(int32), intent(in):: ilayer ! Time instant

        integer(int32)  :: np   ! Number of nodes/particle
        integer(int32)  :: lastp ! Last particle while calculating forces
        integer(int32)  :: master, slave  ! Indices
        real(real64), parameter :: beta = 0.5d0 ! Maximum value of coefficient to be multiplied to resting lengths
        real(real64) :: nbeta ! Value of coefficient to be multiplied to each resting length

        type(vec) :: F

        ! Count the number of nodes/particles in the immersed boundary
        np = size(B%boundary)

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

            ! Calculate the value of coefficient for each resting length
            if ((sin(4*PI*ti)).gt.0) then

                ! Reset the value of switchcounter
                switchcounter = 0

                if ((ilayer == 1)) then
                        nbeta = -beta/(np-1) * (master-1) + 1
                    else
                        nbeta = 1.0d0
                endif

            else

                ! Increment the value of switchcounter
                switchcounter = switchcounter + 1

                if ((ilayer == 2)) then
                        nbeta = -beta/(np-1) * (master-1) + 1;
                    else
                        nbeta = 1.0d0
                end if
            end if

            ! print *, 'switchcounter=', switchcounter 
            ! Calculate the spring force between master and slave nodes
            F = spring_force(B%boundary(master),B%boundary(slave),nbeta*ko,ti,nbeta)

            ! Assign fores to master node
            B%boundary(master)%Fx = B%boundary(master)%Fx + F%x
            B%boundary(master)%Fy = B%boundary(master)%Fy + F%y

            ! Assign forces to slave node
            B%boundary(slave)%Fx = B%boundary(slave)%Fx - F%x
            B%boundary(slave)%Fy = B%boundary(slave)%Fy - F%y
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

    subroutine create_cilia(C,W,dp,origin)
        ! Creates a vertical cilia
        class(cilia), intent(in out) :: C
        class(vec), intent(in) :: origin
        real(real64), intent(in) :: W ! Width of cilia
        real(real64), intent(in) :: dp ! Spacing of particles in a cilia layer

        integer(int32) :: nl
        integer(int32) :: np

        integer(int32) :: il,ip

        nl = C%nl
        np = C%np

        ! Assign locations to each of the layers of the cilia
        do il = 1,C%nl
            call initialize_ib(C%layers(il))
            
            C%layers(il)%boundary(1)%x = origin%x + W*(il-1)
            C%layers(il)%boundary(1)%y = origin%y
            do ip = 2,C%np
                C%layers(il)%boundary(ip)%x = C%layers(il)%boundary(1)%x
                C%layers(il)%boundary(ip)%y = C%layers(il)%boundary(ip-1)%y + dp
            end do
        end do

        ! do il = 1,C%nl
        !    do ip = 1,np
        !        print *, C%layers(il)%boundary(ip)%x
        !        print *, C%layers(il)%boundary(ip)%y
        !    end do
        ! end do
    end subroutine create_cilia

    subroutine calculate_cilia_force(C,ko,kd,Rl,ti)
        class(cilia), intent(in out) :: C       ! Cilia structure
        real(real64), intent(in)  :: ko         ! Horizontal spring stiffness
        real(real64), intent(in)  :: kd         ! Diagonal spring stiffness
        real(real64), intent(in)   :: Rl        ! Resting length
        real(real64),       intent(in)      :: ti       ! Time instant

        character(1)  :: t = 'o'  ! Boundary type (Open for cilia)
        integer(int32) :: il, ip

        ! Initialize the forces to zero on all the nodes at every time step
        do il = 1,C%nl 
            do ip = 1,C%np
                C%layers(il)%boundary(ip)%Fx = 0.0d0
                C%layers(il)%boundary(ip)%Fy = 0.0d0
            end do
        end do

        ! Calculate forces on the layers
        do il = 1,C%nl
            call calculate_spring_force(C%layers(il),ko,kd,Rl,t,ti,il)
        end do



        ! Calculate forces on the horizontal links
        call calculate_horizontal_link_force(C%layers(1),C%layers(2),ko,Rl,ti)

        ! Calculate forces on the diagonal links (negative slope)
        call calculate_diagonal_link_force(C%layers(1),C%layers(2),kd,Rl,ti)
        call calculate_diagonal_link_force_pos(C%layers(2),C%layers(1),kd,Rl,ti)

        ! Calculate forces on the diagonal links (Positive slope)
        ! call calculate_diagonal_link_force(C%layers(1),C%layers(2),ks,Rl)

        ! Fix the last particle in each layer (Specify the forces to be zero)
        do il = 1,C%nl
            C%layers(il)%boundary(1)%Fx = 0.0d0 
            C%layers(il)%boundary(1)%Fy = 0.0d0 
        end do

    end subroutine calculate_cilia_force

    subroutine calculate_horizontal_link_force(masterL,slaveL,ks,Rlo,t)
        class(ib), intent(in out) :: masterL, slaveL    ! Immersed boundary
        real(real64), intent(in)  :: ks   ! Spring stiffness
        real(real64), intent(in)  :: Rlo   ! Resting length
        real(real64), intent(in)  :: t       ! Time instant


        real(real64), parameter :: beta = 0.5d0 ! Maximum value of coefficient to be multiplied to resting lengths
        real(real64) :: nbeta ! Value of coefficient to be multiplied to each resting length
        integer(int32)  :: np   ! Number of nodes/particle
        integer(int32)  :: ip   ! Indices


        type(vec) :: F

        ! Count the number of nodes/particles in the immersed boundary layers
        ! (Assuming both the layers have the same number of particles)
        ! np = size(masterL%boundary)
        np = masterL%np

        do ip = 1,np

            ! Calculate the value of coefficient for each resting length
            ! nbeta = -beta/(np-1) * (ip-1) + 1
            nbeta = 1

            ! Calculate the spring force between master and slave nodes
            F = spring_force(masterL%boundary(ip),slaveL%boundary(ip),ks,t,nbeta)

            ! Assign fores to master node
            masterL%boundary(ip)%Fx = masterL%boundary(ip)%Fx + F%x
            masterL%boundary(ip)%Fy = masterL%boundary(ip)%Fy + F%y

            ! Assign forces to slave node
            slaveL%boundary(ip)%Fx = slaveL%boundary(ip)%Fx - F%x
            slaveL%boundary(ip)%Fy = slaveL%boundary(ip)%Fy - F%y
        end do

    end subroutine calculate_horizontal_link_force

    subroutine calculate_diagonal_link_force_pos(masterL,slaveL,kd,Rlo,t)
        class(ib), intent(in out) :: masterL, slaveL    ! Immersed boundary
        real(real64), intent(in)  :: kd   ! Spring stiffness
        real(real64),intent(in)   :: Rlo   ! Resting length
        real(real64),       intent(in)      :: t       ! Time instant

        integer(int32)  :: np   ! Number of nodes/particle
        integer(int32)  :: ip, master, slave   ! Indices

        type(vec) :: F

        ! Count the number of nodes/particles in the immersed boundary layers
        ! (Assuming both the layers have the same number of particles)
        ! np = size(masterL%boundary)
        np = masterL%np

        ! Calculates for negative slope diagonal links
        do ip = 2,np
            master = ip
            slave = master - 1

            ! Calculate the spring force between master and slave nodes
            F = spring_force(masterL%boundary(master),slaveL%boundary(slave),kd,t)
            
            ! Assign fores to master node
            masterL%boundary(master)%Fx = masterL%boundary(master)%Fx + F%x
            masterL%boundary(master)%Fy = masterL%boundary(master)%Fy + F%y

            ! Assign forces to slave node
            slaveL%boundary(slave)%Fx = slaveL%boundary(slave)%Fx - F%x
            slaveL%boundary(slave)%Fy = slaveL%boundary(slave)%Fy - F%y
        end do

        ! For closed loop cilia
        ! The first node of the outer circle connects to the last node of the inner circle
        if (masterL%t.eq.'c') then
            master = 1
            slave = np

            ! Calculate the spring force between master and slave nodes
            F = spring_force(masterL%boundary(master),slaveL%boundary(slave),kd,t)

            ! Assign fores to master node
            masterL%boundary(master)%Fx = masterL%boundary(master)%Fx + F%x
            masterL%boundary(master)%Fy = masterL%boundary(master)%Fy + F%y

            ! Assign forces to slave node
            slaveL%boundary(slave)%Fx = slaveL%boundary(slave)%Fx - F%x
            slaveL%boundary(slave)%Fy = slaveL%boundary(slave)%Fy - F%y
        end if

    end subroutine calculate_diagonal_link_force_pos

    subroutine calculate_diagonal_link_force(masterL,slaveL,kd,Rlo,t)
        class(ib), intent(in out) :: masterL, slaveL    ! Immersed boundary
        real(real64), intent(in)  :: kd   ! Spring stiffness
        real(real64),intent(in)   :: Rlo   ! Resting length
        real(real64),       intent(in)      :: t       ! Time instant

        integer(int32)  :: np   ! Number of nodes/particle
        integer(int32)  :: ip, master, slave   ! Indices
        real(real64)    :: d    ! Distance between two nodes

        real(real64) :: Fmx, Fslx, Fmy, Fsly ! Forces (Master(m) and Slave(sl) node)
        real(real64) :: xm, xsl, ym, ysl ! Location of master and slave nodes

        ! Original resting length
        real(real64)   :: Rl   ! Resting length
        type(vec) :: F

        ! Count the number of nodes/particles in the immersed boundary layers
        ! (Assuming both the layers have the same number of particles)
        np = size(masterL%boundary)

        ! Calculates for negative slope diagonal links
        do ip = 1,np-1
            master = ip
            slave = master + 1

            ! Calculate the spring force between master and slave nodes
            F = spring_force(masterL%boundary(master),slaveL%boundary(slave),kd,t)

            ! Assign fores to master node
            masterL%boundary(master)%Fx = masterL%boundary(master)%Fx + F%x
            masterL%boundary(master)%Fy = masterL%boundary(master)%Fy + F%y

            ! Assign forces to slave node
            slaveL%boundary(slave)%Fx = slaveL%boundary(slave)%Fx - F%x
            slaveL%boundary(slave)%Fy = slaveL%boundary(slave)%Fy - F%y
        end do

        ! For closed loops of cilia
        ! The last node of the outer circle connects to the first node of the inner circle
        if (masterL%t.eq.'c') then
            master = np
            slave = 1
            
            ! Calculate the spring force between master and slave nodes
            F = spring_force(masterL%boundary(master),slaveL%boundary(slave),kd,t)

            ! Assign fores to master node
            masterL%boundary(master)%Fx = masterL%boundary(master)%Fx + F%x
            masterL%boundary(master)%Fy = masterL%boundary(master)%Fy + F%y

            ! Assign forces to slave node
            slaveL%boundary(slave)%Fx = slaveL%boundary(slave)%Fx - F%x
            slaveL%boundary(slave)%Fy = slaveL%boundary(slave)%Fy - F%y
        end if

    end subroutine calculate_diagonal_link_force

    subroutine create_cilia_array(CA,W,dc,dp,origin)
        class(cilia_array), intent(in out) :: CA
        real(real64),   intent(in) :: W     ! Width of cilia
        real(real64),   intent(in) :: dc    ! Spacing of cilia
        real(real64),   intent(in) :: dp    ! Spacing of particles in a layer
        class(vec),     intent(in) :: origin

        integer(int32) :: ic
        type(vec) :: corigin

        corigin = origin

        ! CA = cilia_array(nc,nl,np)

        do ic = 1,CA%nc
            corigin%x = origin%x + (ic-1)*(dc+W)
            ! if ((ic==1).or.(ic==2).or.(ic==(CA%nc-1)).or.(ic==CA%nc)) then
            if ((ic==1).or.(ic==CA%nc)) then
                call create_cilia(CA%array(ic),W,1.3*dp,corigin)
            else
                call create_cilia(CA%array(ic),W,dp,corigin)
            endif
        end do

    end subroutine create_cilia_array

    subroutine calculate_cilia_array_force(CA,ko,kd,Rl,t)
        class(cilia_array), intent(in out)  :: CA      ! Cilia array
        real(real64),       intent(in)      :: ko      ! Horizontal spring stiffness
        real(real64),       intent(in)      :: kd      ! Diagonal spring stiffness
        real(real64),       intent(in)      :: Rl      ! Resting length
        real(real64),       intent(in)      :: t       ! Time instant

        integer(int32) :: ic

        ! Calculate forces for all the cilia within the array
        do ic = 1,CA%nc
            !if ((ic==1).or.(ic==2).or.(ic==(CA%nc-1)).or.(ic==CA%nc)) then
            if ((ic==1).or.(ic==CA%nc)) then
                ! call calculate_cilia_force(CA%array(ic),4*ko,4*kd,Rl,t)
                call calculate_cilia_force(CA%array(ic),ko,kd,Rl,t)
            else
                call calculate_cilia_force(CA%array(ic),ko,kd,Rl,t)
            endif
        end do
    end subroutine calculate_cilia_array_force

    subroutine update_cilia(C,dt)
        class(cilia), intent(in out)    :: C
        real(real64), intent(in)        :: dt

        integer(int32) :: il

        do il = 1,C%nl
            call update_ib(C%layers(il),dt)
        end do
    end subroutine update_cilia
    
    subroutine update_cilia_array(CA,dt)
        class(cilia_array), intent(in out)    :: CA
        real(real64), intent(in)        :: dt

        integer(int32) :: ic

        do ic = 1,CA%nc 
            call update_cilia(CA%array(ic),dt)
        end do
    end subroutine update_cilia_array

    subroutine spread_force_cilia(M,C,Fx,Fy)
        class(mesh), intent(in) :: M
        class(cilia), intent(in) :: C
        real(real64), intent(in out) :: Fx(:,:), Fy(:,:)

        integer(int32) :: il

        do il = 1,C%nl
            call spread_force_compact(M,C%layers(il),Fx,Fy)
        end do
    end subroutine

    subroutine spread_force_cilia_array(M,CA,Fx,Fy)
        class(mesh), intent(in) :: M
        class(cilia_array), intent(in) :: CA
        real(real64), intent(in out) :: Fx(:,:), Fy(:,:)

        integer(int32) :: ic

        do ic = 1,CA%nc
            call spread_force_cilia(M,CA%array(ic),Fx,Fy)
        end do
    end subroutine spread_force_cilia_array

    subroutine interpolate_velocity_cilia(M,C,u,v)
        class(mesh), intent(in) :: M
        class(cilia), intent(in out) :: C
        real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)

        integer(int32) :: il

        do il = 1,C%nl
            call interpolate_velocity_compact(M,C%layers(il),u,v)
        end do
    end subroutine interpolate_velocity_cilia

    subroutine interpolate_velocity_cilia_array(M,CA,u,v)
        class(mesh), intent(in) :: M
        class(cilia_array), intent(in out) :: CA
        real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)

        integer(int32) :: ic

        do ic = 1,CA%nc
            call interpolate_velocity_cilia(M,CA%array(ic),u,v)
        end do
    end subroutine interpolate_velocity_cilia_array

    subroutine initialize_velocity_cilia_array(CA)
        class(cilia_array), intent(inout) :: CA

        integer(int32) :: ic,il,ip
        
        do ic = 1,CA%nc
            do il = 1,CA%array(ic)%nl
                do ip = 1,size(CA%array(ic)%layers(il)%boundary)
                    CA%array(ic)%layers(il)%boundary(ip)%Ux = 0.0d0
                    CA%array(ic)%layers(il)%boundary(ip)%Uy = 0.0d0
                end do
            end do
        end do

    end subroutine initialize_velocity_cilia_array

    subroutine apply_tip_force_cilia(C,Ftip,t)
        class(cilia), intent(inout) :: C
        real(real64), intent(in) :: Ftip
        real(real64), intent(in) :: t

        integer(int32) :: il, np
        integer(int32) :: pt

        pt = np
        if (C%layers(1)%t.eq.'c') then
           pt = 1 
        end if

        np = size(C%layers(1)%boundary)
        do il = 1,C%nl
            C%layers(il)%boundary(pt)%Fx = Ftip!*cos(2*PI/2.0*t)
        end do

    end subroutine apply_tip_force_cilia

    subroutine apply_tip_force_cilia_array(CA,Ftip,t)
        class(cilia_array), intent(inout) :: CA
        real(real64), intent(in) :: Ftip
        real(real64), intent(in) :: t

        integer(int32) :: ic

        do ic = 1,CA%nc
            call apply_tip_force_cilia(CA%array(ic),Ftip,t)
        end do

    end subroutine apply_tip_force_cilia_array

    subroutine store_original_locations(CA)
        class(cilia_array), intent(inout) :: CA

        integer(int32) :: ic, il, ip

        do ic = 1,CA%nc
            do il = 1,CA%array(ic)%nl
                do ip = 1,CA%array(ic)%layers(il)%np
                    CA%array(ic)%layers(il)%boundary(ip)%xo = CA%array(ic)%layers(il)%boundary(ip)%x
                    CA%array(ic)%layers(il)%boundary(ip)%yo = CA%array(ic)%layers(il)%boundary(ip)%y
                end do
            end do
        end do

    end subroutine store_original_locations

    function spring_force(master,slave,k,t,beta) result(F)
        class(particle), intent(in) :: master, slave
        real(real64), intent(in) :: k
        real(real64), intent(in)      :: t       ! Time instant
        real(real64), intent(in), optional      :: beta    ! Resting length coefficient

        type(vec) :: F
        real(real64) :: xm,ym,xsl,ysl,Rl,d 
        real(real64) :: alpha
        real(real64) :: omga

        ! Calculate the original spacing between the particles
        ! Master node location
        xm = master%xo
        ym = master%yo

        ! Slave node location
        xsl = slave%xo
        ysl = slave%yo

        ! Calculate distance between master and slave nodes
        Rl = norm2([(xsl-xm),(ysl-ym)])



        ! Calculate the resting length factor if the theta is nonzero
        ! For a theta value of zero, alpha is going to be 1
        ! But this implementation nedds to be checked and mofified

        omga = 2*PI/1.0d0
        alpha = 1 + 0.40d0*sin(omga*t) * sin(2*master%theta)! + phirnd*2*PI)

        ! If the alpha value is already provided
        if (present(beta)) then
            alpha = beta
        end if


        Rl = alpha * Rl
        ! Calculate the current spacing between particles
        ! Master node location
        xm = master%x
        ym = master%y

        ! Slave node location
        xsl = slave%x
        ysl = slave%y

        ! Calculate distance between master and slave nodes
        d = norm2([(xsl-xm),(ysl-ym)])
        
        ! Calculate forces (Master node)
        F%x = k*(1.0d0-Rl/d)*(xsl-xm)
        F%y = k*(1.0d0-Rl/d)*(ysl-ym)
    end function spring_force

end module mod_ibm
