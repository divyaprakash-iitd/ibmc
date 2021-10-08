module mod_ibm
    use iso_fortran_env, only: int32, real32, real64
    use mod_mesh
    use mod_ib
    implicit none
   
    real(real64), parameter :: PI = 3.141592653589793

    private
    public :: create_ib, spread_force
contains
   
    subroutine create_ib(B)
        class(ib), intent(in out) :: B

        ! Creates and Initializes an immersed boundary structure
        ! Location
        B%boundary(1)%x = 0.5d0
        B%boundary(1)%y = 0.5d0
        ! Forces
        B%boundary(1)%Fx = -0.1d0
        B%boundary(1)%Fy = 0.0d0
    end subroutine

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
                    phi = 0
                end if
            
                dirac = dirac*phi
            end do        

            dirac = 1/h**2 * dirac

        end function dirac
    end subroutine spread_force


end module mod_ibm