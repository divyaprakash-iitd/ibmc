module mod_mesh
    use iso_fortran_env, only: int32, real32, int64, real64
    use mod_vec
    use mod_dims
    
    implicit none

    private
    public :: mesh

    type :: mesh
        character(len=:), allocatable :: name
        real(real32) :: Lx, Ly, dx, dy
        integer(int32) :: Nx, Ny
        type(dims) :: xu, yu, xv, yv, xp, yp
        type(vec), allocatable :: u_mesh(:,:), v_mesh(:,:), p_mesh(:,:)
    end type mesh

    interface mesh
        module procedure :: mesh_constructor
    end interface mesh

contains

    pure type(mesh) function mesh_constructor(name,Lx,Ly,Nx,Ny) result(self)
            character(len=*), intent(in) :: name
            real(real32), intent(in) :: Lx, Ly
            integer(int32), intent(in) :: Nx, Ny

            integer(int32) :: imin, imax, jmin, jmax

            imin = 1
            imax = imin + Nx - 1
            jmin = 1
            jmax = jmin + Ny - 1 

            self%name = name
            self%Lx = Lx
            self%Ly = Ly

            self%dx = Lx/Nx
            self%dy = Ly/Ny

            self%xu = dims('xu',imin,imax+1)
            self%yu = dims('yu',jmin-1,jmax+1)

            self%xv = dims('xv',imin-1,imax+1)
            self%yv = dims('yv',jmin, jmax+1)

            self%xp = dims('xp',imin,imax)
            self%yp = dims('yp',jmin,jmax)

            ! Allocate mesh storage arrays
            allocate(self%u_mesh(self%xu%lb:self%xu%ub,self%yu%lb:self%yu%ub))
            allocate(self%v_mesh(self%xv%lb:self%xv%ub,self%yv%lb:self%yv%ub))
            allocate(self%p_mesh(self%xp%lb:self%xp%ub,self%yp%lb:self%yp%ub))

        end function mesh_constructor
    
    subroutine generate_mesh(self)
        class(mesh), intent(in out) :: self

        integer(int32) :: i, j

        ! Assign coordinates to mesh points 
        ! u_mesh
        ! x-ccordinates
        do i = 2,size(self%u_mesh,1)
            self%u_mesh(i,:)%x = self%u_mesh(i,:)%x + self%dx 
        end do
        ! y-coordinates
        ! The u cells start from below the boundary in the y-direction
        self%u_mesh(:,1)%y = -self%dy/2.0d0
        do j = 2,size(self%u_mesh,2)
            self%u_mesh(:,j)%y = self%u_mesh(:,j)%y + self%dy
        end do

        ! v_mesh
        ! x-coordinates
        self%v_mesh(1,:)%x = -self%dx/2.0d0
        do i = 2,size(self%v_mesh,1)
            self%v_mesh(i,:)%x = self%v_mesh(i,:)%x + self%dx
        end do
        ! y-coordinates
        ! The v cells start from the left side of the boundary in the x-direction
        do j = 2,size(self%v_mesh),2
            self%v_mesh(:,j)%y = self%v_mesh(:,j)%y + self%dy
        end do
    end subroutine generate_mesh


end module mod_mesh