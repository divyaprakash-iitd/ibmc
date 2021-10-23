module mod_mesh
    use iso_fortran_env, only: int32, real64, int32, real64
    use mod_vec
    use mod_dims
    
    implicit none

    private
    public :: mesh

    type :: mesh
        character(len=:), allocatable :: name
        real(real64) :: Lx, Ly, dx, dy
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
            real(real64), intent(in) :: Lx, Ly
            integer(int32), intent(in) :: Nx, Ny

            integer(int32) :: imin, imax, jmin, jmax
            integer(int32) :: i, j

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

            ! Initialize
            self%u_mesh(:,:)%x = 0.0d0
            self%u_mesh(:,:)%y = 0.0d0
            self%v_mesh(:,:)%x = 0.0d0
            self%v_mesh(:,:)%y = 0.0d0
            self%p_mesh(:,:)%x = 0.0d0
            self%p_mesh(:,:)%y = 0.0d0

            ! Assign coordinates to mesh points 
            ! u_mesh
            ! x-cordinates
            do i = self%xu%lb+1,self%xu%ub ! Since the first one is zero
                self%u_mesh(i,:)%x = self%u_mesh(i-1,:)%x + self%dx 
            end do
            ! y-coordinates
            ! The u cells start from below the boundary in the y-direction
            self%u_mesh(:,self%yu%lb)%y = -self%dy/2.0d0
            do j = self%yu%lb+1,self%yu%ub
                self%u_mesh(:,j)%y = self%u_mesh(:,j-1)%y + self%dy
            end do
            
            ! v_mesh
            ! x-coordinates
            self%v_mesh(self%xv%lb,:)%x = -self%dx/2.0d0
            do i = self%xv%lb+1,self%xv%ub
                self%v_mesh(i,:)%x = self%v_mesh(i-1,:)%x + self%dx
            end do
            ! y-coordinates
            ! The v cells start from the left side of the boundary in the x-direction
            do j = self%yv%lb+1,self%yv%ub
                self%v_mesh(:,j)%y = self%v_mesh(:,j-1)%y + self%dy
            end do

            ! p_mesh
            ! x-coordinates
            self%p_mesh(self%xp%lb,:)%x = self%dx/2.0d0
            do i = self%xp%lb+1,self%xp%ub
                self%p_mesh(i,:)%x = self%p_mesh(i-1,:)%x + self%dx
            end do
            ! y-coordinates
            self%p_mesh(:,self%yp%lb)%y = self%dy/2.0d0
            do j = self%yp%lb+1,self%yp%ub
                self%p_mesh(:,j)%y = self%p_mesh(:,j-1)%y + self%dy
            end do
        end function mesh_constructor
    
end module mod_mesh