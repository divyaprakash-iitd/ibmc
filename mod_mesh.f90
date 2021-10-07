module mod_mesh
    use iso_fortran_env, only: int32, real32, int64, real64
    use mod_dims
    
    implicit none

    private
    public :: mesh

    type :: mesh
        character(len=:), allocatable :: name
        real(real32) :: Lx, Ly, dx, dy
        integer(int32) :: Nx, Ny
        type(dims) :: xu, yu, xv, yv, xp, yp
        real(real64), allocatable :: u_mesh(:,:), v_mesh(:,:), p_mesh(:,:)
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
           
            ! Initialize the arrays
            self%u_mesh = 0.0d0
            self%v_mesh = 0.0d0
            self%p_mesh = 0.0d0

        end function mesh_constructor
    
    ! subroutine generate_mesh(self)
    !     class(mesh), intent(in out) :: self

    !     integer(int32) :: i, j

    !     ! Assign coordinates to mesh points 
    !     ! u_mesh
    !     do concurrent (j = 1:self%Nx, i = 1:self%Ny)
    !         u_mesh(i,j) = u_mesh(i,j)
    !     end do
        
        
    ! end subroutine generate_mesh


end module mod_mesh