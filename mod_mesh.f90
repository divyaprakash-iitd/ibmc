module mod_mesh
    use iso_fortran_env, only: int32, real32, int64, real64
    use mod_dims
    
    implicit none

    Private

    type :: mesh
        character(len=*) :: mesh
        real(real32) :: Lx, Ly, dx, dy
        integer(int32) :: Nx, Ny
        type(dims) :: xu, yu, xv, yv, xp, yp
        ! To-do
        ! real(real64), allocatable :: u_mesh(:,:), v_mesh, p_mesh
    end type mesh

    interface mesh
        module procedure :: mesh_constructor
    end interface mesh

contains

    pure type(mesh) function mesh_constructor(name,Lx,Ly,Nx,Ny)
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

        end function mesh_constructor
    
end module mod_mesh