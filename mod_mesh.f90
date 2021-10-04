module mod_mesh
    use iso_fortran_env, only: int32, real32, int64, real64
    use mod_dims
    
    implicit none

    Private

    type :: mesh
        character(len=*) :: mesh
        type(dims) :: xu, yu, xv, yv, xp, yp
        real(real64), allocatable :: u_mesh(:,:), v_mesh, p_mesh
    end type mesh

contains
    
end module mod_mesh