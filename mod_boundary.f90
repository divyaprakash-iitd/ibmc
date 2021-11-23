module mod_boundary
    use iso_fortran_env, only: int32, real32, int64, real64
    use mod_mesh
    implicit none
   
    private
    public :: apply_boundary

contains

    pure subroutine apply_boundary(M,u,v,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
        class(mesh), intent(in) :: M
        real(real32), intent(in) :: utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright
        real(real64), intent(in out) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)

        ! Apply boundary conditions
        u(:,M%yu%lb) = 2*ubottom - u(:,M%yu%lb+1);
        u(:,M%yu%ub) = 2*utop    - u(:,M%yu%ub-1);
        v(M%xv%lb,:) = 2*vleft   - v(M%xv%lb+1,:);
        v(M%xv%ub,:) = 2*vright  - v(M%xv%ub-1,:);
    

        u(M%xu%lb,:)    = uleft;
        u(M%xu%ub,:)    = uright;
        v(:,M%yv%ub)    = vtop;
        v(:,M%yv%lb)    = vbottom;
    end subroutine apply_boundary
    
end module mod_boundary
