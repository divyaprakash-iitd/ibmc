module mod_boundary
    use iso_fortran_env, only: int32, real64, int32, real64
    use mod_mesh
    implicit none
   
    private
    public :: apply_boundary

contains

    pure subroutine apply_boundary(M,u,v,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
        class(mesh), intent(in) :: M
        real(real64), intent(in) :: utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright
        real(real64), intent(in out) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)

        ! Apply boundary conditions
        u(:,M%yu%lb) = 2*ubottom - u(:,M%yu%lb+1);
        u(:,M%yu%ub) = 2*utop    - u(:,M%yu%ub-1);
        v(M%xv%lb,:) = 2*vleft   - v(M%xv%lb+1,:);
        v(M%yv%ub,:) = 2*vright  - v(M%yv%ub-1,:);
    

        u(M%xu%lb,:)    = uleft;
        u(M%xu%ub,:)    = uright;
        v(:,M%yv%ub)    = vtop;
        v(:,M%yv%lb)    = vbottom;
    end subroutine apply_boundary
    
    pure subroutine apply_boundary_periodic(M,u,v,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
        class(mesh), intent(in)         :: M
        real(real64), intent(in)        :: utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright
        real(real64), intent(in out)    :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)

        ! Create temporary variables for the left and right boundary velocity
        real(real64) :: tempL(size(u,1))
        real(real64) :: tempR(size(u,1))

        ! Store velocity values in temporary variables
        tempL = u(M%xu%lb,:)
        tempR = u(M%xu%ub,:)

        ! Apply periodic boundary conditions to u and v cells
        u(M%xu%lb,:) = tempR ! Left boundary
        u(M%xu%ub,:) = tempL ! Right boundary

        v(M%xv%lb,:) = v(M%xv%ub-1,:) ! Left boundary
        v(M%xv%ub,:) = v(M%xv%lb+1,:) ! Right boundary

        u(:,M%yu%lb) = 2*ubottom - u(:,M%yu%lb+1); ! Top boundary
        u(:,M%yu%ub) = 2*utop    - u(:,M%yu%ub-1); ! Bottom boundary
    
        v(:,M%yv%ub)    = vtop;     ! Top boundary
        v(:,M%yv%lb)    = vbottom;  ! Bottom boundary
    end subroutine apply_boundary_periodic

end module mod_boundary
