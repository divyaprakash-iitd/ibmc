module mod_boundary
    use iso_fortran_env, only: int32, real64, int32, real64
    use mod_mesh
    implicit none
   
    private
    public :: apply_boundary, apply_boundary_channel, &
              apply_parabolic_inlet, apply_parabolic_initialization, &
              apply_pulsating_inlet

contains

    pure subroutine apply_boundary(M,u,v,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
        class(mesh), intent(in) :: M
        real(real64), intent(in) :: utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright
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
    
    pure subroutine apply_boundary_channel(M,u,v,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
        class(mesh), intent(in) :: M
        real(real64), intent(in) :: utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright
        real(real64), intent(in out) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)

        ! Apply boundary conditions
        ! Top
        u(:,M%yu%lb) = 2*ubottom - u(:,M%yu%lb+1);
        u(:,M%yu%ub) = 2*utop    - u(:,M%yu%ub-1);
        ! Bottom
        v(:,M%yv%ub)    = vtop;
        v(:,M%yv%lb)    = vbottom;
        ! Outlet
        u(M%xu%ub,:)    = u(M%xu%ub-1,:);
        v(M%xv%ub,:)    = v(M%xv%ub-1,:);
        ! Inlet
        u(M%xu%lb,:)    = u(M%xu%lb+1,:);
        v(M%xv%lb,:)    = v(M%xv%lb+1,:); 
        ! u(M%xu%lb,:)    = uleft;
        ! v(M%xv%lb,:) = 2*vleft   - v(M%xv%lb+1,:);
    end subroutine apply_boundary_channel

    subroutine apply_parabolic_inlet(M,u,uleft)
        class(mesh), intent(in)         :: M
        real(real64), intent(in out)    :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
        real(real64), intent(in)        :: uleft

        integer(int32) :: i
        real(real64) :: y

        do concurrent (i = M%yu%lb:M%yu%ub)
            y = M%u_mesh(M%xu%lb,i)%y 
            u(M%xu%lb,i) = 6*uleft*y/M%Ly*(1-y/M%Ly) 
        end do 

    end subroutine apply_parabolic_inlet

    subroutine apply_pulsating_inlet(M,u,uleft,t)
        class(mesh), intent(in)         :: M
        real(real64), intent(in out)    :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
        real(real64), intent(in)        :: uleft
        real(real64), intent(in)        :: t

        integer(int32) :: i
        real(real64) :: y

        do concurrent (i = M%yu%lb:M%yu%ub)
            y = M%u_mesh(M%xu%lb,i)%y 
            ! u(M%xu%lb,i) = 6*uleft*y/M%Ly*(1-y/M%Ly) 
            u(M%xu%lb,i) = sin(2*3.14159*t)*uleft 
        end do 

    end subroutine apply_pulsating_inlet
    
    subroutine apply_parabolic_initialization(M,u,uleft)
        class(mesh), intent(in)         :: M
        real(real64), intent(in out)    :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
        real(real64), intent(in)        :: uleft

        integer(int32) :: i,j
        real(real64) :: y

        do concurrent (j = M%yu%lb:M%yu%ub)
            do concurrent (i = M%xu%lb:M%xu%ub)
                y = M%u_mesh(M%xu%lb,j)%y 
                u(i,j) = 6*uleft*y/M%Ly*(1-y/M%Ly) 
        end do
    end do 
    end subroutine apply_parabolic_initialization

end module mod_boundary
