module mod_time
    use iso_fortran_env, only: int32, real64, real64
    use mod_mesh
    implicit none
    
contains

    pure subroutine predictor(M,u,v,us,vs,nu,dt)
        class(mesh), intent(in) :: M
        real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in out) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: nu, dt

        integer(int32) :: i,j
        real(real64) :: ucenter, vcenter, dxi, dyi

        dxi = 1/M%dx
        dyi = 1/M%dy

        ! Perform predictor step
       
        ! us 
        ! (u-velocity cell)
        do concurrent (j = M%yu%lb+1:M%yu%ub-1, i = M%xu%lb+1:M%xu%ub-1)
            vcenter = 0.25*(v(i-1,j) + v(i-1,j+1) + v(i,j) + v(i,j+1))
            us(i,j) = u(i,j) + dt* &
                        ( nu*(u(i-1,j) - 2*u(i,j) + u(i+1,j))*dxi**2 &
                        + nu*(u(i,j-1) -2*u(i,j) + u(i,j+1))*dyi**2 &
                        - u(i,j)*(u(i+1,j) - u(i-1,j))*0.5*dxi &
                        - vcenter*(u(i,j+1)-u(i,j-1))*0.5*dyi)
        end do

        ! vs 
        ! (v-velocity cell)
        do concurrent (j = M%yv%lb+1:M%yv%ub-1, i = M%xv%lb+1:M%xv%ub-1)
            ucenter = 0.25*(u(i,j-1)+u(i,j)+u(i+1,j-1)+u(i+1,j))
            vs(i,j) = v(i,j) + dt* &
                        ( nu*(v(i-1,j) - 2*v(i,j) + v(i+1,j))*dxi**2 &
                        + nu*(v(i,j-1) - 2*v(i,j) + v(i,j+1))*dyi**2 &
                        - ucenter*(v(i+1,j) - v(i-1,j))*0.5*dxi &
                        - v(i,j)*(v(i,j+1)-v(i,j-1))*0.5*dyi)
        end do

    end subroutine predictor
    
    pure subroutine predictor_periodic(M,u,v,us,vs,nu,dt)
        class(mesh), intent(in) :: M
        real(real64), intent(in out) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in out) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: nu, dt

        integer(int32) :: i,j
        real(real64) :: ucenter, vcenter, dxi, dyi

        ! Padded u at the left and right boundary for applying periodic boundary condtions
        real(real64) :: uPad(M%xu%lb-1:M%xu%ub+1,M%yu%lb:M%yu%ub)

        ! Apply boundary conditions to u via padding
        uPad(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub) = u
        uPad(M%xu%lb-1,:) = u(M%xu%ub,:)
        uPad(M%xu%ub+1,:) = u(M%xu%lb,:)

        ! Apply periodic boundary conditions to v
        v(M%xv%lb,:) = v(M%xv%ub-1,:) ! Left boundary
        v(M%xv%ub,:) = v(M%xv%lb+1,:) ! Right boundary

        dxi = 1/M%dx
        dyi = 1/M%dy

        ! Perform predictor step
       
        ! us 
        ! (u-velocity cell)
        do concurrent (j = M%yu%lb+1:M%yu%ub-1, i = M%xu%lb:M%xu%ub)
            vcenter = 0.25*(v(i-1,j) + v(i-1,j+1) + v(i,j) + v(i,j+1))
            us(i,j) = uPad(i,j) + dt* &
                        ( nu*(uPad(i-1,j) - 2*uPad(i,j) + uPad(i+1,j))*dxi**2 &
                        + nu*(uPad(i,j-1) -2*uPad(i,j) + uPad(i,j+1))*dyi**2 &
                        - uPad(i,j)*(uPad(i+1,j) - uPad(i-1,j))*0.5*dxi &
                        - vcenter*(uPad(i,j+1)-uPad(i,j-1))*0.5*dyi)
        end do

        ! vs 
        ! (v-velocity cell)
        ! No changes in v-velocity calculation for periodic boundary condition
        do concurrent (j = M%yv%lb+1:M%yv%ub-1, i = M%xv%lb+1:M%xv%ub-1)
            ucenter = 0.25*(u(i,j-1)+u(i,j)+u(i+1,j-1)+u(i+1,j))
            vs(i,j) = v(i,j) + dt* &
                        ( nu*(v(i-1,j) - 2*v(i,j) + v(i+1,j))*dxi**2 &
                        + nu*(v(i,j-1) - 2*v(i,j) + v(i,j+1))*dyi**2 &
                        - ucenter*(v(i+1,j) - v(i-1,j))*0.5*dxi &
                        - v(i,j)*(v(i,j+1)-v(i,j-1))*0.5*dyi)
        end do

    end subroutine predictor_periodic

    pure subroutine corrector(M,u,v,us,vs,p,rho,dt)
        class(mesh), intent(in) :: M
        real(real64), intent(in out) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: P(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub)
        real(real64), intent(in) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: rho, dt

        integer(int32) :: i,j
        real(real64) :: dxi, dyi

        dxi = 1/M%dx
        dyi = 1/M%dy

        ! Perform the corrector steps
        ! u 
        ! (u-velocity cell)
        do concurrent (j = M%yu%lb+1:M%yu%ub-1, i = M%xu%lb+1:M%xu%ub-1)
            u(i,j) = us(i,j) - dt/rho * (p(i,j) - p(i-1,j)) * dxi
        end do
        ! 
        ! v
        ! (v-velocity cell)
        do concurrent (j = M%yv%lb+1:M%yv%ub-1, i = M%xv%lb+1:M%xv%ub-1)
            v(i,j) = vs(i,j) - dt/rho * (p(i,j) - p(i,j-1)) * dyi
        end do
    end subroutine corrector

    pure subroutine corrector_periodic(M,u,v,us,vs,p,rho,dt)
        class(mesh), intent(in) :: M
        real(real64), intent(in out) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: P(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub)
        real(real64), intent(in out) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: rho, dt

        integer(int32) :: i,j,imin,imax,jmin,jmax
        real(real64), allocatable :: Ppad(:,:)
        real(real64) :: dxi, dyi

        ! Padded u at the left and right boundary for applying periodic boundary condtions
        real(real64) :: uPad(M%xu%lb-1:M%xu%ub+1,M%yu%lb:M%yu%ub)

        ! Apply boundary conditions to us via padding
        uPad(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub) = us
        uPad(M%xu%lb-1,:) = us(M%xu%ub,:)
        uPad(M%xu%ub+1,:) = us(M%xu%lb,:)

        ! Apply periodic boundary conditions to vs
        vs(M%xv%lb,:) = vs(M%xv%ub-1,:) ! Left boundary
        vs(M%xv%ub,:) = vs(M%xv%lb+1,:) ! Right boundary
        
        ! Define limits
        imin = lbound(P,1)
        imax = ubound(P,1)
        jmin = lbound(P,2)
        jmax = ubound(P,2)

        allocate(Ppad(imin-1:imax+1,jmin-1:jmax+1))
        Ppad = 0.0d0
        Ppad(imin:imax,jmin:jmax) = P
        Ppad(imin-1,:) = Ppad(imax,:)
        Ppad(imax+1,:) = Ppad(imin,:)

        dxi = 1/M%dx
        dyi = 1/M%dy

        ! Perform the corrector steps
        ! u 
        ! (u-velocity cell)
        do concurrent (j = M%yu%lb+1:M%yu%ub-1, i = M%xu%lb:M%xu%ub)
            u(i,j) = uPad(i,j) - dt/rho * (Ppad(i,j) - Ppad(i-1,j)) * dxi
        end do
        ! 
        ! v
        ! (v-velocity cell)
        do concurrent (j = M%yv%lb:M%yv%ub-1, i = M%xv%lb+1:M%xv%ub-1)
            v(i,j) = vs(i,j) - dt/rho * (Ppad(i,j) - Ppad(i,j-1)) * dyi
        end do

    end subroutine corrector_periodic

    pure subroutine calculate_rhs(M,us,vs,R,rho,dt)
        class(mesh), intent(in) :: M
        real(real64), intent(in) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in out) :: R(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub)
        real(real64), intent(in) :: rho, dt

        integer(int32) :: i,j
        real(real64) :: dxi, dyi

        dxi = 1/M%dx
        dyi = 1/M%dy

        ! Form the right hand side of the pressure poisson equation
        ! (Pressure cell)
        do concurrent (j = M%yp%lb:M%yp%ub, i = M%xp%lb:M%xp%ub)
            R(i,j) = -rho/dt* &
                    ( (us(i+1,j) - us(i,j))*dxi &
                    + (vs(i,j+1) - vs(i,j))*dyi)
        end do
    end subroutine calculate_rhs

    pure subroutine cdu(M,u,v,us,nu,Fx)
        ! Calculates convection and diffusion terms for u
        class(mesh), intent(in) :: M
        real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: Fx(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
        real(real64), intent(in) :: nu 

        real(real64), intent(in out) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
        integer(int32) :: i,j
        real(real64) :: vcenter, dxi, dyi 

        dxi = 1/M%dx
        dyi = 1/M%dy

        ! (u-velocity cell)
        do concurrent (j = M%yu%lb+1:M%yu%ub-1, i = M%xu%lb+1:M%xu%ub-1)
            vcenter = 0.25*(v(i-1,j) + v(i-1,j+1) + v(i,j) + v(i,j+1))
            us(i,j) =   nu*(u(i-1,j) - 2*u(i,j) + u(i+1,j))*dxi**2 &
                        + nu*(u(i,j-1) -2*u(i,j) + u(i,j+1))*dyi**2 &
                        - u(i,j)*(u(i+1,j) - u(i-1,j))*0.5*dxi &
                        - vcenter*(u(i,j+1)-u(i,j-1))*0.5*dyi + Fx(i,j)
        end do

    end subroutine cdu

    pure subroutine cdu_p(M,u,v,us,nu,Fx)
        ! Calculates convection and diffusion terms for u
        class(mesh), intent(in) :: M
        real(real64), intent(in out) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: Fx(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
        real(real64), intent(in) :: nu 

        real(real64), intent(in out) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
        integer(int32) :: i,j
        real(real64) :: vcenter, dxi, dyi 

        ! Padded u at the left and right boundary for applying periodic boundary condtions
        real(real64) :: uPad(M%xu%lb-1:M%xu%ub+1,M%yu%lb:M%yu%ub)

        ! Apply boundary conditions to u via padding
        uPad(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub) = u
        uPad(M%xu%lb-1,:) = u(M%xu%ub,:)
        uPad(M%xu%ub+1,:) = u(M%xu%lb,:)

        ! Apply periodic boundary conditions to v
        v(M%xv%lb,:) = v(M%xv%ub-1,:) ! Left boundary
        v(M%xv%ub,:) = v(M%xv%lb+1,:) ! Right boundary

        dxi = 1/M%dx
        dyi = 1/M%dy

        ! Perform predictor step
       
        ! us 
        ! (u-velocity cell)
        do concurrent (j = M%yu%lb+1:M%yu%ub-1, i = M%xu%lb:M%xu%ub)
            vcenter = 0.25*(v(i-1,j) + v(i-1,j+1) + v(i,j) + v(i,j+1))
            us(i,j) =   nu*(uPad(i-1,j) - 2*uPad(i,j) + uPad(i+1,j))*dxi**2 &
                        + nu*(uPad(i,j-1) -2*uPad(i,j) + uPad(i,j+1))*dyi**2 &
                        - uPad(i,j)*(uPad(i+1,j) - uPad(i-1,j))*0.5*dxi &
                        - vcenter*(uPad(i,j+1)-uPad(i,j-1))*0.5*dyi + Fx(i,j)
        end do

    end subroutine cdu_p

    pure subroutine cdv_p(M,u,v,vs,nu,Fy)
        ! Calculates convection and diffusion terms for u
        class(mesh), intent(in) :: M
        real(real64), intent(in out) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub) 
        real(real64), intent(in) :: Fy(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: nu 

        real(real64), intent(in out) :: vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        integer(int32) :: i,j
        real(real64) :: ucenter, dxi, dyi

        ! Apply periodic boundary conditions to v
        v(M%xv%lb,:) = v(M%xv%ub-1,:) ! Left boundary
        v(M%xv%ub,:) = v(M%xv%lb+1,:) ! Right boundary
        dxi = 1/M%dx
        dyi = 1/M%dy

        ! vs 
        ! (v-velocity cell)
        do concurrent (j = M%yv%lb+1:M%yv%ub-1, i = M%xv%lb+1:M%xv%ub-1)
            ucenter = 0.25*(u(i,j-1)+u(i,j)+u(i+1,j-1)+u(i+1,j))
            vs(i,j) =   nu*(v(i-1,j) - 2*v(i,j) + v(i+1,j))*dxi**2 &
                        + nu*(v(i,j-1) - 2*v(i,j) + v(i,j+1))*dyi**2 &
                        - ucenter*(v(i+1,j) - v(i-1,j))*0.5*dxi &
                        - v(i,j)*(v(i,j+1)-v(i,j-1))*0.5*dyi + Fy(i,j)
        end do

    end subroutine cdv_p

    pure subroutine cdv(M,u,v,vs,nu,Fy)
        ! Calculates convection and diffusion terms for u
        class(mesh), intent(in) :: M
        real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub) 
        real(real64), intent(in) :: Fy(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: nu 

        real(real64), intent(in out) :: vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        integer(int32) :: i,j
        real(real64) :: ucenter, dxi, dyi

        dxi = 1/M%dx
        dyi = 1/M%dy

        ! vs 
        ! (v-velocity cell)
        do concurrent (j = M%yv%lb+1:M%yv%ub-1, i = M%xv%lb+1:M%xv%ub-1)
            ucenter = 0.25*(u(i,j-1)+u(i,j)+u(i+1,j-1)+u(i+1,j))
            vs(i,j) =   nu*(v(i-1,j) - 2*v(i,j) + v(i+1,j))*dxi**2 &
                        + nu*(v(i,j-1) - 2*v(i,j) + v(i,j+1))*dyi**2 &
                        - ucenter*(v(i+1,j) - v(i-1,j))*0.5*dxi &
                        - v(i,j)*(v(i,j+1)-v(i,j-1))*0.5*dyi + Fy(i,j)
        end do

    end subroutine cdv

    pure function cdu_f(M,u,v,nu,Fx)
        ! Calculates convection and diffusion terms for u

        class(mesh), intent(in) :: M
        real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: Fx(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
        real(real64), intent(in) :: nu 

        real(real64) :: cdu_f(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
        integer(int32) :: i,j
        real(real64) :: vcenter, dxi, dyi 

        dxi = 1/M%dx
        dyi = 1/M%dy

        ! (u-velocity cell)
        do concurrent (j = M%yu%lb+1:M%yu%ub-1, i = M%xu%lb+1:M%xu%ub-1)
            vcenter = 0.25*(v(i-1,j) + v(i-1,j+1) + v(i,j) + v(i,j+1))
            cdu_f(i,j) =   nu*(u(i-1,j) - 2*u(i,j) + u(i+1,j))*dxi**2 &
                        + nu*(u(i,j-1) -2*u(i,j) + u(i,j+1))*dyi**2 &
                        - u(i,j)*(u(i+1,j) - u(i-1,j))*0.5*dxi &
                        - vcenter*(u(i,j+1)-u(i,j-1))*0.5*dyi + Fx(i,j)
        end do

    end function cdu_f

    pure function cdv_f(M,u,v,nu,Fy)
        ! Calculates convection and diffusion terms for u
        class(mesh), intent(in) :: M
        real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub) 
        real(real64), intent(in) :: Fy(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: nu 

        real(real64)  :: cdv_f(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        integer(int32) :: i,j
        real(real64) :: ucenter, dxi, dyi

        dxi = 1/M%dx
        dyi = 1/M%dy

        ! vs 
        ! (v-velocity cell)
        do concurrent (j = M%yv%lb+1:M%yv%ub-1, i = M%xv%lb+1:M%xv%ub-1)
            ucenter = 0.25*(u(i,j-1)+u(i,j)+u(i+1,j-1)+u(i+1,j))
            cdv_f(i,j) =   nu*(v(i-1,j) - 2*v(i,j) + v(i+1,j))*dxi**2 &
                        + nu*(v(i,j-1) - 2*v(i,j) + v(i,j+1))*dyi**2 &
                        - ucenter*(v(i+1,j) - v(i-1,j))*0.5*dxi &
                        - v(i,j)*(v(i,j+1)-v(i,j-1))*0.5*dyi + Fy(i,j)
        end do

    end function cdv_f
end module mod_time
