module mod_time
    use iso_fortran_env, only: int32, real32, real64
    use mod_mesh
    implicit none
    
contains

    pure subroutine predictor(M,u,v,us,vs,nu,dt)
        class(mesh), intent(in) :: M
        real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in out) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real32), intent(in) :: nu, dt

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
    
    pure subroutine corrector(M,u,v,us,vs,p,rho,dt)
        class(mesh), intent(in) :: M
        real(real64), intent(in out) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: P(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub)
        real(real64), intent(in) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real32), intent(in) :: rho, dt

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

    pure subroutine calculate_rhs(M,us,vs,R,rho,dt)
        class(mesh), intent(in) :: M
        real(real64), intent(in) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in out) :: R(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub)
        real(real32), intent(in) :: rho, dt

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
        real(real64), intent(in) :: Fx(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub)
        real(real32), intent(in) :: nu 

        real(real64), intent(in out) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
        integer(int32) :: i,j
        real(real64) :: vcenter, dxi, dyi, Fxavg

        dxi = 1/M%dx
        dyi = 1/M%dy

        ! (u-velocity cell)
        do concurrent (j = M%yu%lb+1:M%yu%ub-1, i = M%xu%lb+1:M%xu%ub-1)
            Fxavg = 0.5d0*(Fx(i-1,j) + Fx(i+1,j))
            vcenter = 0.25*(v(i-1,j) + v(i-1,j+1) + v(i,j) + v(i,j+1))
            us(i,j) =   nu*(u(i-1,j) - 2*u(i,j) + u(i+1,j))*dxi**2 &
                        + nu*(u(i,j-1) -2*u(i,j) + u(i,j+1))*dyi**2 &
                        - u(i,j)*(u(i+1,j) - u(i-1,j))*0.5*dxi &
                        - vcenter*(u(i,j+1)-u(i,j-1))*0.5*dyi + Fxavg
        end do

    end subroutine cdu

    pure subroutine cdv(M,u,v,vs,nu,Fy)
        ! Calculates convection and diffusion terms for u
        class(mesh), intent(in) :: M
        real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub) 
        real(real64), intent(in) :: Fy(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub)
        real(real32), intent(in) :: nu 

        real(real64), intent(in out) :: vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        integer(int32) :: i,j
        real(real64) :: ucenter, dxi, dyi, Fyavg

        dxi = 1/M%dx
        dyi = 1/M%dy

        ! vs 
        ! (v-velocity cell)
        do concurrent (j = M%yv%lb+1:M%yv%ub-1, i = M%xv%lb+1:M%xv%ub-1)
            Fyavg = 0.5d0*(Fy(i,j-1) + Fy(i,j+1))
            ucenter = 0.25*(u(i,j-1)+u(i,j)+u(i+1,j-1)+u(i+1,j))
            vs(i,j) =   nu*(v(i-1,j) - 2*v(i,j) + v(i+1,j))*dxi**2 &
                        + nu*(v(i,j-1) - 2*v(i,j) + v(i,j+1))*dyi**2 &
                        - ucenter*(v(i+1,j) - v(i-1,j))*0.5*dxi &
                        - v(i,j)*(v(i,j+1)-v(i,j-1))*0.5*dyi + Fyavg
        end do

    end subroutine cdv

end module mod_time