module mod_time_stepping
    use iso_fortran_env, only: int32, real32, real64
    use mod_time
    implicit none
    
contains

    pure subroutine euler(M,u,v,us,vs,nu,dt,Fx,Fy)
        class(mesh), intent(in) :: M
        real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in out) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: Fx(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
        real(real64), intent(in) :: Fy(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real32), intent(in) :: nu, dt
    
        
        call cdu(M,u,v,us,nu,Fx)
        call cdv(M,u,v,vs,nu,Fy) 

        us = u + dt*us
        vs = v + dt*vs

    end subroutine euler
    
    ! subroutine RK2(M,u,v,us,vs,nu,dt)
    !     class(mesh), intent(in) :: M
    !     real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
    !     real(real64), intent(in out) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
    !     real(real32), intent(in) :: nu, dt

    !     real(real64) :: du1(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), dv1(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
    !     real(real64) :: du2(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), dv2(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
    !     real(real64) :: k1u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), k1v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
    !     real(real64) :: k2u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), k2v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        
    !     ! Step 1
    !     call cdu(M,u,v,k1u,nu)
    !     du1 = k1u*dt
    !     call cdv(M,u,v,k1v,nu)
    !     dv1 = k1v*dt
        
    !     ! Step 2
    !     call cdu(M,(u+du1),(v+dv1),k2u,nu)
    !     call cdv(M,(u+du1),(v+dv1),k2v,nu)

    !     us = u + dt/2*(k1u + k2u)
    !     vs = v + dt/2*(k1v + k2v)

    ! end subroutine RK2

    ! subroutine RK4(M,u,v,us,vs,nu,dt)
    !     class(mesh), intent(in) :: M
    !     real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
    !     real(real64), intent(in out) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
    !     real(real32), intent(in) :: nu, dt

    !     real(real64) :: du1(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), dv1(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
    !     real(real64) :: du2(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), dv2(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
    !     real(real64) :: du3(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), dv3(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
    !     real(real64) :: k1u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), k1v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
    !     real(real64) :: k2u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), k2v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
    !     real(real64) :: k3u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), k3v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
    !     real(real64) :: k4u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), k4v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)

    !     ! Step 1
    !     call cdu(M,u,v,k1u,nu)
    !     du1 = k1u*dt/2
    !     call cdv(M,u,v,k1v,nu)
    !     dv1 = k1v*dt/2
        
    !     ! Step 2
    !     call cdu(M,(u+du1),(v+dv1),k2u,nu)
    !     du2 = k2u*dt/2
    !     call cdv(M,(u+du1),(v+dv1),k2v,nu)
    !     dv2 = k2v*dt/2

    !     ! Step 3
    !     call cdu(M,(u+du2),(v+dv2),k3u,nu)
    !     du3 = dt*k3u
    !     call cdv(M,(u+du2),(v+dv2),k3v,nu)
    !     dv3 = dt*k3v

    !     ! Step 4
    !     call cdu(M,(u+du3),(v+dv3),k4u,nu)
    !     call cdv(M,(u+du3),(v+dv3),k4v,nu)

    !     us = u + 1.0d0/6*dt*(k1u + 2*k2u + 2*k3u + k4u)
    !     vs = v + 1.0d0/6*dt*(k1v + 2*k2v + 2*k3v + k4v)

    !     ! us = us/2
    !     ! vs = vs/2

    !     ! call cdu(M,u,v,us,nu)
    !     ! call cdv(M,u,v,vs,nu) 

    !     ! us = u + dt*us
    !     ! vs = v + dt*vs
    ! end subroutine RK4

end module mod_time_stepping
