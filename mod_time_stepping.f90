module mod_time_stepping
    use iso_fortran_env, only: int32, real32, real64
    use mod_time
    implicit none
    
contains

    subroutine euler(M,u,v,us,vs,nu,dt)
        class(mesh), intent(in) :: M
        real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in out) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real32), intent(in) :: nu, dt
    
        ! print *, lbound(fu(M,u,v,nu),1)
        
        call fu(M,u,v,us,nu)
        call fv(M,u,v,vs,nu) 

        us = u + dt*us
        vs = v + dt*vs

    end subroutine euler
    
    subroutine RK2(M,u,v,us,vs,nu,dt)
        class(mesh), intent(in) :: M
        real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in out) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real32), intent(in) :: nu, dt

        real(real64) :: du1(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), dv1(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64) :: du2(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), dv2(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64) :: tmpu(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), tmpv(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        
        call fu(M,u,v,tmpu,nu)
        du1 = real(dt,real64)*tmpu

        call fv(M,u,v,tmpv,nu)
        dv1 = real(dt,real64)*tmpv
        
        call fu(M,(u+du1),(v+dv1),tmpu,nu)
        du2 = real(dt,real64)*tmpu
        
        call fv(M,(u+du1),(v+dv1),tmpv,nu)
        dv2 = real(dt,real64)*tmpv


        us = u + (du1+du2)/2
        vs = v + (dv1+dv2)/2
    end subroutine RK2

    subroutine RK4(M,u,v,us,vs,nu,dt)
        class(mesh), intent(in) :: M
        real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in out) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real32), intent(in) :: nu, dt

        real(real64) :: h
        real(real64) :: du1(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), dv1(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64) :: du2(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), dv2(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64) :: du3(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), dv3(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64) :: du4(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), dv4(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64) :: tmpu(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), tmpv(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
       
        h = real(dt,real64)

        ! Step 1
        call fu(M,u,v,tmpu,nu)
        du1 = tmpu*h/2
        call fv(M,u,v,tmpv,nu)
        dv1 = tmpv*h/2
        
        ! Step 2
        call fu(M,(u+du1),(v+dv1),tmpu,nu)
        du2 = h*tmpu/2
        call fv(M,(u+du1),(v+dv1),tmpv,nu)
        dv2 = h*tmpv/2

        ! Step 3
        call fu(M,(u+du2),(v+dv2),tmpu,nu)
        du3 = h*tmpu/2
        call fv(M,(u+du2),(v+dv2),tmpv,nu)
        dv3 = h*tmpv/2

        ! Step 4
        call fu(M,(u+du3),(v+dv3),tmpu,nu)
        du4 = tmpu
        call fv(M,(u+du3),(v+dv3),tmpv,nu)
        dv4 = tmpv

        us = u + 1.0d0/6*h*(du1*2/h+2*du2+2*du3+du4)
        vs = v + 1.0d0/6*h*(dv1*2/h+2*dv2+2*dv3+dv4)

    end subroutine RK4

end module mod_time_stepping
