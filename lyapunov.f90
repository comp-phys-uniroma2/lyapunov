program lyapunov 
 use precision
 use functions
 use solvers
 use gram_schmidt
 implicit none

 real(dp) :: h, t, tfin, abserr, error, s 
 real(dp), dimension(3) :: u, unext, lambda, cum
 real(dp), dimension(3,3) :: VO, VA
 integer :: i, k, Nstep, Niters 
 character(10) :: arg, fname
 
 if(iargc()<6) then
  print*,'lyapunov h Nstep Niters ss rr bb'  
  stop  
 endif

 call getarg(1,arg)
 read(arg,*) h
 
 call getarg(2,arg)
 read(arg,*) Nstep
 
 call getarg(3,arg)
 read(arg,*) Niters 

 call getarg(4,arg)
 read(arg,*) ss
 
 call getarg(5,arg)
 read(arg,*) rr
   
 call getarg(6,arg)
 read(arg,*) bb
 
 ! ------------------------------------------------------------------
 open(194, file='values.dat')

 t=0.0_dp
 u(1) = 1.0_dp
 u(2) = 0.0_dp
 u(3) = 0.0_dp
 s = 1.0_dp

 tfin = 0.0_dp
 VO = 0.0_dp
 VO(1,1) = 1.0_dp 
 VO(2,2) = 1.0_dp 
 VO(3,3) = 1.0_dp
 cum = 0.0_dp

 ! Qualsiasi base iniziale funziona, e.g., anche:
 ![ 1 0 0 ]
 ![ 0 1 0 ]
 ![ 0 0 1 ]
 
 ! remove transient
 do k = 1, Niters*Nstep 
   call dopri54(lorenz, t, h, u, unext, abserr)
   t = t + h
   u = unext
 end do

    do k = 1, Niters ! numero di iterazioni tra una call di gs e l'altra 
        !print*,"iter:",k
        VA=VO
        ! Propago per un po' di step
        do i = 1, Nstep
            call dopri54(lorenz, t, h, u, unext, abserr)
            u = unext
            !print*,t,y,abserr
            u1=u(1)
            u2=u(2)
            u3=u(3)
            call dopri54(lorenz_linear, t, h, VA(:,1), VA(:,1), abserr)
            call dopri54(lorenz_linear, t, h, VA(:,2), VA(:,2), abserr)
            call dopri54(lorenz_linear, t, h, VA(:,3), VA(:,3), abserr)
            t = t + h
        enddo
    call gs(VA,VO,lambda)
    tfin = tfin + Nstep*h
    if (any(lambda<0)) then
        stop 'lambda<0: reduce Nstep'
    elseif (any(isnan(lambda))) then
        stop 'lambda=NaN: reduce Nstep'
    else
        cum = cum + log(lambda)
    end if
    write(194,*) tfin, cum/tfin
    enddo  

 close(194)

 ! Faccio la media sugli ultimo 1/3 del campione (convergenza)
 open(194,file='values.dat')
 do i = 1, int(2*Niters/3)
   read(194,*) tfin, cum
 end do  
 lambda=0.0_dp
 do i = int(2*Niters/3)+1, Niters
   read(194,*) tfin, cum
   lambda=lambda+cum
 end do
 lambda(1) = lambda(1)/(real(Niters,dp)/3.0_dp)
 lambda(2) = lambda(2)/(real(Niters,dp)/3.0_dp)
 lambda(3) = lambda(3)/(real(Niters,dp)/3.0_dp)
 print*,'lyapunov exponents: '
 print*, lambda(1)
 print*, lambda(2)
 print*, lambda(3)
 print*, "D_kap_yor=", 2.d0+(lambda(3)+lambda(3))/abs(lambda(1))
  
 close(194)

end program lyapunov 
