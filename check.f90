subroutine check_time
    use par
    implicit none
    real(kind=dp)::vp
    real(kind=dp)::Rmit,Zmit,Gast,capt
    nout=nout+1
    
    !time step restrict
    dt=mindt
    
    Cdmax=0._dp
    Rmit=mindt
    Zmit=mindt
    
    RS=1
    ZS=1
    RE=NR-1
    ZE=NZ-1	
    do i=RS,RE
            do k=ZS,ZE
                Rmit=min(Rmit,abs(dr/(1.e-30_dp+Uo(i,k))))
                Zmit=min(Zmit,abs(dz/(1.e-30_dp+Wo(i,k))))
            end do
    end do
    Gast=Adic/Avis*min(dr,dz)**2
    capt=sqrt((Wdic+Adic)/2.0*min(dr,dz)**3/pi/epiron)
    
    dt=0.5_dp*min(Rmit,Zmit,mindt,Gast,capt)
    
    RS=1
    ZS=1
    RE=NR-1
    ZE=NZ-1	
    do i=RS,RE
            do k=ZS,ZE
                vp=RC(i)*dr*dz
                Cdmax=max(Cdmax,max(-U(i,k)*RP(i)*dz/vp,0._dp)+max(U(i+1,k)*RP(i+1)*dz/vp,0._dp)+&
                    max(-W(i,k)*RC(i)*dr/vp,0._dp)+max(W(i,k+1)*RC(i)*dr/vp,0._dp))
            end do
    end do
    

    dt=min(dt,0.1_dp/(Cdmax+1.e-10_dp))

    Cdmax=Cdmax*dt
    time=time+dt
    Cyc=Cyc+1

    return
    end subroutine check_time

subroutine check_ns
    use par
    implicit none
    real(kind=dp)::residual,res,vp,residualo
    integer::io,ko,iter
    residual=0._dp
    RS=1
    ZS=1
    RE=NR-1
    ZE=NZ-1 
    do i=RS,RE
            do k=ZS,ZE
                vp=RC(i)*dr*dz
                res=U(i+1,k)*RP(i+1)*dz-U(i,k)*RP(i)*dz+&
                        W(i,k+1)*RC(i)*dr-W(i,k)*RC(i)*dr
                res=abs(res)/vp
                if(residual<res) then
                    residual=res
                    io=i
                    ko=k
                end if
            end do
    end do
    iter=1
    do while(residual>mass.and.iter<=itermaxzhu)!when residual<mass stop and.iter<=5

        do i=RS,RE
                do k=ZS,ZE
                    Um(i,k)=U(i,k)
                    Wm(i,k)=W(i,k)
                end do
        end do
        
        call boundary_Um
        
        call vel
        
        residual=1e-10_dp
        RS=1
        ZS=1
        RE=NR-1
        ZE=NZ-1 
        do i=RS,RE
                do k=ZS,ZE
                    vp=RC(i)*dr*dz
                    res=U(i+1,k)*RP(i+1)*dz-U(i,k)*RP(i)*dz+&
                        W(i,k+1)*RC(i)*dr-W(i,k)*RC(i)*dr
                    res=abs(res)/vp

                    if(residual<res) then
                        residual=res
                        io=i
                        ko=k
                    end if
                end do
        end do
        iter=iter+1
    end do
       
    write(*,20) nout,dt,time,Cdmax,residual

 20   FORMAT(I12,6X,4(F15.10, 4X))
    if(residual>mass) then
        write(*,*)'this procession is not converge'
        write(*,*)'io,ko is',io,ko
    end if
    return
    end subroutine check_ns
