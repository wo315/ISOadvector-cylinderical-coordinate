subroutine com_advection
    use par
    implicit none
    real(kind=dp)::UF,Phi(4),PhiHR,flux(4)
    real(kind=dp)::Aface
    real(kind=dp),pointer::Ua(:,:)
    real(kind=dp),pointer,dimension(:,:)::Xfluxr,Xfluxz
    real(kind=dp),pointer,dimension(:,:)::Zfluxr,Zfluxz
    
    RS=1
    RE=NR-1
    
    ZS=1
    ZE=NZ-1
    
    allocate(Xfluxr(RS:RE,ZS:ZE),Xfluxz(RS:RE,ZS:ZE))
    allocate(Zfluxr(RS:RE,ZS:ZE),Zfluxz(RS:RE,ZS:ZE))
    
    do i=RS,RE
        do k=ZS,ZE
            Xfluxr(i,k)=0.
            Zfluxr(i,k)=0.
            Xfluxz(i,k)=0.
            Zfluxz(i,k)=0.
        end do
    end do
    
     !..............................
    Ua=>Uo
    
    RS=2
    RE=NR-1
    
    ZS=1
    ZE=NZ-1
     !....................

    do i=RS,RE
        do k=ZS,ZE
            Aface=RC(i-1)*dz
            if(i==RS) then
                Uf=0.5_dp*(Uo(i,k)+Uo(i-1,k))
                Phi(2)=Ua(i-1,k)
                Phi(3)=Ua(i,k)
                call upwind(Uf,Phi,PhiHR)
                Xfluxr(i,k)=Aface*Uf*PhiHR
            else
                Uf=0.5_dp*(Uo(i,k)+Uo(i-1,k))
                Phi(1)=Ua(i-2,k)
                Phi(2)=Ua(i-1,k)
                Phi(3)=Ua(i,k)
                Phi(4)=Ua(i+1,k)
                call HRscheme(Uf,Phi,PhiHR)
                Xfluxr(i,k)=Aface*Uf*PhiHR
            end if
          !..................X-monment Z-direction
            Aface=RP(i)*dr
            if(k==ZS) then
                Uf=0.5_dp*(Wo(i,k)+Wo(i-1,k))
                Xfluxz(i,k)=Aface*Uf*(Ua(i,k)+Ua(i,k-1))*0.5_dp  !(U(i,k)+U(i,k-1))*0.5_dp 边界上速度值
            else if(k==ZS+1) then
                Uf=0.5_dp*(Wo(i,k)+Wo(i-1,k))
                Phi(2)=Ua(i,k-1)
                Phi(3)=Ua(i,k)
                call upwind(Uf,Phi,PhiHR)
                Xfluxz(i,k)=Aface*Uf*PhiHR
            else
                Uf=0.5_dp*(Wo(i,k)+Wo(i-1,k))
                Phi(1)=Ua(i,k-2)
                Phi(2)=Ua(i,k-1)
                Phi(3)=Ua(i,k)
                Phi(4)=Ua(i,k+1)
                call HRscheme(Uf,Phi,PhiHR)
                Xfluxz(i,k)=Aface*Uf*PhiHR
            end if
        end do
    end do 

    !...................compute the advection of X-monmuent
    do i=RS,RE
        do k=ZS,ZE
            flux(1)=Xfluxr(i,k)
            if(i==RE) then
                Aface=RC(i)*dz
                Uf=0.5_dp*(Uo(i,k)+Uo(i+1,k))
                Phi(2)=Ua(i,k)
                Phi(3)=Ua(i+1,k)
                call upwind(Uf,Phi,PhiHR)
                flux(2)=Aface*Uf*PhiHR
            else
                flux(2)=Xfluxr(i+1,k)
            end if
            
            flux(3)=Xfluxz(i,k)
            if(k==ZE) then
                Aface=RP(i)*dr
                Uf=0.5_dp*(Wo(i,k+1)+Wo(i-1,k+1))
                flux(4)=Aface*Uf*(Ua(i,k)+Ua(i,k+1))*0.5_dp
            else
                flux(4)=Xfluxz(i,k+1)
            end if
            radv(i,k)=flux(4)-flux(3)+flux(2)-flux(1)
        end do
    end do

    !....................................................................    
    Ua=>Wo
    
    RS=1
    RE=NR-1

    
    ZS=2
    ZE=NZ-1          
       !................................................
    do i=RS,RE
        do k=ZS,ZE
            Aface=RP(i)*dz
            if(i==RS) then
                Uf=0.5_dp*(Uo(i,k)+Uo(i,k-1))
                Zfluxr(i,k)=Aface*Uf*(Ua(i,k)+Ua(i-1,k))*0.5_dp
            else if(i==RS+1) then
                Uf=0.5_dp*(Uo(i,k)+Uo(i,k-1))
                Phi(2)=Ua(i-1,k)
                Phi(3)=Ua(i,k)
                call upwind(Uf,Phi,PhiHR)
                Zfluxr(i,k)=Aface*Uf*PhiHR
            else
                Uf=0.5_dp*(Uo(i,k)+Uo(i,k-1))
                Phi(1)=Ua(i-2,k)
                Phi(2)=Ua(i-1,k)
                Phi(3)=Ua(i,k)
                Phi(4)=Ua(i+1,k)
                call HRscheme(Uf,Phi,PhiHR)
                Zfluxr(i,k)=Aface*Uf*PhiHR
            end if

             !..................Z-monment Z-direction
            Aface=RC(i)*dr
            if(k==ZS) then
                Uf=0.5_dp*(Wo(i,k)+Wo(i,k-1))
                Phi(2)=Ua(i,k-1)
                Phi(3)=Ua(i,k)
                call upwind(Uf,Phi,PhiHR)
                Zfluxz(i,k)=Aface*Uf*PhiHR
            else
                Uf=0.5_dp*(Wo(i,k)+Wo(i,k-1))
                Phi(1)=Ua(i,k-2)
                Phi(2)=Ua(i,k-1)
                Phi(3)=Ua(i,k)
                Phi(4)=Ua(i,k+1)
                call HRscheme(Uf,Phi,PhiHR)
                Zfluxz(i,k)=Aface*Uf*PhiHR
            end if
        end do
    end do

    !...............................................
    do i=RS,RE
        do k=ZS,ZE
            flux(1)=Zfluxr(i,k)
            if(i==RE) then
                Aface=RP(i+1)*dz
                Uf=0.5_dp*(Uo(i+1,k)+Uo(i+1,k-1))
                flux(2)=Aface*Uf*(Ua(i+1,k)+Ua(i,k))*0.5_dp
            else
                flux(2)=Zfluxr(i+1,k)
            end if
            
            flux(3)=Zfluxz(i,k)
            if(k==ZE) then
                Aface=RC(i)*dr
                Uf=0.5_dp*(Wo(i,k+1)+Wo(i,k))
                Phi(2)=Ua(i,k)
                Phi(3)=Ua(i,k+1)
                call upwind(Uf,Phi,PhiHR)
                flux(4)=Aface*Uf*PhiHR
            else
                flux(4)=Zfluxz(i,k+1)
            end if
            zadv(i,k)=flux(4)-flux(3)+flux(2)-flux(1)
        end do
    end do

    deallocate(Xfluxr,Xfluxz)
    deallocate(Zfluxr,Zfluxz)
    return
    end subroutine com_advection  

subroutine com_diffusion
    use par 
    implicit none
    
    real(kind=dp)::dif(4),gradu,denm,vism,Aface
    real(kind=dp),pointer::Ua(:,:),Den(:,:),Vis(:,:)
    
    real(kind=dp),pointer,dimension(:,:)::Xdifr,Xdifz
    real(kind=dp),pointer,dimension(:,:)::Zdifr,Zdifz
    
    RS=1
    RE=NR-1
    
    ZS=1
    ZE=NZ-1
    allocate(Xdifr(RS:RE,ZS:ZE),Xdifz(RS:RE,ZS:ZE))
    allocate(Zdifr(RS:RE,ZS:ZE),Zdifz(RS:RE,ZS:ZE))
    do i=RS,RE
            do k=ZS,ZE
                Xdifr(i,k)=0_dp
                Xdifz(i,k)=0_dp
                Zdifr(i,k)=0_dp
                Zdifz(i,k)=0_dp
            end do
    end do


   
 !remember that in the X and Z direction,the i,do not match with the physical property and Qua,so when treat with the density and viscosity and the Gradu
 !we should be careful that those values' tag is ahead the i 1
    !......................................................  
 !density and viscosity has additional cells,so don't need to consider the cycle at S_direction
 !Quo,Qvo,Qwo don't has additional cells,so need to sondier the cycle at S-Direction
 !Gradu fallow the normal direction of the face
    Ua=>Uo
    Den=>Denc
    Vis=>Visc
    
    RS=2
    RE=NR-1
    
    ZS=1
    ZE=NZ-1
 
    do i=RS,RE
            do k=ZS,ZE
             !..................X-monment R-direction
                Aface=RC(i-1)*dz
                Vism=0.125_dp*(Vis(i,k)+Vis(i,k)+Vis(i,k+1)+Vis(i,k+1)+&
                        Vis(i-1,k)+Vis(i-1,k)+Vis(i-1,k+1)+Vis(i-1,k+1))
                Gradu=(Ua(i,k)-Ua(i-1,k))/dr
                Xdifr(i,k)=Vism*Gradu*Aface
                 !..................X-monment Z-direction
                Aface=RP(i)*dr
                Vism=0.5_dp*(Vis(i,k)+Vis(i,k))
                Gradu=(Ua(i,k)-Ua(i,k-1))/dz
                Xdifz(i,k)=Vism*Gradu*Aface
            end do
    end do

    do i=RS,RE
        do k=ZS,ZE
            dif(1)=Xdifr(i,k)
            if(i==RE) then
                Aface=RC(i)*dz
                Vism=0.125_dp*(Vis(i+1,k)+Vis(i+1,k)+Vis(i+1,k+1)+Vis(i+1,k+1)+&
                        Vis(i,k)+Vis(i,k)+Vis(i,k+1)+Vis(i,k+1))
                Gradu=(Ua(i+1,k)-Ua(i,k))/dr
                dif(2)=Vism*Gradu*Aface
            else
                dif(2)=Xdifr(i+1,k)
            end if
               
            dif(3)=Xdifz(i,k)
            if(k==ZE) then
                Aface=RP(i)*dr
                Vism=0.5_dp*(Vis(i,k+1)+Vis(i,k+1))
                Gradu=(Ua(i,k+1)-Ua(i,k))/dz
                dif(4)=Vism*Gradu*Aface
            else
                dif(4)=Xdifz(i,k+1)
            end if
            rdif(i,k)=dif(4)-dif(3)+dif(2)-dif(1)
        end do
    end do

     !...................................................................................................................          
    !compute the diffusion term of Z-monment
    Ua=>Wo
    Den=>Denc
    Vis=>Visc
    
    RS=1
    RE=NR-1
    
    ZS=2
    ZE=NZ-1
        
    do i=RS,RE
            do k=ZS,ZE
                 !..................Z-monment R-direction
                Aface=RP(i)*dz
                Vism=0.5_dp*(Vis(i,k)+Vis(i,k))
                Gradu=(Ua(i,k)-Ua(i-1,k))/dr
                Zdifr(i,k)=Vism*Gradu*Aface
                !..................Z-monment Z-direction
                Aface=RC(i)*dr
                Vism=0.125_dp*(Vis(i,k)+Vis(i,k)+Vis(i,k-1)+Vis(i,k-1)+&
                        Vis(i+1,k)+Vis(i+1,k-1)+Vis(i+1,k)+Vis(i+1,k-1))
                Gradu=(Ua(i,k)-Ua(i,k-1))/dz
                Zdifz(i,k)=Vism*Gradu*Aface
            end do
    end do

    do i=RS,RE
        do k=ZS,ZE
            dif(1)=Zdifr(i,k)
            if(i==RE) then
                Aface=RP(i+1)*dz
                Vism=0.5_dp*(Vis(i+1,k)+Vis(i+1,k))
                Gradu=(Ua(i+1,k)-Ua(i,k))/dr
                dif(2)=Vism*Gradu*Aface
            else
                dif(2)=Zdifr(i+1,k)
            end if
            dif(3)=Zdifz(i,k)
            if(k==ZE) then
                Aface=RC(i)*dr
                Vism=0.125_dp*(Vis(i,k)+Vis(i,k)+Vis(i,k+1)+Vis(i,k+1)+&
                        Vis(i+1,k)+Vis(i+1,k+1)+Vis(i+1,k)+Vis(i+1,k+1))
                Gradu=(Ua(i,k+1)-Ua(i,k))/dz
                dif(4)=Vism*Gradu*Aface
            else
                dif(4)=Zdifz(i,k+1)
            end if
            zdif(i,k)=dif(4)-dif(3)+dif(2)-dif(1)
        end do
    end do
    deallocate(Xdifr,Xdifz)
    deallocate(Zdifr,Zdifz)
    return
    end subroutine com_diffusion


subroutine com_source
    use par 
    implicit none
    real(kind=dp)::DpDn,DuDz,DuDr,DwDz,DwDr
    real(kind=dp),pointer::Ua(:,:),Den(:,:),Vis(:,:)
    real(kind=dp)::Vism,Denm,forceif,visaff

    
    call acceleration(time,Agx,Agz)
    
    Ua=>Uo

    Den=>Denc
    Vis=>Visc
    
    RS=2
    RE=NR-1
    
    
    ZS=1
    ZE=NZ-1
     do i=RS,RE
            do k=ZS,ZE
                Denm=0.25_dp*(Den(i,k)+Den(i,k)+Den(i,k+1)+Den(i,k+1))
                Vism=0.25_dp*(Vis(i,k)+Vis(i,k)+Vis(i,k+1)+Vis(i,k+1))
                DpDn=(Po(i,k)-Po(i-1,k))/dr
                Uap(i,k)=-Vism/RP(i)/RP(i)*Uo(i,k)
                forceif=Rcsf(i,k)-Dpdn
                !viscosity influence of pressure jump at the interface
                DwDr=(Wo(i,k+1)+Wo(i,k)-Wo(i-1,k+1)-Wo(i-1,k))/2.0/dr
                DwDz=(Wo(i-1,k+1)+Wo(i,k+1)-Wo(i-1,k)-Wo(i-1,k+1))/2.0/dz
                DuDr=-Uo(i,k)/RP(i)-DwDz
                visaff=(Wvis-Avis)*(F(i,k+1)+F(i-1,k+1)-F(i,k-1)-F(i-1,k-1))/4.0/dz*DwDr
                visaff=visaff+(Wvis-Avis)*(F(i,k)-F(i-1,k))/dr*DuDr
                if(abs(Rcsf(i,k))<1.e-12) then
                    visaff=0.
                end if 
                !if(abs(Rcsf(i,k))>1.e-12) then
                    !forceif=forceif-max(min(forceif,0.01*Rcsfavge),-0.01*Rcsfavge)
                !end if 
               
                Rsource(i,k)=forceif/Denm+visaff/Denm+Uap(i,k)+Agx
            end do
     end do
    
    
    
    RS=1
    RE=NR-1
    
    
    ZS=2
    ZE=NZ-1
    
    do i=RS,RE
            do k=ZS,ZE
                 !..................Z-monment
                Denm=0.25_dp*(Den(i,k)+Den(i+1,k)+Den(i,k)+Den(i+1,k))
                DpDn=(Po(i,k)-Po(i,k-1))/dz
                forceif=Zcsf(i,k)-DpDn
                !viscosity jump at the interface influence 
                DuDz=(Uo(i,k)+Uo(i+1,k)-Uo(i,k-1)-Uo(i+1,k-1))/dz/2.0
                DuDr=(Uo(i+1,k)+Uo(i+1,k-1)-Uo(i,k)-Uo(i,k-1))/dr/2.0
                DwDz=-(Uo(i,k)+Uo(i+1,k)+Uo(i,k-1)+Uo(i+1,k-1))/4.0/RC(i)-DuDr
                visaff=(Wvis-Avis)*(F(i,k)-F(i,k-1))/dz*DwDz
                Visaff=visaff+(Wvis-Avis)*(F(i+1,k)+F(i+1,k-1)-F(i-1,k)-F(i-1,k-1))/dr/4.0*DuDz
                if(abs(Zcsf(i,k))<1.e-12) then
                    visaff=0.
                end if 
                !..........................................
                !if(abs(Zcsf(i,k))>1.e-12) then
                    !forceif=forceif-max(min(forceif,0.01*Zcsfavge),-0.01*Zcsfavge)
                !end if 
                
                Zsource(i,k)=forceif/Denm+visaff/Denm+Agz
            end do
    end do

    
    return
    
    end subroutine com_source