subroutine cal_Um
    use par
    implicit none
    real(kind=dp)::vp
    !x-direction
    !excpt volecity tag is alone,the adv,diff,Source,Uap,Vap term's tag is along with the cell tag
    call property
    call csf
    call com_advection
    call com_diffusion
    call com_source
  

    RS=2
    RE=NR-1
    
    
    ZS=1
    ZE=NZ-1
    !...........
     do i=RS,RE
            do k=ZS,ZE
                vp=RP(i)*dr*dz
                if(cyc<3) then
                    Um(i,k)=Uo(i,k)+dt*(Rdif(i,k)-Radv(i,k))/vp+Rsource(i,k)*dt
                else
                    Um(i,k)=Uo(i,k)+dt*(1.5_dp*Rdif(i,k)-0.5_dp*Rdifo(i,k)-(1.5_dp*Radv(i,k)-0.5_dp*Radvo(i,k)))/vp+&
                   Rsource(i,k)*dt
                end if
                Radvo(i,k)=Radv(i,k)
                Rdifo(i,k)=Rdif(i,k)
                Rsourceo(i,k)=Rsource(i,k)
            end do
     end do
     
    
    
    RS=1
    RE=NR-1
    
    
    ZS=2
    ZE=NZ-1
    do i=RS,RE
            do k=ZS,ZE
                vp=RC(i)*dr*dz
                if(cyc<3) then
                    Wm(i,k)=Wo(i,k)+dt*(Zdif(i,k)-Zadv(i,k))/vp+Zsource(i,k)*dt
                else
                    Wm(i,k)=Wo(i,k)+dt*(1.5_dp*Zdif(i,k)-0.5_dp*Zdifo(i,k)-(1.5_dp*Zadv(i,k)-0.5_dp*Zadvo(i,k)))/vp+&
                    Zsource(i,k)*dt
                end if
                Zadvo(i,k)=Zadv(i,k)
                Zdifo(i,k)=Zdif(i,k)
                Zsourceo(i,k)=Zsource(i,k)
            end do
    end do
    
    call boundary_Um
    
    return 
    end subroutine cal_Um     

subroutine vel
   use par
   implicit none
   real(kind=dp)::Denm,Poo
   real(kind=dp),pointer::Den(:,:)
   
   call possion !cofficient
   
   call solver_possion
   
   Den=>Denc

   !.....................
   RS=2
   RE=NR-1
   ZS=1
   ZE=NZ-1
   

   
   do i=RS,RE
          do k=ZS,ZE
              Denm=0.25_dp*(Den(i,k)+Den(i,k)+Den(i,k+1)+Den(i,k+1))
              U(i,k)=Um(i,k)-dt*(Pm(i,k)-Pm(i-1,k))/dr/Denm
          end do
   end do
   

   
   RS=1
   RE=NR-1

   ZS=2
   ZE=NZ-1
    do i=RS,RE
          do k=ZS,ZE
              Denm=(Den(i,k)+Den(i,k)+Den(i+1,k)+Den(i+1,k))*0.25_dp
              W(i,k)=Wm(i,k)-dt*(Pm(i,k)-Pm(i,k-1))/dz/Denm
          end do
    end do
   
       
   RS=1
   RE=NR-1

   ZS=1
   ZE=NZ-1


   do i=RS,RE
          do k=ZS,ZE
              P(i,k)=Po(i,k)+Pm(i,k)!
          end do
   end do
   
   
   call boundary_velocity
   
   
   return 
   
   end subroutine vel 