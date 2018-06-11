subroutine property
    use par 
    implicit none
    real(kind=dp),pointer,dimension(:,:)::Dich,Visko
    allocate(Dich(0:NR,0:NZ),Visko(0:NR,0:NZ)) !ÑØ×Å¸÷¸ö·½ÏòÀ©ÕÅ1¸öÍø¸ñ
    do i=0,NR
            do k=0,NZ
                Dich(i,k)=0.
                Visko(i,k)=0.
            end do
    end do
    
    do i=0,NR
            do k=0,NZ
                !»ñµÃËùÓÐÖÐÐÄ´¦ÎïÐÔ²ÎÊý£¬some additional cells
                !to get the physical properties of central of the cell,I add additional cell
                Dich(i,k)=Wdic*F(i,k)+Adic*(1.-F(i,k))! Dichte means density
                Visko(i,k)=Wvis/Wdic*F(i,k)+Avis/Adic*(1.-F(i,k))!visko means viscosity
            end do
    end do
    
    !»ñÈ¡½ÚµãÏà½»´¦µÄÎïÐÔ²ÎÊý
    do i=1,NR
            do k=1,NZ
                Denc(i,k)=0.125_dp*(Dich(i,k)+Dich(i,k)+Dich(i,k-1)+Dich(i,k-1)+&
                Dich(i-1,k)+Dich(i-1,k)+Dich(i-1,k-1)+Dich(i-1,k-1))!ÏòÐ¡µÄ·½ÏòÇÐÈë45¶È(i,k) and (i-1,k-1)¶Ô½ÇÏßÁ¬³ÉµÄÕý·½Ìå
                Visc(i,k)=0.125_dp*(Visko(i,k)+Visko(i,k)+Visko(i,k-1)+&
                Visko(i,k-1)+Visko(i-1,k)+Visko(i-1,k)+Visko(i-1,k-1)+Visko(i-1,k-1))
            end do
    end do
    
    deallocate(Dich,Visko)
    return
    end subroutine property