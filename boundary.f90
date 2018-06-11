subroutine  boundary_Um
    use par
    implicit none
    RS=1
    RE=NR-1
    
    ZS=1
    ZE=NZ-1

    do k=ZS,ZE
        i=RS   
        if(BCXLEF==1.or.BCXLEF==2.or.BCXLEF==3) then
            Um(i,k)=0._dp
        else if(BCXLEF==4) then 
            Um(i,k)=Um(i+1,k)
        end if
        i=RE+1
        if(BCXRIG==1.or.BCXRIG==2) then
            Um(i,k)=0._dp
        else if(BCXRIG==4) then
            Um(i,k)=Um(i-1,k)*RP(i-1)/RP(i)
        end if
    end do

    RS=1
    RE=NR-1
    
    ZS=1
    ZE=NZ-1

    do i=RS,RE
            k=ZS
            if(BCYBOT==1.or.BCYBOT==2) then
                Wm(i,k)=0._dp
            else if(BCYBOT==4) then
                Wm(i,k)=Wm(i,k+1)
            end if
            k=ZE+1
            if(BCYTOP==1.or.BCYTOP==2) then
                Wm(i,k)=0._dp
            else if(BCYTOP==4) then
                Wm(i,k)=Wm(i,k-1)
            end if
    end do
    return
    end subroutine boundary_Um
    
subroutine boundary_f  !½ü±ÚÃæ0ÌÝ¶È¼ÙÉè
    use par
    implicit none
    real(kind=dp),pointer,dimension(:,:)::FI
    FI=>F
    RS=1
    RE=NR-1

    ZS=1
    ZE=NZ-1
    
    do k=ZS,ZE
    do i=-2,0
        if(BCXLEF==3) then
            FI(i,k)=FI(RS-i,k)
        else
            FI(i,k)=FI(RS,k)
        end if 
    enddo
    enddo
    
    do k=ZS,ZE
        do i=RE+1,RE+3
            if(BCXRIG==4) then
                if(U(NR,k)>0.) then
                    FI(i,k)=FI(RE,k)
                else
                    FI(i,k)=0.
                end if 
            else 
                 FI(i,k)=FI(RE,k)
            end if
        end do
    end do
    
    do i=-2,RE+3  !full extend is used 
            do k=-2,0
                FI(i,k)=FI(i,ZS)
            end do
            do k=ZE+1,ZE+3
                if(BCYTOP==4) then
                    if(W(i,NZ)>0.) then
                        FI(i,k)=FI(i,ZE)
                    else
                        FI(i,k)=0.
                    end if 
                else 
                    FI(i,k)=FI(i,ZE)
                end if
            end do
    end do
    return
    end subroutine boundary_f

subroutine boundary_FN
    use par 
    implicit none

    RS=1
    RE=NR-1

    ZS=1
    ZE=NZ-1

    do k=ZS,ZE
            !two additional values
        do i=-3,-1
            if(BCXLEF==3) then
                FN(RS+i,k)=FN(RS-i-1,k)
            else
                FN(RS+i,k)=FN(RS,k)
            end if

            if(BCXRIG==4) then
                if(U(NR,k)>0.) then
                    FN(RE-i,k)=FN(RE,k)
                else
                    FN(RE-i,k)=0.
                end if 
            else 
                 FN(RE-i,k)=FN(RE,k)
            end if
        end do
    end do
    
    do i=RS-3,RE+3
        do k=-3,-1
            FN(i,ZS+k)=FN(i,ZS)

            if(BCYTOP==4) then
                    if(W(i,NZ)>0.) then
                        FN(i,ZE-k)=FN(i,ZE)
                    else
                        FN(i,ZE-k)=0.
                    end if 
                else 
                   FN(i,ZE-k)=FN(i,ZE)
                end if
        end do
    end do
    return
    end subroutine boundary_FN


subroutine boundary_cofficient_FN  !use the zero gradient of F on the boundary
    use par 
    implicit none 
    RS=1
    RE=NR-1

    
    ZS=1
    ZE=NZ-1
    do k=ZS,ZE
        i=RS
        FAP(i,k)=FAP(i,k)+FAW(i,k)
        FAW(i,k)=0._dp

        i=RE
        FAP(i,k)=FAP(i,k)+FAE(i,k)
        FAE(i,k)=0._dp
    end do
    
    do i=RS,RE
        k=ZS
        FAP(i,k)=FAP(i,k)+FAB(i,k)
        FAB(i,k)=0._dp

        k=ZE
        FAP(i,k)=FAP(i,k)+FAT(i,k)   !fix boundary conditions
        FAT(i,k)=0._dp
    end do 

    return 

    end subroutine boundary_cofficient_FN
 
subroutine boundary_possion  
  use par
  implicit none
    RS=1
    RE=NR-1

    
    ZS=1
    ZE=NZ-1
    do k=ZS,ZE
        i=RS
        if(BCXLEF==1.or.BCXLEF==2.or.BCXLEF==3) then
            PAP(i,k)=PAP(i,k)+PAW(i,k)
            PAW(i,k)=0._dp
        else if(BCXLEF==4) then
            PAP(i,k)=PAP(i,k)+PAW(i,k)+PAE(i,k)
            PAW(i,k)=0._dp
            PAE(i,k)=0._dp
        end if
        i=RE
        if(BCXRIG==1.or.BCXRIG==2) then
            PAP(i,k)=PAP(i,k)+PAE(i,k)
            PAE(i,k)=0._dp
        else if(BCXRIG==4) then
            PAP(i,k)=PAP(i,k)+PAE(i,k)+PAW(i,k)
            PAE(i,k)=0._dp  
            PAW(i,k)=0._dp
        end if
    end do
    
    do i=RS,RE
        k=ZS
        if(BCYBOT==1.or.BCYBOT==2) then
            PAP(i,k)=PAP(i,k)+PAB(i,k)
            PAB(i,k)=0._dp
        else if(BCYBOT==4) then
            PAP(i,k)=PAP(i,k)+PAB(i,k)+PAT(i,k)
            PAB(i,k)=0._dp
            PAT(i,k)=0._dp
        end if
        k=ZE
        if(BCYTOP==1.or.BCYTOP==2) then
            PAP(i,k)=PAP(i,k)+PAT(i,k)   !fix boundary conditions
            PAT(i,k)=0._dp
        else if(BCYTOP==4) then
            PAP(i,k)=PAP(i,k)+PAB(i,k)+PAT(i,k)
            PAB(i,k)=0._dp
            PAT(i,k)=0._dp
        end if
    end do
    
    return
    end subroutine boundary_possion


    
    
!ËÙ¶È±ß½çÌõ¼þÉèÖÃ 
subroutine boundary_velocity    !no slip condition½ü±ÚÃænormalÎª0£¬ÇÐÏß·½ÏòÎªÌÝ¶ÈÎª0
    use par
    implicit none

    RS=1
    RE=NR-1

    ZS=1
    ZE=NZ-1
    !ghost cell method
    if(BCXLEF==1) then        !no slip
        do k=1,NZ
            i=RS
            U(i,k)=0._dp 

            U(i-1,k)=-U(i+1,k)    !no slip
            W(i-1,k)=-W(i  ,k)    !no slip,W'=-W
        enddo
    else if(BCXLEF==2) then  !slip
        do k=1,NZ
            i=RS
            U(i,k)=0._dp   

            U(i-1,k)=-U(i+1,k)  !slip,du/dr=0 
            W(i-1,k)= W(i  ,k) !slip W'=W
        enddo
    else if(BCXLEF==3) then !axisymetric
        do k=1,NZ
            i=RS
            U(i,k)=0._dp

            U(i-1,k)=-U(i+1,k)  !du/dr=0 require by cyclinder coordinate
            W(i-1,k)=W(i,k)    !dw/dr=0
        end do
    else if(BCXLEF==4) then   !open boundary
        do k=1,NZ
            i=RS
            U(i,k)=U(i+1,k)
            W(i-1,k)=W(i,k)   !dw/dr=0
        end do
    endif

    if(BCXRIG==1) then     !no slip
      do k=1,NZ   
            i=RE+1
            U(i,k)=0._dp
            W(i,k)=-W(i-1,k)

            U(i+1,k)=-U(i-1,k)
     !out boundary 
     enddo
    elseif(BCXRIG==2) then   !slip
      do k=1,NZ   
           i=RE+1
            U(i,k)=0._dp
            W(i,k)=W(i-1,k)

            U(i+1,k)=-U(i-1,k)
     enddo
    else if(BCXRIG==4) then  !zero_gradient
         do k=1,NZ   
            i=RE+1
            U(i,k)=U(i-1,k)*RP(i-1)/RP(i) !use zero flux boundary to replace the zero gradient boundary condition
            W(i,k)=W(i-1,k)

           ! U(i+1,k)=U(i-1,k) do not need the value at NR+1
        enddo
       
    endif

    if(BCYBOT==1) then     !no slip condition
        do i=1,NR
            k=ZS
            W(i,k)=0._dp

            W(i,k-1)=-W(i,k+1) 
            U(i,k-1)=-U(i,k)
        end do
    else if(BCYBOT==2) then
        do i=1,NR
            k=ZS
            W(i,k)=0._dp
            
            W(i,k-1)=-W(i,k+1) 
            U(i,k-1)=U(i,k)
        end do
    else if(BCYBOT==4) then
            do i=1,NR
                k=ZS
                W(i,k)=W(i,k+1)
                U(i,k-1)=U(i,k)
            end do
    endif


    if(BCYTOP==1) then  !no slip condition
      do i=1,NR
        k=ZE+1
        W(i,k)=0._dp
        U(i,k)=-U(i,k-1)

        W(i,k+1)=-W(i,k-1)
      end do
    else    if(BCYTOP==2) then 
      do i=1,NR
        k=ZE+1
        W(i,k)=0._dp
        U(i,k)=U(i,k-1)

        W(i,k+1)=-W(i,k-1)
      end do

    else   if(BCYTOP==4)  then
        do i=1,NR
            k=ZE+1
            W(i,k)=W(i,k-1)
            U(i,k)=U(i,k-1)
      end do
    endif

   return
    end subroutine boundary_velocity
