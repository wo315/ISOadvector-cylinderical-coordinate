subroutine vof_explicit
	use par
	implicit none
	call vof_cfl
	call vof_no
	call vof_co
	call vof_update
	return
	end subroutine vof_explicit

subroutine vof_cfl   !cf
    use par
    implicit none
    
    real(kind=dp)::vp
    RS=1
    RE=NR-1
    ZS=1
    ZE=NZ-1
    !....................................................................
    do i=0,RE
            do k=ZS,ZE
                vp=RC(i)*dr*dz
                Cfl(i,k)=max(U(i,k)*RP(i)*dz/vp*dt,0._dp)+max(-U(i+1,k)*RP(i+1)*dz/vp*dt,0._dp)+&
                    max(W(i,k)*RC(i)*dr/vp*dt,0._dp)+max(-W(i,k+1)*RC(i)*dr/vp*dt,0._dp)
            end do
    end do


    return
    end subroutine vof_cfl

subroutine vof_no   !F½ü±ÚÃæ0ÌÝ¶È¼ÙÉè
    use par
    implicit none
    real(kind=dp)::Afw,Afe,Afb,Aft,Ffw,Ffe,Ffb,Fft
    real(kind=dp),pointer::Fnode(:,:)
    real(kind=dp)::vp,mold,Rcomponent,Zcomponent
    
    allocate(Fnode(1:NR,1:NZ))
    do i=1,NR
            do k=1,NZ
                Fnode(i,k)=0.
            end do
    end do
    
    do i=-1,NR+1
        do k=-1,NZ+1
            FN(i,k)=F(i,k)
        end do 
    end do 

    call smooth_f

    do i=1,NR 
        do k=1,NZ 
            Fnode(i,k)=0.25*(FN(i,k)+FN(i-1,k)+FN(i,k-1)+FN(i-1,k-1))
        end do 
    end do 

    RS=1
    ZS=1 
    RE=NR-1
    ZE=NZ-1

    do i=RS,RE 
        do k=ZS,ZE
            vp=RC(i)*dr*dz
            Afw=RP(i)*dz
            Afe=RP(i+1)*dz
            Afb=RC(i)*dr
            Aft=RC(i)*dr

            Ffw=0.5_dp*(Fnode(i,k)+Fnode(i,k+1))
            Ffe=0.5_dp*(Fnode(i+1,k)+Fnode(i+1,k+1))
            Ffb=0.5_dp*(Fnode(i,k)+Fnode(i+1,k))
            Fft=0.5_dp*(Fnode(i,k+1)+Fnode(i+1,k+1))
            Rcomponent=(Ffe-Ffw)/dr
            Zcomponent=(FFt-Ffb)/dz

            mold=sqrt(Rcomponent**2+Zcomponent**2)
            Nor(i,k)=Rcomponent/(mold+1.0e-20)
            Noz(i,k)=Zcomponent/(mold+1.0e-20)
        end do 
    end do 
             
    !boundary   Nor,Nos£¬Noz in these additional cells will be use
     do k=ZS,ZE
            i=RS-1 !symetric
            Nor(i,k)=-Nor(RS,k)
            i=RE+1
            Nor(i,k)=Nor(RE,k)
    end do
    
    do i=RS-1,RE+1
            k=ZS-1
            Noz(i,k)=Noz(i,ZS)!Noz(i,ZS)
            k=ZE+1
            Noz(i,k)=Noz(i,ZE)
    end do
    
    
    deallocate(Fnode)
    
    return
    end  subroutine vof_no

subroutine vof_co
    use par
    implicit none
    real(kind=dp)::Uf,Fi(4),No(2),Cf(2),Fcoe
    RS=1
    RE=NR-1

    ZS=1
    ZE=NZ-1
    do k=ZS,ZE  !i=1ÒÔ¼°i=NR´¦µÄÃæ¿ÉÒÔ²»Ëã
            do i=1,RE   !i=1 for symetric
                Uf=U(i,k)!UfÅÐ¶ÎË³·çÄæ·ç£¬¿ÉÒÔÓÃQu´úÈ¡
                Fi(1)=F(i-2,k)
                Fi(2)=F(i-1,k)
                Fi(3)=F(i,k)
                Fi(4)=F(i+1,k)
                No(1)=Nor(i-1,k)
                No(2)=Nor(i,k)
                Cf(1)=Cfl(i-1,k)
                Cf(2)=Cfl(i,k)
                call vof_HRmo(Uf,Fi,Cf,No,Fcoe)
                Coer(i,k)=Fcoe
            end do
    end do
    
    
    do i=RS,RE
            do k=2,ZE !k==1,and k==NZ
                Uf=W(i,k)
                Fi(1)=F(i,k-2)
                Fi(2)=F(i,k-1)
                Fi(3)=F(i,k)
                Fi(4)=F(i,k+1)
                No(1)=Noz(i,k-1)
                No(2)=Noz(i,k)
                Cf(1)=Cfl(i,k-1)
                Cf(2)=Cfl(i,k)
                call vof_HRmo(Uf,Fi,Cf,No,Fcoe)
                Coez(i,k)=Fcoe
            end do
    end do
    
    return
    end subroutine vof_co

subroutine vof_update  
    use par 
    implicit none
    integer::vofstep,it
    real(kind=dp)::err,vp
    real(kind=dp)::Ff,Ffo
    real(kind=dp)::Fluw,Flue,Flut,Flub
    real(kind=dp)::Ffe,Ffw,Fft,Ffb,Ffet,Ffwt,Ffbt,Fftt
    real(kind=dp)::Artiw,Artie,Artib,Artit
    real(kind=dp)::FD,FA,FDO,FAO
    real(kind=dp)::divv,temp,dtt
    real(kind=dp)::deltaF,express,co
    real(kind=dp)::DfnDr,DfnDz,mo
    real(kind=dp),pointer,dimension(:,:)::Nar,Naz
    allocate(Nar(1:NR,1:NZ),Naz(1:NR,1:NZ))
    vofstep=5
    
    dtt=2._dp/3._dp*dt
    RS=1
    RE=NR-1
    
    ZS=1
    ZE=NZ-1
    do i=RS,RE
            do k=ZS,ZE
                Coero(i,k)=Coer(i,k)
                Coezo(i,k)=Coez(i,k)
            end do
    end do

    do i=1,NR
            do k=1,NZ
                Nar(i,k)=0.
                Naz(i,k)=0.
            end do
    end do

     !.........................dual time step
    it=1
    err=1._dp
    do while(it<vofstep.and.err>1.e-6)
        !prepare FN at pesudao time
        do i=-1,NR+1  
            do k=-1,NZ+1
                FN(i,k)=F(i,k)
            end do
        end do
        call smooth_f
     !.........................
     !compute DFNDR/|SXX|
         do i=RS,NR
            do k=ZS,ZE
                DfnDr=(FN(i,k)-FN(i-1,k))/dr
                DfnDz=0.25_dp*(FN(i,k+1)-FN(i,k-1)+&
                    FN(i-1,k+1)-FN(i-1,k-1))/dz
                mo=sqrt(DfnDr*DfnDr+DfnDz*DfnDz+1.e-20)
                Nar(i,k)=DFnDr/mo
            end do
        end do
        do k=ZS,NZ
            do i=RS,RE
                DfnDz=(FN(i,k)-FN(i,k-1))/dz
                DfnDr=(FN(i+1,k)-FN(i-1,k)+&
                    FN(i+1,k-1)-FN(i-1,k-1))/dr*0.25_dp
                mo=sqrt(DfnDr*DfnDr+DfnDz*DfnDz+1.e-20)
                Naz(i,k)=DFnDz/mo
            end do
        end do

        !.........................
        it=it+1
        err=0._dp
        RS=1
        RE=NR-1
        ZS=1
        ZE=NZ-1
        !correction processÃ¿¸öÑ­»·ÐÞÕýFD,FA,bf
        !check F's value
        !compute F,and correct F on the surface
        do i=RS,RE
            do k=ZS,ZE
                !.......................................w face
                Fluw=-U(i,k)*RP(i)*dz
                if(-Fluw>=0._dp) then
                    FD=F(i-1,k)
                    FA=F(i,k)
                    FDO=FO(i-1,k)
                    FAO=FO(i,k)
                else
                    FA=F(i-1,k)
                    FD=F(i,k)
                    FAO=FO(i-1,k)
                    FDO=FO(i,k)
                end if
                !symetric in include
                Ff=(1._dp-Coer(i,k))*FD+Coer(i,k)*FA
                Ffo=(1._dp-Coero(i,k))*FDO+Coero(i,k)*FAO
                Ffw=0.5_dp*(Fluw*Ff+Fluw*Ffo)
                Ffwt=Ff*(1.0_dp-Ff)*abs(Fluw)
                !.......................................e face
                Flue=U(i+1,k)*RP(i+1)*dz
                if(Flue>=0._dp) then
                    FD=F(i,k)
                    FA=F(i+1,k)
                    FDO=FO(i,k)
                    FAO=FO(i+1,k)
                else
                    FA=F(i,k)
                    FD=F(i+1,k)
                    FAO=FO(i,k)
                    FDO=FO(i+1,k)
                end if
                     
                if(i==RE) then
                    Ff=0.5_dp*(F(i,k)+F(i+1,k))
                    Ffo=0.5_dp*(Fo(i,k)+Fo(i+1,k))
                    Ffe=0.5_dp*(Flue*Ff+Flue*Ffo)
                    Ffet=Ff*(1.0_dp-Ff)*abs(Flue)
                else
                    Ff=(1._dp-Coer(i+1,k))*FD+Coer(i+1,k)*FA
                    Ffo=(1._dp-Coero(i+1,k))*FDO+Coero(i+1,k)*FAO
                    Ffe=0.5_dp*(Flue*Ff+Flue*Ffo)
                    Ffet=Ff*(1.0_dp-Ff)*abs(Flue)
                end if
                !....................................... b face
                Flub=-W(i,k)*RC(i)*dr
                if(-Flub>0._dp) then
                    FD=F(i,k-1)
                    FA=F(i,k)
                    FDO=FO(i,k-1)
                    FAO=FO(i,k)
                else
                    FD=F(i,k)
                    FA=F(i,k-1)
                    FDO=FO(i,k)
                    FAO=FO(i,k-1)
                end if
                     
                if(k==ZS) then
                    Ff=0.5_dp*(F(i,k)+F(i,k-1))
                    Ffo=0.5_dp*(Fo(i,k)+Fo(i,k-1))
                    Ffb=0.5_dp*(Flub*Ff+Flub*Ffo)
                    Ffbt=Ff*(1.0_dp-Ff)*abs(Flub)
                else
                    Ff=(1._dp-Coez(i,k))*FD+Coez(i,k)*FA
                    Ffo=(1._dp-Coezo(i,k))*FDO+Coezo(i,k)*FAO
                    Ffb=0.5_dp*(Flub*Ff+Flub*Ffo)
                    Ffbt=Ff*(1.0_dp-Ff)*abs(Flub)
                end if
                !....................................... t face
                Flut=W(i,k+1)*RC(i)*dr
                if(Flut>0._dp) then
                    FD=F(i,k)
                    FA=F(i,k+1)
                    FDO=FO(i,k)
                    FAO=FO(i,k+1)
                else
                    FD=F(i,k+1)
                    FA=F(i,k)
                    FDO=FO(i,k+1)
                    FAO=FO(i,k)
                end if
                        
                if(k==ZE) then
                    Ff=0.5_dp*(F(i,k)+F(i,k+1))
                    Ffo=0.5_dp*(F(i,k)+F(i,k+1))
                    Fft=0.5_dp*(Flut*Ff+Flut*Ffo)
                    Fftt=Ff*(1.0_dp-Ff)*abs(Flut)
                else 
                    Ff=(1._dp-Coez(i,k+1))*FD+Coez(i,k+1)*FA
                    Ffo=(1._dp-coezo(i,k+1))*FDO+Coezo(i,k+1)*FAO
                    Fft=0.5_dp*(Flut*Ff+Flut*Ffo)
                    Fftt=Ff*(1.0_dp-Ff)*abs(Flut)
                end if

                !solve the vof equationµÚÒ»´Î¸üÐÂ
                vp=Rc(i)*dr*dz
                Artiw=0.1_dp*Ffwt*Nar(i,k)
                Artie=0.1_dp*Ffet*Nar(i+1,k)
                Artib=0.1_dp*Ffbt*Naz(i,k)
                Artit=0.1_dp*Fftt*Naz(i,k+1)
                Temp=F(i,k)-(Ffw+Ffe+Ffb+Fft)*dtt/vp-(F(i,k)-FO(i,k))/dt*dtt-(Artie-Artiw+Artit-Artib)*dtt/vp
                TEMP=MIN(1.,MAX(TEMP,0.))
                err=max(err,abs(Temp-F(i,k)))
                !until the err is less than 1e-6,go to the next real time
                if(abs(temp-F(i,k))>1.e-6) then
                    F(i,k)=Temp
                else
                    F(i,k)=F(i,k)
                end if
                !Íê³ÉÖµµÄ¸üÐÂ£¬½øÐÐÐÞÕý
                !µÚ¶þ´Î¸üÐÂ(Êµ¼ÊÉÏÊÇÐÞÕý)
                do while(F(i,k)<0._dp.or.F(i,k)>1.0_dp)
                    !.......................................w face
                    Fluw=-U(i,k)*RP(i)*dz
                    if(-Fluw>=0._dp) then
                        FD=F(i-1,k)
                        FA=F(i,k)
                        FDO=FO(i-1,k)
                        FAO=FO(i,k)
                    else
                        FA=F(i-1,k)
                        FD=F(i,k)
                        FAO=FO(i-1,k)
                        FDO=FO(i,k)
                    end if
                        
                    deltaF=0.5_dp*(FA+FAO-FD-FDO)
                    if(F(i,k)<0._dp) then
                        express=max(-FD,0._dp)
                        FD=FD+express
                        FA=FA-express!FDO,FAO²»ÐÞÕý
                    else if(F(i,k)>1._dp) then
                        express=max(FD-1._dp,0._dp)
                        FD=FD-express
                        FA=FA+express
                    end if
                     
                    if(i==RS) then
                        Ffw=Fluw*0.5_dp*(F(i,k)+F(i-1,k))
                    else
                        if(F(i,k)<0._dp.or.F(i,k)>1.0_dp) then
                            call correct_co(deltaF,express,coer(i,k),Cfl(i,k),F(i,k),co)!bf=coe
                            Coer(i,k)=Coer(i,k)-co
                            Coero(i,k)=Coer(i,k)
                        else
                            Coer(i,k)=Coer(i,k)
                            Coero(i,k)=Coer(i,k)
                        end if
                        Ffw=0.5_dp*((1._dp-Coer(i,k))*FD+Coer(i,k)*FA+&
                        (1._dp-Coero(i,k))*FDO+Coero(i,k)*FAO)*Fluw
                    end if
                    !.......................................e face
                    Flue=U(i+1,k)*RP(i+1)*dz
                    if(Flue>=0._dp) then
                        FD=F(i,k)
                        FA=F(i+1,k)
                        FDO=FO(i,k)
                        FAO=FO(i+1,k)
                    else
                        FA=F(i,k)
                        FD=F(i+1,k)
                        FAO=FO(i,k)
                        FDO=FO(i+1,k)
                    end if
                     
                    deltaF=0.5_dp*(FA+FAO-FD-FDO)
                    if(F(i,k)<0._dp) then
                        express=max(-FD,0._dp)
                        FD=FD+express
                        FA=FA-express
                    else if(F(i,k)>1._dp) then
                        express=max(FD-1._dp,0._dp)
                        FD=FD-express
                        FA=FA+express
                    end if
                     
                    if(i==RE) then
                        Ffe=Flue*0.5_dp*(F(i,k)+F(i+1,k))
                    else
                        if(F(i,k)<0._dp.or.F(i,k)>1.0_dp) then
                            call correct_co(deltaF,express,Coer(i+1,k),Cfl(i,k),F(i,k),co)
                            Coer(i+1,k)=Coer(i+1,k)-co
                            Coero(i+1,k)=Coer(i+1,k)
                        else
                            Coer(i+1,k)=Coer(i+1,k)
                            Coero(i+1,k)=Coer(i+1,k)
                        end if
                        Ffe=0.5_dp*((1._dp-Coer(i+1,k))*FD+Coer(i+1,k)*FA+&
                        (1._dp-Coero(i+1,k))*FDO+Coero(i+1,k)*FAO)*Flue
                    end if
                    !....................................... b face
                    Flub=-W(i,k)*RC(i)*dr
                    if(-Flub>0._dp) then
                        FD=F(i,k-1)
                        FA=F(i,k)
                        FDO=FO(i,k-1)
                        FAO=FO(i,k)
                    else
                        FD=F(i,k)
                        FA=F(i,k-1)
                        FDO=FO(i,k)
                        FAO=FO(i,k-1)
                    end if
                     
                    deltaF=0.5_dp*(FA+FAO-FD-FDO)
                    if(F(i,k)<0._dp) then
                        express=max(-FD,0._dp)
                        FD=FD+express
                        FA=FA-express
                    else if(F(i,k)>1._dp) then
                        express=max(FD-1._dp,0._dp)
                        FD=FD-express
                        FA=FA+express
                    end if
                     
                    if(k==ZS) then
                        Ffb=Flub*0.5_dp*(F(i,k)+F(i,k-1))
                    else
                        if(F(i,k)<0._dp.or.F(i,k)>1.0_dp) then
                            call correct_co(deltaF,express,coez(i,k),Cfl(i,k),F(i,k),co)
                            Coez(i,k)=Coez(i,k)-co
                            Coezo(i,k)=Coez(i,k)
                        else
                            Coez(i,k)=Coez(i,k)
                            Coezo(i,k)=Coez(i,k)
                        end if
                        Ffb=0.5_dp*((1._dp-Coez(i,k))*FD+Coez(i,k)*FA+&
                        (1._dp-Coezo(i,k))*FDO+Coezo(i,k)*FAO)*Flub
                    end if
                    !....................................... t face
                    Flut=W(i,k+1)*RC(i)*dr
                    if(Flut>0._dp) then
                        FD=F(i,k)
                        FA=F(i,k+1)
                        FDO=FO(i,k)
                        FAO=FO(i,k+1)
                    else
                        FD=F(i,k+1)
                        FA=F(i,k)
                        FDO=FO(i,k+1)
                        FAO=FO(i,k)
                    end if
                        
                    deltaF=0.5_dp*(FA+FAO-FD-FDO)
                    if(F(i,k)<0._dp) then
                        express=max(-FD,0._dp)
                        FD=FD+express
                        FA=FA-express
                    else if(F(i,k)>1._dp) then
                        express=max(FD-1._dp,0._dp)
                        FD=FD-express
                        FA=FA+express
                    end if
                        
                    if(k==ZE) then
                        Fft=Flut*0.5_dp*(F(i,k)+F(i,k+1))  
                    else 
                        if(F(i,k)<0._dp.or.F(i,k)>1.0_dp) then
                            call correct_co(deltaF,express,coez(i,k+1),Cfl(i,k),F(i,k),co)
                            Coez(i,k+1)=Coez(i,k+1)-co
                            Coezo(i,k+1)=Coez(i,k+1)
                        else
                            Coez(i,k+1)=Coez(i,k+1)
                            Coezo(i,k+1)=Coez(i,k+1)
                        end if
                        Fft=0.5_dp*((1._dp-Coez(i,k+1))*FD+Coez(i,k+1)*FA+&
                        (1._dp-coezo(i,k+1))*FDO+Coezo(i,k+1)*FAO)*Flut
                    end if
                    !ÔÙ´Î¸üÐÂFÖµ
                    vp=RC(i)*dr*dz
                    Temp=F(i,k)-(Ffw+Ffe+Ffb+Fft)*dtt/vp-(F(i,k)-FO(i,k))/dt*dtt
                    !until the err is less than 1e-6,go to the next real time
                    err=max(err,abs(Temp-F(i,k)))
                    if(abs(temp-F(i,k))>1.e-6) then
                        F(i,k)=Temp
                    else
                        F(i,k)=F(i,k)
                    end if
                end do!½áÊøÐÞÕý
            end do
        end do
        call boundary_f
        call vof_no
        call vof_co
         
    end do!end of do while

    deallocate(Nar,Naz)
    return
    
    end subroutine vof_update

subroutine vof_HR(Ufin,Fiin,Cfin,Noin,Fcoeout)
    implicit none
    integer,parameter::dp=selected_real_kind(p=15)
    real(kind=dp),intent(in)::Ufin,Fiin(4),Noin(2),Cfin(2)
    real(kind=dp),intent(out)::Fcoeout
    real(kind=dp)::FUin,FDin,FAin,anglein,courantin
    real(kind=dp)::FDNin,FCBCin,FUQin
    if(Ufin>0._dp) then
        FUin=Fiin(1)
        FDin=Fiin(2)
        FAin=Fiin(3)
        anglein=Noin(1)
        courantin=Cfin(1)
    else
        FUin=Fiin(4)
        FDin=Fiin(3)
        FAin=Fiin(2)
        anglein=Noin(2)
        courantin=Cfin(2)
    end if
    
    if(abs(FAin-FUin)>1.E-20) then
        FDNin=(FDin-FUin)/(FAin-FUin)
    else
        FDNin=0._dp !´ËÊ±FfÈ¡FU£¬FD£¬FAÃ»ÓÐÇø±ð
    end if
    !CBC value
    if(FDNin>=0._dp.and.FDNin<=1._dp.and.abs(courantin)>1.e-20_dp) then
        FCBCin=min(FDNin/(courantin+1.e-30_dp),1._dp)
    else
        FCBCin=FDNin
    end if
    !UQ value
    if(FDNin>=0._dp.and.FDNin<=1._dp) then
        FUQin=min((8._dp*courantin*FDNin+(1._dp-courantin)*(6._dp*FDNin+3._dp))/8._dp,FCBCin)
    else
        FUQin=FDNin
    end if
    anglein=abs(anglein)
    anglein=acos(anglein)
    anglein=min((cos(2._dp*anglein)+1._dp)*0.5_dp,1._dp)
    Fcoeout=anglein*FCBCin+(1._dp-anglein)*FUQin
    
    if(abs(1._dp-FDNin)>1.e-20_dp) then
        Fcoeout=(Fcoeout-FDNin)/(1._dp-FDNin)
    else
        Fcoeout=0._dp
    end if
    
    Fcoeout=min(max(Fcoeout,0._dp),1._dp)
    return
    end subroutine vof_HR

subroutine vof_HRmo(Ufin,Fiin,Cfin,Noin,Fcoeout)
    implicit none
    integer,parameter::dp=selected_real_kind(p=15)
    real(kind=dp),intent(in)::Ufin,Fiin(4),Noin(2),Cfin(2)
    real(kind=dp),intent(out)::Fcoeout
    real(kind=dp)::FUin,FDin,FAin,anglein,courantin
    real(kind=dp)::FDNin,FCBCin,FUQin
    if(Ufin>0._dp) then
        FUin=Fiin(1)
        FDin=Fiin(2)
        FAin=Fiin(3)
        anglein=Noin(1)
        courantin=Cfin(1)
    else
        FUin=Fiin(4)
        FDin=Fiin(3)
        FAin=Fiin(2)
        anglein=Noin(2)
        courantin=Cfin(2)
    end if
    
    if(abs(FAin-FUin)>1.E-20) then
        FDNin=(FDin-FUin)/(FAin-FUin)
    else
        FDNin=0._dp !´ËÊ±FfÈ¡FU£¬FD£¬FAÃ»ÓÐÇø±ð
    end if
    !CBC value
    if(FDNin>=0._dp.and.FDNin<=1._dp.and.abs(courantin)>1.e-20_dp) then
        FCBCin=min(2.0*FDNin,1._dp)
    else
        FCBCin=FDNin
    end if
    !UQ value
    if(FDNin>=0._dp.and.FDNin<=1._dp) then
        FUQin=min(0.25_dp+FDNin,FCBCin)
    else
        FUQin=FDNin
    end if
    anglein=abs(anglein)
    anglein=anglein**0.25
    Fcoeout=anglein*FCBCin+(1._dp-anglein)*FUQin
    
    if(abs(1._dp-FDNin)>1.e-20_dp) then
        Fcoeout=(Fcoeout-FDNin)/(1._dp-FDNin)
    else
        Fcoeout=0._dp
    end if
    
    Fcoeout=min(max(Fcoeout,0._dp),1._dp)
    return
end subroutine vof_HRmo
                         
    
subroutine correct_co(apha,ep,bf,cf,Fr,coff)
    implicit none
    integer,parameter::dp=selected_real_kind(p=15)
    real(kind=dp),intent(in)::apha,ep,bf,cf,Fr
    real(kind=dp),intent(out)::coff
    if(Fr<0._dp) then
        if(apha>ep) then
            coff=min(bf,ep*(2.0_dp+cf-2.0_dp*cf*bf)*0.5_dp/(cf+1.0e-30)/(apha-ep+1.0e-30))
        else
            coff=0._dp
        end if
    else if(Fr>1._dp)then
        if(apha<-ep) then
            coff=min(bf,ep*(2.0_dp+cf-2.0_dp*cf*bf)*0.5_dp/(cf+1.0e-30)/(-apha-ep+1.0e-30))
        else
            coff=0._dp
        end if
    end if
    end subroutine correct_co


    subroutine smooth_f
    use par
    implicit none
    integer::it,itm
    real(kind=dp)::Ffe,Ffw,Ffb,Fft
    real(kind=dp)::Afe,Afw,Afb,Aft
    real(kind=dp),pointer,dimension(:,:)::Ftem
    allocate(Ftem(1:NR-1,1:NZ-1))
    
    it=1
    itm=3
    
    do i=1,NR-1
            do k=1,NZ-1
                Ftem(i,k)=0_dp
            end do
    end do

    !......................................
    do while(it<itm)
        it=it+1 
        RS=1
        RE=NR-1
        ZS=1
        ZE=NZ-1 
        do i=RS,RE
                do k=ZS,ZE
                    !the i,k respent the cell number
                    !get the smooth F
                    Ffw=0.5_dp*(FN(i,k)+FN(i-1,k))
                    Ffe=0.5_dp*(FN(i,k)+FN(i+1,k))
                    Ffb=0.5_dp*(FN(i,k)+FN(i,k-1))
                    Fft=0.5_dp*(FN(i,k)+FN(i,k+1))
                    
                    Afw=dz*RP(i)
                    Afe=dz*RP(i+1)
                    Afb=dr*RC(i)
                    Aft=dr*RC(i)

                    Ftem(i,k)=(Ffe*Afe+Ffw*Afw+Ffb*Afb+Fft*Aft)
                    Ftem(i,k)=Ftem(i,k)/(Afe+Afw+Afb+Aft)
                end do
        end do
        
        RS=1
        RE=NR-1
        ZS=1
        ZE=NZ-1 
        do i=RS,RE
                do k=ZS,ZE      
                    FN(i,k)=Ftem(i,k)
                end do
        end do
        
    end do
    
    !............................................................................ 
    !boundary condition use to compute the additional csf 
    call boundary_FN
    
    deallocate(Ftem)
    
    return
    end subroutine smooth_f